#!/usr/bin/env python3
"""
tune_params.py
==============
**Fine tuning** de los parametros internos de CIRCUST-Python sobre datos
sinteticos de calibracion.

Se ejecuta *antes* que los experimentos comparativos del Capitulo 3: el
objetivo es fijar de forma justificada los parametros del algoritmo, para
que las diferencias que despues se atribuyan a las aportaciones
metodologicas no procedan en realidad de un parametro mal elegido.

Que se ajusta y que no
----------------------
Se distingue entre **parametros** y **ejes de comparacion**:

- Se ajustan aqui los parametros internos que no forman parte de ninguna
  aportacion: los umbrales que definen el conjunto TOP (etapa 3) y los de
  la estimacion robusta (etapa 4). Ver :data:`SWEEPS`.

- **No** se ajustan ``fmm_method``, ``robust_method``, ``aggregation`` ni
  ``consensus_scale``. Esos son precisamente las aportaciones que los
  experimentos E1 y E2 contrastan: "optimizarlos" aqui vaciaria de
  contenido esa comparacion. Durante el tuning se mantienen fijos en la
  configuracion mas proxima al procedimiento original (FMM base, best_k
  con mediana de R2), declarada en :data:`HELD_FIXED`.

  La unica excepcion es ``selection_r2_threshold`` (tau), que es un
  parametro *del* criterio ``n_well_fitted`` y no un eje de comparacion,
  de modo que se calibra aqui. Su barrido fuerza dicho criterio.

Separacion calibracion / evaluacion
-----------------------------------
Los escenarios y semillas de este modulo (:data:`CALIBRATION_SCENARIOS`,
:data:`CALIBRATION_SEEDS`) son **disjuntos** de los que usan los
experimentos comparativos. Ajustar y evaluar sobre los mismos datos
producira mejoras optimistamente sesgadas, exactamente igual que evaluar
un modelo sobre su conjunto de entrenamiento.

Estrategia de barrido
---------------------
OFAT (*one factor at a time*): se parte del punto central
:data:`DEFAULT_CENTER`, que son los valores heredados del articulo
original, y se mueve **un parametro cada vez** dejando el resto en el
centro. No se escalariza el resultado en un unico numero ni se elige un
optimo de forma automatica: se emite el perfil multi-metrica completo
(MedAE, SDAE, % dentro de 2 h, correlacion circular) y la eleccion se
hace leyendo la frontera. Para una segunda pasada alrededor de un centro
distinto basta con ``--center r2_fmm_threshold=0.6,n_reps=10``.

Coste
-----
El diseno explota la estructura de dependencias del pipeline:

- Las etapas 0 a 2 (preprocesado, genes core, CPCA, orden preliminar) son
  **invariantes** a todos los parametros ajustados aqui, asi que se
  calculan una sola vez por dataset y se reutilizan.
- Un cambio en un parametro TOP obliga a rehacer las etapas 3 y 4.
- Un cambio en un parametro robusto solo obliga a rehacer la etapa 4, de
  modo que todo el bloque robusto comparte un unico resultado de etapa 3.

Salidas
-------
``tuning_profile.csv``
    Formato largo, una fila por (bloque, parametro, valor, escenario,
    replica) con todas las metricas y los tiempos por etapa.

``tuning_summary.csv``
    Agregado por (bloque, parametro, valor): mediana y rango
    intercuartilico de cada metrica sobre los escenarios de calibracion.
    Es la tabla que se lee para decidir.

Uso
---
    # perfil completo (todos los bloques)
    python validation/tune_params.py -o output/tuning

    # solo el bloque TOP, 2 replicas, prueba rapida
    python validation/tune_params.py -o output/tuning_top \\
        --blocks top --reps 2 --quick

    # segunda pasada con un centro desplazado
    python validation/tune_params.py -o output/tuning_pass2 \\
        --center r2_fmm_threshold=0.6,n_reps=10
"""
from __future__ import annotations

import argparse
import json
import sys
import time
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))
sys.path.insert(0, str(PROJECT_ROOT / "validation"))

from circust.preprocessing import Preprocessor
from circust.cpca          import CPCA
from circust.synchronizer  import CircularSynchronizer
from circust.top_genes     import TopGeneSelector
from circust.robust_order  import RobustOrderEstimator
from circust.core_genes    import CoreGeneSelector

from generate_synthetic import generate, PERIOD_HOURS
from metrics            import evaluate_all


# ===========================================================================
# Configuracion invariante
# ===========================================================================

#: Parametros de las etapas 0 a 2. No se ajustan: el preprocesado y la
#: deteccion de atipicas son decisiones de formato y de control de calidad
#: del metodo original, y los genes ancla (ARNTL / DBP / CRY1) son una
#: hipotesis biologica, no un hiperparametro. Se replican tal cual de
#: ``scripts/run_pipeline.py`` para que el tuning corra el mismo pipeline.
FIXED_CONFIG: dict[str, Any] = {
    "sparse_threshold":       0.3,
    "n_outlier_candidates":   8,
    "tight_radius":           0.10,
    "loose_radius":           0.15,
    "outlier_fmm_threshold":  3.0,
    "max_outlier_fraction":   0.05,
    "fmm_alpha_grid":         48,
    "fmm_omega_grid":         24,
    "fmm_reps":               3,
    "anchor_gene":            "ARNTL",
    "direction_gene":         "DBP",
}

#: Ejes de comparacion, congelados durante el tuning. Se fijan en la
#: configuracion mas proxima al procedimiento original de Larriba et al.,
#: de modo que los parametros se calibren sobre el algoritmo de partida y
#: no sobre una de las variantes que luego hay que evaluar.
HELD_FIXED: dict[str, Any] = {
    "fmm_method":          "base",
    "robust_method":       "best_k",
    "selection_criterion": "r2_median",
}

#: Punto central del barrido OFAT: los valores actuales del algoritmo,
#: heredados del articulo original.
DEFAULT_CENTER: dict[str, Any] = {
    # --- bloque TOP (etapa 3) ---
    "r2_ori_threshold":       0.5,
    "r2_fmm_threshold":       0.5,
    "omega_min":              0.1,
    "n_sectors":              4,
    "top_fraction":           None,
    # --- bloque robusto (etapa 4) ---
    "n_reps":                 5,
    "sample_size_fraction":   2.0 / 3.0,
    "r2_min":                 0.5,
    "selection_r2_threshold": 0.65,
}

#: A que bloque pertenece cada parametro. Determina que etapas hay que
#: recalcular al moverlo, y por tanto el coste del barrido.
PARAM_BLOCK: dict[str, str] = {
    "r2_ori_threshold":       "top",
    "r2_fmm_threshold":       "top",
    "omega_min":              "top",
    "n_sectors":              "top",
    "top_fraction":           "top",
    "n_reps":                 "robust",
    "sample_size_fraction":   "robust",
    "r2_min":                 "robust",
    "selection_r2_threshold": "robust",
}

#: Rejilla de valores a barrer por parametro. Cada rejilla incluye el
#: valor por defecto, de forma que el centro aparezca siempre en el perfil
#: y sirva de referencia visual.
SWEEPS: dict[str, tuple] = {
    # --- bloque TOP ---
    # Filtro de entrada no parametrico. Bajarlo deja pasar mas genes al
    # ajuste FMM (mas coste, mas ruido), subirlo restringe el pool.
    "r2_ori_threshold":     (0.3, 0.4, 0.5, 0.6, 0.7),
    # Umbral que define el conjunto TOP propiamente dicho. Es el que mas
    # margen tiene: gobierna cuantos genes alimentan toda la etapa 4.
    "r2_fmm_threshold":     (0.3, 0.4, 0.5, 0.6, 0.7),
    # Descarte de perfiles 'spiked'. omega_min = 0 desactiva el filtro.
    "omega_min":            (0.0, 0.05, 0.1, 0.15, 0.2),
    # Sectores en que se reparte la cobertura circular de los picos.
    "n_sectors":            (2, 4, 6, 8),
    # Cuota de genes retenidos. None = sin cuota (todos los que pasan).
    "top_fraction":         (None, 0.25, 0.50, 0.75, 1.0),
    # --- bloque robusto ---
    # Numero K de repeticiones. El coste crece linealmente con K, asi que
    # su lectura es tanto de precision como de presupuesto.
    "n_reps":               (3, 5, 10, 15, 20),
    # Fraccion del conjunto TOP que entra en cada submuestra.
    "sample_size_fraction": (0.50, 2.0 / 3.0, 0.80, 0.90),
    # R2 minimo para entrar en el pool de candidatos de cada repeticion.
    "r2_min":               (0.3, 0.4, 0.5, 0.6, 0.7),
    # tau. Se incluye 0,5 solo para evidenciar la saturacion por
    # construccion descrita en la Seccion 2.2.2 de la memoria.
    "selection_r2_threshold": (0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90),
}

#: Rejillas reducidas para ``--quick`` (humo, no para resultados).
QUICK_SWEEPS: dict[str, tuple] = {
    "r2_ori_threshold":       (0.4, 0.5, 0.6),
    "r2_fmm_threshold":       (0.4, 0.5, 0.6),
    "omega_min":              (0.05, 0.1, 0.2),
    "n_sectors":              (2, 4, 8),
    "top_fraction":           (None, 0.5),
    "n_reps":                 (3, 5, 10),
    "sample_size_fraction":   (0.5, 2.0 / 3.0),
    "r2_min":                 (0.4, 0.5, 0.6),
    "selection_r2_threshold": (0.55, 0.65, 0.75, 0.85),
}


# ===========================================================================
# Escenarios de calibracion
# ===========================================================================

#: Numero total de genes de la matriz simulada.
#:
#: Se fija a escala transcriptomica realista y **no** al valor pequeno con
#: el que se prototipo el generador. El motivo no es de fidelidad estetica
#: sino de potencia del experimento: ``n_total`` acota por construccion el
#: tamano del conjunto TOP, y con un TOP pequeno las K repeticiones de la
#: etapa 4 dejan de ser distintas entre si. Medido sobre este pipeline
#: (muestreo equiespaciado, sigma = 1, n = 48, prop = 0,05):
#:
#:     n_total    genes ritmicos    nTOP    reps distintas de 10    MedAE
#:       1200            60           14           4                 2,13
#:       5000           250           20          10                 1,22
#:      12000           600           32          10                 0,96
#:      20000          1000           44          10                 1,04
#:
#: Con 1200 genes solo 4 de cada 10 repeticiones son distintas, de modo que
#: la coleccion de K ordenaciones degenera y el procedimiento de agregacion
#: de la Seccion 2.2.2 no tiene variabilidad sobre la que operar: el
#: experimento seria incapaz de detectar la aportacion que pretende medir.
#: El coste apenas cambia (2,8 s por dataset con 1200 genes frente a 3,0 s
#: con 20000), asi que no hay contrapartida.
N_TOTAL_REALISTIC: int = 20000

#: Escenarios sobre los que se calibra. Se eligen para **cubrir los ejes**
#: del espacio de simulacion sin ser el grid de evaluacion: un escenario
#: comodo, uno de muestreo irregular, uno de onda apuntada y uno de senal
#: escasa. Asi los valores elegidos no quedan ajustados a un unico regimen.
CALIBRATION_SCENARIOS: dict[str, dict[str, Any]] = {
    "base": dict(sampling="equispaced", omega=0.5, sigma=1.0,
                 prop_rhythmic=0.05, n_samples=48, n_total=N_TOTAL_REALISTIC),
    "irregular": dict(sampling="clustered_gap", omega=0.5, sigma=1.0,
                      prop_rhythmic=0.05, n_samples=48, n_total=N_TOTAL_REALISTIC),
    "apuntado": dict(sampling="equispaced", omega=0.1, sigma=1.0,
                     prop_rhythmic=0.05, n_samples=48, n_total=N_TOTAL_REALISTIC),
    "senal_debil": dict(sampling="equispaced", omega=0.5, sigma=0.5,
                        prop_rhythmic=0.01, n_samples=24, n_total=N_TOTAL_REALISTIC),
}

#: Semillas de calibracion. **Disjuntas** de las de evaluacion (que parten
#: de 2000 en el runner de experimentos), para que el tuning no contamine
#: la medida posterior de las aportaciones.
CALIBRATION_SEEDS: tuple[int, ...] = (1000, 1001, 1002, 1003, 1004)

#: Metricas que lideran la lectura del perfil, en orden de importancia.
#:
#: La concordancia se mide con ``rho`` (regresion circular-circular) y no
#: con el CCC de Jammalamadaka. El CCC degenera cuando las fases
#: verdaderas son uniformes sobre el ciclo, que es exactamente lo que
#: produce el muestreo equiespaciado: al anularse la resultante media, la
#: direccion que centra el coeficiente queda indeterminada. Medido sobre
#: este mismo pipeline, el CCC devuelve valores en torno a -0,9 sobre
#: reconstrucciones con error mediano de 30 minutos. ``rho`` depende solo
#: de los residuos del ajuste de rotacion y reflexion y se mantiene valida.
PROFILE_METRICS: tuple[str, ...] = (
    "MedAE", "SDAE", "pct_within_2h", "rho",
)


# ===========================================================================
# Contexto del pipeline: etapas 0 a 2, calculadas una sola vez
# ===========================================================================

@dataclass
class PipelineContext:
    """
    Resultado de las etapas 0 a 2 para un dataset concreto.

    Estas etapas no dependen de **ningun** parametro ajustado en este
    modulo (los umbrales del conjunto TOP entran en la etapa 3 y los de la
    estimacion robusta en la 4), de modo que se calculan una vez por
    dataset y se reutilizan a lo largo de todo el barrido. Es lo que hace
    viable el tuning.
    """
    scenario:     str
    seed:         int
    order_result: Any
    cpca_result:  Any
    core_found:   list[str]
    true_phase_h: pd.Series
    t_setup_s:    float


def build_context(scenario: str, cfg: dict, seed: int,
                  verbose: bool = False) -> PipelineContext:
    """
    Genera el dataset sintetico y ejecuta sobre el las etapas 0 a 2.

    El ``fmm_method`` usado en CPCA es el de :data:`HELD_FIXED`, no un
    parametro libre.
    """
    t0 = time.perf_counter()

    ds = generate(
        name          = f"{scenario}__seed{seed}",
        n_samples     = cfg["n_samples"],
        sampling      = cfg["sampling"],
        omega         = cfg["omega"],
        sigma         = cfg["sigma"],
        prop_rhythmic = cfg["prop_rhythmic"],
        n_total       = cfg["n_total"],
        seed          = seed,
    )

    prep = Preprocessor(
        sparse_threshold = FIXED_CONFIG["sparse_threshold"],
        verbose          = verbose,
    ).run(ds.expr)

    core_result = CoreGeneSelector.from_string(None).select(prep.expr_norm)

    cpca_result = CPCA(
        core_genes            = core_result.genes,
        n_outlier_candidates  = FIXED_CONFIG["n_outlier_candidates"],
        tight_radius          = FIXED_CONFIG["tight_radius"],
        loose_radius          = FIXED_CONFIG["loose_radius"],
        fmm_threshold         = FIXED_CONFIG["outlier_fmm_threshold"],
        max_outlier_fraction  = FIXED_CONFIG["max_outlier_fraction"],
        fmm_length_alpha_grid = FIXED_CONFIG["fmm_alpha_grid"],
        fmm_length_omega_grid = FIXED_CONFIG["fmm_omega_grid"],
        fmm_num_reps          = FIXED_CONFIG["fmm_reps"],
        fmm_method            = HELD_FIXED["fmm_method"],
        verbose               = verbose,
    ).run(prep.expr_norm)

    core_found = list(cpca_result.core_genes_found)

    order_result = CircularSynchronizer(
        anchor_gene    = FIXED_CONFIG["anchor_gene"],
        direction_gene = FIXED_CONFIG["direction_gene"],
        verbose        = verbose,
    ).run(cpca_result, core_found)

    # Ground truth plegado a [0, 24), indexado por nombre de muestra. Se
    # indexa por nombre y no por posicion porque la deteccion de atipicas
    # puede haber eliminado muestras.
    true_phase = pd.Series(
        np.asarray(ds.true_times.values, dtype=np.float64) % PERIOD_HOURS,
        index = ds.true_times.index,
    )

    return PipelineContext(
        scenario     = scenario,
        seed         = seed,
        order_result = order_result,
        cpca_result  = cpca_result,
        core_found   = core_found,
        true_phase_h = true_phase,
        t_setup_s    = time.perf_counter() - t0,
    )


# ===========================================================================
# Etapas 3 y 4 parametrizadas
# ===========================================================================

def run_top_stage(ctx: PipelineContext, params: dict, verbose: bool = False):
    """Etapa 3 con los parametros TOP dados. Devuelve ``(top_result, segundos)``."""
    t0 = time.perf_counter()
    top_result = TopGeneSelector(
        r2_ori_threshold      = params["r2_ori_threshold"],
        r2_fmm_threshold      = params["r2_fmm_threshold"],
        omega_min             = params["omega_min"],
        n_sectors             = params["n_sectors"],
        top_fraction          = params["top_fraction"],
        fmm_length_alpha_grid = FIXED_CONFIG["fmm_alpha_grid"],
        fmm_length_omega_grid = FIXED_CONFIG["fmm_omega_grid"],
        fmm_num_reps          = FIXED_CONFIG["fmm_reps"],
        fmm_method            = HELD_FIXED["fmm_method"],
        lambda_star           = ctx.cpca_result.lambda_star,
        n_jobs                = -1,
        verbose               = verbose,
    ).run(
        expr_norm      = ctx.order_result.expr_ordered,
        circular_scale = ctx.order_result.circular_scale,
        seed_genes     = ctx.core_found,
        sample_order   = ctx.order_result.sample_order,
    )
    return top_result, time.perf_counter() - t0


def run_robust_stage(ctx: PipelineContext, top_result, params: dict,
                     seed: int, verbose: bool = False):
    """
    Etapa 4 con los parametros robustos dados. Devuelve
    ``(phases_h_por_muestra, segundos)``, donde el primer elemento es una
    Serie indexada por nombre de muestra.

    El criterio de seleccion se toma de :data:`HELD_FIXED` salvo cuando se
    esta barriendo tau, en cuyo caso hay que forzar ``n_well_fitted``: tau
    es un parametro *de ese criterio* y con ``r2_median`` no tendria ningun
    efecto, de modo que el barrido saldria plano por construccion.
    """
    t0 = time.perf_counter()
    criterion = params.get("_selection_criterion", HELD_FIXED["selection_criterion"])

    robust_result = RobustOrderEstimator(
        n_reps                 = params["n_reps"],
        sample_size_fraction   = params["sample_size_fraction"],
        method                 = HELD_FIXED["robust_method"],
        selection_criterion    = criterion,
        selection_r2_threshold = params["selection_r2_threshold"],
        r2_min                 = params["r2_min"],
        anchor_gene            = FIXED_CONFIG["anchor_gene"],
        direction_gene         = FIXED_CONFIG["direction_gene"],
        fmm_length_alpha_grid  = FIXED_CONFIG["fmm_alpha_grid"],
        fmm_length_omega_grid  = FIXED_CONFIG["fmm_omega_grid"],
        fmm_num_reps           = FIXED_CONFIG["fmm_reps"],
        fmm_method             = HELD_FIXED["fmm_method"],
        lambda_star            = ctx.cpca_result.lambda_star,
        n_jobs                 = -1,
        seed                   = seed,
        verbose                = verbose,
    ).run(
        top_result     = top_result,
        expr_full_norm = ctx.order_result.expr_ordered,
        core_genes     = ctx.core_found,
    )

    # Misma convencion que scripts/run_pipeline.py: circular_scale[j] es la
    # fase de la muestra situada en la posicion j del orden estimado.
    names = [ctx.order_result.expr_ordered.columns[i]
             for i in robust_result.sample_order]
    phases_h = (np.asarray(robust_result.circular_scale, dtype=np.float64)
                * PERIOD_HOURS / (2.0 * np.pi))
    return pd.Series(phases_h, index=names), time.perf_counter() - t0


def score(ctx: PipelineContext, phases_h: pd.Series) -> dict:
    """
    Metricas de la estimacion frente al ground truth. Empareja por nombre
    de muestra, no por posicion, porque la etapa 1 puede haber eliminado
    atipicas.
    """
    common = [s for s in phases_h.index if s in ctx.true_phase_h.index]
    if len(common) < 3:
        return {"status": "muestras_insuficientes", "n": len(common)}
    res = evaluate_all(
        predicted = phases_h.loc[common].to_numpy(dtype=np.float64),
        true      = ctx.true_phase_h.loc[common].to_numpy(dtype=np.float64),
        period    = PERIOD_HOURS,
        align     = True,
    )
    res["status"] = "ok"
    return res


# ===========================================================================
# Barrido OFAT
# ===========================================================================

def _resolve_center(overrides: str | None) -> dict:
    """
    Punto central del barrido, con posibles sobreescrituras de la forma
    ``"r2_fmm_threshold=0.6,n_reps=10"``. Sirve para hacer una segunda
    pasada alrededor de un centro desplazado.
    """
    center = dict(DEFAULT_CENTER)
    if not overrides:
        return center
    for item in overrides.split(","):
        if "=" not in item:
            raise SystemExit(f"--center mal formado: {item!r}. Formato: clave=valor")
        key, raw = item.split("=", 1)
        key = key.strip()
        if key not in DEFAULT_CENTER:
            raise SystemExit(
                f"--center: parametro desconocido {key!r}. "
                f"Disponibles: {sorted(DEFAULT_CENTER)}"
            )
        raw = raw.strip()
        if raw.lower() in ("none", "null"):
            center[key] = None
        elif isinstance(DEFAULT_CENTER[key], int) and not isinstance(DEFAULT_CENTER[key], bool):
            center[key] = int(raw)
        else:
            center[key] = float(raw)
    return center


def _fmt(value: Any) -> str:
    """Representacion compacta de un valor de parametro para las tablas."""
    if value is None:
        return "None"
    if isinstance(value, float):
        return f"{value:.4g}"
    return str(value)


def sweep_context(
    ctx:      PipelineContext,
    center:   dict,
    sweeps:   dict[str, tuple],
    blocks:   set[str],
    verbose:  bool = False,
) -> list[dict]:
    """
    Barre todos los parametros solicitados sobre **un** dataset.

    Explota la estructura de dependencias: el bloque robusto comparte un
    unico resultado de etapa 3 (el del centro), mientras que cada valor de
    un parametro TOP obliga a recalcular las etapas 3 y 4. Los resultados
    de etapa 3 se cachean por configuracion, de modo que el valor central
    de cada parametro TOP se calcula una sola vez aunque aparezca en todas
    las rejillas.
    """
    rows: list[dict] = []
    top_cache: dict[tuple, tuple] = {}

    def _top_key(params: dict) -> tuple:
        return tuple(params[k] for k in
                     ("r2_ori_threshold", "r2_fmm_threshold",
                      "omega_min", "n_sectors", "top_fraction"))

    def _get_top(params: dict):
        key = _top_key(params)
        if key not in top_cache:
            top_cache[key] = run_top_stage(ctx, params, verbose=verbose)
        return top_cache[key]

    for param, values in sweeps.items():
        block = PARAM_BLOCK[param]
        if block not in blocks:
            continue

        for value in values:
            params = dict(center)
            params[param] = value
            # tau solo actua bajo el criterio al que pertenece.
            if param == "selection_r2_threshold":
                params["_selection_criterion"] = "n_well_fitted"

            row = {
                "block":      block,
                "param":      param,
                "value":      _fmt(value),
                "value_num":  (np.nan if value is None else float(value)),
                "is_center":  bool(value == center[param]),
                "scenario":   ctx.scenario,
                "seed":       ctx.seed,
                "t_setup_s":  ctx.t_setup_s,
            }

            try:
                top_result, t_top = _get_top(params)
                n_top = int(top_result.candidate_matrix.shape[0])
                if n_top < 3:
                    row.update(status="top_vacio", n_top=n_top, t_top_s=t_top)
                    rows.append(row)
                    continue

                phases, t_rob = run_robust_stage(
                    ctx, top_result, params, seed=ctx.seed, verbose=verbose)
                row.update(n_top=n_top, t_top_s=t_top, t_robust_s=t_rob)
                row.update(score(ctx, phases))

            except Exception as exc:                      # pragma: no cover
                # Una configuracion puede ser degenerada (por ejemplo un
                # umbral que vacia el conjunto TOP). Se registra y se sigue:
                # que una configuracion falle es informacion util del perfil.
                row.update(status=f"error: {type(exc).__name__}: {exc}")

            rows.append(row)

    return rows


# ===========================================================================
# Agregacion del perfil
# ===========================================================================

def summarise(profile: pd.DataFrame) -> pd.DataFrame:
    """
    Agrega el perfil por (bloque, parametro, valor) sobre los escenarios y
    replicas de calibracion.

    Se usan **mediana e IQR** y no media y desviacion tipica porque basta
    un escenario patologico para desplazar la media, y lo que interesa es
    el comportamiento tipico del parametro. Se reporta ademas el numero de
    ejecuciones validas, que delata las configuraciones degeneradas.
    """
    ok = profile[profile["status"] == "ok"]
    if ok.empty:
        return pd.DataFrame()

    metrics = [m for m in PROFILE_METRICS if m in ok.columns]
    extra   = [c for c in ("n_top", "t_top_s", "t_robust_s") if c in ok.columns]

    grouped = ok.groupby(["block", "param", "value", "value_num", "is_center"],
                         dropna=False)
    out = grouped.agg(
        n_runs = ("status", "size"),
        **{f"{m}_med": (m, "median") for m in metrics},
        **{f"{m}_iqr": (m, lambda s: s.quantile(0.75) - s.quantile(0.25))
           for m in metrics},
        **{f"{c}_med": (c, "median") for c in extra},
    ).reset_index()

    # Numero de ejecuciones que fallaron o degeneraron, por celda.
    total = (profile.groupby(["block", "param", "value"], dropna=False)
                    .size().rename("n_total").reset_index())
    out = out.merge(total, on=["block", "param", "value"], how="left")
    out["n_failed"] = out["n_total"] - out["n_runs"]

    return out.sort_values(["block", "param", "value_num"]).reset_index(drop=True)


def print_profile(summary: pd.DataFrame) -> None:
    """Vuelca el perfil por pantalla, un bloque de texto por parametro."""
    if summary.empty:
        print("\n  (sin ejecuciones validas)")
        return

    for (block, param), sub in summary.groupby(["block", "param"], sort=False):
        print()
        print(f"  {'=' * 74}")
        print(f"  [{block}] {param}")
        print(f"  {'-' * 74}")
        header = (f"  {'valor':>10} {'MedAE':>9} {'(IQR)':>8} {'SDAE':>9} "
                  f"{'%<2h':>7} {'rho':>7} {'nTOP':>6} {'t4 (s)':>8}  ")
        print(header)
        for _, r in sub.iterrows():
            mark = " <- centro" if r["is_center"] else ""
            print(f"  {r['value']:>10} "
                  f"{r.get('MedAE_med', float('nan')):>9.3f} "
                  f"{r.get('MedAE_iqr', float('nan')):>8.3f} "
                  f"{r.get('SDAE_med', float('nan')):>9.3f} "
                  f"{r.get('pct_within_2h_med', float('nan')):>7.1f} "
                  f"{r.get('rho_med', float('nan')):>7.3f} "
                  f"{r.get('n_top_med', float('nan')):>6.0f} "
                  f"{r.get('t_robust_s_med', float('nan')):>8.1f}"
                  f"{mark}")
        n_failed = int(sub["n_failed"].sum())
        if n_failed:
            print(f"  ({n_failed} ejecuciones fallidas o degeneradas en este parametro)")


# ===========================================================================
# CLI
# ===========================================================================

def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=("Fine tuning OFAT de los parametros internos de CIRCUST "
                     "sobre escenarios sinteticos de calibracion."),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("-o", "--output", required=True, type=str,
                   help="Directorio de salida.")
    p.add_argument("--blocks", default="top,robust",
                   help="Bloques a barrer, separados por coma: top, robust. "
                        "Defecto: ambos.")
    p.add_argument("--params", default=None,
                   help="Subconjunto de parametros a barrer (coma). Por "
                        "defecto todos los de los bloques elegidos.")
    p.add_argument("--scenarios", default=None,
                   help=f"Escenarios de calibracion (coma). Disponibles: "
                        f"{list(CALIBRATION_SCENARIOS)}. Defecto: todos.")
    p.add_argument("--reps", type=int, default=3,
                   help="Replicas (semillas) por escenario. Defecto: 3.")
    p.add_argument("--center", default=None,
                   help="Sobreescribe el centro del barrido, formato "
                        "'clave=valor,clave=valor'. Para segundas pasadas.")
    p.add_argument("--quick", action="store_true",
                   help="Rejillas reducidas. Prueba de humo, no para resultados.")
    p.add_argument("--verbose-stages", action="store_true",
                   help="Deja hablar a los modulos del pipeline (muy ruidoso).")
    return p.parse_args()


def main() -> None:
    args = _parse_args()
    out  = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)

    blocks = {b.strip() for b in args.blocks.split(",") if b.strip()}
    unknown = blocks - {"top", "robust"}
    if unknown:
        raise SystemExit(f"Bloques desconocidos: {sorted(unknown)}. Usa top y/o robust.")

    sweeps = dict(QUICK_SWEEPS if args.quick else SWEEPS)
    if args.params:
        wanted = {p.strip() for p in args.params.split(",") if p.strip()}
        unknown_p = wanted - set(sweeps)
        if unknown_p:
            raise SystemExit(f"Parametros desconocidos: {sorted(unknown_p)}")
        sweeps = {k: v for k, v in sweeps.items() if k in wanted}

    scen_names = ([s.strip() for s in args.scenarios.split(",")]
                  if args.scenarios else list(CALIBRATION_SCENARIOS))
    unknown_s = set(scen_names) - set(CALIBRATION_SCENARIOS)
    if unknown_s:
        raise SystemExit(f"Escenarios desconocidos: {sorted(unknown_s)}")

    if args.reps > len(CALIBRATION_SEEDS):
        raise SystemExit(
            f"--reps={args.reps} supera las semillas de calibracion "
            f"disponibles ({len(CALIBRATION_SEEDS)}). Amplia CALIBRATION_SEEDS."
        )
    seeds = CALIBRATION_SEEDS[:args.reps]

    center = _resolve_center(args.center)

    n_datasets = len(scen_names) * len(seeds)
    n_configs  = sum(len(v) for k, v in sweeps.items()
                     if PARAM_BLOCK[k] in blocks)

    print("=" * 78)
    print("  Fine tuning de parametros CIRCUST (OFAT sobre calibracion)")
    print("=" * 78)
    print(f"  Escenarios      : {scen_names}")
    print(f"  Semillas        : {list(seeds)}   (disjuntas de evaluacion)")
    print(f"  Datasets        : {n_datasets}")
    print(f"  Bloques         : {sorted(blocks)}")
    print(f"  Parametros      : {[k for k in sweeps if PARAM_BLOCK[k] in blocks]}")
    print(f"  Configuraciones : {n_configs} por dataset")
    print(f"  Ejes congelados : {HELD_FIXED}")
    print(f"  Centro          : {{{', '.join(f'{k}={_fmt(v)}' for k, v in center.items())}}}")
    print(f"  Salida          : {out.resolve()}")
    print("=" * 78)

    all_rows: list[dict] = []
    t_start = time.time()
    done = 0

    # numba y sklearn son ruidosos en este bucle y sus avisos no aportan.
    warnings.filterwarnings("ignore", category=RuntimeWarning)

    for scen in scen_names:
        cfg = CALIBRATION_SCENARIOS[scen]
        for seed in seeds:
            done += 1
            tag = f"[{done}/{n_datasets}] {scen} seed={seed}"
            print(f"\n{tag}  preparando etapas 0-2 ...", flush=True)

            try:
                ctx = build_context(scen, cfg, seed,
                                    verbose=args.verbose_stages)
            except Exception as exc:
                print(f"      FALLO en etapas 0-2: {type(exc).__name__}: {exc}")
                all_rows.append({"scenario": scen, "seed": seed,
                                 "status": f"setup_error: {exc}"})
                continue

            print(f"      etapas 0-2 en {ctx.t_setup_s:.1f}s, barriendo ...",
                  flush=True)
            rows = sweep_context(ctx, center, sweeps, blocks,
                                 verbose=args.verbose_stages)
            all_rows.extend(rows)

            n_ok = sum(1 for r in rows if r.get("status") == "ok")
            elapsed = time.time() - t_start
            eta = elapsed / done * (n_datasets - done)
            print(f"      {n_ok}/{len(rows)} configuraciones OK  "
                  f"| ETA {eta / 60:.1f} min")

            # Persistencia incremental: un barrido completo es largo y no
            # conviene perderlo por una interrupcion.
            pd.DataFrame(all_rows).to_csv(out / "tuning_profile.csv", index=False)

    profile = pd.DataFrame(all_rows)
    profile.to_csv(out / "tuning_profile.csv", index=False)

    summary = summarise(profile)
    summary.to_csv(out / "tuning_summary.csv", index=False)

    manifest = {
        "scenarios":       scen_names,
        "seeds":           list(seeds),
        "blocks":          sorted(blocks),
        "center":          {k: (None if v is None else v) for k, v in center.items()},
        "held_fixed":      HELD_FIXED,
        "fixed_config":    FIXED_CONFIG,
        "sweeps":          {k: list(v) for k, v in sweeps.items()
                            if PARAM_BLOCK[k] in blocks},
        "quick":           args.quick,
        "elapsed_min":     (time.time() - t_start) / 60.0,
    }
    with open(out / "manifest.json", "w") as f:
        json.dump(manifest, f, indent=2, default=str)

    print()
    print("=" * 78)
    print(f"  Perfil de calibracion  ({(time.time() - t_start) / 60:.1f} min)")
    print("=" * 78)
    print_profile(summary)
    print()
    print(f"  tuning_profile.csv : {len(profile)} filas (formato largo)")
    print(f"  tuning_summary.csv : {len(summary)} filas (mediana + IQR)")
    print(f"  Guardado en        : {out.resolve()}")
    print()
    print("  Lee el perfil y fija los valores a mano. No se elige optimo de")
    print("  forma automatica: con varias metricas sin escalarizar la eleccion")
    print("  es una decision, no un argmin.")


if __name__ == "__main__":
    main()
