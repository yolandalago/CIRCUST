"""
Ejecucion de CIRCUST sobre conjuntos sinteticos.

Es una copia de ``scripts/run_pipeline.py`` adaptada: sin CLI, sin escritura
a disco, sin impresion por pantalla, y partida en dos entradas para que el
bloque de agregacion pueda reutilizar las etapas comunes.

La configuracion del algoritmo **no se redeclara**, se importa de
``run_pipeline`` para que no pueda separarse en silencio de lo que se ejecuta
con datos reales.

Diferencias deliberadas frente a ``run_pipeline.py``
----------------------------------------------------
1. ``anchor_gene`` y ``direction_gene`` son los sinteticos, no ARNTL y DBP.
2. Los ``consistency_pairs`` se derivan de la verdad y se pasan **tambien** a
   la Etapa 4, que ``run_pipeline.py`` deja con los pares por defecto. Sin
   esto la Etapa 4 buscaria ARNTL, DBP y CRY1, que no existen aqui.
3. ``n_jobs=1`` en todas partes: se paraleliza a nivel de celda.
4. Se recogen cuatro indicadores observables sin la verdad, ver mas abajo.

El atlas
--------
Se construye con el mismo ``build_atlas`` de ``robust_order`` que usa
``run_pipeline``, no con una copia, porque el bloque de regularizacion lo
necesita: CircaFMM **no puede** cambiar el orden de muestras, que sale de una SVD sobre los
genes nucleares donde el FMM no interviene en ningun paso. Donde si cambia
las cosas es en los parametros por gen, y eso es lo que el atlas publica.
Comparar base contra CircaFMM sobre el orden seria compararlos donde son
identicos por construccion.

Indicadores observables
-----------------------
``error_orden_h`` y ``tau_circular`` necesitan los tiempos verdaderos, asi que
sirven para caracterizar el metodo en simulacion pero no para aplicarlo a un
tejido nuevo. Estos cuatro si se calculan sin la verdad, y el estudio de
simulacion sirve para calibrar su relacion con el error real:

``r2_anchor``     lo bien determinado que esta el pico que fija el origen.
``pares_ok``      fraccion de pares a priori que se cumplen en el marco final.
``evidencia``     suma de |sin(delta)| sobre esos pares, o sea cuanta
                  informacion sobre la orientacion habia disponible.
``desacuerdo_K``  dispersion entre las K repeticiones de la Etapa 4.
``n_top``         tamano del conjunto TOP, o sea cuanta senal ritmica habia.
"""
from __future__ import annotations

import sys
import time
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
for _extra in (PROJECT_ROOT / "src", PROJECT_ROOT / "scripts"):
    if str(_extra) not in sys.path:
        sys.path.insert(0, str(_extra))

from circust.core_genes    import CoreGeneSelector
from circust.cpca          import CPCA
from circust.preprocessing import Preprocessor
from circust.robust_order  import RobustOrderEstimator, RobustRepetitions, build_atlas
from circust.synchronizer  import CircularSynchronizer
from circust.top_genes     import TopGeneSelector

from run_pipeline import DEFAULT_CONFIG   # no se redeclara

from simulation.salidas import fases_por_muestra
from simulation.truth import TWO_PI, SimulationTruth

# Mismos valores que usa scripts/run_pipeline.py en la Etapa 3. No estan en
# DEFAULT_CONFIG porque alli son argumentos de la CLI.
R2_ORI_THRESHOLD       = 0.5
R2_FMM_THRESHOLD       = 0.5
OMEGA_MIN              = 0.1
N_SECTORS              = 8
TOP_FRACTION           = 0.005
SELECTION_R2_THRESHOLD = 0.5
SAMPLE_SIZE_FRACTION   = 2.0 / 3.0
# R2 minimo que exige el Paso 2.2 para que un par de consistencia cuente.
CONSISTENCY_R2_MIN     = 0.5


@dataclass(frozen=True)
class RobustConfig:
    """Configuracion de la Etapa 4, lo unico que cambia entre metodos."""

    robust_method:       str = "best_k"      # best_k|aggregate
    selection_criterion: str = "r2_median"   # r2_median|n_well_fitted|medoid
    aggregation:         str = "circmean"    # circmean|tsp|hod
    n_reps:              int = 5
    seed:                int = 42

    @property
    def consensus_scale(self) -> str:
        """Como se obtienen las fases del consenso, que ya no se elige aparte.

        Se deduce del metodo: la media circular da sus propias fases y TSP y
        Hodge las reconstruyen por CIRE bajo su orden. Se conserva como
        columna del CSV para que las salidas antiguas, en las que se elegia,
        sigan siendo comparables.
        """
        return "circmean" if self.aggregation == "circmean" else "cire"

    @property
    def label(self) -> str:
        """Etiqueta corta para la columna del CSV."""
        if self.robust_method == "best_k":
            return f"best_k:{self.selection_criterion}"
        if self.robust_method == "aggregate":
            return f"aggregate:{self.aggregation}:{self.consensus_scale}"
        return self.robust_method


@dataclass
class PreRobustState:
    """Salida de las Etapas 0 a 3, entrada de la Etapa 4."""

    top_matrix:     pd.DataFrame
    theta_pre:      np.ndarray
    expr_full_norm: pd.DataFrame
    core_genes:     list[str]
    lambda_star:    float | None
    sample_names:   list[str]
    n_outliers:     int
    n_top:          int
    r2_anchor:      float
    pares_ok:       float
    evidencia:      float
    timing:         dict = field(default_factory=dict)


@dataclass
class Stage4Result:
    """Salida de la Etapa 4."""

    theta_hat:       pd.Series     # indexada por nombre de muestra
    desacuerdo_K:    float
    direction_flips: int
    atlas:           pd.DataFrame | None = None
    selected_rep:    int | None = None      # solo best_k
    msce:            float | None = None    # solo aggregate


@dataclass
class RunResult:
    """Salida completa de una ejecucion."""

    theta_hat:       pd.Series
    timing:          dict
    n_top:           int
    n_outliers:      int
    r2_anchor:       float
    pares_ok:        float
    evidencia:       float
    desacuerdo_K:    float
    direction_flips: int
    atlas:           pd.DataFrame | None = None


def _diagnostico_pares(
    core_genes: list[str],
    peak_times: np.ndarray,
    r2_fmm:     np.ndarray,
    pares:      list[tuple[str, str]],
) -> tuple[float, float]:
    """
    Dos indicadores observables sobre los pares a priori en el marco final.

    ``pares_ok``
        Fraccion de pares que se cumplen, con la misma regla de signo que usa
        el Paso 2.2, ``sin(peak(g2) - peak(g1)) > 0``.
    ``evidencia``
        Suma de ``|sin(delta)|``, o sea cuanta informacion sobre la
        orientacion habia realmente disponible. Su maximo es el numero de
        pares evaluables y se alcanza cuando todos estan separados 6 h. Cae a
        cero cuando los pares son antipodales, que no distinguen el sentido, o
        cuando los genes que los portan estan mal ajustados y no superan el
        filtro de R2.

    Es la variable candidata a decidir si se puede confiar en el sentido de
    giro de un tejido sin conocer la verdad. El umbral no se fija aqui: lo
    calibra el estudio, comparando esta cantidad contra el signo de ``tau_circular``.

    Ambos son ``nan`` si ningun par es evaluable, que es informacion distinta
    de que ningun par se cumpla.
    """
    indice = {g: i for i, g in enumerate(core_genes)}
    n_valid = n_agree = 0
    evidencia = 0.0
    for g_antes, g_despues in pares:
        i_a, i_d = indice.get(g_antes), indice.get(g_despues)
        if i_a is None or i_d is None:
            continue
        if r2_fmm[i_a] < CONSISTENCY_R2_MIN or r2_fmm[i_d] < CONSISTENCY_R2_MIN:
            continue
        n_valid += 1
        aporte = np.sin(peak_times[i_d] - peak_times[i_a])
        n_agree += int(aporte > 0.0)
        evidencia += abs(aporte)
    if not n_valid:
        return float("nan"), float("nan")
    return float(n_agree / n_valid), float(evidencia)


def _desacuerdo_entre_repeticiones(circular_scales: np.ndarray) -> float:
    """
    Dispersion circular entre las K repeticiones, en ``[0, 1]``.

    Para cada muestra se toma la resultante de sus K fases estimadas y se
    devuelve ``1 - R`` promediado. Vale 0 si las K coinciden y se acerca a 1
    si se reparten por el circulo. Es el indicador de estabilidad que en el
    babuino resulto estar gobernado por la cantidad de senal y no por el
    ruido.
    """
    escalas = np.asarray(circular_scales, dtype=float)
    if escalas.ndim != 2 or escalas.shape[0] < 2:
        return float("nan")
    c = np.cos(escalas).mean(axis=0)
    s = np.sin(escalas).mean(axis=0)
    return float(np.mean(1.0 - np.hypot(c, s)))


def run_stages_0_3(
    matrix:          pd.DataFrame,
    truth:           SimulationTruth,
    n_core_declared: int,
    fmm_method:      str = "base",
) -> PreRobustState:
    """
    Etapas 0 a 3: preprocesado, CPCA, sincronizacion y seleccion del TOP.

    Ninguna de ellas depende de ``robust_method``, asi que su salida se
    reutiliza entre los metodos de resumen del bloque de agregacion.
    """
    cfg    = DEFAULT_CONFIG
    timing: dict = {}

    core_solicitados  = truth.core_genes(n_core_declared)
    consistency_pairs = truth.consistency_pairs()

    t = time.perf_counter()
    prep = Preprocessor(
        sparse_threshold = cfg["sparse_threshold"],
        verbose          = False,
    ).run(matrix)
    timing["preprocess"] = time.perf_counter() - t

    t = time.perf_counter()
    core_result = CoreGeneSelector.from_string(
        ",".join(core_solicitados), verbose=False,
    ).select(prep.expr_norm)
    timing["core_genes"] = time.perf_counter() - t

    t = time.perf_counter()
    cpca_result = CPCA(
        core_genes              = core_result.genes,
        tight_radius            = cfg["tight_radius"],
        loose_radius            = cfg["loose_radius"],
        spline_gene_threshold   = cfg["spline_gene_threshold"],
        spline_outlier_fraction = cfg["spline_outlier_fraction"],
        max_outlier_fraction    = cfg["max_outlier_fraction"],
        fmm_length_alpha_grid   = cfg["fmm_alpha_grid"],
        fmm_length_omega_grid   = cfg["fmm_omega_grid"],
        fmm_num_reps            = cfg["fmm_reps"],
        fmm_method              = fmm_method,
        verbose                 = False,
    ).run(prep.expr_norm)
    timing["cpca"] = time.perf_counter() - t

    core_encontrados = list(cpca_result.core_genes_found)

    t = time.perf_counter()
    order_result = CircularSynchronizer(
        anchor_gene       = truth.anchor,
        direction_gene    = truth.direction,
        consistency_pairs = consistency_pairs,
        verbose           = False,
    ).run(
        expr_ordered   = cpca_result.expr_norm_final,
        circular_scale = cpca_result.circular_scale,
        core_genes     = core_encontrados,
        fmm_fits       = cpca_result.fmm_fits_final,
    )
    timing["preliminary_order"] = time.perf_counter() - t

    # Observables, sobre el marco final que es el que se reporta
    r2_fmm    = np.asarray(order_result.r2_fmm, dtype=float)
    i_ancla   = (core_encontrados.index(truth.anchor)
                 if truth.anchor in core_encontrados else None)
    r2_anchor = float(r2_fmm[i_ancla]) if i_ancla is not None else float("nan")
    pares_ok, evidencia = _diagnostico_pares(
        core_encontrados, np.asarray(order_result.peak_times, dtype=float),
        r2_fmm, consistency_pairs,
    )

    t = time.perf_counter()
    top_result = TopGeneSelector(
        r2_ori_threshold      = R2_ORI_THRESHOLD,
        r2_fmm_threshold      = R2_FMM_THRESHOLD,
        omega_min             = OMEGA_MIN,
        n_sectors             = N_SECTORS,
        top_fraction          = TOP_FRACTION,
        fmm_length_alpha_grid = cfg["fmm_alpha_grid"],
        fmm_length_omega_grid = cfg["fmm_omega_grid"],
        fmm_num_reps          = cfg["fmm_reps"],
        fmm_method            = fmm_method,
        lambda_star           = cpca_result.lambda_star,
        n_jobs                = 1,
        verbose               = False,
    ).run(
        expr_norm      = order_result.expr_ordered,
        circular_scale = order_result.circular_scale,
        seed_genes     = core_encontrados,
        sample_order   = order_result.sample_order,
    )
    timing["top_genes"] = time.perf_counter() - t

    return PreRobustState(
        top_matrix     = top_result.candidate_matrix,
        theta_pre      = top_result.circular_scale,
        expr_full_norm = order_result.expr_ordered,
        core_genes     = core_encontrados,
        lambda_star    = cpca_result.lambda_star,
        sample_names   = list(order_result.expr_ordered.columns),
        n_outliers     = len(cpca_result.samples_dropped),
        n_top          = int(top_result.candidate_matrix.shape[0]),
        r2_anchor      = r2_anchor,
        pares_ok       = pares_ok,
        evidencia      = evidencia,
        timing         = timing,
    )


def _estimador(
    state:      PreRobustState,
    truth:      SimulationTruth,
    robust_cfg: RobustConfig,
    fmm_method: str,
) -> RobustOrderEstimator:
    """El estimador de la Etapa 4 con la configuracion de ``run_pipeline``."""
    cfg = DEFAULT_CONFIG
    return RobustOrderEstimator(
        n_reps                 = robust_cfg.n_reps,
        sample_size_fraction   = SAMPLE_SIZE_FRACTION,
        method                 = robust_cfg.robust_method,
        selection_criterion    = robust_cfg.selection_criterion,
        selection_r2_threshold = SELECTION_R2_THRESHOLD,
        aggregation_method     = robust_cfg.aggregation,
        anchor_gene            = truth.anchor,
        direction_gene         = truth.direction,
        consistency_pairs      = truth.consistency_pairs(),
        fmm_length_alpha_grid  = cfg["fmm_alpha_grid"],
        fmm_length_omega_grid  = cfg["fmm_omega_grid"],
        fmm_num_reps           = cfg["fmm_reps"],
        fmm_method             = fmm_method,
        lambda_star            = state.lambda_star,
        n_jobs                 = 1,
        seed                   = robust_cfg.seed,
        verbose                = False,
    )


def run_stage_4_repetitions(
    state:      PreRobustState,
    truth:      SimulationTruth,
    fmm_method: str = "base",
    robust_cfg: RobustConfig | None = None,
) -> RobustRepetitions:
    """
    Las K repeticiones de la Etapa 4, sin resumir.

    Es donde termina la Etapa 4 del CIRCUST original, asi que es el punto en
    que se para el cronometro al comparar con R. Solo dependen de ``n_reps`` y
    de la semilla, no del metodo de resumen, y se comparten entre los metodos
    de un mismo grupo.
    """
    robust_cfg = robust_cfg or RobustConfig()
    return _estimador(state, truth, robust_cfg, fmm_method).run_repetitions(
        top_matrix     = state.top_matrix,
        theta_pre      = state.theta_pre,
        expr_full_norm = state.expr_full_norm,
        core_genes     = state.core_genes,
    )


def run_stage_4(
    state:      PreRobustState,
    truth:      SimulationTruth,
    robust_cfg: RobustConfig,
    fmm_method: str = "base",
    con_atlas:  bool = True,
    reps:       RobustRepetitions | None = None,
) -> Stage4Result:
    """
    Etapa 4: estimacion robusta del orden a partir del estado intermedio.

    Con ``reps`` se resumen unas repeticiones ya calculadas por
    :func:`run_stage_4_repetitions` (con el mismo ``n_reps`` y la misma
    semilla) en vez de volver a calcularlas. El resultado es identico.

    Devuelve las fases estimadas **indexadas por nombre de muestra**, no por
    posicion: el pipeline elimina outliers y reordena columnas.
    """
    cfg = DEFAULT_CONFIG
    if reps is None:
        reps = run_stage_4_repetitions(state, truth, fmm_method, robust_cfg)
    resultado = _estimador(state, truth, robust_cfg, fmm_method).summarize(reps)

    if resultado.sample_order.size == 0:
        raise ValueError(
            f"El metodo {robust_cfg.robust_method!r} no produce orden de "
            f"muestras, no se puede evaluar sobre metricas de orden"
        )

    nombres = [state.sample_names[i] for i in resultado.sample_order]
    fases   = np.asarray(resultado.circular_scale, dtype=float) % TWO_PI

    atlas = None
    if con_atlas:
        # n_jobs=1: el estudio paraleliza a nivel de celda, y con aggregate
        # esta llamada re-ajusta los genes TOP sobre el orden consenso.
        atlas = build_atlas(
            method_used           = resultado.method_used,
            selected_fmm          = resultado.selected_fmm,
            top_matrix            = state.top_matrix,
            sample_order          = resultado.sample_order,
            circular_scale        = resultado.circular_scale,
            fmm_method            = fmm_method,
            fmm_length_alpha_grid = cfg["fmm_alpha_grid"],
            fmm_length_omega_grid = cfg["fmm_omega_grid"],
            fmm_num_reps          = cfg["fmm_reps"],
            lambda_star           = state.lambda_star,
            n_jobs                = 1,
        )

    flips = np.asarray(resultado.direction_flipped)
    return Stage4Result(
        theta_hat       = pd.Series(fases, index=nombres, name="theta_hat"),
        # por muestra: circular_scales[k] viene en el orden de la repeticion k
        desacuerdo_K    = _desacuerdo_entre_repeticiones(
            fases_por_muestra(resultado.sample_orders, resultado.circular_scales)),
        direction_flips = int(flips.sum()) if flips.size else 0,
        atlas           = atlas,
        selected_rep    = resultado.selected_rep,
        msce            = resultado.aggregation_msce,
    )


def run_full(
    matrix:          pd.DataFrame,
    truth:           SimulationTruth,
    n_core_declared: int,
    robust_cfg:      RobustConfig,
    fmm_method:      str = "base",
    con_atlas:       bool = True,
) -> RunResult:
    """Pipeline completo sobre una matriz sintetica."""
    estado = run_stages_0_3(matrix, truth, n_core_declared, fmm_method)

    t = time.perf_counter()
    etapa4 = run_stage_4(estado, truth, robust_cfg, fmm_method, con_atlas)
    estado.timing["robust_order"] = time.perf_counter() - t
    estado.timing["total"] = sum(estado.timing.values())

    return RunResult(
        theta_hat       = etapa4.theta_hat,
        timing          = dict(estado.timing),
        n_top           = estado.n_top,
        n_outliers      = estado.n_outliers,
        r2_anchor       = estado.r2_anchor,
        pares_ok        = estado.pares_ok,
        evidencia       = estado.evidencia,
        desacuerdo_K    = etapa4.desacuerdo_K,
        direction_flips = etapa4.direction_flips,
        atlas           = etapa4.atlas,
    )
