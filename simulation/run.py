"""
Ejecutor de los bloques del estudio de simulacion.

Vuelca un CSV largo, una fila por celda, replica y metodo de resumen. Es
reanudable: relee lo ya escrito y salta esas filas, asi que se puede
interrumpir y continuar sin perder trabajo.

Junto al CSV guarda las salidas del pipeline de cada grupo y metodo
(``simulation/salidas.py``), para poder recalcular o anadir medidas sin
volver a ejecutar nada.

Las K repeticiones de la Etapa 4 no dependen del metodo de resumen, asi que se
calculan una vez por grupo y se resumen con cada metodo. ``t_reps_s`` es su
coste, el punto en que termina la Etapa 4 del CIRCUST original, y
``t_resumen_s`` el del resumen y el atlas de cada metodo. ``t_robust_s`` es la
suma, lo que costaria la Etapa 4 de ese metodo ejecutada sola.

Uso
---
    python -m simulation.run --bloque sin_ruido --reps 3
    python -m simulation.run --bloque principal --reps 20 --jobs 8
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from simulation.atlas import COLUMNAS_ATLAS, atlas_metrics
from simulation.generate import generate_dataset
from simulation.grid import BLOQUES, Cell, build_grid
from simulation.metrics import ajuste_fmm, cre, evaluate, genes_ritmicos
from simulation.runner import (
    RobustConfig, run_stage_4, run_stage_4_repetitions, run_stages_0_3,
)
from simulation.salidas import Salidas, grupo_id, salidas_grupo, salidas_metodo
from simulation.truth import SimConfig

COLUMNAS: list[str] = [
    # identidad
    "bloque", "row_id", "rep", "seed",
    # factores de generacion
    "n_samples", "sampling", "omega", "A", "sigma", "snr",
    "n_genes", "prop_rhythmic", "n_core",
    # configuracion del algoritmo
    "fmm_method", "robust_method", "selection_criterion",
    "aggregation", "consensus_scale",
    # estado
    "status", "error",
    # metricas en el marco tal cual, es lo que se reporta
    "error_fase_medio_h", "error_fase_h", "rho_circ", "rho_fl", "tau",
    "p1", "p2", "p4",
    # metricas invariantes al marco, es el diagnostico
    "tau_circular", "error_orden_h", "cre",
    # observables sin la verdad, para calibrar el aviso de fase absoluta
    "r2_anchor", "pares_ok", "evidencia", "desacuerdo_K",
    "n_top", "n_outliers", "n_evaluadas", "direction_flips",
    # que eligio el metodo: la repeticion (best_k) o el MSCE del consenso
    "rep_elegida", "msce",
    # calidad del atlas: es lo que separa el FMM base de CircaFMM, que
    # sobre el orden de muestras son identicos por construccion
    *COLUMNAS_ATLAS,
    # tiempos
    "t_preprocess_s", "t_cpca_s", "t_sync_s", "t_top_s",
    "t_reps_s", "t_resumen_s", "t_robust_s", "t_total_s",
]

# Filas que se acumulan antes de escribir. Un bloque tarda horas, asi que el
# CSV se va volcando por lotes en vez de al final: lo que se pierde en una
# interrupcion es a lo sumo un lote, y ``existing_ids`` reanuda solo.
LOTE_ESCRITURA = 50

_METRICAS = ("error_fase_medio_h", "error_fase_h", "rho_circ", "rho_fl", "tau", "tau_circular",
             "p1", "p2", "p4", "error_orden_h", "cre")


class _ReferenciaCRE:
    """
    Ajustes FMM de referencia del CRE, con los tiempos verdaderos.

    Las Etapas 0 a 3 son comunes a todos los metodos de un grupo, asi que la
    referencia se ajusta una sola vez y la comparten. Se guarda por conjunto
    de muestras evaluadas por si algun metodo devolviera menos.

    El ajustador es el modelo FMM del grupo, y con CircaFMM la ``lambda_star``
    que calibro el pipeline, igual para la referencia y para las estimadas.
    """

    def __init__(self, expr: pd.DataFrame, truth, fmm_method: str = "base",
                 lambda_star: float | None = None) -> None:
        self.expr  = expr.loc[genes_ritmicos(truth.gene_params, expr.index)]
        self.truth = truth
        self.modelo = dict(fmm_method=fmm_method, lambda_star=lambda_star)
        self._cache: dict[tuple, tuple[np.ndarray, np.ndarray]] = {}

    def cre(self, theta_hat: pd.Series) -> float:
        if self.expr.empty:
            return float("nan")
        nombres = tuple(sorted(theta_hat.index))
        if nombres not in self._cache:
            Y = self.expr[list(nombres)].to_numpy()
            self._cache[nombres] = (
                Y, ajuste_fmm(Y, self.truth.theta_for(list(nombres)), **self.modelo),
            )
        Y, referencia = self._cache[nombres]
        return cre(referencia, ajuste_fmm(Y, theta_hat.loc[list(nombres)].to_numpy(),
                                          **self.modelo))


def _fila_base(
    bloque:     str,
    cfg:        SimConfig,
    n_core:     int,
    fmm_method: str,
    robust_cfg: RobustConfig,
    rep:        int,
) -> dict:
    """Fila con los factores rellenos y todo lo demas a vacio."""
    fila = {c: np.nan for c in COLUMNAS}
    fila.update(
        bloque              = bloque,
        rep                 = rep,
        seed                = cfg.seed,
        n_samples           = cfg.n_samples,
        sampling            = cfg.sampling,
        omega               = str(cfg.omega),
        A                   = cfg.amplitude,
        sigma               = cfg.sigma,
        snr                 = cfg.snr,
        n_genes             = cfg.n_genes,
        prop_rhythmic       = cfg.prop_rhythmic,
        n_core              = n_core,
        fmm_method          = fmm_method,
        robust_method       = robust_cfg.robust_method,
        selection_criterion = robust_cfg.selection_criterion,
        aggregation         = robust_cfg.aggregation,
        consensus_scale     = robust_cfg.consensus_scale,
        status              = "ok",
        error               = "",
    )
    fila["row_id"] = Cell(bloque, cfg, n_core, fmm_method, robust_cfg, rep).row_id
    return fila


def run_cell_group(
    cfg:        SimConfig,
    n_core:     int,
    fmm_method: str,
    metodos:    list[RobustConfig],
    bloque:     str,
    rep:        int,
) -> tuple[list[dict], Salidas]:
    """
    Genera una matriz, corre las Etapas 0 a 3 y las K repeticiones una vez, y
    el resumen de la Etapa 4 por metodo.

    Devuelve una fila por metodo y las salidas del grupo. Un fallo no
    interrumpe nada: se registra en la fila con ``status='fallo'`` porque la
    tasa de fallo por celda es un resultado del estudio y no algo que
    esconder.
    """
    t0 = time.perf_counter()
    filas = [
        _fila_base(bloque, cfg, n_core, fmm_method, m, rep) for m in metodos
    ]
    gid = grupo_id(bloque, cfg, n_core, fmm_method, rep)
    salidas = Salidas()
    if len({(m.n_reps, m.seed) for m in metodos}) > 1:
        raise ValueError(
            "Los metodos de un grupo comparten las K repeticiones, asi que "
            "tienen que tener el mismo n_reps y la misma semilla"
        )

    try:
        matriz, truth = generate_dataset(cfg)
        estado = run_stages_0_3(matriz, truth, n_core, fmm_method)
        t = time.perf_counter()
        reps = run_stage_4_repetitions(estado, truth, fmm_method, metodos[0])
        t_reps = time.perf_counter() - t
        salidas.unir(salidas_grupo(gid, estado, reps, t_reps))
    except Exception as exc:                      # noqa: BLE001
        for fila in filas:
            fila["status"] = "fallo"
            fila["error"] = f"{type(exc).__name__}: {exc}"
            fila["t_total_s"] = time.perf_counter() - t0
        return filas, salidas

    tiempos = estado.timing
    referencia_cre = _ReferenciaCRE(estado.expr_full_norm, truth,
                                    fmm_method, estado.lambda_star)
    for fila, metodo in zip(filas, metodos):
        t_metodo = time.perf_counter()
        theta_hat = None
        try:
            etapa4     = run_stage_4(estado, truth, metodo, fmm_method, reps=reps)
            t_resumen  = time.perf_counter() - t_metodo
            salidas.unir(salidas_metodo(fila["row_id"], gid, etapa4))
            theta_hat  = etapa4.theta_hat
            fila["rep_elegida"] = etapa4.selected_rep
            fila["msce"]        = etapa4.msce
            theta_true = truth.theta_for(list(theta_hat.index))
            fila.update(evaluate(theta_true, theta_hat.values))
            fila["n_evaluadas"]     = len(theta_hat)
            fila["desacuerdo_K"]    = etapa4.desacuerdo_K
            fila["direction_flips"] = etapa4.direction_flips
            # quitar_giro: el giro global de la escala es error de anclaje, ya
            # medido por error_fase_medio_h, y dejarlo dentro taparia la calidad del ajuste
            fila.update(atlas_metrics(etapa4.atlas, truth.gene_params,
                                      quitar_giro=True))
        except Exception as exc:                  # noqa: BLE001
            fila["status"] = "fallo"
            fila["error"] = f"{type(exc).__name__}: {exc}"
            theta_hat = None
            t_resumen = time.perf_counter() - t_metodo
            for clave in (*_METRICAS, *COLUMNAS_ATLAS):
                fila[clave] = np.nan
        fila["n_top"]          = estado.n_top
        fila["n_outliers"]     = estado.n_outliers
        fila["r2_anchor"]      = estado.r2_anchor
        fila["pares_ok"]       = estado.pares_ok
        fila["evidencia"]      = estado.evidencia
        fila["t_preprocess_s"] = tiempos.get("preprocess", np.nan)
        fila["t_cpca_s"]       = tiempos.get("cpca", np.nan)
        fila["t_sync_s"]       = tiempos.get("preliminary_order", np.nan)
        fila["t_top_s"]        = tiempos.get("top_genes", np.nan)
        fila["t_reps_s"]       = t_reps
        fila["t_resumen_s"]    = t_resumen
        fila["t_robust_s"]     = t_reps + t_resumen
        fila["t_total_s"]      = (
            sum(v for v in tiempos.values()) + fila["t_robust_s"]
        )
        # Despues de cronometrar: el CRE ajusta un FMM por gen ritmico y es
        # coste de la evaluacion, no del metodo. Un fallo aqui no invalida la
        # fila, solo deja el CRE vacio.
        if theta_hat is not None:
            try:
                fila["cre"] = referencia_cre.cre(theta_hat)
            except Exception as exc:              # noqa: BLE001
                fila["error"] = f"cre: {type(exc).__name__}: {exc}"
    return filas, salidas


def existing_ids(path: Path) -> set[str]:
    """``row_id`` ya presentes en el CSV, para poder reanudar."""
    path = Path(path)
    if not path.exists():
        return set()
    previo = pd.read_csv(path, usecols=["row_id"])
    return set(previo["row_id"].astype(str))


def _agrupar(celdas: list[Cell]) -> list[tuple]:
    """
    Agrupa las celdas que comparten matriz y Etapas 0 a 3.

    Es lo que hace barato el bloque de agregacion: las seis variantes de
    ``robust_method`` de una misma celda se resuelven con un solo pipeline
    hasta la Etapa 3.
    """
    grupos: dict[tuple, list[RobustConfig]] = {}
    for c in celdas:
        clave = (c.bloque, c.cfg, c.n_core_declared, c.fmm_method, c.rep)
        grupos.setdefault(clave, []).append(c.robust_cfg)
    return [(*clave, metodos) for clave, metodos in grupos.items()]


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Ejecuta un bloque del estudio de simulacion de CIRCUST.",
    )
    parser.add_argument("--bloque", required=True, choices=list(BLOQUES))
    parser.add_argument("--reps", type=int, default=20,
                        help="Replicas por celda. Por defecto: 20.")
    parser.add_argument("--out", type=str, default="output/simulacion_2026-10",
                        help="Directorio de salida. Por defecto: output/simulacion_2026-10")
    parser.add_argument("--jobs", type=int, default=-1,
                        help="Procesos en paralelo a nivel de celda. Por defecto: todos.")
    parser.add_argument("--limit", type=int, default=None,
                        help="Corta la rejilla a los primeros N grupos (pruebas).")
    args = parser.parse_args()

    salida = Path(args.out)
    salida.mkdir(parents=True, exist_ok=True)
    destino = salida / "resultados.csv"

    celdas    = build_grid(args.bloque, args.reps)
    ya_hechos = existing_ids(destino)
    celdas    = [c for c in celdas if c.row_id not in ya_hechos]

    grupos = _agrupar(celdas)
    if args.limit is not None:
        grupos = grupos[: args.limit]

    print(f"Bloque {args.bloque}: {len(celdas)} filas pendientes "
          f"en {len(grupos)} grupos ({len(ya_hechos)} ya hechas)")
    if not grupos:
        return

    t0 = time.perf_counter()
    # Generador, no lista: las filas se escriben segun llegan. Un bloque tarda
    # horas y antes el CSV no se tocaba hasta el final, de modo que cualquier
    # interrupcion tiraba todo el trabajo. Ahora solo se pierde el ultimo lote
    # y ``existing_ids`` reanuda por donde iba.
    resultados = Parallel(n_jobs=args.jobs, verbose=10,
                          return_as="generator_unordered")(
        delayed(run_cell_group)(cfg, n_core, fmm, metodos, bloque, rep)
        for bloque, cfg, n_core, fmm, rep, metodos in grupos
    )

    n_filas = n_fallos = 0
    cpu = 0.0
    lote: list[dict] = []
    salidas_lote = Salidas()

    def volcar() -> None:
        nonlocal lote
        if not lote:
            return
        # Primero las salidas y despues el CSV, que es el registro con el que
        # se reanuda: si se corta entre los dos, el grupo se repite y
        # ``salidas.leer`` descarta la copia vieja.
        salidas_lote.escribir(salida)
        df = pd.DataFrame(lote, columns=COLUMNAS)
        df.to_csv(destino, mode="a", header=not destino.exists(), index=False)
        lote = []

    for grupo, salidas_grupo_ in resultados:
        salidas_lote.unir(salidas_grupo_)
        lote.extend(grupo)
        n_filas += len(grupo)
        n_fallos += sum(1 for f in grupo if f["status"] == "fallo")
        cpu += sum(f["t_total_s"] for f in grupo
                   if f["t_total_s"] == f["t_total_s"])   # descarta nan
        if len(lote) >= LOTE_ESCRITURA:
            volcar()
    volcar()
    elapsed = time.perf_counter() - t0

    print(f"\nEscritas {n_filas} filas en {destino}")
    print(f"Fallos: {n_fallos} de {n_filas}")
    print(f"Tiempo de pared: {elapsed:.1f} s "
          f"({elapsed / max(len(grupos), 1):.2f} s por grupo)")
    print(f"Tiempo de CPU acumulado: {cpu:.1f} s "
          f"({cpu / max(n_filas, 1):.2f} s por fila)")


if __name__ == "__main__":
    main()
