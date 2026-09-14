"""
Ejecutor de los bloques del estudio de simulacion.

Vuelca un CSV largo, una fila por celda, replica y metodo de resumen. Es
reanudable: relee lo ya escrito y salta esas filas, asi que se puede
interrumpir y continuar sin perder trabajo.

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
from simulation.metrics import evaluate
from simulation.runner import RobustConfig, run_stage_4, run_stages_0_3
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
    "mae_h", "medae_h", "rho_circ", "rho_fl", "tau",
    "p1", "p2", "p4",
    # metricas invariantes al marco, es el diagnostico
    "tau_circ", "orientacion_ok", "mae_h_oracle",
    # observables sin la verdad, para calibrar el aviso de fase absoluta
    "r2_anchor", "pares_ok", "evidencia", "desacuerdo_K",
    "n_top", "n_outliers", "n_evaluadas", "direction_flips",
    # calidad del atlas: es lo que separa el FMM base del estabilizado, que
    # sobre el orden de muestras son identicos por construccion
    *COLUMNAS_ATLAS,
    # tiempos
    "t_preprocess_s", "t_cpca_s", "t_sync_s", "t_top_s", "t_robust_s",
    "t_total_s",
]

# Filas que se acumulan antes de escribir. Un bloque tarda horas, asi que el
# CSV se va volcando por lotes en vez de al final: lo que se pierde en una
# interrupcion es a lo sumo un lote, y ``existing_ids`` reanuda solo.
LOTE_ESCRITURA = 50

_METRICAS = ("mae_h", "medae_h", "rho_circ", "rho_fl", "tau", "tau_circ",
             "orientacion_ok", "p1", "p2", "p4", "mae_h_oracle")


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
) -> list[dict]:
    """
    Genera una matriz, corre las Etapas 0 a 3 una vez y la Etapa 4 por metodo.

    Devuelve una fila por metodo. Un fallo no interrumpe nada: se registra en
    la fila con ``status='fallo'`` porque la tasa de fallo por celda es un
    resultado del estudio y no algo que esconder.
    """
    t0 = time.perf_counter()
    filas = [
        _fila_base(bloque, cfg, n_core, fmm_method, m, rep) for m in metodos
    ]

    try:
        matriz, truth = generate_dataset(cfg)
        estado = run_stages_0_3(matriz, truth, n_core, fmm_method)
    except Exception as exc:                      # noqa: BLE001
        for fila in filas:
            fila["status"] = "fallo"
            fila["error"] = f"{type(exc).__name__}: {exc}"
            fila["t_total_s"] = time.perf_counter() - t0
        return filas

    tiempos = estado.timing
    for fila, metodo in zip(filas, metodos):
        t_metodo = time.perf_counter()
        try:
            etapa4     = run_stage_4(estado, truth, metodo, fmm_method)
            theta_hat  = etapa4.theta_hat
            theta_true = truth.theta_for(list(theta_hat.index))
            fila.update(evaluate(theta_true, theta_hat.values))
            fila["n_evaluadas"]     = len(theta_hat)
            fila["desacuerdo_K"]    = etapa4.desacuerdo_K
            fila["direction_flips"] = etapa4.direction_flips
            # quitar_giro: el giro global de la escala es error de anclaje, ya
            # medido por mae_h, y dejarlo dentro taparia la calidad del ajuste
            fila.update(atlas_metrics(etapa4.atlas, truth.gene_params,
                                      quitar_giro=True))
        except Exception as exc:                  # noqa: BLE001
            fila["status"] = "fallo"
            fila["error"] = f"{type(exc).__name__}: {exc}"
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
        fila["t_robust_s"]     = time.perf_counter() - t_metodo
        fila["t_total_s"]      = (
            sum(v for v in tiempos.values()) + fila["t_robust_s"]
        )
    return filas


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

    Es lo que hace barato el bloque de agregacion: las cinco variantes de
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
    parser.add_argument("--out", type=str, default="output/simulacion",
                        help="Directorio de salida. Por defecto: output/simulacion")
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

    def volcar() -> None:
        nonlocal lote
        if not lote:
            return
        df = pd.DataFrame(lote, columns=COLUMNAS)
        df.to_csv(destino, mode="a", header=not destino.exists(), index=False)
        lote = []

    for grupo in resultados:
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
