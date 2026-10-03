#!/usr/bin/env python3
"""
run_circust.py
=======================
Ejecuta CIRCUST (Larriba et al., 2023) sobre una matriz de expresion.

Etapas
------
    1    Preprocesado: filtro del 30 % de ceros y normalizacion a [-1, 1]
    2    CPCA con los genes core, deteccion de outliers y orden preliminar
         sincronizado (ARNTL en pi y sentido)
    3    Genes TOP (ORI, intervalos de la circunferencia, FMM)
    4    K subconjuntos aleatorios, K ordenes sincronizados, ajuste FMM de los
         genes TOP en cada orden y atlas con la mediana de los K parametros

Ejemplos
--------
    python scripts/run_circust.py --data data/matrixIn_mini.parquet -o resultados
    python scripts/run_circust.py --data datos.csv --core-genes zhang --n-reps 10
    python scripts/run_circust.py --data datos.csv --values original --matrix-format parquet
"""
from __future__ import annotations

import argparse
import time
from pathlib import Path

from circust.core_genes import CoreGeneSelector
from circust.cpca import CPCA
from circust.outputs import write_results
from circust.preprocessing import Preprocessor, load_expression_matrix
from circust.robust_order import RobustOrderEstimator
from circust.synchronizer import CircularSynchronizer
from circust.top_genes import TopGeneSelector

_GENE_COL = {".parquet": "gene_id"}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Ejecuta CIRCUST (etapas 1 a 4).",
        formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    p.add_argument("--data", required=True, help="Matriz de expresion (CSV, TSV, Parquet, Excel).")
    p.add_argument("--gene-column", default=None,
                   help="Columna con los simbolos de gen. Por defecto gene_id en Parquet y la primera en el resto.")
    p.add_argument("--core-genes", default=None,
                   help="Preset (circust, zhang, ruben) o lista separada por comas. Por defecto circust.")
    p.add_argument("-o", "--output", default="output_circust", help="Carpeta de salida.")
    p.add_argument("--n-reps", type=int, default=5, help="Numero K de ordenes. Por defecto 5.")
    p.add_argument("--seed", type=int, default=42, help="Semilla del sorteo de subconjuntos. Por defecto 42.")
    p.add_argument("--values", choices=["normalized", "original"], default="normalized",
                   help="Valores de ordered_matrix_k<k>: normalizados a [-1,1] u originales.")
    p.add_argument("--matrix-format", choices=["csv", "parquet"], default="csv",
                   help="Formato de las matrices completas ordenadas.")
    p.add_argument("--n-intervals", type=int, default=8,
                   help="Intervalos de la circunferencia en la Etapa 3. Por defecto 8.")
    p.add_argument("--genes-per-interval", type=int, default=50,
                   help="Genes admitidos por intervalo en la Etapa 3. Por defecto 50.")
    p.add_argument("--jobs", type=int, default=-1, help="Procesos en paralelo. -1 usa todos.")
    p.add_argument("--quiet", action="store_true", help="No imprimir el progreso.")
    return p.parse_args()


def main() -> None:
    a = parse_args()
    v = not a.quiet
    t0 = time.perf_counter()
    timing: dict[str, float] = {}

    data = Path(a.data)
    gene_col = a.gene_column or _GENE_COL.get(data.suffix.lower())
    selector = CoreGeneSelector.from_string(a.core_genes, verbose=v)

    tic = time.perf_counter()
    raw = load_expression_matrix(str(data), gene_column=gene_col)
    prep = Preprocessor(sparse_threshold=0.3, verbose=v).run(raw)
    core = selector.select(prep.expr_norm).genes
    timing["preprocess"] = time.perf_counter() - tic

    tic = time.perf_counter()
    cpca = CPCA(core_genes=core, verbose=v).run(prep.expr_norm)
    sync = CircularSynchronizer(verbose=v).run(
        cpca.expr_norm_final, cpca.circular_scale,
        list(cpca.core_genes_found), cpca.fmm_fits_final)
    timing["preliminary_order"] = time.perf_counter() - tic

    tic = time.perf_counter()
    top = TopGeneSelector(n_intervals=a.n_intervals, genes_per_interval=a.genes_per_interval,
                          n_jobs=a.jobs, verbose=v).run(
        sync.expr_ordered, sync.circular_scale, list(cpca.core_genes_found))
    timing["top_genes"] = time.perf_counter() - tic

    tic = time.perf_counter()
    robust = RobustOrderEstimator(n_reps=a.n_reps, seed=a.seed, n_jobs=a.jobs, verbose=v).run(
        top.candidate_matrix, top.circular_scale, list(cpca.core_genes_found),
        forced_genes=top.forced_genes)
    timing["robust_estimation"] = time.perf_counter() - tic
    timing["total"] = time.perf_counter() - t0

    manifest = {
        "input": {"data": str(data), "gene_column": gene_col, "core_genes": list(core),
                  "genes_in": prep.n_genes_in, "genes_after_filter": prep.n_genes_out,
                  "samples_in": raw.shape[1], "samples_after_outliers": len(cpca.sample_order)},
        "algorithm": {"n_reps": a.n_reps, "seed": a.seed, "n_intervals": a.n_intervals,
                      "genes_per_interval": a.genes_per_interval,
                      "top_genes": len(top.gene_names), "forced_core": top.forced_genes,
                      "outliers": cpca.samples_dropped_names,
                      "direction_flipped_per_k": robust.direction_flipped.tolist(),
                      "random_draws": robust.n_attempts, "relaxed_selections": robust.n_relaxed},
        "timing_s": {k: round(x, 3) for k, x in timing.items()},
    }
    files = write_results(
        Path(a.output), robust, sync.expr_ordered, cpca, top, manifest,
        values=a.values, matrix_format=a.matrix_format,
        raw_matrix=raw if a.values == "original" else None)

    print(f"\nCIRCUST completado en {timing['total']:.1f} s")
    print(f"  Genes TOP: {len(top.gene_names)}   Outliers: {len(cpca.samples_dropped)}   K = {a.n_reps}")
    print(f"  {len(files)} ficheros en {Path(a.output).resolve()}")


if __name__ == "__main__":
    main()
