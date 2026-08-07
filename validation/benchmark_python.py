#!/usr/bin/env python3
"""
benchmark_python.py
===================
Lado **Python** del Experimento 1 (R vs Python) de la validacion del TFG.

Corre el pipeline CIRCUST en Python (``scripts/run_pipeline.py``) sobre la
matriz canonica ``data/matrixIn.parquet`` --- la misma que el benchmark R
(``rscripts/benchmark_timing.R``) usa desde ``data/matrixIn.RData`` --- N
veces, y recoge el cronometraje **por etapas** que el pipeline vuelca en
``manifest.json`` (bloque ``timing``).

Salida
------
    ``validation/reference/timing_py.csv``   columnas: rep, stage, seconds

Las etapas Python (load, preprocess, core_genes, cpca, preliminary_order,
top_genes, robust_order, build_atlas, elapsed_total) se alinean con las del
benchmark R en ``analyze_results.py`` (el crosswalk de etapas vive alli, no
aqui, para mantener cada benchmark fiel a su propia implementacion).

Uso (desde la raiz del repo)
----------------------------
    python validation/benchmark_python.py --reps 10
    python validation/benchmark_python.py --reps 1 --data data/matrixIn_sub.parquet  # smoke
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
RUN_PIPELINE = PROJECT_ROOT / "scripts" / "run_pipeline.py"
DEFAULT_DATA = PROJECT_ROOT / "data" / "matrixIn.parquet"
OUT_DIR      = PROJECT_ROOT / "validation" / "reference"


def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=("Benchmark de tiempos del pipeline CIRCUST en Python sobre "
                     "la matriz canonica. Lado Python del Exp 1."),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--reps", type=int, default=10,
                   help="Numero de repeticiones. Defecto: 10.")
    p.add_argument("--data", type=str, default=str(DEFAULT_DATA),
                   help="Matriz de entrada. Defecto: data/matrixIn.parquet.")
    p.add_argument("--fmm-method", choices=["base", "stabilized"], default="base",
                   help="Modelo FMM. Defecto: base (el comparable al R antiguo).")
    p.add_argument("--n-reps", type=int, default=5,
                   help="K repeticiones de RobustOrder (como runCIRCUST K=5). Defecto: 5.")
    p.add_argument("--seed", type=int, default=42, help="Semilla. Defecto: 42.")
    p.add_argument("--keep-runs", action="store_true",
                   help="No borrar los directorios de salida de cada repeticion.")
    return p.parse_args()


def main() -> None:
    args = _parse_args()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    runs_dir = OUT_DIR / "_py_runs"
    runs_dir.mkdir(exist_ok=True)

    data_path = Path(args.data).resolve()
    if not data_path.exists():
        raise SystemExit(f"No existe la matriz de entrada: {data_path}")

    print(f"╔══════════════════════════════════════════════════════════════════╗")
    print(f"  Benchmark Python (Exp 1)")
    print(f"    data       : {data_path}")
    print(f"    reps       : {args.reps}")
    print(f"    fmm-method : {args.fmm_method}")
    print(f"    K (n-reps) : {args.n_reps}")
    print(f"╚══════════════════════════════════════════════════════════════════╝")

    rows: list[dict] = []
    for r in range(1, args.reps + 1):
        run_d = runs_dir / f"rep{r:02d}"
        log_p = runs_dir / f"rep{r:02d}.log"
        print(f"  [{r}/{args.reps}] corriendo pipeline ...", flush=True)

        t0 = time.time()
        with open(log_p, "w") as log:
            proc = subprocess.run(
                [sys.executable, str(RUN_PIPELINE),
                 "--data", str(data_path),
                 "--label", "matrixIn_benchmark",
                 "--fmm-method", args.fmm_method,
                 "--robust-method", "best_k",
                 "--n-reps", str(args.n_reps),
                 "--seed", str(args.seed),
                 "-o", str(run_d)],
                stdout=log, stderr=subprocess.STDOUT,
            )
        wall = time.time() - t0
        if proc.returncode != 0:
            print(f"      ERROR (ver {log_p}); se omite esta repeticion.")
            continue

        manifest = json.loads((run_d / "manifest.json").read_text())
        timing = manifest.get("timing", {})
        for stage, secs in timing.items():
            name = stage[:-2] if stage.endswith("_s") else stage   # elapsed_total_s -> elapsed_total
            rows.append({"rep": r, "stage": name, "seconds": secs})
        # tiempo de pared del proceso (incluye arranque del interprete)
        rows.append({"rep": r, "stage": "wall_clock_process", "seconds": round(wall, 4)})
        print(f"      ok  total_algoritmo={timing.get('elapsed_total_s', float('nan')):.1f}s  "
              f"wall={wall:.1f}s")

    if not rows:
        raise SystemExit("Ninguna repeticion termino correctamente.")

    df = pd.DataFrame(rows)
    out_csv = OUT_DIR / "timing_py.csv"
    df.to_csv(out_csv, index=False)

    if not args.keep_runs:
        import shutil
        shutil.rmtree(runs_dir, ignore_errors=True)

    print()
    print(f"  ✓ Guardado: {out_csv}  ({len(df)} filas)")


if __name__ == "__main__":
    main()
