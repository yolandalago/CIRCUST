#!/usr/bin/env python3
"""
run_baboon_validation.py
========================
Orquesta la validacion de CIRCUST sobre el atlas de babuino (Mure et al. 2018),
**un tejido a la vez**, contra el ZT real del muestreo. Es el espejo de
``run_simulation_grid.py`` pero iterando por tejido en lugar de por celda
factorial. Por cada tejido encadena:

    prepare_baboon (una vez)  ->  run_pipeline  ->  validate

y agrega los ``summary.csv`` en un unico ``results.csv`` (una fila por tejido),
listo para ``analyze_baboon.py``. Cada configuracion (FMM base/stabilized,
robust best_k/aggregate) se vuelca a su propio ``-o``, de modo que luego se
pueden contrastar pareando por tejido.

Estructura de salida
--------------------
    <out>/                      (p.ej. output/baboon/base_bestk)
      <TEJIDO>/
        run/       sample_order.csv, circust_atlas.csv, manifest.json
        validate/  summary.csv, per_sample.csv, manifest.json
      results.csv               agregado (tejido + config + metricas + tiempos)

Los datasets por tejido (matrix.csv + true_times.csv) se comparten entre
configuraciones en ``<out>/../_tissues/`` (se extraen una sola vez).

Ejemplos
--------
    # Canonico: FMM base + best_k sobre los 64 tejidos
    python validation/run_baboon_validation.py -o output/baboon/base_bestk

    # Contraste posterior: FMM stabilized
    python validation/run_baboon_validation.py -o output/baboon/base_stab \\
        --fmm-method stabilized

    # Contraste posterior: agregacion TSP3 de muestras
    python validation/run_baboon_validation.py -o output/baboon/base_tsp3 \\
        --robust-method aggregate --aggregation tsp3

    # Smoke test de un tejido
    python validation/run_baboon_validation.py -o /tmp/baboon_smoke --tissues LIV

    # Solo enumerar (no ejecuta nada)
    python validation/run_baboon_validation.py -o /tmp/x --dry-run
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "validation"))

from prepare_baboon import extract_all, DEFAULT_CSV

RUN_PIPELINE = PROJECT_ROOT / "scripts" / "run_pipeline.py"
VALIDATE     = PROJECT_ROOT / "validation" / "validate.py"


# ===========================================================================
# Ejecucion de subprocesos
# ===========================================================================

def _run_subprocess(cmd: list[str], log_path: Path) -> tuple[bool, float]:
    """Lanza un comando, vuelca stdout+stderr a log_path. Devuelve (ok, segundos)."""
    log_path.parent.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    with open(log_path, "w") as log:
        proc = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT)
    return proc.returncode == 0, time.time() - t0


def run_tissue(
    tissue:        str,
    tissue_dir:    Path,
    out_root:      Path,
    fmm_method:    str,
    robust_method: str,
    aggregation:   str,
    n_reps_pipe:   int,
    seed:          int,
    with_cre:      bool,
) -> dict:
    """Corre pipeline + validate de un tejido. Devuelve la fila (dict) con
    tejido, config, metricas y tiempos, o con ``status`` de error."""
    matrix_csv  = tissue_dir / "matrix.csv"
    true_csv    = tissue_dir / "true_times.csv"
    tissue_out  = out_root / tissue
    run_d       = tissue_out / "run"
    val_d       = tissue_out / "validate"

    row = {
        "tissue":        tissue,
        "fmm_method":    fmm_method,
        "robust_method": robust_method,
        "aggregation":   aggregation if robust_method == "aggregate" else None,
        "seed":          seed,
        "status":        "ok",
        "t_pipeline_s":  np.nan,
        "t_validate_s":  np.nan,
    }

    if not matrix_csv.exists() or not true_csv.exists():
        row["status"] = "missing_tissue_data"
        return row

    # 1) Pipeline CIRCUST (gen column por defecto = primera col = Symbol)
    pipe_cmd = [sys.executable, str(RUN_PIPELINE),
                "--data", str(matrix_csv),
                "--label", tissue,
                "--fmm-method", fmm_method,
                "--robust-method", robust_method,
                "--aggregation", aggregation,
                "--n-reps", str(n_reps_pipe),
                "--seed", str(seed),
                "-o", str(run_d)]
    ok, secs = _run_subprocess(pipe_cmd, tissue_out / "pipeline.log")
    row["t_pipeline_s"] = secs
    if not ok:
        row["status"] = "pipeline_error"
        return row

    # 1b) Timing por etapa desde el manifest del pipeline (perf_counter)
    run_manifest = run_d / "manifest.json"
    if run_manifest.exists():
        try:
            tdict = json.loads(run_manifest.read_text()).get("timing", {})
            for stage, secs_stage in tdict.items():
                name = stage[:-2] if stage.endswith("_s") else stage
                row[f"t_pipe_{name}_s"] = secs_stage
        except (json.JSONDecodeError, OSError):
            pass

    # 2) Validacion contra el ZT real
    val_cmd = [sys.executable, str(VALIDATE),
               "--predicted", str(run_d / "sample_order.csv"),
               "--true",      str(true_csv),
               "--label",     tissue,
               "-o",          str(val_d)]
    if with_cre:
        val_cmd += ["--expr",  str(matrix_csv),
                    "--atlas", str(run_d / "circust_atlas.csv"),
                    "--fmm-method", fmm_method]
    ok, secs = _run_subprocess(val_cmd, tissue_out / "validate.log")
    row["t_validate_s"] = secs
    if not ok:
        row["status"] = "validate_error"
        return row

    # 3) Recoger metricas del summary.csv (sin pisar tejido/config)
    summary_path = val_d / "summary.csv"
    if summary_path.exists():
        summ = pd.read_csv(summary_path).iloc[0].to_dict()
        for k, v in summ.items():
            if k != "label":
                row[k] = v
    else:
        row["status"] = "no_summary"

    return row


# ===========================================================================
# CLI
# ===========================================================================

def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=("Orquesta prepare_baboon->run_pipeline->validate por "
                     "tejido y agrega los resultados en results.csv."),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("-o", "--output", required=True, type=str,
                   help="Directorio de salida de ESTA configuracion "
                        "(p.ej. output/baboon/base_bestk).")
    p.add_argument("--data", default=str(DEFAULT_CSV), type=str,
                   help="CSV del atlas (solo se usa si hay que extraer).")
    p.add_argument("--tissues-dir", default=None, type=str,
                   help="Carpeta _tissues/ con los datasets por tejido. "
                        "Defecto: <output>/../_tissues. Si no existe, se extrae.")
    p.add_argument("--tissues", default=None, type=str,
                   help="Subconjunto de tejidos (coma). Sin flag: todos.")
    # Paso a traves al pipeline
    p.add_argument("--fmm-method", choices=["base", "stabilized"], default="base",
                   help="Modelo FMM del pipeline. Defecto: base.")
    p.add_argument("--robust-method", choices=["best_k", "aggregate"], default="best_k",
                   help="Eleccion del orden final: best_k o aggregate. Defecto: best_k.")
    p.add_argument("--aggregation", choices=["auto", "tsp3", "hods"], default="auto",
                   help="Metodo de agregacion (solo si robust-method=aggregate). Defecto: auto.")
    p.add_argument("--n-reps-pipeline", type=int, default=5,
                   help="K repeticiones de RobustOrder en el pipeline. Defecto: 5.")
    p.add_argument("--seed", type=int, default=42,
                   help="Semilla del pipeline. Defecto: 42.")
    p.add_argument("--no-cre", action="store_true",
                   help="Desactiva el calculo de CRE (re-fit FMM en validate).")
    p.add_argument("--limit", type=int, default=None,
                   help="Limita el numero de tejidos (tras enumerar).")
    p.add_argument("--resume", action="store_true",
                   help="Salta tejidos cuya validate/summary.csv ya exista.")
    p.add_argument("--dry-run", action="store_true",
                   help="Enumera tejidos, no ejecuta nada.")
    return p.parse_args()


def main() -> None:
    args = _parse_args()
    out_root = Path(args.output)

    only = ([t.strip() for t in args.tissues.split(",") if t.strip()]
            if args.tissues else None)

    # ── Resolver / construir los datasets por tejido ──────────────────────
    tissues_dir = (Path(args.tissues_dir) if args.tissues_dir
                   else out_root.parent / "_tissues")

    def _has_tissue_data(d: Path) -> bool:
        return d.exists() and any(
            sub.is_dir() and (sub / "matrix.csv").exists()
            and (sub / "true_times.csv").exists()
            for sub in d.iterdir()
        )

    if not _has_tissue_data(tissues_dir):
        print(f"  Datasets por tejido no encontrados en {tissues_dir} — extrayendo...")
        extract_all(Path(args.data), tissues_dir.parent, only=only, verbose=True)

    # Enumerar tejidos disponibles (carpetas con matrix.csv + true_times.csv)
    available = sorted(
        d.name for d in tissues_dir.iterdir()
        if d.is_dir() and (d / "matrix.csv").exists()
        and (d / "true_times.csv").exists()
    )
    if only:
        available = [t for t in available if t in set(only)]
    if args.limit is not None:
        available = available[:args.limit]

    n_tissues = len(available)
    tag = (f"{args.fmm_method}_{args.robust_method}"
           + (f"_{args.aggregation}" if args.robust_method == "aggregate" else ""))

    print("╔══════════════════════════════════════════════════════════════════╗")
    print("  Validacion de CIRCUST sobre babuino (Mure 2018) vs ZT")
    print(f"    tejidos        : {n_tissues}")
    print(f"    config         : {tag}")
    print(f"    fmm-method     : {args.fmm_method}")
    print(f"    robust-method  : {args.robust_method}"
          + (f" ({args.aggregation})" if args.robust_method == "aggregate" else ""))
    print(f"    n-reps pipeline: {args.n_reps_pipeline}")
    print(f"    CRE            : {'no' if args.no_cre else 'si'}")
    print(f"    _tissues       : {tissues_dir.resolve()}")
    print(f"    output         : {out_root.resolve()}")
    print("╚══════════════════════════════════════════════════════════════════╝")

    if args.dry_run:
        print(f"\n  [dry-run] tejidos: {available}")
        return

    out_root.mkdir(parents=True, exist_ok=True)
    rows: list[dict] = []
    t_start = time.time()

    for i, tissue in enumerate(available, 1):
        val_summary = out_root / tissue / "validate" / "summary.csv"
        head = f"[{i}/{n_tissues}] {tissue}"

        if args.resume and val_summary.exists():
            print(f"  {head}  (skip: ya existe)")
            summ = pd.read_csv(val_summary).iloc[0].to_dict()
            row = {"tissue": tissue, "fmm_method": args.fmm_method,
                   "robust_method": args.robust_method,
                   "aggregation": (args.aggregation
                                   if args.robust_method == "aggregate" else None),
                   "seed": args.seed, "status": "cached"}
            row.update({k: v for k, v in summ.items() if k != "label"})
            rows.append(row)
            pd.DataFrame(rows).to_csv(out_root / "results.csv", index=False)
            continue

        print(f"  {head} ...", flush=True)
        row = run_tissue(
            tissue=tissue, tissue_dir=tissues_dir / tissue, out_root=out_root,
            fmm_method=args.fmm_method, robust_method=args.robust_method,
            aggregation=args.aggregation, n_reps_pipe=args.n_reps_pipeline,
            seed=args.seed, with_cre=not args.no_cre,
        )
        rows.append(row)
        elapsed = time.time() - t_start
        eta = elapsed / i * (n_tissues - i)
        print(f"      status={row['status']}  "
              f"MedAE={row.get('MedAE', float('nan')):.2f}h  "
              f"%2h={row.get('pct_within_2h', float('nan')):.0f}  "
              f"t_pipe={row.get('t_pipeline_s', float('nan')):.1f}s  "
              f"| ETA {eta/60:.1f} min")

        # Persistir parcial cada vez (robustez ante cortes)
        pd.DataFrame(rows).to_csv(out_root / "results.csv", index=False)

    results = pd.DataFrame(rows)
    results.to_csv(out_root / "results.csv", index=False)

    elapsed = time.time() - t_start
    n_ok = int((results["status"] == "ok").sum()) if "status" in results else 0
    print()
    print(f"  ✓ Validacion completada en {elapsed/60:.1f} min")
    print(f"    {n_ok}/{len(results)} tejidos OK")
    if "MedAE" in results.columns and n_ok:
        ok = results[results["status"] == "ok"]
        print(f"    MedAE mediana={ok['MedAE'].median():.2f}h  "
              f"%within2h mediana={ok['pct_within_2h'].median():.0f}%")
    print(f"    Agregado: {out_root / 'results.csv'}")


if __name__ == "__main__":
    main()
