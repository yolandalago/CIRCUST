#!/usr/bin/env python3
"""
run_simulation_grid.py
======================
Orquesta el **estudio de simulacion factorial** de validacion de CIRCUST.

Recorre el producto cartesiano completo de los factores definidos en
``generate_synthetic.FACTOR_LEVELS`` (muestreo x omega x sigma x prop
ritmica x n) por ``B`` replicas, y para cada celda encadena:

    generate_synthetic  ->  run_pipeline  ->  validate

Cada ``summary.csv`` de ``validate.py`` se anota con los factores de su
celda y se agregan todos en un unico ``results.csv`` *tidy* (una fila por
celda x replica), listo para las tablas y figuras del TFG.

Estructura de salida
--------------------
    <output>/
      <cell_id>/
        rep00/
          synth/     matrix.csv, true_times.csv, gene_params.csv, manifest.json
          run/       sample_order.csv, circust_atlas.csv, ...
          validate/  summary.csv, per_sample.csv, manifest.json
        rep01/ ...
      results.csv    agregado (factores + metricas + tiempos)

Grid completo vs submuestreo
----------------------------
Sin flags de factor se ejecuta el **grid completo**. Para el exploratorio
de tiempos se puede restringir cada factor a un subconjunto::

    # 3 celdas para cronometrar (n=24,48,96 en escenario base)
    python validation/run_simulation_grid.py -o output/grid_timing \\
        --sampling equispaced --omega 0.5 --sigma 1.0 \\
        --prop 0.05 --n 24,48,96 --reps 1

    # grid completo, B=5
    python validation/run_simulation_grid.py -o output/grid_full --reps 5

    # solo contar celdas / estimar tamano (no ejecuta nada)
    python validation/run_simulation_grid.py -o /tmp/x --dry-run
"""
from __future__ import annotations

import argparse
import itertools
import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "validation"))

from generate_synthetic import (
    FACTOR_LEVELS,
    SAMPLING_SCHEMES,
    generate,
    write_dataset,
)

RUN_PIPELINE = PROJECT_ROOT / "scripts" / "run_pipeline.py"
VALIDATE     = PROJECT_ROOT / "validation" / "validate.py"


# ===========================================================================
# Enumeracion del grid
# ===========================================================================

def _cell_id(sampling: str, omega: float, sigma: float,
             prop: float, n: int) -> str:
    """Identificador de carpeta para una celda (factores legibles)."""
    return (f"samp-{sampling}__om{omega:g}__sig{sigma:g}"
            f"__prop{prop:g}__n{n}")


def enumerate_cells(levels: dict) -> list[dict]:
    """Producto cartesiano de los 5 factores -> lista de dicts de celda."""
    keys = ["sampling", "omega", "sigma", "prop_rhythmic", "n_samples"]
    cells = []
    for combo in itertools.product(*(levels[k] for k in keys)):
        cell = dict(zip(keys, combo))
        cell["cell_id"] = _cell_id(
            cell["sampling"], cell["omega"], cell["sigma"],
            cell["prop_rhythmic"], cell["n_samples"],
        )
        cells.append(cell)
    return cells


# ===========================================================================
# Ejecucion de una celda x replica
# ===========================================================================

def _run_subprocess(cmd: list[str], log_path: Path) -> tuple[bool, float]:
    """Lanza un comando, vuelca stdout+stderr a log_path. Devuelve (ok, segundos)."""
    t0 = time.time()
    with open(log_path, "w") as log:
        proc = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT)
    return proc.returncode == 0, time.time() - t0


def run_cell_rep(
    cell:          dict,
    rep:           int,
    out_root:      Path,
    n_total:       int,
    seed:          int,
    fmm_method:    str,
    robust_method: str,
    aggregation:   str,
    n_reps_pipe:   int,
    skip_pipeline: bool,
    with_cre:      bool,
) -> dict:
    """
    Genera, corre el pipeline y valida una celda x replica. Devuelve una
    fila (dict) con factores, metricas y tiempos, o con ``status`` de error.
    """
    rep_dir = out_root / cell["cell_id"] / f"rep{rep:02d}"
    synth_d = rep_dir / "synth"
    run_d   = rep_dir / "run"
    val_d   = rep_dir / "validate"

    row = {
        "cell_id":       cell["cell_id"],
        "rep":           rep,
        "sampling":      cell["sampling"],
        "omega":         cell["omega"],
        "sigma":         cell["sigma"],
        "prop_rhythmic": cell["prop_rhythmic"],
        "n_samples":     cell["n_samples"],
        "n_total":       n_total,
        "seed":          seed,
        "fmm_method":    fmm_method,
        "robust_method": robust_method,
        "aggregation":   aggregation if robust_method == "aggregate" else None,
        "status":        "ok",
        "t_generate_s":  np.nan,
        "t_pipeline_s":  np.nan,
        "t_validate_s":  np.nan,
    }

    # 1) Generar (in-process)
    t0 = time.time()
    try:
        ds = generate(
            name          = f"{cell['cell_id']}__rep{rep:02d}",
            n_samples     = cell["n_samples"],
            sampling      = cell["sampling"],
            omega         = cell["omega"],
            sigma         = cell["sigma"],
            prop_rhythmic = cell["prop_rhythmic"],
            n_total       = n_total,
            seed          = seed,
        )
        write_dataset(ds, synth_d)
    except Exception as e:  # pragma: no cover - defensivo
        row["status"] = f"generate_error: {e}"
        return row
    row["t_generate_s"] = time.time() - t0

    if skip_pipeline:
        row["status"] = "generate_only"
        return row

    # 2) Pipeline CIRCUST
    pipe_cmd = [sys.executable, str(RUN_PIPELINE),
                "--data", str(synth_d / "matrix.csv"),
                "--label", cell["cell_id"],
                "--fmm-method", fmm_method,
                "--robust-method", robust_method,
                "--aggregation", aggregation,
                "--n-reps", str(n_reps_pipe),
                "--seed", str(seed),
                "-o", str(run_d)]
    ok, secs = _run_subprocess(
        pipe_cmd, run_d.parent / f"rep{rep:02d}_pipeline.log")
    row["t_pipeline_s"] = secs
    if not ok:
        row["status"] = "pipeline_error"
        return row

    # 2b) Timing por etapa desde el manifest del pipeline (perf_counter).
    run_manifest = run_d / "manifest.json"
    if run_manifest.exists():
        try:
            tdict = json.loads(run_manifest.read_text()).get("timing", {})
            for stage, secs_stage in tdict.items():
                name = stage[:-2] if stage.endswith("_s") else stage
                row[f"t_pipe_{name}_s"] = secs_stage
        except (json.JSONDecodeError, OSError):
            pass

    # 3) Validacion
    val_cmd = [sys.executable, str(VALIDATE),
               "--predicted", str(run_d / "sample_order.csv"),
               "--true",      str(synth_d / "true_times.csv"),
               "--label",     cell["cell_id"],
               "-o",          str(val_d)]
    if with_cre:
        val_cmd += ["--expr",  str(synth_d / "matrix.csv"),
                    "--atlas", str(run_d / "circust_atlas.csv")]
    ok, secs = _run_subprocess(val_cmd, run_d.parent / f"rep{rep:02d}_validate.log")
    row["t_validate_s"] = secs
    if not ok:
        row["status"] = "validate_error"
        return row

    # 4) Recoger metricas del summary.csv
    summary_path = val_d / "summary.csv"
    if summary_path.exists():
        summ = pd.read_csv(summary_path).iloc[0].to_dict()
        for k, v in summ.items():
            if k not in ("label", "n_samples"):  # no pisar factores
                row[k] = v
    else:
        row["status"] = "no_summary"

    return row


# ===========================================================================
# CLI
# ===========================================================================

def _parse_list(arg: str | None, cast, default: tuple):
    """Convierte 'a,b,c' en una tupla casteada; None => niveles por defecto."""
    if arg is None:
        return default
    return tuple(cast(x) for x in arg.split(","))


def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=("Orquesta el grid factorial generate->pipeline->validate "
                     "y agrega los resultados en results.csv."),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("-o", "--output", required=True, type=str,
                   help="Directorio raiz de salida.")
    p.add_argument("--reps", type=int, default=5,
                   help="Numero B de replicas por celda. Defecto: 5.")
    p.add_argument("--seed-base", type=int, default=1000,
                   help="Semilla base; la replica r usa seed-base + r. Defecto: 1000.")
    p.add_argument("--n-total", type=int, default=1200,
                   help="Genes totales por dataset. Defecto: 1200.")
    # Submuestreo de factores (para exploratorio). None => nivel completo.
    p.add_argument("--sampling", default=None,
                   help=f"Subconjunto de muestreos (coma). Disp.: {list(SAMPLING_SCHEMES)}.")
    p.add_argument("--omega", default=None,
                   help="Subconjunto de omega (coma).")
    p.add_argument("--sigma", default=None,
                   help="Subconjunto de sigma (coma).")
    p.add_argument("--prop", default=None,
                   help="Subconjunto de prop_rhythmic (coma).")
    p.add_argument("--n", default=None,
                   help="Subconjunto de n_samples (coma).")
    p.add_argument("--limit", type=int, default=None,
                   help="Limita el numero total de celdas (tras enumerar).")
    # Paso a traves al pipeline / validate
    p.add_argument("--fmm-method", choices=["base", "stabilized"], default="base",
                   help="Modelo FMM del pipeline. Defecto: base.")
    p.add_argument("--robust-method", choices=["best_k", "aggregate"], default="best_k",
                   help="Eleccion del orden final: best_k o aggregate. Defecto: best_k.")
    p.add_argument("--aggregation", choices=["auto", "tsp3", "hods"], default="auto",
                   help="Metodo de agregacion (solo si robust-method=aggregate). Defecto: auto.")
    p.add_argument("--n-reps-pipeline", type=int, default=5,
                   help="K repeticiones de RobustOrder en el pipeline. Defecto: 5.")
    p.add_argument("--skip-pipeline", action="store_true",
                   help="Solo genera datasets (mide coste de generacion).")
    p.add_argument("--resume", action="store_true",
                   help="Salta celdas cuya validate/summary.csv ya exista.")
    p.add_argument("--dry-run", action="store_true",
                   help="Enumera y cuenta celdas, no ejecuta nada.")
    return p.parse_args()


def main() -> None:
    args = _parse_args()
    out_root = Path(args.output)

    levels = {
        "sampling":      _parse_list(args.sampling, str,   FACTOR_LEVELS["sampling"]),
        "omega":         _parse_list(args.omega,    float, FACTOR_LEVELS["omega"]),
        "sigma":         _parse_list(args.sigma,    float, FACTOR_LEVELS["sigma"]),
        "prop_rhythmic": _parse_list(args.prop,     float, FACTOR_LEVELS["prop_rhythmic"]),
        "n_samples":     _parse_list(args.n,        int,   FACTOR_LEVELS["n_samples"]),
    }
    cells = enumerate_cells(levels)
    if args.limit is not None:
        cells = cells[:args.limit]

    n_cells = len(cells)
    n_runs  = n_cells * args.reps

    print("╔══════════════════════════════════════════════════════════════════╗")
    print("  Grid de simulacion CIRCUST")
    for k, v in levels.items():
        print(f"    {k:14s}: {list(v)}")
    print(f"    {'reps (B)':14s}: {args.reps}")
    print(f"    {'celdas':14s}: {n_cells}")
    print(f"    {'runs totales':14s}: {n_runs}  (celdas x reps)")
    print(f"    {'n_total genes':14s}: {args.n_total}")
    print(f"    {'fmm-method':14s}: {args.fmm_method}")
    print(f"    {'robust-method':14s}: {args.robust_method}"
          + (f" ({args.aggregation})" if args.robust_method == "aggregate" else ""))
    print(f"    {'output':14s}: {out_root.resolve()}")
    print("╚══════════════════════════════════════════════════════════════════╝")

    if args.dry_run:
        print("\n  [dry-run] No se ejecuta nada.")
        return

    out_root.mkdir(parents=True, exist_ok=True)
    rows: list[dict] = []
    t_start = time.time()
    done = 0

    for cell in cells:
        for rep in range(args.reps):
            done += 1
            seed = args.seed_base + rep
            val_summary = (out_root / cell["cell_id"] / f"rep{rep:02d}"
                           / "validate" / "summary.csv")
            tag = f"[{done}/{n_runs}] {cell['cell_id']} rep{rep:02d}"

            if args.resume and val_summary.exists():
                print(f"  {tag}  (skip: ya existe)")
                summ = pd.read_csv(val_summary).iloc[0].to_dict()
                row = {**cell, "rep": rep, "n_total": args.n_total,
                       "seed": seed,
                       "fmm_method": args.fmm_method,
                       "robust_method": args.robust_method,
                       "aggregation": (args.aggregation
                                       if args.robust_method == "aggregate" else None),
                       "status": "cached"}
                row.update({k: v for k, v in summ.items()
                            if k not in ("label", "n_samples")})
                rows.append(row)
                continue

            print(f"  {tag} ...", flush=True)
            row = run_cell_rep(
                cell=cell, rep=rep, out_root=out_root,
                n_total=args.n_total, seed=seed,
                fmm_method=args.fmm_method,
                robust_method=args.robust_method,
                aggregation=args.aggregation,
                n_reps_pipe=args.n_reps_pipeline,
                skip_pipeline=args.skip_pipeline,
                with_cre=True,   # CRE siempre activo
            )
            rows.append(row)
            elapsed = time.time() - t_start
            eta = elapsed / done * (n_runs - done)
            print(f"      status={row['status']}  "
                  f"t_gen={row.get('t_generate_s', float('nan')):.1f}s "
                  f"t_pipe={row.get('t_pipeline_s', float('nan')):.1f}s "
                  f"t_val={row.get('t_validate_s', float('nan')):.1f}s  "
                  f"| ETA {eta/60:.1f} min")

            # Persistir parcial cada vez (robustez ante cortes).
            pd.DataFrame(rows).to_csv(out_root / "results.csv", index=False)

    results = pd.DataFrame(rows)
    results.to_csv(out_root / "results.csv", index=False)

    elapsed = time.time() - t_start
    n_ok = int((results["status"] == "ok").sum())
    print()
    print(f"  ✓ Grid completado en {elapsed/60:.1f} min")
    print(f"    {n_ok}/{len(results)} runs OK")
    print(f"    Agregado: {out_root / 'results.csv'}")


if __name__ == "__main__":
    main()
