#!/usr/bin/env python3
"""
analyze_results.py
==================
Capa de **reporte** de la validacion comparativa del nuevo CIRCUST. Produce
las tablas y figuras de los tres experimentos del TFG a partir de los
artefactos generados por el grid (``run_simulation_grid.py``) y los
benchmarks de tiempo (``benchmark_timing.R`` / ``benchmark_python.py``):

- **Exp 1 (R vs Python)**  : tiempos por etapa y total sobre la matriz canonica.
- **Exp 2 (FMM base vs robusto)** : G1(stabilized) vs G0(base), pareado por celda.
- **Exp 3 (best-K vs TSP-3)**     : G2(tsp3)       vs G0(best_k), pareado por celda.

Las comparaciones 2 y 3 son **pareadas**: como el grid usa la misma semilla por
celda/replica, cada (cell_id, rep) de G0/G1/G2 comparte el dataset, asi que el
delta de metricas aisla el efecto del cambio algoritmico.

Uso (desde la raiz del repo)
----------------------------
    python validation/analyze_results.py
    python validation/analyze_results.py \
        --base output/grid/base_bestk/results.csv \
        --fmm  output/grid/stab_bestk/results.csv \
        --agg  output/grid/base_tsp3/results.csv \
        -o output/analysis
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PROJECT_ROOT = Path(__file__).resolve().parent.parent

# Factores del grid (para desgloses)
FACTORS = ["sampling", "omega", "sigma", "prop_rhythmic", "n_samples"]

# Metricas y su sentido (True = mejor cuanto MAYOR)
METRICS = {
    "MedAE":         False,
    "SDAE":          False,
    "pct_within_1h": True,
    "pct_within_2h": True,
    "CCC":           True,
    "rTrim":         True,
    "CRE":           False,
    "t_pipeline_s":  False,
}

# Crosswalk de etapas: cada etapa R agrega 1+ etapas Python equivalentes.
STAGE_CROSSWALK = {
    "preprocess":        ["load", "preprocess", "core_genes", "cpca"],
    "preliminary_order": ["preliminary_order"],
    "top_genes":         ["top_genes"],
    "robust_order":      ["robust_order", "build_atlas"],
    "elapsed_total":     ["elapsed_total"],
}


# ===========================================================================
# Carga
# ===========================================================================

def _load(path: Path, label: str) -> pd.DataFrame | None:
    if not path.exists():
        print(f"  [aviso] no existe {label}: {path} — se omite.")
        return None
    df = pd.read_csv(path)
    print(f"  cargado {label}: {len(df)} filas ({path})")
    return df


# ===========================================================================
# Exp 1 — R vs Python (tiempos)
# ===========================================================================

def exp1_timing(timing_r: pd.DataFrame, timing_py: pd.DataFrame,
                out: Path) -> None:
    print("\n── Exp 1: R vs Python (tiempos por etapa) ─────────────")

    # Python: sumar las etapas mapeadas por cada etapa R, por replica
    py_wide = timing_py.pivot_table(index="rep", columns="stage",
                                    values="seconds", aggfunc="sum")
    py_stage_means, py_stage_sds = {}, {}
    for r_stage, py_stages in STAGE_CROSSWALK.items():
        present = [s for s in py_stages if s in py_wide.columns]
        agg = py_wide[present].sum(axis=1) if present else pd.Series(dtype=float)
        py_stage_means[r_stage] = agg.mean()
        py_stage_sds[r_stage]   = agg.std(ddof=0)

    # R: media/sd por etapa
    r_grp = timing_r.groupby("stage")["seconds"]
    r_means, r_sds = r_grp.mean().to_dict(), r_grp.std(ddof=0).to_dict()

    rows = []
    for stage in STAGE_CROSSWALK:
        rm, pm = r_means.get(stage, np.nan), py_stage_means.get(stage, np.nan)
        rows.append({
            "stage":      stage,
            "R_mean_s":   rm, "R_sd_s":  r_sds.get(stage, np.nan),
            "Py_mean_s":  pm, "Py_sd_s": py_stage_sds.get(stage, np.nan),
            "speedup_RxPy": (rm / pm) if (pm and pm > 0) else np.nan,
        })
    tab = pd.DataFrame(rows)
    tab.to_csv(out / "exp1_timing.csv", index=False)
    print(tab.to_string(index=False))

    # Figura: barras R vs Python por etapa (sin total)
    try:
        st = tab[tab["stage"] != "elapsed_total"]
        x = np.arange(len(st)); w = 0.38
        fig, ax = plt.subplots(figsize=(9, 5))
        ax.bar(x - w/2, st["R_mean_s"],  w, yerr=st["R_sd_s"],  label="R (antiguo)",  capsize=3)
        ax.bar(x + w/2, st["Py_mean_s"], w, yerr=st["Py_sd_s"], label="Python (nuevo)", capsize=3)
        ax.set_xticks(x); ax.set_xticklabels(st["stage"], rotation=20, ha="right")
        ax.set_ylabel("segundos"); ax.set_title("Exp 1: tiempo por etapa, R vs Python")
        ax.legend()
        fig.tight_layout(); fig.savefig(out / "exp1_timing.png", dpi=130)
        plt.close(fig)
        print(f"  figura: {out / 'exp1_timing.png'}")
    except Exception as e:  # pragma: no cover
        print(f"  [aviso] figura Exp 1 fallo: {e}")


# ===========================================================================
# Exp 2 / 3 — comparaciones pareadas
# ===========================================================================

def paired_compare(df_base: pd.DataFrame, df_var: pd.DataFrame,
                   var_label: str, prefix: str, out: Path) -> None:
    print(f"\n── {prefix}: {var_label} vs base (pareado) ───────────")

    keys = ["cell_id", "rep"]
    metric_cols = [m for m in METRICS if m in df_base.columns and m in df_var.columns]
    cols = keys + FACTORS + metric_cols

    merged = df_base[cols].merge(
        df_var[keys + metric_cols], on=keys, suffixes=("_base", "_var"))
    if merged.empty:
        print(f"  [aviso] sin celdas pareadas para {prefix}; se omite.")
        return

    # Deltas (var - base) por metrica
    for m in metric_cols:
        merged[f"d_{m}"] = merged[f"{m}_var"] - merged[f"{m}_base"]
    merged.to_csv(out / f"{prefix}_paired.csv", index=False)

    # Resumen global: delta medio + tasa de victorias del variante
    summary = []
    for m in metric_cols:
        higher_better = METRICS[m]
        d = merged[f"d_{m}"]
        win = (d > 0) if higher_better else (d < 0)
        summary.append({
            "metric":       m,
            "base_mean":    merged[f"{m}_base"].mean(),
            "var_mean":     merged[f"{m}_var"].mean(),
            "delta_mean":   d.mean(),
            "delta_median": d.median(),
            "win_rate_var": win.mean(),
            "higher_better": higher_better,
        })
    summ = pd.DataFrame(summary)
    summ.to_csv(out / f"{prefix}_summary.csv", index=False)
    print(summ.to_string(index=False))

    # Desglose de MedAE (precision) y t_pipeline_s (tiempo) por factor
    for m in ("MedAE", "t_pipeline_s"):
        if f"d_{m}" not in merged.columns:
            continue
        for f in FACTORS:
            by = merged.groupby(f)[f"d_{m}"].agg(["mean", "std", "count"])
            by.to_csv(out / f"{prefix}_{m}_by_{f}.csv")

    # Figura: delta medio de MedAE por nivel de cada factor
    try:
        if "d_MedAE" in merged.columns:
            fig, axes = plt.subplots(1, len(FACTORS), figsize=(4 * len(FACTORS), 4),
                                     sharey=True)
            for ax, f in zip(np.atleast_1d(axes), FACTORS):
                g = merged.groupby(f)["d_MedAE"].mean()
                ax.bar([str(i) for i in g.index], g.values)
                ax.axhline(0, color="k", lw=0.8)
                ax.set_title(f); ax.set_xlabel(f)
                ax.tick_params(axis="x", rotation=30)
            axes_list = np.atleast_1d(axes)
            axes_list[0].set_ylabel(f"Δ MedAE (h)  [{var_label} − base]")
            fig.suptitle(f"{prefix}: Δ MedAE por factor  (negativo = {var_label} mejora)")
            fig.tight_layout(); fig.savefig(out / f"{prefix}_dMedAE_by_factor.png", dpi=130)
            plt.close(fig)
            print(f"  figura: {out / f'{prefix}_dMedAE_by_factor.png'}")
    except Exception as e:  # pragma: no cover
        print(f"  [aviso] figura {prefix} fallo: {e}")


# ===========================================================================
# CLI
# ===========================================================================

def _parse_args() -> argparse.Namespace:
    g = PROJECT_ROOT / "output" / "grid"
    ref = PROJECT_ROOT / "validation" / "reference"
    p = argparse.ArgumentParser(
        description="Tablas y figuras de los 3 experimentos de validacion.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--base", default=str(g / "base_bestk" / "results.csv"),
                   help="results.csv de G0 (base + best_k).")
    p.add_argument("--fmm",  default=str(g / "stab_bestk" / "results.csv"),
                   help="results.csv de G1 (stabilized + best_k).")
    p.add_argument("--agg",  default=str(g / "base_tsp3" / "results.csv"),
                   help="results.csv de G2 (base + aggregate/tsp3).")
    p.add_argument("--timing-r",  default=str(ref / "timing_R.csv"),
                   help="timing_R.csv del benchmark R.")
    p.add_argument("--timing-py", default=str(ref / "timing_py.csv"),
                   help="timing_py.csv del benchmark Python.")
    p.add_argument("-o", "--output", default=str(PROJECT_ROOT / "output" / "analysis"),
                   help="Directorio de salida. Defecto: output/analysis.")
    return p.parse_args()


def main() -> None:
    args = _parse_args()
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)

    print("╔══════════════════════════════════════════════════════════════════╗")
    print("  Analisis de resultados de validacion CIRCUST")
    print("╚══════════════════════════════════════════════════════════════════╝")

    df_base = _load(Path(args.base),     "G0 base+best_k")
    df_fmm  = _load(Path(args.fmm),      "G1 stabilized+best_k")
    df_agg  = _load(Path(args.agg),      "G2 base+tsp3")
    t_r     = _load(Path(args.timing_r), "timing_R")
    t_py    = _load(Path(args.timing_py),"timing_py")

    if t_r is not None and t_py is not None:
        exp1_timing(t_r, t_py, out)

    if df_base is not None and df_fmm is not None:
        paired_compare(df_base, df_fmm, "stabilized", "exp2_fmm", out)

    if df_base is not None and df_agg is not None:
        paired_compare(df_base, df_agg, "tsp3", "exp3_agg", out)

    print(f"\n  ✓ Analisis en {out.resolve()}")


if __name__ == "__main__":
    main()
