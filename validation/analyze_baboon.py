#!/usr/bin/env python3
"""
analyze_baboon.py
=================
Capa de reporte de la validacion de CIRCUST sobre el atlas de babuino
(Mure 2018), a partir del ``results.csv`` que produce
``run_baboon_validation.py``. Reproduce el tipo de analisis que hizo el paper
original (Larriba et al. 2023, Fig. 7A): la **distribucion entre tejidos** de
la calidad del orden temporal inferido contra el ZT real.

Dos modos
---------
* **Agregado (1 results.csv):** distribucion entre tejidos de la correlacion
  circular estimado-vs-ZT (``CCC``), ``MedAE``, ``% within 2h``, ``AUC`` y
  ``CRE``. Escribe ``baboon_summary.csv`` (mediana/IQR + recuentos por umbral)
  y ``baboon_by_tissue.csv`` (ranking). Figuras: scatter fase estimada vs ZT
  (estilo Fig 7A) para los tejidos mejor ordenados + boxplot de CCC/MedAE.

* **Contraste (--compare A.csv B.csv):** comparacion **pareada por tejido**
  (mismo ZT => comparacion justa) con deltas y win-rate, para contrastar
  configuraciones (base vs stabilized, best_k vs aggregate). Calca la logica de
  ``analyze_results.paired_compare`` pero con clave = tejido. Ademas de los CSV
  (``_paired.csv``, ``_summary.csv``) genera tres figuras autoexplicativas para
  la seccion de validacion del TFG:
    - ``_scatter.png``        : scatter pareado base-vs-variante (MedAE y rho),
                                un punto por tejido, zona de mejora sombreada.
    - ``_delta_by_tissue.png``: barras de Delta por tejido (verde mejora / rojo
                                empeora), ordenadas por mejora en MedAE.
    - ``_summary.png``        : win-rate por metrica con la mediana del Delta
                                anotada (vista de una ojeada del summary).

Nota metodologica (para la defensa)
-----------------------------------
El ZT del babuino es equiespaciado por diseno (12 puntos cada 2 h), asi que
``R_bar(true) ≈ 0`` y el heuristico ``CCC_reliability`` de validate.py marca los
tejidos como "inestable". Eso mide concentracion de los tiempos reales, no
fiabilidad de la metrica: la correlacion circular de Jammalamadaka (``CCC``) y
el ``MedAE`` son validos y son justo los que reporto el paper. Aqui se titula
con ``CCC``/``MedAE``/%within y se ignora esa etiqueta.

Uso
---
    # Reporte agregado (config canonica)
    python validation/analyze_baboon.py output/baboon/base_bestk/results.csv \\
        -o output/baboon/analysis

    # Contraste pareado entre dos configuraciones
    python validation/analyze_baboon.py --compare \\
        output/baboon/base_bestk/results.csv output/baboon/base_stab/results.csv \\
        -o output/baboon/analysis
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "validation"))

from metrics import circular_regression_rho

# Metricas del babuino y su sentido (True = mejor cuanto MAYOR).
#
# IMPORTANTE: el ZT del babuino es equiespaciado por diseno (12 puntos cada 2 h)
# => R_bar(true) ≈ 0. En ese regimen la correlacion circular de Jammalamadaka
# (CCC) es NUMERICAMENTE INESTABLE y puede dar valores negativos incluso para
# tejidos casi perfectamente ordenados (la media circular de los tiempos reales
# queda mal definida; ver docstring de metrics.circular_correlation). Por eso el
# reporte LIDERA con metricas de error absoluto (MedAE, SDAE, %within), que son
# robustas a este diseno, y deja CCC/rTrim como secundarias informativas.
PRIMARY = "MedAE"   # metrica de ordenacion y seleccion (menor = mejor)
METRICS = {
    "MedAE":         False,
    "rho_cc":        True,    # GOF regresion circular-circular = METRICA DEL PAPER
    "SDAE":          False,
    "pct_within_2h": True,
    "pct_within_1h": True,
    "AUC":           True,
    "rTrim":         True,    # correlacion circular robusta (Mahmood 2022)
    "CCC":           True,    # Jammalamadaka — INESTABLE con ZT uniforme
    "CRE":           False,
    "t_pipeline_s":  False,
}

MEDAE_THRESHOLDS = (1.0, 2.0, 3.0)


# ===========================================================================
# Carga
# ===========================================================================

def _load_results(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise SystemExit(f"No existe results.csv: {path}")
    df = pd.read_csv(path)
    print(f"  cargado {path}: {len(df)} tejidos "
          f"({int((df.get('status') == 'ok').sum())} OK)")
    return df


def _ok(df: pd.DataFrame) -> pd.DataFrame:
    """Filtra a los tejidos con status ok/cached y MedAE finito (la metrica
    primaria; CCC se ignora aqui por ser inestable con ZT uniforme)."""
    if "status" in df.columns:
        df = df[df["status"].isin(["ok", "cached"])]
    return df[np.isfinite(pd.to_numeric(df.get(PRIMARY), errors="coerce"))].copy()


def _ensure_rho(df: pd.DataFrame, config_dir: Path) -> pd.DataFrame:
    """Garantiza la columna ``rho_cc`` (GOF de regresion circular-circular).

    Las corridas previas a esta metrica no la traen en ``results.csv``, asi que
    se recalcula por tejido desde su ``per_sample.csv`` (columnas ``true_h`` y
    ``predicted_h``) sin necesidad de re-ejecutar el pipeline."""
    if "rho_cc" not in df.columns:
        df["rho_cc"] = np.nan
    for idx, row in df.iterrows():
        if np.isfinite(pd.to_numeric(row.get("rho_cc"), errors="coerce")):
            continue
        ps = config_dir / str(row.get("tissue")) / "validate" / "per_sample.csv"
        if not ps.exists():
            continue
        d = pd.read_csv(ps)
        if {"true_h", "predicted_h"}.issubset(d.columns):
            df.at[idx, "rho_cc"] = circular_regression_rho(
                d["predicted_h"].to_numpy(dtype=float),
                d["true_h"].to_numpy(dtype=float),
            )
    return df


# ===========================================================================
# Modo agregado
# ===========================================================================

def aggregate_report(results_path: Path, out: Path, top_n: int = 9) -> None:
    df_all = _ensure_rho(_load_results(results_path), results_path.parent)
    df = _ok(df_all)
    n_total = len(df_all)
    n_ok    = len(df)
    print(f"\n── Distribucion entre tejidos ({n_ok}/{n_total} OK) ───────")

    # 1) Resumen por metrica (distribucion entre tejidos)
    metric_cols = [m for m in METRICS if m in df.columns]
    rows = []
    for m in metric_cols:
        s = pd.to_numeric(df[m], errors="coerce").dropna()
        if s.empty:
            continue
        rows.append({
            "metric":  m,
            "n":       int(s.size),
            "median":  s.median(),
            "q25":     s.quantile(0.25),
            "q75":     s.quantile(0.75),
            "mean":    s.mean(),
            "min":     s.min(),
            "max":     s.max(),
        })
    summary = pd.DataFrame(rows)

    # 2) Recuentos por umbral de calidad del orden (metrica primaria: MedAE)
    counts = {"n_ok_tissues": n_ok, "n_total_tissues": n_total}
    if "MedAE" in df.columns:
        med = pd.to_numeric(df["MedAE"], errors="coerce")
        for thr in MEDAE_THRESHOLDS:
            counts[f"n_MedAE_le_{thr:g}h"] = int((med <= thr).sum())
    if "pct_within_2h" in df.columns:
        pw = pd.to_numeric(df["pct_within_2h"], errors="coerce")
        counts["n_pct_within_2h_ge_50"] = int((pw >= 50.0).sum())

    summary.to_csv(out / "baboon_summary.csv", index=False)
    pd.DataFrame([counts]).to_csv(out / "baboon_counts.csv", index=False)
    print(summary.to_string(index=False))
    print("  " + "  ".join(f"{k}={v}" for k, v in counts.items()))
    print("  [nota] rho_cc = GOF de regresion circular-circular (metrica del "
          "paper), robusta al ZT equiespaciado. CCC (Jammalamadaka) es "
          "inestable en este diseno; se lidera con MedAE/rho_cc/%within.")

    # 3) Ranking por tejido (mejor recuperacion = MedAE menor)
    keep = ["tissue"] + [c for c in metric_cols if c in df.columns] + \
           [c for c in ("n_samples", "status") if c in df.columns]
    by_tissue = (df[keep].sort_values(PRIMARY, ascending=True)
                 .reset_index(drop=True))
    by_tissue.to_csv(out / "baboon_by_tissue.csv", index=False)

    # 4) Figura A — boxplot MedAE y rho_cc entre tejidos (metricas robustas).
    #    rho_cc es la GOF de regresion circular-circular del paper (Tabla A).
    try:
        fig, axes = plt.subplots(1, 2, figsize=(8, 4))
        if "MedAE" in df.columns:
            axes[0].boxplot(pd.to_numeric(df["MedAE"], errors="coerce").dropna(),
                            vert=True, showmeans=True)
            axes[0].set_title("Error mediano absoluto vs ZT")
            axes[0].set_ylabel("MedAE (h)")
            axes[0].axhline(6.0, color="r", lw=0.6, ls=":")  # azar ~ 6 h
        if "rho_cc" in df.columns:
            axes[1].boxplot(pd.to_numeric(df["rho_cc"], errors="coerce").dropna(),
                            vert=True, showmeans=True)
            axes[1].set_title("GOF regresion circular-circular")
            axes[1].set_ylabel(r"$\rho$ (paper, Tabla A)")
            axes[1].set_ylim(0, 1)
            axes[1].axhline(0, color="k", lw=0.6)
        fig.suptitle(f"CIRCUST sobre babuino — distribucion entre {n_ok} tejidos")
        fig.tight_layout()
        fig.savefig(out / "baboon_distribution.png", dpi=130)
        plt.close(fig)
        print(f"  figura: {out / 'baboon_distribution.png'}")
    except Exception as e:  # pragma: no cover
        print(f"  [aviso] figura distribucion fallo: {e}")

    # 5) Figura B — scatter fase estimada vs ZT (estilo Fig 7A), top tejidos
    try:
        _fig7a_scatter(by_tissue, results_path.parent, out, top_n=top_n)
    except Exception as e:  # pragma: no cover
        print(f"  [aviso] figura Fig7A fallo: {e}")

    print(f"\n  ✓ Reporte agregado en {out.resolve()}")


def _fig7a_scatter(by_tissue: pd.DataFrame, config_dir: Path,
                   out: Path, top_n: int) -> None:
    """Scatter fase estimada (alineada) vs ZT real para los top_n tejidos mejor
    ordenados (MedAE menor), leyendo cada per_sample.csv. La diagonal indica
    recuperacion correcta del orden temporal (admite reflexion, ya alineada por
    validate.py). by_tissue llega ordenado por MedAE ascendente."""
    tops = by_tissue.head(top_n)
    n = len(tops)
    if n == 0:
        return
    ncol = min(3, n)
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.2 * ncol, 3.2 * nrow),
                             squeeze=False)
    for ax in axes.flat:
        ax.set_visible(False)

    for ax, (_, r) in zip(axes.flat, tops.iterrows()):
        ps = config_dir / r["tissue"] / "validate" / "per_sample.csv"
        if not ps.exists():
            continue
        d = pd.read_csv(ps)
        ax.set_visible(True)
        ax.scatter(d["true_h"], d["aligned_pred_h"], s=28, alpha=0.85)
        ax.plot([0, 24], [0, 24], color="k", lw=0.7, ls="--")
        ax.set_xlim(-1, 25); ax.set_ylim(-1, 25)
        ax.set_xticks([0, 6, 12, 18, 24]); ax.set_yticks([0, 6, 12, 18, 24])
        rho = r["rho_cc"] if "rho_cc" in r and np.isfinite(r["rho_cc"]) else float("nan")
        ax.set_title(f"{r['tissue']}  (MedAE={r['MedAE']:.1f}h, "
                     rf"$\rho$={rho:.2f})", fontsize=9)

    # Etiquetas comunes
    for ax in axes[-1, :]:
        if ax.get_visible():
            ax.set_xlabel("ZT real (h)")
    for ax in axes[:, 0]:
        if ax.get_visible():
            ax.set_ylabel("fase estimada (h)")
    fig.suptitle("CIRCUST sobre babuino — fase estimada vs ZT (top tejidos)")
    fig.tight_layout()
    fig.savefig(out / "baboon_fig7a_scatter.png", dpi=130)
    plt.close(fig)
    print(f"  figura: {out / 'baboon_fig7a_scatter.png'}")


# ===========================================================================
# Modo contraste (pareado por tejido)
# ===========================================================================

def _config_tag(df: pd.DataFrame) -> str:
    if df.empty:
        return "config"
    r = df.iloc[0]
    tag = f"{r.get('fmm_method', '?')}_{r.get('robust_method', '?')}"
    if r.get("robust_method") == "aggregate":
        tag += f"_{r.get('aggregation', '?')}"
    return tag


# ---------------------------------------------------------------------------
# Etiquetas legibles de metricas para las figuras del TFG
# ---------------------------------------------------------------------------
_METRIC_LABEL = {
    "MedAE":         "MedAE (h)",
    "rho_cc":        r"$\rho$ (regresion circular-circular)",
    "SDAE":          "SDAE (h)",
    "pct_within_2h": "% muestras < 2 h",
    "pct_within_1h": "% muestras < 1 h",
    "AUC":           "AUC del CDF de error",
    "rTrim":         "rTrim",
    "CCC":           "CCC (Jammalamadaka)*",
    "CRE":           "CRE",
    "t_pipeline_s":  "tiempo pipeline (s)",
}
# Verde = la variante mejora; rojo = empeora; gris = empate.
_C_BETTER, _C_WORSE, _C_TIE = "#2ca02c", "#d62728", "#9e9e9e"


def _improves(delta: np.ndarray, higher_better: bool) -> np.ndarray:
    """Mascara booleana: True donde la variante mejora a la base."""
    return (delta > 0) if higher_better else (delta < 0)


def _fig_paired_scatter(merged, base_lab, var_lab, out, prefix) -> None:
    """Scatter pareado base-vs-variante (un punto por tejido) para las dos
    metricas titulares. Punto bajo la diagonal (MedAE) o sobre la diagonal
    (rho) = la variante mejora. La region de mejora va sombreada en verde."""
    headline = [m for m in ("MedAE", "rho_cc")
                if f"{m}_base" in merged.columns and f"{m}_var" in merged.columns]
    if not headline:
        return
    fig, axes = plt.subplots(1, len(headline), figsize=(5.0 * len(headline), 4.8),
                             squeeze=False)
    for ax, m in zip(axes[0], headline):
        hb = METRICS[m]
        b  = pd.to_numeric(merged[f"{m}_base"], errors="coerce").to_numpy(float)
        v  = pd.to_numeric(merged[f"{m}_var"],  errors="coerce").to_numpy(float)
        ok = np.isfinite(b) & np.isfinite(v)
        b, v = b[ok], v[ok]
        if b.size == 0:
            continue
        better = _improves(v - b, hb)
        lo = float(min(b.min(), v.min())); hi = float(max(b.max(), v.max()))
        pad = 0.05 * (hi - lo + 1e-9); lo -= pad; hi += pad

        # Region sombreada donde la variante es mejor
        if hb:   # mejor por encima de la diagonal
            ax.fill_between([lo, hi], [lo, hi], [hi, hi], color=_C_BETTER,
                            alpha=0.12, zorder=0)
        else:    # mejor por debajo de la diagonal
            ax.fill_between([lo, hi], [lo, lo], [lo, hi], color=_C_BETTER,
                            alpha=0.12, zorder=0)
        ax.plot([lo, hi], [lo, hi], "k--", lw=0.9, zorder=1)
        ax.scatter(b[better],  v[better],  s=30, c=_C_BETTER, edgecolor="k",
                   lw=0.3, zorder=3, label=f"mejora ({int(better.sum())})")
        ax.scatter(b[~better], v[~better], s=30, c=_C_WORSE, edgecolor="k",
                   lw=0.3, zorder=3, label=f"empeora/iguala ({int((~better).sum())})")
        ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel(f"{_METRIC_LABEL.get(m, m)}\n{base_lab}", fontsize=9)
        ax.set_ylabel(f"{var_lab}\n{_METRIC_LABEL.get(m, m)}", fontsize=9)
        dmed = float(np.median(v - b))
        sentido = "mayor = mejor" if hb else "menor = mejor"
        ax.set_title(f"{m}  ({sentido})\n"
                     f"gana en {int(better.sum())}/{b.size} tejidos · "
                     f"$\\Delta$mediana = {dmed:+.3f}", fontsize=9.5)
        ax.legend(fontsize=7.5, loc="upper left" if hb else "lower right",
                  framealpha=0.9)
    fig.suptitle(f"Contraste pareado por tejido: {var_lab} vs {base_lab}\n"
                 f"cada punto es un tejido (mismo ZT real); zona verde = la "
                 f"variante mejora", fontsize=11)
    fig.tight_layout()
    fig.savefig(out / f"{prefix}_scatter.png", dpi=150)
    plt.close(fig)


def _fig_delta_by_tissue(merged, base_lab, var_lab, out, prefix) -> None:
    """Barras horizontales de Delta por tejido (variante - base) para MedAE y
    rho. Tejidos ordenados por mejora en MedAE (mismo orden en ambos paneles).
    Verde = la variante mejora ese tejido, rojo = lo empeora."""
    headline = [m for m in ("MedAE", "rho_cc") if f"d_{m}" in merged.columns]
    if "MedAE" not in headline:
        return
    # Orden compartido: mejora en MedAE (mas negativo = mejor) arriba
    order = merged.assign(_k=merged["d_MedAE"]).sort_values(
        "_k", ascending=False).reset_index(drop=True)
    tissues = order["tissue"].to_numpy()
    y = np.arange(len(order))

    fig, axes = plt.subplots(1, len(headline),
                             figsize=(4.6 * len(headline), max(6.0, 0.20 * len(order))),
                             squeeze=False, sharey=True)
    for ax, m in zip(axes[0], headline):
        hb = METRICS[m]
        d  = pd.to_numeric(order[f"d_{m}"], errors="coerce").to_numpy(float)
        imp = _improves(d, hb)
        colors = [(_C_BETTER if i else (_C_TIE if dd == 0 else _C_WORSE))
                  for i, dd in zip(imp, d)]
        ax.barh(y, d, color=colors, edgecolor="k", lw=0.2)
        ax.axvline(0, color="k", lw=0.9)
        ax.set_title(f"$\\Delta${m} = {var_lab} - {base_lab}\n"
                     f"({'mayor' if hb else 'menor'} = mejora · "
                     f"{int(imp.sum())}/{len(d)} mejoran)", fontsize=9.5)
        ax.set_xlabel(f"$\\Delta${_METRIC_LABEL.get(m, m)}", fontsize=9)
    axes[0][0].set_yticks(y)
    axes[0][0].set_yticklabels(tissues, fontsize=5.5)
    axes[0][0].set_ylim(-1, len(order))
    axes[0][0].invert_yaxis()  # mejor (arriba)
    # Titulo arriba del todo y leyenda justo debajo (sin solaparse)
    handles = [plt.Rectangle((0, 0), 1, 1, color=_C_BETTER),
               plt.Rectangle((0, 0), 1, 1, color=_C_WORSE)]
    top = 1.0 - 0.55 / fig.get_figheight()   # reserva ~0.55in para cabecera
    fig.suptitle(f"Efecto por tejido: {var_lab} vs {base_lab}",
                 y=0.997, fontsize=11)
    fig.legend(handles, ["la variante mejora el tejido", "la variante lo empeora"],
               loc="upper center", bbox_to_anchor=(0.5, top + 0.012),
               ncol=2, fontsize=8, frameon=False)
    fig.tight_layout(rect=(0, 0, 1, top))
    fig.savefig(out / f"{prefix}_delta_by_tissue.png", dpi=150)
    plt.close(fig)


def _fig_summary_overview(summ, base_lab, var_lab, out, prefix) -> None:
    """Vista de una ojeada del summary: para cada metrica de calidad, fraccion
    de tejidos donde la variante gana (win-rate) con la mediana del Delta
    anotada. Resume visualmente compare_*_summary.csv."""
    drop = {"t_pipeline_s", "CRE"}   # coste y CRE (vacio con --no-cre) aparte
    q = summ[~summ["metric"].isin(drop)].dropna(subset=["win_rate_var"]).copy()
    if q.empty:
        return
    q = q.iloc[::-1]   # primera metrica arriba
    y = np.arange(len(q))
    wr = q["win_rate_var"].to_numpy(float)
    colors = [_C_BETTER if w > 0.5 else (_C_TIE if w == 0.5 else _C_WORSE)
              for w in wr]
    fig, ax = plt.subplots(figsize=(8.2, 0.6 * len(q) + 1.6))
    ax.barh(y, wr, color=colors, edgecolor="k", lw=0.3)
    ax.axvline(0.5, color="k", ls="--", lw=1.0)
    ax.text(0.5, len(q) - 0.3, "sin diferencia (0.5)", rotation=90,
            va="top", ha="right", fontsize=7.5, color="0.3")
    ax.set_xlim(0, 1)
    ax.set_yticks(y)
    ax.set_yticklabels([_METRIC_LABEL.get(m, m) for m in q["metric"]], fontsize=9)
    ax.set_xlabel(f"fraccion de tejidos donde {var_lab} supera a {base_lab}")
    # Anotar la mediana del delta de cada metrica
    for yi, (_, row) in zip(y, q.iterrows()):
        ax.text(min(row["win_rate_var"] + 0.02, 0.98), yi,
                f"$\\Delta$med={row['delta_median']:+.3f}",
                va="center", ha="left", fontsize=7.5)
    ax.set_title(f"Resumen del contraste: {var_lab} vs {base_lab}\n"
                 f"verde = favorece a la variante · rojo = favorece a la base",
                 fontsize=10.5)
    ax.annotate("* CCC es inestable con el ZT equiespaciado del babuino "
                "(no usar como titular)", xy=(0, -0.5), xycoords=("axes fraction",
                "data"), fontsize=7, color="0.35")
    fig.tight_layout()
    fig.savefig(out / f"{prefix}_summary.png", dpi=150)
    plt.close(fig)


def compare_report(base_path: Path, var_path: Path, out: Path) -> None:
    df_base = _ok(_ensure_rho(_load_results(base_path), base_path.parent))
    df_var  = _ok(_ensure_rho(_load_results(var_path),  var_path.parent))
    base_lab, var_lab = _config_tag(df_base), _config_tag(df_var)
    prefix = f"compare_{var_lab}_vs_{base_lab}"
    print(f"\n── Contraste pareado por tejido: {var_lab} vs {base_lab} ──")

    metric_cols = [m for m in METRICS
                   if m in df_base.columns and m in df_var.columns]
    merged = df_base[["tissue"] + metric_cols].merge(
        df_var[["tissue"] + metric_cols], on="tissue",
        suffixes=("_base", "_var"))
    if merged.empty:
        print("  [aviso] sin tejidos pareados; se omite.")
        return

    for m in metric_cols:
        merged[f"d_{m}"] = merged[f"{m}_var"] - merged[f"{m}_base"]
    merged.to_csv(out / f"{prefix}_paired.csv", index=False)

    rows = []
    for m in metric_cols:
        higher_better = METRICS[m]
        d = merged[f"d_{m}"]
        win = (d > 0) if higher_better else (d < 0)
        rows.append({
            "metric":        m,
            "base_mean":     merged[f"{m}_base"].mean(),
            "var_mean":      merged[f"{m}_var"].mean(),
            "base_median":   merged[f"{m}_base"].median(),
            "var_median":    merged[f"{m}_var"].median(),
            "delta_mean":    d.mean(),
            "delta_median":  d.median(),
            "win_rate_var":  win.mean(),
            "higher_better": higher_better,
        })
    summ = pd.DataFrame(rows)
    summ.to_csv(out / f"{prefix}_summary.csv", index=False)
    print(f"  tejidos pareados: {len(merged)}")
    print(summ.to_string(index=False))

    # Figuras autoexplicativas para la seccion de validacion del TFG
    try:
        _fig_paired_scatter (merged, base_lab, var_lab, out, prefix)
        _fig_delta_by_tissue(merged, base_lab, var_lab, out, prefix)
        _fig_summary_overview(summ, base_lab, var_lab, out, prefix)
        print(f"  figuras: {prefix}_scatter.png, _delta_by_tissue.png, _summary.png")
    except Exception as e:  # pragma: no cover
        print(f"  [aviso] figuras de contraste fallaron: {e}")

    print(f"\n  ✓ Contraste en {out.resolve()}")


# ===========================================================================
# CLI
# ===========================================================================

def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Reporte de la validacion de CIRCUST sobre babuino.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("results", nargs="+", type=str,
                   help="results.csv (1 = agregado; con --compare se esperan 2).")
    p.add_argument("--compare", action="store_true",
                   help="Modo contraste pareado: pasa 2 results.csv (base var).")
    p.add_argument("--top-n", type=int, default=9,
                   help="Tejidos en la figura tipo Fig 7A. Defecto: 9.")
    p.add_argument("-o", "--output",
                   default=str(PROJECT_ROOT / "output" / "baboon" / "analysis"),
                   help="Directorio de salida. Defecto: output/baboon/analysis.")
    return p.parse_args()


def main() -> None:
    args = _parse_args()
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)

    print("╔══════════════════════════════════════════════════════════════════╗")
    print("  Analisis CIRCUST sobre babuino (Mure 2018) vs ZT")
    print("╚══════════════════════════════════════════════════════════════════╝")

    if args.compare:
        if len(args.results) != 2:
            raise SystemExit("--compare requiere exactamente 2 results.csv "
                             "(base var).")
        compare_report(Path(args.results[0]), Path(args.results[1]), out)
    else:
        if len(args.results) != 1:
            raise SystemExit("El modo agregado espera 1 results.csv "
                             "(usa --compare para 2).")
        aggregate_report(Path(args.results[0]), out, top_n=args.top_n)


if __name__ == "__main__":
    main()
