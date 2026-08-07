#!/usr/bin/env python3
"""
validate.py
===========
Compara la salida de ``run_pipeline.py`` (``sample_order.csv``) contra un
ground truth temporal (``true_times.csv``) y calcula el panel completo de
metricas definidas en :mod:`validation.metrics`.

Este script es **generico**: funciona tanto para datos sinteticos
generados con ``generate_synthetic.py`` como para datos reales con tiempo
conocido (skeletal muscle de Larriba, baboons de Mure, etc.) — basta con
que el ground truth se entregue en formato ``true_times.csv``.

Entradas requeridas
-------------------
``--predicted PATH/sample_order.csv``
    CSV con columnas ``sample_name``, ``circular_phase`` (rad) y
    ``phase_h`` (horas). Salida directa de ``run_pipeline.py``.

``--true PATH/true_times.csv``
    CSV con columnas ``sample_name`` y al menos ``true_phase_h``
    (horas, plegado a [0, 24)). Salida directa de ``generate_synthetic.py``
    o equivalente para datos reales.

Entradas opcionales (para CRE)
------------------------------
``--expr PATH/matrix.csv``
    Matriz de expresion **cruda** (la misma que se alimento a
    ``run_pipeline.py``). El script aplica el mismo preprocesado del
    pipeline (drop sparse + min-max normalize) y refittea FMM por gen
    TOP bajo ambos ordenes (predicho vs real). Necesaria para CRE.

``--atlas PATH/circust_atlas.csv``
    Atlas del pipeline. De aqui se extrae la lista de genes TOP a
    refittear. Tambien necesario para CRE.

``--fmm-method {base,stabilized}``
    Modelo FMM a usar para el re-fit. Defecto: ``base`` (mas simple y
    no requiere un dataset preparado por cada conjunto de fases).

Salidas
-------
    ``<output>/summary.csv``      Una fila con todas las metricas
                                  (incluye ``CRE`` si se proporciono
                                  ``--expr`` y ``--atlas``).
    ``<output>/per_sample.csv``   Por muestra: true_h, pred_h, error_h.
    ``<output>/manifest.json``    Trazabilidad de archivos usados.

Workflow tipico
---------------
    python validation/generate_synthetic.py --dataset SynDST4 --seed 42 \\
            -o output/synth/SynDST4
    python scripts/run_pipeline.py \\
            --data output/synth/SynDST4/matrix.csv \\
            -o output/run/SynDST4
    # solo metricas de fases
    python validation/validate.py \\
            --predicted output/run/SynDST4/sample_order.csv \\
            --true      output/synth/SynDST4/true_times.csv \\
            -o          output/validate/SynDST4
    # con CRE
    python validation/validate.py \\
            --predicted output/run/SynDST4/sample_order.csv \\
            --true      output/synth/SynDST4/true_times.csv \\
            --expr      output/synth/SynDST4/matrix.csv \\
            --atlas     output/run/SynDST4/circust_atlas.csv \\
            -o          output/validate/SynDST4
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "validation"))
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from metrics import (
    PERIOD_HOURS,
    align_phases_by_mad,
    auc_error_cdf,
    circular_abs_diff,
    circular_correlation,
    circular_correlation_trimmed,
    circular_regression_rho,
    cre,
    mean_resultant_length,
    median_absolute_error,
    pct_within,
    std_absolute_error,
)


# ===========================================================================
# Carga de entradas
# ===========================================================================

def _load_predicted(path: Path) -> pd.Series:
    """
    Lee ``sample_order.csv`` y devuelve una Serie ``phase_h`` indexada
    por ``sample_name``.
    """
    df = pd.read_csv(path)
    required = {"sample_name", "phase_h"}
    missing = required - set(df.columns)
    if missing:
        raise SystemExit(
            f"'{path}' no tiene las columnas requeridas {required}. "
            f"Faltan: {missing}. Columnas presentes: {list(df.columns)}"
        )
    return pd.Series(
        df["phase_h"].values,
        index=df["sample_name"].values,
        name="predicted_h",
    )


def _load_true(path: Path) -> pd.Series:
    """
    Lee ``true_times.csv`` y devuelve una Serie ``true_phase_h``
    (plegada a [0, 24)) indexada por ``sample_name``.

    Acepta dos esquemas:
      - Si existe ``true_phase_h`` la usa directamente.
      - Si solo existe ``true_time_h`` la pliega modulo 24.
    """
    df = pd.read_csv(path)
    if "sample_name" not in df.columns:
        raise SystemExit(
            f"'{path}' debe tener una columna 'sample_name'. "
            f"Columnas presentes: {list(df.columns)}"
        )
    if "true_phase_h" in df.columns:
        true_h = df["true_phase_h"].values
    elif "true_time_h" in df.columns:
        true_h = df["true_time_h"].values % PERIOD_HOURS
    else:
        raise SystemExit(
            f"'{path}' debe tener al menos 'true_phase_h' o 'true_time_h'. "
            f"Columnas presentes: {list(df.columns)}"
        )
    return pd.Series(true_h, index=df["sample_name"].values, name="true_h")


# ===========================================================================
# CRE: re-fit FMM bajo orden predicho y orden real
# ===========================================================================

def _refit_fmm_both_orders(
    expr_path:    Path,
    atlas_path:   Path,
    sample_names: list,
    aligned_h:    np.ndarray,
    true_h:       np.ndarray,
    fmm_method:   str = "base",
) -> tuple[np.ndarray, np.ndarray, list, int]:
    """
    Aplica el preprocesado de CIRCUST a la matriz cruda, restringe a los
    genes TOP del atlas y refittea FMM por gen dos veces: con las fases
    predichas alineadas y con las fases reales. Devuelve dos matrices
    ``(G, n)`` listas para alimentar :func:`validation.metrics.cre`.

    Notas
    -----
    - Aligned phases are used (rotacion + reflexion respecto al ground
      truth). La rotacion no afecta los valores ajustados FMM por muestra
      (el curve solo se desplaza), pero la reflexion si, asi que aqui
      pasamos las fases ya alineadas para que el orden flippeado se
      compare con su mejor interpretacion biologica.
    - Para ``stabilized``, se preparan dos datasets distintos (uno por
      conjunto de fases) ya que el grid/cache depende del conjunto.

    Returns
    -------
    fit_pred : np.ndarray, forma ``(G, n)``
    fit_real : np.ndarray, forma ``(G, n)``
    top_used : list[str]
        Nombres de los genes TOP efectivamente refittados (los que
        sobreviven al preprocesado).
    n_dropped_genes : int
        Numero de genes TOP del atlas que el preprocesado descarto.
    """
    from circust.preprocessing      import load_expression_matrix, Preprocessor
    from circust.fitting.fmm        import FMMModel
    from circust.fitting.fmm_stabilized import (
        StabilizedFMMModel, prepare_dataset,
    )

    # 1) Cargar matriz cruda y aplicar el MISMO preprocesado del pipeline
    expr_raw = load_expression_matrix(str(expr_path))
    prep     = Preprocessor(sparse_threshold=0.3, verbose=False).run(expr_raw)
    expr_norm = prep.expr_norm

    # 2) Filtrar a los genes TOP del atlas que sobrevivieron
    atlas       = pd.read_csv(atlas_path)
    top_names   = atlas["gene_name"].tolist()
    top_present = [g for g in top_names if g in expr_norm.index]
    n_dropped   = len(top_names) - len(top_present)

    # 3) Restringir a los samples del cruce
    cols_present = [s for s in sample_names if s in expr_norm.columns]
    if len(cols_present) != len(sample_names):
        missing = set(sample_names) - set(cols_present)
        raise SystemExit(
            f"{len(missing)} muestras del cruce no estan en la matriz "
            f"de expresion (probablemente filtradas). Ejemplos: "
            f"{list(missing)[:3]}"
        )
    top_expr = expr_norm.loc[top_present, cols_present].values.astype(np.float64)

    # 4) Convertir fases a radianes (lo que esperan los FMM)
    two_pi_over_p = 2.0 * np.pi / PERIOD_HOURS
    pred_rad = (aligned_h % PERIOD_HOURS) * two_pi_over_p
    true_rad = (true_h    % PERIOD_HOURS) * two_pi_over_p

    # 5) Construir los modelos. Para stabilized hay dos datasets distintos.
    if fmm_method == "stabilized":
        prep_pred = prepare_dataset(pred_rad)
        prep_real = prepare_dataset(true_rad)
        fmm_pred  = StabilizedFMMModel(prepared=prep_pred)
        fmm_real  = StabilizedFMMModel(prepared=prep_real)
    else:
        fmm_pred = FMMModel(
            length_alpha_grid = 48,
            length_omega_grid = 24,
            num_reps          = 3,
        )
        fmm_real = fmm_pred  # mismo modelo, distintas fases

    # 6) Re-fit por gen
    G = len(top_present)
    n = len(cols_present)
    fit_pred = np.zeros((G, n), dtype=np.float64)
    fit_real = np.zeros((G, n), dtype=np.float64)
    for k in range(G):
        y = top_expr[k]
        fit_pred[k] = fmm_pred.fit(y, pred_rad).fitted
        fit_real[k] = fmm_real.fit(y, true_rad).fitted

    return fit_pred, fit_real, top_present, n_dropped


# ===========================================================================
# Diagnostico
# ===========================================================================

def _ccc_reliability(r_bar: float) -> str:
    """Etiqueta cualitativa de fiabilidad del CCC segun R_bar de los true."""
    if   r_bar >= 0.30: return "fiable"
    elif r_bar >= 0.10: return "frontera"
    else:               return "inestable"


# ===========================================================================
# CLI
# ===========================================================================

def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Compara sample_order.csv (output de run_pipeline.py) contra "
            "true_times.csv y reporta MedAE, SDAE, AUC, CCC, rTrim, %within."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--predicted", required=True, type=str,
                   help="Ruta a sample_order.csv (output del pipeline).")
    p.add_argument("--true",      required=True, type=str,
                   help="Ruta a true_times.csv (ground truth).")
    p.add_argument("--expr",      default=None, type=str,
                   help="Matriz de expresion cruda (matrix.csv). Necesaria "
                        "para CRE.")
    p.add_argument("--atlas",     default=None, type=str,
                   help="circust_atlas.csv. Necesario para CRE (lista de "
                        "genes TOP).")
    p.add_argument("--fmm-method", default="base",
                   choices=["base", "stabilized"],
                   help="Modelo FMM a usar para el re-fit de CRE. "
                        "Defecto: base.")
    p.add_argument("--no-reflection", action="store_true",
                   help="Desactiva la prueba de reflexion en la alineacion. "
                        "Util si se quiere preservar el sentido del circulo.")
    p.add_argument("--trim-fraction", type=float, default=0.10,
                   help="Fraccion de poda para rTrim (Mahmood 2022). Defecto: 0.10.")
    p.add_argument("--label",     default=None, type=str,
                   help="Etiqueta del dataset en summary.csv. Si se omite, "
                        "se infiere del nombre del directorio padre.")
    p.add_argument("-o", "--output", required=True, type=str,
                   help="Directorio de salida (se crea si no existe).")
    return p.parse_args()


# ===========================================================================
# Main
# ===========================================================================

def main() -> None:
    args = _parse_args()
    out  = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)

    pred_path = Path(args.predicted).resolve()
    true_path = Path(args.true).resolve()
    label = args.label or pred_path.parent.name

    print(f"╔══════════════════════════════════════════════════════════════════╗")
    print(f"  Validacion : {label}")
    print(f"  predicted  : {pred_path}")
    print(f"  true       : {true_path}")
    print(f"  output     : {out.resolve()}")
    print(f"╚══════════════════════════════════════════════════════════════════╝")

    t_start = time.time()

    # 1) Cargar y cruzar por sample_name
    pred = _load_predicted(pred_path)
    true = _load_true(true_path)

    common = pred.index.intersection(true.index)
    n_common = len(common)
    if n_common == 0:
        raise SystemExit(
            "Los archivos predicted y true no comparten ningun sample_name.\n"
            f"  predicted: {list(pred.index[:3])}... ({len(pred)} en total)\n"
            f"  true     : {list(true.index[:3])}... ({len(true)} en total)"
        )

    n_only_pred = len(pred.index.difference(true.index))
    n_only_true = len(true.index.difference(pred.index))
    if n_only_pred or n_only_true:
        print(f"  Aviso: {n_only_pred} muestras solo en predicted, "
              f"{n_only_true} solo en true. Se ignoran.")

    pred_arr = pred.loc[common].values.astype(np.float64)
    true_arr = true.loc[common].values.astype(np.float64) % PERIOD_HOURS

    # 2) Alinear (rotacion + reflexion opcional)
    alig = align_phases_by_mad(
        pred_arr, true_arr,
        period         = PERIOD_HOURS,
        try_reflection = not args.no_reflection,
    )

    # 3) Metricas sobre fases
    aligned = alig.aligned_pred

    medae      = median_absolute_error      (aligned, true_arr, PERIOD_HOURS)
    sdae       = std_absolute_error         (aligned, true_arr, PERIOD_HOURS)
    auc        = auc_error_cdf              (aligned, true_arr, PERIOD_HOURS)
    ccc        = circular_correlation       (aligned, true_arr, PERIOD_HOURS)
    rtrim      = circular_correlation_trimmed(
        aligned, true_arr, PERIOD_HOURS, trim_fraction=args.trim_fraction,
    )
    # rho de regresion circular-circular (metrica del paper para babuino):
    # robusta al muestreo equiespaciado. Se ajusta sobre la prediccion cruda
    # porque la regresion incorpora su propia rotacion + reflexion.
    rho_cc     = circular_regression_rho    (pred_arr, true_arr, PERIOD_HOURS)
    pct1       = pct_within                 (aligned, true_arr, 1.0, PERIOD_HOURS)
    pct2       = pct_within                 (aligned, true_arr, 2.0, PERIOD_HOURS)
    r_bar_true = mean_resultant_length      (true_arr, PERIOD_HOURS)
    r_bar_pred = mean_resultant_length      (aligned,  PERIOD_HOURS)

    # 4) per_sample
    per_sample = pd.DataFrame({
        "sample_name":    common,
        "true_h":         true_arr,
        "predicted_h":    pred_arr,
        "aligned_pred_h": aligned,
        "abs_error_h":    circular_abs_diff(aligned, true_arr, PERIOD_HOURS),
    })
    per_sample.to_csv(out / "per_sample.csv", index=False)

    # 5) CRE (opcional, requiere --expr y --atlas)
    cre_value:   float | None = None
    n_top_used:  int   | None = None
    n_top_lost:  int   | None = None
    if args.expr and args.atlas:
        print()
        print(f"  ── CRE: re-fit FMM ({args.fmm_method}) ─────────────────")
        t_cre = time.time()
        fit_pred, fit_real, top_used, n_top_lost = _refit_fmm_both_orders(
            expr_path    = Path(args.expr).resolve(),
            atlas_path   = Path(args.atlas).resolve(),
            sample_names = list(common),
            aligned_h    = aligned,
            true_h       = true_arr,
            fmm_method   = args.fmm_method,
        )
        cre_value  = cre(fit_real, fit_pred)
        n_top_used = len(top_used)
        cre_secs   = time.time() - t_cre
        print(f"    Genes TOP refittados : {n_top_used}  "
              f"(+ {n_top_lost} descartados por preprocesado)")
        print(f"    Tiempo               : {cre_secs:.1f}s")
        print(f"    CRE                  : {cre_value:+.4f}")
    elif args.expr or args.atlas:
        print()
        print(f"  Aviso: --expr y --atlas se necesitan AMBOS para CRE; "
              f"se ignora la entrada parcial.")

    # 6) summary
    summary_row = {
        "label":             label,
        "n_samples":         n_common,
        "MedAE":             medae,
        "SDAE":              sdae,
        "AUC":               auc,
        "CCC":               ccc,
        "rTrim":             rtrim,
        "rho_cc":            rho_cc,
        "CRE":               cre_value,
        "n_top_genes_CRE":   n_top_used,
        "pct_within_1h":     pct1,
        "pct_within_2h":     pct2,
        "R_bar_true":        r_bar_true,
        "R_bar_pred":        r_bar_pred,
        "CCC_reliability":   _ccc_reliability(r_bar_true),
        "alignment_delta_h": alig.delta,
        "alignment_sign":    alig.sign,
        "alignment_mad_h":   alig.mad,
    }
    pd.DataFrame([summary_row]).to_csv(out / "summary.csv", index=False)

    # 7) manifest
    manifest = {
        "label":          label,
        "predicted_file": str(pred_path),
        "true_file":      str(true_path),
        "expr_file":      str(Path(args.expr).resolve())  if args.expr  else None,
        "atlas_file":     str(Path(args.atlas).resolve()) if args.atlas else None,
        "fmm_method":     args.fmm_method if (args.expr and args.atlas) else None,
        "n_samples":      n_common,
        "try_reflection": not args.no_reflection,
        "trim_fraction":  args.trim_fraction,
        "period_h":       PERIOD_HOURS,
    }
    with open(out / "manifest.json", "w") as f:
        json.dump(manifest, f, indent=2, default=str)

    # 7) Reporte
    elapsed = time.time() - t_start
    print()
    print(f"  ── Metricas (fases) ──────────────────────────────────")
    print(f"    n samples       = {n_common}")
    print(f"    MedAE           = {medae:7.3f} h")
    print(f"    SDAE            = {sdae:7.3f} h")
    print(f"    AUC             = {auc:7.3f}")
    print(f"    CCC             = {ccc:+7.3f}  "
          f"[R_bar(true)={r_bar_true:.3f} → {_ccc_reliability(r_bar_true)}]")
    print(f"    rTrim (δ={args.trim_fraction:.2f}) = {rtrim:+7.3f}")
    print(f"    rho (reg c-c)   = {rho_cc:+7.3f}  (robusta a ZT equiespaciado)")
    print(f"    % within 1h     = {pct1:6.1f}%")
    print(f"    % within 2h     = {pct2:6.1f}%")
    print(f"  ── Alineacion ────────────────────────────────────────")
    print(f"    delta           = {alig.delta:+.3f} h")
    print(f"    sign            = {alig.sign:+d}")
    print(f"    MAD (post-align)= {alig.mad:.3f} h")

    if cre_value is not None:
        print(f"  ── CRE ───────────────────────────────────────────────")
        print(f"    CRE             = {cre_value:+.4f}  "
              f"(genes={n_top_used})")

    print()
    print(f"  ✓ Validacion en {elapsed:.1f}s")
    print(f"    Archivos en {out.resolve()}:")
    print(f"      - summary.csv")
    print(f"      - per_sample.csv")
    print(f"      - manifest.json")


if __name__ == "__main__":
    main()
