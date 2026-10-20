"""
circust/outputs.py
===========================
Escritura de los resultados de CIRCUST.

CIRCUST no entrega un unico orden de muestras sino K ordenes sincronizados,
asi que las salidas por orden se escriben K veces, una por cada k.

Ficheros (``k`` = 1..K)
-----------------------
sample_order_k<k>.csv     muestra, posicion, fase del orden k (radianes en [0, 2pi))
ordered_matrix_k<k>.csv   matriz COMPLETA (todos los genes tras el preprocesado)
                          reordenada segun el orden k. Con ``values="original"``
                          lleva los valores de entrada sin normalizar.
fitted_top_k<k>.csv       valores ajustados del FMM de los genes TOP segun el orden k
parameters_by_k.csv       K filas por gen TOP: M, A, alpha, beta, omega, t_U, t_L, R2
circust_atlas.csv         atlas [M^TOP]: mediana de los K parametros por gen TOP
outliers.csv              muestras eliminadas por la deteccion de outliers
top_genes.csv             genes TOP, con su intervalo y sus valores en el orden preliminar
manifest.json             configuracion, semilla y tiempos
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from circust.cpca import CPCAResult
from circust.robust_order import RobustOrderResult
from circust.top_genes import TopGeneResult


def _write_matrix(df: pd.DataFrame, path: Path, fmt: str, index_label: str) -> None:
    if fmt == "parquet":
        out = df.copy()
        out.index.name = index_label
        out.reset_index().to_parquet(path.with_suffix(".parquet"), index=False)
    else:
        df.to_csv(path.with_suffix(".csv"), index_label=index_label)


def write_results(
    out_dir:        Path,
    robust:         RobustOrderResult,
    full_matrix:    pd.DataFrame,
    cpca:           CPCAResult,
    top:            TopGeneResult,
    manifest:       dict,
    values:         str = "normalized",
    matrix_format:  str = "csv",
    raw_matrix:     pd.DataFrame | None = None,
) -> list[str]:
    """Escribe todas las salidas y devuelve la lista de ficheros creados.

    Parameters
    ----------
    full_matrix : matriz completa normalizada en el orden preliminar, con las
        mismas columnas que la matriz TOP.
    values : ``"normalized"`` o ``"original"``. Que valores lleva
        ``ordered_matrix_k<k>``. Con ``"original"`` hay que pasar ``raw_matrix``,
        la matriz de entrada sin normalizar.
    matrix_format : ``"csv"`` o ``"parquet"`` para las matrices completas.
    """
    if values not in ("normalized", "original"):
        raise ValueError("values debe ser 'normalized' u 'original'.")
    if matrix_format not in ("csv", "parquet"):
        raise ValueError("matrix_format debe ser 'csv' o 'parquet'.")
    if values == "original" and raw_matrix is None:
        raise ValueError("values='original' necesita raw_matrix.")

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    written: list[str] = []

    def done(name: str) -> None:
        written.append(name)

    source = full_matrix if values == "normalized" else raw_matrix
    cols = list(full_matrix.columns)
    for k in range(robust.n_reps):
        names = robust.ordered_names(k)
        pd.DataFrame({
            "sample_name":    names,
            "position":       np.arange(1, len(names) + 1),
            "circular_phase": robust.circular_scales[k],
        }).to_csv(out_dir / f"sample_order_k{k + 1}.csv", index=False)
        done(f"sample_order_k{k + 1}.csv")

        ordered = source.loc[:, names] if values == "original" else \
            full_matrix.iloc[:, robust.sample_orders[k]]
        ordered.columns = names
        _write_matrix(ordered, out_dir / f"ordered_matrix_k{k + 1}", matrix_format, "gene_id")
        done(f"ordered_matrix_k{k + 1}.{matrix_format}")

        robust.fitted[k].to_csv(out_dir / f"fitted_top_k{k + 1}.csv", index_label="gene_id")
        done(f"fitted_top_k{k + 1}.csv")

    robust.parameters.to_csv(out_dir / "parameters_by_k.csv", index=False)
    done("parameters_by_k.csv")
    robust.atlas.sort_values("r2_fmm", ascending=False).to_csv(
        out_dir / "circust_atlas.csv", index=False)
    done("circust_atlas.csv")

    pd.DataFrame({"sample_name": cpca.samples_dropped_names,
                  "reason": ["|Z|>4" if i in cpca.outliers_uni else "|Z|>3 + radial"
                             for i in cpca.samples_dropped]}
                 ).to_csv(out_dir / "outliers.csv", index=False)
    done("outliers.csv")

    pd.DataFrame({
        "gene_name": top.gene_names, "interval": top.interval,
        "forced_core": [g in set(top.forced_genes) for g in top.gene_names],
        "peak_prelim_rad": top.fmm_peaks, "omega_prelim": top.omega,
        "r2_fmm_prelim": top.r2_fmm, "r2_ori": top.r2_ori,
    }).to_csv(out_dir / "top_genes.csv", index=False)
    done("top_genes.csv")

    (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2, default=str),
                                           encoding="utf-8")
    done("manifest.json")
    return written
