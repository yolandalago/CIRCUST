"""
Recalcula las medidas de validacion a partir de las salidas guardadas.

``run.py`` guarda, junto al CSV, lo que el pipeline produce en cada grupo y
metodo (``simulation/salidas.py``). Este script lee esas salidas, regenera la
verdad y la matriz de cada grupo a partir de la rejilla, que es determinista,
y vuelve a calcular las medidas sin ejecutar el pipeline. Sirve para corregir
una medida o anadir una nueva: basta con tocar :func:`medidas_fila` y volver
a correrlo.

Uso
---
    python -m simulation.reevaluar --bloque rejilla --reps 10 \\
        --out output/simulacion_2026-10

Escribe ``medidas.csv`` en ``--out``, una fila por ``row_id``. Las columnas
de factores y de configuracion siguen en ``resultados.csv`` y se cruzan por
``row_id``.
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
for _extra in (PROJECT_ROOT, PROJECT_ROOT / "src", PROJECT_ROOT / "scripts"):
    if str(_extra) not in sys.path:
        sys.path.insert(0, str(_extra))

from simulation.atlas import atlas_metrics
from simulation.generate import generate_dataset
from simulation.grid import BLOQUES, build_grid
from simulation.metrics import ajuste_fmm, cre, evaluate, genes_ritmicos
from simulation.salidas import grupo_id, leer, matriz_evaluacion


def medidas_fila(truth, theta_hat: pd.Series, atlas: pd.DataFrame | None,
                 Y: np.ndarray, referencia: np.ndarray, muestras_T: list[str],
                 modelo: dict | None = None) -> dict:
    """Todas las medidas de una fila. Es el sitio donde anadir una nueva.

    ``modelo`` = ``{"fmm_method", "lambda_star"}`` del grupo, el ajustador
    del CRE (ver ``metrics.ajuste_fmm``).
    """
    modelo = modelo or {}
    salida = evaluate(truth.theta_for(list(theta_hat.index)), theta_hat.to_numpy(float))
    if atlas is not None and not atlas.empty:
        salida.update(atlas_metrics(atlas, truth.gene_params, quitar_giro=True))
    if Y.size:
        salida["cre"] = cre(referencia, ajuste_fmm(
            Y, theta_hat.loc[muestras_T].to_numpy(float), **modelo))
    return salida


def reevaluar_grupo(cfg, fases: pd.DataFrame, atlas: pd.DataFrame,
                    muestras: list[str], fmm_method: str = "base",
                    lambda_star: float | None = None) -> list[dict]:
    """Medidas de todas las filas de un grupo, con su modelo FMM y su lambda."""
    modelo = dict(fmm_method=fmm_method,
                  lambda_star=None if lambda_star is None or np.isnan(lambda_star)
                  else float(lambda_star))
    matriz, truth = generate_dataset(cfg)
    expr = matriz_evaluacion(matriz, muestras)
    T = genes_ritmicos(truth.gene_params, expr.index)
    Y = expr.loc[T, muestras].to_numpy()
    referencia = (ajuste_fmm(Y, truth.theta_for(muestras), **modelo)
                  if T else np.empty((0, 0)))

    filas = []
    for row_id, f in fases.groupby("row_id", sort=False):
        theta_hat = pd.Series(f["theta_hat"].to_numpy(float), index=f["sample"].tolist())
        a = atlas[atlas["row_id"] == row_id].drop(columns=["row_id", "grupo_id"])
        filas.append({"row_id": row_id,
                      **medidas_fila(truth, theta_hat, a, Y, referencia, muestras, modelo)})
    return filas


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument("--bloque", required=True, choices=list(BLOQUES))
    parser.add_argument("--reps", type=int, required=True,
                        help="Las replicas con las que se ejecuto el bloque.")
    parser.add_argument("--out", type=str, default="output/simulacion_2026-10")
    parser.add_argument("--jobs", type=int, default=-1)
    args = parser.parse_args()
    raiz = Path(args.out)

    # La configuracion de cada grupo sale de la rejilla, no del CSV: la
    # semilla y los parametros por defecto del generador no estan todos ahi.
    celda_de = {grupo_id(c.bloque, c.cfg, c.n_core_declared, c.fmm_method, c.rep): c
                for c in build_grid(args.bloque, args.reps)}

    fases = leer(raiz, "fases")
    atlas = leer(raiz, "atlas")
    prelim = leer(raiz, "preliminar")
    # la lambda que calibro el pipeline en cada grupo, para el CRE con CircaFMM
    lambda_de = leer(raiz, "grupos").set_index("grupo_id")["lambda_star"].to_dict()
    gids = [g for g in fases["grupo_id"].unique() if g in celda_de]
    print(f"{len(gids)} grupos del bloque {args.bloque} con salidas guardadas")

    t0 = time.perf_counter()
    trozos = Parallel(n_jobs=args.jobs, verbose=5)(
        delayed(reevaluar_grupo)(
            celda_de[g].cfg,
            fases[fases["grupo_id"] == g],
            atlas[atlas["grupo_id"] == g],
            prelim.loc[prelim["grupo_id"] == g, "sample"].tolist(),
            celda_de[g].fmm_method,
            lambda_de.get(g),
        )
        for g in gids
    )
    df = pd.DataFrame([f for t in trozos for f in t])
    df.to_csv(raiz / "medidas.csv", index=False)
    print(f"{len(df)} filas en {raiz / 'medidas.csv'} ({time.perf_counter() - t0:.0f} s)")


if __name__ == "__main__":
    main()
