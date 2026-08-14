#!/usr/bin/env python3
"""
simulation/cells.py
===================
Volcado de celdas sinteticas a disco, en la forma de unidad que espera la capa
de evaluacion de ``validation/``.

Por que hace falta
------------------
``generate()`` devuelve un objeto en memoria, pero ``comparar_unidad`` recibe
la ruta de un fichero y ``descubrir`` recorre directorios. Una celda sintetica
escrita con este modulo es indistinguible de un tejido real para esa capa, que
es lo que permite reutilizarla sin tocarla y que las tablas del capitulo de
simulaciones y las de los datos reales salgan del mismo codigo.

Forma de una celda
------------------
    <destino>/
        matriz.parquet   genes x muestras, la entrada del pipeline
        tiempos.csv      sample_name, true_time_h
        params.csv       parametros FMM reales de los genes ritmicos
        meta.json        todo lo necesario para reconstruir la celda

``tiempos.csv`` y ``params.csv`` son dos verdades distintas y complementarias.
El primero permite medir el error del **orden de las muestras**, el segundo el
error de la **caracterizacion de cada gen**. Esta segunda no puede medirse
sobre datos reales, porque alli no se conocen los parametros verdaderos, y es
una de las carencias del estudio de simulacion original de CIRCUST.

Formato y ciclo de vida
-----------------------
Por defecto la matriz se escribe en Parquet, que sobre la celda mayor de la
rejilla (56.200 genes por 480 muestras) ocupa 270 MB y tarda 1 segundo, frente
a los 547 MB y 20 segundos del CSV. Se escribe ademas el CSV solo en las celdas
que tiene que leer el CIRCUST original en R, para no obligar a la instalacion
de R a depender del paquete ``arrow``.

Una celda es una funcion pura de sus parametros y su semilla, de modo que la
matriz es material regenerable y no hace falta conservarla. El ciclo normal de
la rejilla es escribir, evaluar y borrar la matriz con :func:`matrix_on_disk`,
conservando siempre el ground truth y el manifiesto, que pesan poco. Si mas
adelante hace falta volver sobre una celda concreta, :func:`regenerate` la
reconstruye desde su manifiesto.
"""
from __future__ import annotations

import hashlib
import json
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator

import numpy as np
import pandas as pd

from simulation.generate import SyntheticDataset, generate

# Nombre de la columna que lleva los identificadores de gen. Es el mismo para
# CSV y para Parquet, de modo que la regla de carga es unica:
# ``load_expression_matrix(ruta, gene_column=GENE_COLUMN)``. Parquet, ademas,
# exige indicarla de forma explicita.
GENE_COLUMN: str = "gene"

FICHERO_TIEMPOS: str = "tiempos.csv"
FICHERO_PARAMS:  str = "params.csv"
FICHERO_META:    str = "meta.json"

FORMATOS: tuple[str, ...] = ("parquet", "csv")


@dataclass
class Truth:
    """El ground truth de una celda, leido de disco.

    Attributes
    ----------
    tiempos
        Serie indexada por nombre de muestra con el instante real en horas. Es
        la forma que ``comparar_unidad`` espera para los tiempos.
    params
        Una fila por gen ritmico, con sus parametros FMM verdaderos. Los genes
        de fondo no aparecen: son deducibles, porque no son ritmicos y su
        numero esta en ``meta["n_genes"] - meta["n_rhythmic"]``.
    meta
        Manifiesto de generacion.
    """
    tiempos: pd.Series
    params:  pd.DataFrame
    meta:    dict


def _huella(matriz: pd.DataFrame) -> str:
    """SHA-1 de los valores de la matriz.

    Sirve para verificar que una celda regenerada es la misma que se evaluo,
    que es lo que sostiene la afirmacion de reproducibilidad de la memoria.
    """
    valores = np.ascontiguousarray(matriz.to_numpy(dtype=np.float64))
    return hashlib.sha1(valores.tobytes()).hexdigest()


def _escribir_matriz(matriz: pd.DataFrame, destino: Path, formato: str) -> Path:
    if formato not in FORMATOS:
        raise ValueError(
            f"formato {formato!r} desconocido, los disponibles son {FORMATOS}.")
    tabla = matriz.rename_axis(GENE_COLUMN).reset_index()
    ruta = destino / f"matriz.{formato}"
    if formato == "parquet":
        tabla.to_parquet(ruta, index=False)
    else:
        tabla.to_csv(ruta, index=False)
    return ruta


def write_cell(
    ds:          SyntheticDataset,
    destino:     Path | str,
    *,
    formato:     str = "parquet",
    tambien_csv: bool = False,
) -> Path:
    """Escribe una celda completa y devuelve la ruta de la matriz.

    Parameters
    ----------
    ds
        La celda generada.
    destino
        Directorio de la celda. Se crea junto con sus padres si no existe.
    formato
        ``"parquet"`` o ``"csv"``.
    tambien_csv
        Escribe ademas una copia en CSV. Se usa en las celdas que tiene que
        leer el CIRCUST original en R.
    """
    destino = Path(destino)
    destino.mkdir(parents=True, exist_ok=True)

    ruta = _escribir_matriz(ds.matrix, destino, formato)
    if tambien_csv and formato != "csv":
        _escribir_matriz(ds.matrix, destino, "csv")

    pd.DataFrame({
        "sample_name": list(ds.matrix.columns),
        "true_time_h": np.asarray(ds.true_times_h, dtype=float),
    }).to_csv(destino / FICHERO_TIEMPOS, index=False)

    # Solo los ritmicos. Los de fondo no tienen parametros que registrar y
    # multiplicarian por veinte el tamano del ground truth de la rejilla.
    ritmicos = ds.gene_params[ds.gene_params["ritmico"]]
    ritmicos.to_csv(destino / FICHERO_PARAMS, index=False)

    meta = dict(ds.meta)
    meta["matriz_sha1"] = _huella(ds.matrix)
    meta["formato"] = formato
    (destino / FICHERO_META).write_text(
        json.dumps(meta, indent=2, ensure_ascii=False, sort_keys=True))

    return ruta


def matrix_path(destino: Path | str) -> Path | None:
    """Ruta de la matriz de una celda, o ``None`` si ya se borro."""
    destino = Path(destino)
    for formato in FORMATOS:
        ruta = destino / f"matriz.{formato}"
        if ruta.exists():
            return ruta
    return None


def drop_matrix(destino: Path | str) -> None:
    """Borra las matrices de una celda y conserva el ground truth.

    Es idempotente a proposito, porque la rejilla puede reanudarse sobre
    celdas ya limpiadas.
    """
    destino = Path(destino)
    for formato in FORMATOS:
        (destino / f"matriz.{formato}").unlink(missing_ok=True)


@contextmanager
def matrix_on_disk(
    ds:        SyntheticDataset,
    destino:   Path | str,
    *,
    conservar: bool = False,
    **kwargs,
) -> Iterator[Path]:
    """Escribe la celda, cede la ruta de la matriz y la borra al salir.

    Es el ciclo normal de la rejilla. El ground truth y el manifiesto se
    conservan siempre. El borrado ocurre tambien si la evaluacion falla, de
    modo que una rejilla larga no se queda sin disco por el camino.

    Con ``conservar=True`` la matriz se mantiene, que es lo que se necesita en
    las celdas sobre las que despues corre el CIRCUST original en R.
    """
    ruta = write_cell(ds, destino, **kwargs)
    try:
        yield ruta
    finally:
        if not conservar:
            drop_matrix(destino)


def read_truth(destino: Path | str) -> Truth:
    """Lee el ground truth de una celda, este o no la matriz."""
    destino = Path(destino)
    tiempos = pd.read_csv(destino / FICHERO_TIEMPOS)
    return Truth(
        tiempos=tiempos.set_index("sample_name")["true_time_h"],
        params=pd.read_csv(destino / FICHERO_PARAMS),
        meta=json.loads((destino / FICHERO_META).read_text()),
    )


def regenerate(destino: Path | str) -> SyntheticDataset:
    """Reconstruye una celda desde su manifiesto.

    Devuelve exactamente la misma celda que se escribio, incluida la matriz,
    porque ``generate`` es determinista dada la semilla. Permite volver sobre
    una celda cuya matriz ya se borro.
    """
    meta = json.loads((Path(destino) / FICHERO_META).read_text())
    # Las claves que no son argumentos de generate se descartan: son
    # descriptivas (``ruido``) o derivadas (``n_rhythmic``, ``matriz_sha1``).
    firma = ("n_samples", "n_genes", "prop_rhythmic", "sampling",
             "waveform_mix", "omega", "omega_jitter", "beta_mode", "snr",
             "core_phase_sd_h", "seed", "period")
    kwargs = {k: meta[k] for k in firma if k in meta}
    kwargs.update(meta.get("sampling_kwargs") or {})
    return generate(meta.get("name", "celda"), **kwargs)
