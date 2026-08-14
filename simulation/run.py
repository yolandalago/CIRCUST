#!/usr/bin/env python3
"""
simulation/run.py
=================
Lanzador de la rejilla del estudio de simulacion.

Recorre las celdas que se le pidan, genera cada replica, la evalua con el
pipeline y anota los resultados, incluidas las metricas contra los parametros
verdaderos con que se genero.

Alcance
-------
Cubre los bloques que miden precision: ``fase0``, ``tau``, ``R`` y ``A``. El
bloque ``P`` no pasa por aqui porque necesita cronometraje por etapas, y eso lo
da ``scripts/run_pipeline.py`` a traves de su manifiesto. Ese bloque se lanza
junto con la linea base en R, que es la otra mitad de la misma comparacion.

Como sobrevive a una ejecucion larga
------------------------------------
La rejilla completa son varias horas, de modo que el lanzador esta pensado para
que un corte no obligue a repetir nada:

1. Los resultados se anaden al CSV del bloque **replica a replica**, no al
   final, luego un corte pierde como mucho una replica.
2. Al relanzar se leen los pares celda-replica ya presentes y se saltan.
3. Una celda que falla queda anotada con su error y la ejecucion sigue. En una
   rejilla de horas, abortar por una celda rota seria peor que registrarla.

La matriz de cada celda se borra tras evaluarla, porque es material regenerable
desde su manifiesto. El ground truth se conserva siempre. La excepcion es el
bloque P, cuyas celdas se conservan para que R lea exactamente los mismos
datos.

Uso
---
    python -m simulation.run --bloques fase0
    python -m simulation.run --bloques R A
    python -m simulation.run --bloques R --replicas 1     # prueba rapida
"""
from __future__ import annotations

import argparse
import json
import time
import traceback
from pathlib import Path

import pandas as pd

from simulation.cells import GENE_COLUMN, matrix_on_disk, read_truth
from simulation.grid import (
    BLOQUES, RAIZ_POR_DEFECTO, REFERENCIA, Celda, materializar, ruta_de,
)
from simulation.truth_metrics import evaluar_contra_verdad
from validation.compare_methods import comparar_unidad

# Bloques que este lanzador sabe correr. El P queda fuera a proposito, ver el
# docstring del modulo.
BLOQUES_DE_METRICAS: tuple[str, ...] = ("fase0", "tau", "R", "A")

# Parametros de generacion que se copian a cada fila. Son los que varian en
# algun bloque, de modo que el analisis puede agrupar por ellos sin volver a
# leer los manifiestos.
FACTORES: tuple[str, ...] = (
    "n_genes", "n_samples", "snr", "prop_rhythmic", "sampling", "waveform_mix",
)

FICHERO_REFERENCIA: str = "referencia.json"


def columnas_de_factores() -> tuple[str, ...]:
    return FACTORES


# ===========================================================================
# Rutas y reanudacion
# ===========================================================================

def ruta_resultados(bloque: str, raiz: Path | str = RAIZ_POR_DEFECTO) -> Path:
    return Path(raiz) / "resultados" / f"{bloque}.csv"


def ya_hechas(ruta: Path) -> set[tuple[str, int]]:
    """Pares ``(celda, replica)`` que ya figuran en el fichero de resultados."""
    if not ruta.exists():
        return set()
    hecho = pd.read_csv(ruta, usecols=["celda", "replica"])
    return set(zip(hecho["celda"], hecho["replica"].astype(int)))


def _anadir(ruta: Path, filas: list[dict]) -> None:
    """Anade filas al CSV del bloque respetando la cabecera que ya tenga.

    Escribir en modo anadir sin mirar la cabecera desalinea las columnas en
    cuanto una fila trae un juego distinto, que es exactamente lo que pasa con
    las filas de error, mucho mas cortas que las de exito. El resultado seria
    un fichero corrupto en silencio, y en una rejilla de horas nadie se daria
    cuenta hasta el analisis.

    Por eso las filas nuevas se reordenan contra la cabecera existente, y si
    aparece alguna columna que no estaba, el fichero se reescribe con la union
    en vez de forzarla a encajar.
    """
    if not filas:
        return
    ruta.parent.mkdir(parents=True, exist_ok=True)
    nuevo = pd.DataFrame(filas)

    if not ruta.exists():
        nuevo.to_csv(ruta, index=False)
        return

    cabecera = list(pd.read_csv(ruta, nrows=0).columns)
    desconocidas = [c for c in nuevo.columns if c not in cabecera]
    if desconocidas:
        viejo = pd.read_csv(ruta)
        (pd.concat([viejo, nuevo], ignore_index=True)
           .reindex(columns=cabecera + desconocidas)
           .to_csv(ruta, index=False))
        return

    nuevo.reindex(columns=cabecera).to_csv(ruta, mode="a", index=False,
                                           header=False)


# ===========================================================================
# Celda de referencia calibrada
# ===========================================================================

def cargar_referencia(raiz: Path | str = RAIZ_POR_DEFECTO) -> dict:
    """La referencia que fijo la fase 0, o la provisional si aun no existe.

    Es la puerta entre calibracion y evaluacion. Mientras la fase 0 no haya
    corrido y ``analyze`` no haya escrito su veredicto, los bloques R y A usan
    los valores provisionales de :data:`simulation.grid.REFERENCIA`, que sirven
    para probar la maquinaria pero no para producir las cifras de la memoria.
    """
    ruta = Path(raiz) / FICHERO_REFERENCIA
    if not ruta.exists():
        return dict(REFERENCIA)
    return {**REFERENCIA, **json.loads(ruta.read_text())}


def celdas_de(bloques: list[str], referencia: dict | None = None
              ) -> list[Celda]:
    """Las celdas de los bloques pedidos, en orden."""
    ref = referencia or REFERENCIA
    celdas: list[Celda] = []
    for nombre in bloques:
        if nombre not in BLOQUES_DE_METRICAS:
            raise ValueError(
                f"bloque {nombre!r} desconocido o no lanzable desde aqui, los "
                f"disponibles son {BLOQUES_DE_METRICAS}.")
        constructor = BLOQUES[nombre]
        celdas += constructor() if nombre == "fase0" else constructor(ref)
    return celdas


# ===========================================================================
# Evaluacion de una replica
# ===========================================================================

def _contra_la_verdad(params: pd.DataFrame):
    """Cierre que enchufa las metricas del ground truth al evaluador.

    Existe para que ``validation/compare_methods`` no tenga que importar nada
    de ``simulation``, lo que invertiria la direccion de dependencia entre la
    capa de evaluacion compartida y el estudio de simulacion.
    """
    def _extra(rob, top, atlas, fmm_method) -> dict:
        return evaluar_contra_verdad(atlas, params,
                                     genes_top=list(rob.top_gene_names))
    return _extra


def evaluar_replica(celda: Celda, replica: int,
                    raiz: Path | str = RAIZ_POR_DEFECTO) -> list[dict]:
    """Genera, evalua y limpia una replica. Devuelve sus filas."""
    destino = ruta_de(celda, replica, raiz)
    ds = materializar(celda, replica)

    with matrix_on_disk(ds, destino, conservar=celda.conservar_matriz,
                        tambien_csv=celda.tambien_csv) as ruta_matriz:
        verdad = read_truth(destino)
        filas = comparar_unidad(
            celda.nombre, ruta_matriz, verdad.tiempos.to_dict(),
            target_sizes=[100], fmm_methods=list(celda.fmm_methods),
            metodos=list(celda.metodos), gene_column=GENE_COLUMN,
            robust_kwargs=dict(celda.robust_kwargs),
            con_referencia=celda.con_referencia,
            extra=_contra_la_verdad(verdad.params),
        )

    cabecera = {"bloque": celda.bloque, "celda": celda.nombre,
                "replica": replica, "seed": celda.seed_base + replica}
    cabecera.update({k: ds.meta.get(k) for k in FACTORES})
    return [{**cabecera, **fila, "error": ""} for fila in filas]


def _fila_de_error(celda: Celda, replica: int, exc: Exception) -> dict:
    """Deja constancia de una celda que no pudo evaluarse."""
    return {"bloque": celda.bloque, "celda": celda.nombre, "replica": replica,
            "seed": celda.seed_base + replica,
            **{k: celda.generate.get(k) for k in FACTORES},
            "rho": float("nan"),
            "error": f"{type(exc).__name__}: {exc}"}


# ===========================================================================
# Recorrido de la rejilla
# ===========================================================================

def lanzar(celdas: list[Celda], raiz: Path | str = RAIZ_POR_DEFECTO, *,
           replicas: int | None = None, forzar: bool = False,
           verbose: bool = True) -> pd.DataFrame:
    """Recorre las celdas dadas y devuelve solo lo evaluado en esta llamada.

    Parameters
    ----------
    replicas
        Sobrescribe el numero de replicas de cada celda. Util para una prueba
        rapida de la maquinaria antes de lanzar la rejilla de verdad.
    forzar
        Reevalua tambien lo que ya estuviese hecho.
    """
    raiz = Path(raiz)
    nuevas: list[dict] = []
    fallos = 0

    for celda in celdas:
        ruta = ruta_resultados(celda.bloque, raiz)
        hechas = set() if forzar else ya_hechas(ruta)
        n = replicas if replicas is not None else celda.replicas

        for replica in range(n):
            if (celda.nombre, replica) in hechas:
                continue
            t0 = time.perf_counter()
            try:
                filas = evaluar_replica(celda, replica, raiz)
            except Exception as exc:                    # noqa: BLE001
                filas = [_fila_de_error(celda, replica, exc)]
                fallos += 1
                if verbose:
                    print(f"  FALLO {celda.bloque}/{celda.nombre} "
                          f"rep{replica:02d}: {type(exc).__name__}: {exc}")
                    traceback.print_exc()
            else:
                if verbose:
                    rho = filas[0].get("rho")
                    print(f"  {celda.bloque}/{celda.nombre} rep{replica:02d}  "
                          f"{time.perf_counter() - t0:6.1f} s  "
                          f"rho={rho if rho is None else f'{rho:.3f}'}")
            _anadir(ruta, filas)
            nuevas += filas

    if verbose:
        print(f"\n{len(nuevas)} filas nuevas, {fallos} celdas con error.")
    return pd.DataFrame(nuevas)


# ===========================================================================
# CLI
# ===========================================================================

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--bloques", nargs="+", required=True,
                    choices=BLOQUES_DE_METRICAS,
                    help="bloques a lanzar, en orden.")
    ap.add_argument("--raiz", type=Path, default=RAIZ_POR_DEFECTO,
                    help="directorio de trabajo de la rejilla.")
    ap.add_argument("--replicas", type=int, default=None,
                    help="sobrescribe las replicas de cada celda.")
    ap.add_argument("--forzar", action="store_true",
                    help="reevalua tambien lo que ya estuviese hecho.")
    args = ap.parse_args()

    referencia = cargar_referencia(args.raiz)
    if referencia == REFERENCIA and set(args.bloques) - {"fase0"}:
        print("AVISO: no hay referencia calibrada en "
              f"{Path(args.raiz) / FICHERO_REFERENCIA}, se usan los valores "
              "provisionales. Corre antes la fase 0 y su analisis.")

    celdas = celdas_de(args.bloques, referencia)
    print(f"{len(celdas)} celdas en {args.bloques}, raiz {args.raiz}")
    t0 = time.perf_counter()
    lanzar(celdas, args.raiz, replicas=args.replicas, forzar=args.forzar)
    print(f"Total {(time.perf_counter() - t0) / 60:.1f} min")


if __name__ == "__main__":
    main()
