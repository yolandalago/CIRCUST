#!/usr/bin/env python3
"""
simulation/analyze.py
=====================
Analisis de la fase de calibracion: elige la celda de referencia y resume el
barrido del umbral tau.

Por que hace falta calibrar
---------------------------
Una celda con SNR 3 satura. Todos los metodos dan rho por encima de 0,99 y la
desviacion tipica de rho entre las K repeticiones baja a 0,001, de modo que
ninguna aportacion puede distinguirse de ninguna otra. Sobre babuino real el
rho es 0,703. Montar la rejilla en el regimen saturado produciria un capitulo
plano por construccion, y ademas es donde trabaja el estudio de simulacion del
material suplementario de CIRCUST.

El criterio, fijado antes de mirar resultados
---------------------------------------------
1. ``rho`` mediano en el entorno de 0,8, que es el regimen de los datos reales
   etiquetados disponibles.
2. MedAE del orden entre una y dos horas, comparable a la mediana de 1,94 h
   medida sobre babuino.
3. Desviacion tipica de ``rho`` entre las K repeticiones por encima de 0,03.

El tercero es el que mas se olvida y el mas importante para el bloque A. En los
47 tejidos de babuino vale 0,012 y la agregacion no aporta nada, con un
mecanismo identificado: si las K repeticiones coinciden, no hay desacuerdo que
resolver. En GSE71620, con 146 muestras, esta entre 0,03 y 0,15, que es el
regimen donde si tiene algo que hacer. Calibrar por debajo de 0,03 seria situar
el estudio donde ya se sabe que no informa, y concluir que agregar no sirve
cuando lo que ocurre es que no se le ha dado nada que agregar.

A los tres se anade una **condicion previa**, que no es un criterio sino un
requisito: la celda tiene que correr en todas sus replicas. En el pilotaje,
varias celdas de SNR bajo con muchas muestras no llegan a terminar porque el
conjunto TOP colapsa, y una celda asi no puede ser la referencia.

Uso
---
    python -m simulation.analyze fase0
    python -m simulation.analyze fase0 --escribir
    python -m simulation.analyze tau
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from simulation.grid import RAIZ_POR_DEFECTO, REFERENCIA

# Configuracion sobre la que se evalua el criterio. Es la que representa el
# resumen robusto habitual, y se fija para que la eleccion no dependa de cual
# de los metodos se mire.
CONFIG_CRITERIO: tuple[str, str] = ("base", "best_k:r2_median")

VENTANA_RHO: tuple[float, float] = (0.70, 0.90)
VENTANA_MEDAE_H: tuple[float, float] = (1.0, 2.0)
SD_RHO_MINIMA: float = 0.03

# Objetivo dentro de la ventana. Solo se usa para desempatar si varias celdas
# cumplen, y esta fijado de antemano igual que las ventanas.
RHO_OBJETIVO: float = 0.80

FICHERO_REFERENCIA: str = "referencia.json"

# Parametros de la celda que se trasladan a la referencia. El resto se hereda
# de la referencia provisional.
EJES_CALIBRADOS: tuple[str, ...] = ("snr", "n_samples")


# ===========================================================================
# Lectura
# ===========================================================================

def leer(bloque: str, raiz: Path | str = RAIZ_POR_DEFECTO) -> pd.DataFrame:
    ruta = Path(raiz) / "resultados" / f"{bloque}.csv"
    if not ruta.exists():
        raise SystemExit(f"No existe {ruta}. Lanza antes ese bloque.")
    return pd.read_csv(ruta, keep_default_na=False, na_values=[""])


def _sin_errores(d: pd.DataFrame) -> pd.DataFrame:
    if "error" not in d.columns:
        return d
    return d[d["error"].isna() | (d["error"] == "")]


# ===========================================================================
# Criterio de la fase 0
# ===========================================================================

def criterio_fase0(d: pd.DataFrame) -> pd.DataFrame:
    """Una fila por celda, con las tres condiciones evaluadas.

    ``tasa_exito`` es la fraccion de replicas que llegaron a terminar. Se
    calcula sobre todas las replicas de la celda, no solo sobre las que
    produjeron numeros, porque si no una celda que falla la mitad de las veces
    pareceria tan buena como una que nunca falla.
    """
    fmm, metodo = CONFIG_CRITERIO
    replicas = d.groupby("celda")["replica"].nunique()

    util = _sin_errores(d)
    util = util[(util["fmm_method"] == fmm) & (util["metodo"] == metodo)]

    agrupado = util.groupby("celda")
    t = agrupado.agg(
        snr=("snr", "first"),
        n_samples=("n_samples", "first"),
        n_ok=("replica", "nunique"),
        rho=("rho", "median"),
        mad_h=("mad_h", "median"),
        sd_rho_reps=("sd_rho_reps", "median"),
        n_top=("n_top", "median"),
    )
    # Las celdas que fallaron en TODAS sus replicas no aparecen en util, y
    # tienen que salir igualmente en el resumen.
    t = t.reindex(replicas.index)
    t["n_ok"] = t["n_ok"].fillna(0).astype(int)
    t["tasa_exito"] = t["n_ok"] / replicas

    t["corre"] = t["tasa_exito"] >= 1.0
    t["cumple_rho"] = t["rho"].between(*VENTANA_RHO)
    t["cumple_medae"] = t["mad_h"].between(*VENTANA_MEDAE_H)
    t["cumple_sd"] = t["sd_rho_reps"] > SD_RHO_MINIMA
    t["cumple"] = (t["corre"] & t["cumple_rho"] & t["cumple_medae"]
                   & t["cumple_sd"])
    return t.reset_index()


def elegir_referencia(d: pd.DataFrame,
                      referencia_base: dict | None = None) -> dict:
    """La celda de referencia, a partir de los resultados de la fase 0.

    Si varias celdas cumplen, gana la de ``rho`` mas cercano a
    :data:`RHO_OBJETIVO`. El desempate esta fijado de antemano igual que las
    ventanas, de modo que no se elige mirando cual favorece a que metodo.
    """
    t = criterio_fase0(d)
    candidatas = t[t["cumple"]]
    if candidatas.empty:
        raise ValueError(
            "ninguna celda de la fase 0 cumple el criterio. Hay que ampliar la "
            "rejilla de calibracion, no relajar el criterio a posteriori."
        )
    ganadora = candidatas.iloc[
        (candidatas["rho"] - RHO_OBJETIVO).abs().to_numpy().argmin()]
    base = dict(referencia_base or REFERENCIA)
    for eje in EJES_CALIBRADOS:
        valor = ganadora[eje]
        base[eje] = int(valor) if eje == "n_samples" else float(valor)
    return base


def escribir_referencia(referencia: dict,
                        raiz: Path | str = RAIZ_POR_DEFECTO) -> Path:
    """Deja el fichero que ``simulation.run`` lee para congelar la referencia."""
    raiz = Path(raiz)
    raiz.mkdir(parents=True, exist_ok=True)
    ruta = raiz / FICHERO_REFERENCIA
    ruta.write_text(json.dumps(referencia, indent=2, ensure_ascii=False,
                               sort_keys=True))
    return ruta


# ===========================================================================
# Barrido de tau
# ===========================================================================

def _tau_de(metodo: str) -> float:
    return float(str(metodo).rsplit("@", 1)[-1])


def resumen_tau(d: pd.DataFrame) -> pd.DataFrame:
    """Una fila por umbral, con lo que hace falta para decidir sobre tau.

    ``degenerado`` marca los umbrales en los que ningun gen del TOP supera el
    corte. Cuando eso pasa todas las repeticiones empatan a cero y el criterio
    desempata por mediana de R2, es decir, **colapsa exactamente en
    ``best_k:r2_median``**. Un umbral degenerado no es que rinda mal, es que ha
    dejado de ser un criterio distinto.

    ``variedad_reps`` cuenta cuantas repeticiones distintas llega a elegir el
    criterio a lo largo de las replicas. Si vale uno en todo el rango util, la
    conclusion legitima no es un valor de tau sino que el criterio no aporta y
    procede retirarlo junto con su umbral.
    """
    util = _sin_errores(d).copy()
    util["tau"] = util["metodo"].map(_tau_de)

    t = util.groupby("tau").agg(
        rho=("rho", "median"),
        mad_h=("mad_h", "median"),
        score=("score_seleccion", "median"),
        variedad_reps=("rep_elegida", "nunique"),
        n=("replica", "size"),
    )
    t["degenerado"] = t["score"] <= 0
    return t.reset_index()


# ===========================================================================
# CLI
# ===========================================================================

def _informe_fase0(d: pd.DataFrame, raiz: Path, escribir: bool) -> None:
    t = criterio_fase0(d)
    columnas = ["celda", "snr", "n_samples", "tasa_exito", "rho", "mad_h",
                "sd_rho_reps", "n_top", "cumple"]
    print("Fase 0, criterio fijado de antemano")
    print(f"  rho en {VENTANA_RHO}, MedAE en {VENTANA_MEDAE_H} h, "
          f"sd de rho > {SD_RHO_MINIMA}, y la celda tiene que correr siempre")
    print(f"  medido sobre {CONFIG_CRITERIO[0]} + {CONFIG_CRITERIO[1]}\n")
    print(t[columnas].round(3).to_string(index=False))

    rotas = t[t["tasa_exito"] < 1.0]
    if not rotas.empty:
        print("\nCeldas que no corren siempre, con su tasa de exito:")
        print(rotas[["celda", "tasa_exito"]].round(2).to_string(index=False))

    try:
        ref = elegir_referencia(d)
    except ValueError as exc:
        print(f"\n{exc}")
        return
    print("\nCelda de referencia elegida:")
    print("  " + ", ".join(f"{k}={ref[k]}" for k in EJES_CALIBRADOS))
    if escribir:
        print(f"  escrita en {escribir_referencia(ref, raiz)}")
    else:
        print("  (usa --escribir para congelarla y que run.py la lea)")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("bloque", choices=["fase0", "tau"])
    ap.add_argument("--raiz", type=Path, default=RAIZ_POR_DEFECTO)
    ap.add_argument("--escribir", action="store_true",
                    help="congela la referencia elegida en referencia.json.")
    args = ap.parse_args()

    d = leer(args.bloque, args.raiz)
    if args.bloque == "fase0":
        _informe_fase0(d, args.raiz, args.escribir)
        return

    t = resumen_tau(d)
    print("Barrido de tau sobre best_k:n_well_fitted\n")
    print(t.round(3).to_string(index=False))
    utiles = t[~t["degenerado"]]
    if utiles.empty:
        print("\nTodos los umbrales degeneran: el criterio no discrimina.")
    elif (utiles["variedad_reps"] <= 1).all():
        print("\nEn todo el rango util el criterio elige siempre la misma "
              "repeticion, luego no aporta frente a r2_median y procede "
              "retirarlo junto con su umbral.")


if __name__ == "__main__":
    main()
