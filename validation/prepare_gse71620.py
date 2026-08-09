#!/usr/bin/env python3
"""Prepara GSE71620 (corteza prefrontal humana post mortem) para validacion.

Es el segundo conjunto etiquetado del trabajo, junto con el babuino, y cubre
un regimen muestral distinto: frente a las 12 muestras por tejido del babuino,
aqui hay ~100 por region cerebral. Eso es lo que permite comprobar si el
tamano objetivo del conjunto TOP debe escalar con el numero de muestras o si
un valor absoluto vale en todo el rango.

Es el mismo conjunto que usa el articulo CIRCUST para contrastar contra TOD
(Fig 6 del articulo principal).

Entradas (ya en data/GSE71620/)
-------------------------------
    GSE71620_series_matrix.txt.gz   expresion por sonda + metadatos de muestra
    GPL11532.annot.gz               anotacion sonda -> simbolo de gen

Particularidades resueltas
--------------------------
1. La matriz viene por SONDA (Affymetrix Human Gene 1.1 ST), no por gen. Se
   mapea con la anotacion de la plataforma y las sondas sin simbolo se tiran.
2. Varias sondas apuntan al mismo gen. Se resuelve quedandose con la de mayor
   MAD, que es el criterio que ya usa el preprocesado del pipeline para
   duplicados.
3. Un mismo simbolo puede venir como lista separada por "///" cuando el
   conjunto de sondas cae en una region que la anotacion asigna a varios genes
   a la vez. El valor medido es entonces una mezcla que no se puede separar,
   asi que esas sondas se descartan (2 354 de las 33 297 de la plataforma).

   CONSECUENCIA ASUMIDA: se pierden PER3 (UTS2///PER3) y NR1D1
   (NR1D1///THRA), porque esas son sus unicas sondas en GPL11532. NR1D1 y
   THRA estan solapados en hebras opuestas del cromosoma 17. El pipeline
   trabaja entonces con 10 de los 12 genes core. Se prefiere perderlos a
   introducirlos contaminados: NR1D1 es uno de los genes que fijan la
   direccion del orden, y un pico mal estimado ahi puede voltear la
   sincronizacion del conjunto entero. Los genes core ausentes se registran
   en info.json.
4. Las 420 muestras son 210 de BA47 y 210 de BA11, dos regiones distintas del
   mismo sujeto. Se separan: son experimentos independientes.
5. El TOD viene en horas sobre una escala que va de -5,96 a 17,40, es decir,
   la hora del reloj representada con el corte en las 18:00 en lugar de en
   medianoche. Se lleva a [0, 24) con modulo.

   El muestreo NO es uniforme: 105 de los 146 sujetos fallecieron entre las
   4 y las 16, y solo 21 entre las 20 y las 4. Los genes que pican de noche
   se estiman por tanto con muchas menos muestras.
6. No todas las muestras tienen TOD. Las que traen "NA" se descartan, porque
   sin hora conocida no sirven para validar.

Salida
------
    data/GSE71620/BA11/{matriz.csv,tiempos.csv,info.json}
    data/GSE71620/BA47/{matriz.csv,tiempos.csv,info.json}
    data/GSE71620/regiones.csv
"""
from __future__ import annotations

import gzip
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

ROOT   = Path(__file__).resolve().parent.parent
import sys
sys.path.insert(0, str(ROOT / "src"))
from circust.core_genes import CoreGeneSelector          # noqa: E402

CRUDO  = ROOT / "data" / "GSE71620"
MATRIZ = CRUDO / "GSE71620_series_matrix.txt.gz"
ANOTA  = CRUDO / "GPL11532.annot.gz"


def leer_metadatos() -> pd.DataFrame:
    """Extrae del series matrix una fila por muestra: GSM, region y TOD."""
    campos: dict[str, list[str]] = {}
    with gzip.open(MATRIZ, "rt", encoding="utf-8", errors="replace") as fh:
        for linea in fh:
            if linea.startswith("!series_matrix_table_begin"):
                break
            if not linea.startswith("!Sample_"):
                continue
            partes = [p.strip('"') for p in linea.rstrip("\n").split("\t")]
            clave, valores = partes[0], partes[1:]
            # characteristics_ch1 se repite; se distingue por su prefijo
            if clave == "!Sample_characteristics_ch1":
                etiqueta = valores[0].split(":")[0].strip().lower()
                campos[f"car_{etiqueta}"] = valores
            else:
                campos.setdefault(clave, valores)

    gsm    = campos["!Sample_geo_accession"]
    titulo = campos["!Sample_title"]
    tod    = campos["car_tod"]

    meta = pd.DataFrame({
        "gsm":    gsm,
        "titulo": titulo,
        "region": [t.split()[0] for t in titulo],
        "tod_bruto": [v.split(":", 1)[1].strip() for v in tod],
    })
    meta["tod_h"] = pd.to_numeric(meta["tod_bruto"], errors="coerce") % 24.0
    return meta


def leer_expresion() -> pd.DataFrame:
    """Tabla de expresion por sonda (indice = ID de sonda)."""
    with gzip.open(MATRIZ, "rt", encoding="utf-8", errors="replace") as fh:
        for n, linea in enumerate(fh):
            if linea.startswith("!series_matrix_table_begin"):
                salto = n + 1
                break
    expr = pd.read_csv(MATRIZ, sep="\t", skiprows=salto, comment="!",
                       index_col=0, quotechar='"')
    expr.index = expr.index.astype(str)
    expr.index.name = "sonda"
    return expr


def leer_anotacion() -> pd.Series:
    """Mapa sonda -> simbolo de gen, ya limpio de vacios y ambiguos."""
    with gzip.open(ANOTA, "rt", encoding="utf-8", errors="replace") as fh:
        for n, linea in enumerate(fh):
            if linea.startswith("!platform_table_begin"):
                salto = n + 1
                break
    tabla = pd.read_csv(ANOTA, sep="\t", skiprows=salto, comment="!",
                        low_memory=False)
    tabla = tabla[["ID", "Gene symbol"]].dropna()
    tabla["ID"] = tabla["ID"].astype(str)
    simbolo = tabla.set_index("ID")["Gene symbol"].astype(str).str.strip()
    # Sondas que cubren varios genes: ambiguas, fuera.
    simbolo = simbolo[~simbolo.str.contains("///", regex=False)]
    simbolo = simbolo[simbolo.ne("") & simbolo.ne("nan")]
    return simbolo


def sondas_a_genes(expr: pd.DataFrame, simbolo: pd.Series) -> pd.DataFrame:
    """Colapsa sondas a genes quedandose, por gen, con la sonda de mayor MAD.

    El MAD es el mismo criterio con el que el preprocesado del pipeline
    resuelve genes duplicados, de modo que el conjunto llega en las mismas
    condiciones que el resto.
    """
    comunes = expr.index.intersection(simbolo.index)
    expr = expr.loc[comunes]
    genes = simbolo.loc[comunes]

    valores = expr.to_numpy(dtype=np.float64)
    mad = np.median(np.abs(valores - np.median(valores, axis=1, keepdims=True)),
                    axis=1)
    orden = pd.DataFrame({"gen": genes.to_numpy(), "mad": mad},
                         index=expr.index)
    elegidas = (orden.sort_values("mad", ascending=False)
                     .groupby("gen", sort=False).head(1))

    salida = expr.loc[elegidas.index]
    salida.index = elegidas["gen"].to_numpy()
    salida.index.name = "gene_id"
    return salida.sort_index()


def main() -> None:
    for f in (MATRIZ, ANOTA):
        if not f.exists():
            raise SystemExit(f"Falta {f}")

    print("Leyendo metadatos ...")
    meta = leer_metadatos()
    print(f"  {len(meta)} muestras; con TOD: {meta['tod_h'].notna().sum()}")

    print("Leyendo expresion por sonda ...")
    expr = leer_expresion()
    print(f"  {expr.shape[0]} sondas x {expr.shape[1]} muestras")

    print("Leyendo anotacion de la plataforma ...")
    simbolo = leer_anotacion()
    print(f"  {len(simbolo)} sondas con simbolo univoco")

    print("Colapsando sondas a genes ...")
    expr = sondas_a_genes(expr, simbolo)
    print(f"  {expr.shape[0]} genes")

    core = list(CoreGeneSelector.PRESETS["circust"])
    core_ausentes = [g for g in core if g not in expr.index]
    if core_ausentes:
        print(f"  AVISO: genes core sin sonda univoca: {core_ausentes} "
              f"({len(core) - len(core_ausentes)}/{len(core)} disponibles)")

    filas = []
    for region in sorted(meta["region"].unique()):
        sel = meta[(meta["region"] == region) & meta["tod_h"].notna()]
        sub = expr[sel["gsm"].tolist()].copy()
        # Nombre de muestra legible y estable: BA11_GSM1841558
        sub.columns = [f"{region}_{g}" for g in sel["gsm"]]

        destino = CRUDO / region
        destino.mkdir(parents=True, exist_ok=True)
        sub.to_csv(destino / "matriz.csv")
        pd.DataFrame({
            "sample_name": sub.columns,
            "true_time_h": sel["tod_h"].to_numpy(),
        }).to_csv(destino / "tiempos.csv", index=False)
        (destino / "info.json").write_text(json.dumps({
            "gse":        "GSE71620",
            "plataforma": "GPL11532",
            "region":     region,
            "descripcion": "Corteza prefrontal humana post mortem, area de Brodmann "
                           + region.replace("BA", ""),
            "n_genes":    int(sub.shape[0]),
            "n_muestras": int(sub.shape[1]),
            "n_descartadas_sin_tod": int((meta["region"] == region).sum() - len(sel)),
            "tiempo": "TOD (time of death) en horas, llevado a [0,24) con modulo. "
                      "El valor original va de -5,96 a 17,40: es la hora del "
                      "reloj con el corte en las 18:00. Muestreo no uniforme, "
                      "concentrado en horas diurnas.",
            "core_disponibles": [g for g in core if g in sub.index],
            "core_ausentes": core_ausentes,
            "core_ausentes_motivo":
                "Sus unicas sondas en GPL11532 son multigen (UTS2///PER3, "
                "NR1D1///THRA) y se descartan por no poder separar la senal. "
                "Ver el docstring de validation/prepare_gse71620.py.",
        }, indent=2, ensure_ascii=False), encoding="utf-8")

        filas.append({"region": region, "n_genes": sub.shape[0],
                      "n_muestras": sub.shape[1],
                      "n_sin_tod": int((meta["region"] == region).sum() - len(sel))})
        print(f"  [{region}] {sub.shape[0]} genes x {sub.shape[1]} muestras "
              f"({filas[-1]['n_sin_tod']} descartadas sin TOD)")

    pd.DataFrame(filas).to_csv(CRUDO / "regiones.csv", index=False)
    print(f"\nGuardado en {CRUDO}")


if __name__ == "__main__":
    main()
