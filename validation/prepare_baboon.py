#!/usr/bin/env python3
"""
validation/prepare_baboon.py
============================
Adaptador del atlas de babuino de Mure et al. (2018), GEO GSE98965, al formato
que consume el pipeline de CIRCUST.

Por que hace falta un adaptador
-------------------------------
El fichero original es una unica matriz con las 767 muestras de los 64 tejidos
mezcladas, identificadas por el nombre de columna. CIRCUST estima el orden
temporal **dentro** de un tejido, de modo que hay que partirlo en 64 conjuntos
independientes. Ademas la hora de recogida esta codificada en el propio nombre
de la columna y hay que extraerla para disponer del ground truth.

Que produce
-----------
Un directorio por tejido dentro de ``data/babuinos/``, cada uno con:

``matriz.csv``       genes en filas (indice = simbolo), muestras en columnas.
                     Es lo que espera ``load_expression_matrix``.
``info.json``        procedencia, numero de genes y muestras, e incidencias.

Y en la raiz de ``data/babuinos/`` dos ficheros comunes:

``tiempos.csv``      columnas ``sample_name`` y ``true_time_h``, con las doce
                     horas ZT. Es el ground truth del orden, y sirve para todos
                     los tejidos: las muestras se nombran ``ZT00`` a ``ZT22`` en
                     todos ellos, de modo que las horas no dependen del tejido.
                     El unico caso particular es IRI, que no tiene ZT22 y por
                     tanto usa once de las doce filas.
``tejidos.csv``      indice con el listado y los datos de la Tabla A.

Particularidades del fichero original que este adaptador resuelve
-----------------------------------------------------------------
1. La columna ``Symbol`` es texto y no una muestra. Al forzar el tipo numerico
   se convierte en una columna entera de valores ausentes, y como afecta a
   todos los genes a la vez con una fraccion minima, ningun gen supera el
   filtro de dispersion y la columna acabaria colandose en el analisis como
   una muestra ficticia. Aqui se usa como indice de genes y se retira.
2. El separador entre tejido y hora no es homogeneo. La mayoria usa un punto,
   pero COR, IRI, RPE y SKI usan un guion bajo.
3. Ocho muestras llevan el sufijo ``.approx``, que en el original marca una
   hora de recogida aproximada. Se conservan, porque descartarlas dejaria esos
   tejidos con once muestras, pero se registran en ``info.json``.
4. El tejido IRI tiene once muestras en lugar de doce, le falta ZT22.
"""
from __future__ import annotations

import argparse
import json
import re
from collections import defaultdict
from pathlib import Path

import pandas as pd

PERIOD_HOURS = 24.0

# Patron de nombre de muestra: <TEJIDO><separador>ZT<hora>[.approx]
PATRON = re.compile(r"^([A-Za-z0-9]+)[._]ZT(\d+)(\.approx)?$")

# Tabla A del material suplementario de CIRCUST (Larriba et al., 2023).
#
# Contiene los 47 tejidos que ese trabajo retuvo para la validacion sobre el
# babuino, con tres datos por tejido:
#
#   nombre    nombre del tejido, transcrito literalmente de la tabla.
#   r2_ave    medida media de ritmicidad, el R2 del ajuste FMM de los 12 genes
#             semilla sobre los tiempos estimados. Es el criterio con el que
#             seleccionaron los tejidos: solo conservaron los que superan 0,7,
#             "to guarantee the consistency of the results".
#   rho       bondad de ajuste de la regresion circular-circular entre los
#             tiempos reales y los estimados por CIRCUST. Es su resultado, y
#             sirve de referencia contra la que comparar.
#
# Los 17 tejidos restantes del atlas no figuran en la tabla porque no superaron
# el filtro. Quedan sin nombre en el indice, y no se inventa ninguno.
#
# Nota sobre una probable errata del original: WAP aparece con el mismo nombre
# que WAM, "White Adipose Mesenteric". Se transcribe tal cual, sin corregir.
#
# Sobre rho y la variante CIRCUSTprior. Esas cifras se obtuvieron con
# CIRCUSTprior, que para el babuino sustituye el supuesto del gen de direccion
# por relaciones de orden especificas de la especie, {DBP} < {CRY1, CRY2} <
# {ARNTL} o {NR1D1} < {PER1, PER2, PER3} < {ARNTL}, escogiendo entre ambas la
# que deje mas genes semilla dentro del periodo activo. El pipeline de este
# trabajo implementa la variante estandar.
#
# Aun asi la comparacion con esta columna es legitima, porque esa variante
# afecta unicamente a la eleccion del sentido de recorrido y no al ajuste ni a
# ninguna otra etapa. Un cambio de sentido equivale a reflejar el orden
# completo, y las medidas de ``validation.metrics`` alinean por rotacion y por
# reflexion antes de evaluar, de modo que un sentido u otro producen el mismo
# valor. La unica salvedad es que se desconoce si la regresion circular con la
# que ellos calcularon rho contempla tambien la reflexion.
TABLA_A: dict[str, tuple[str, float, float]] = {
    "IRI": ("Iris", 0.811, 0.800),
    "COR": ("Cornea", 0.846, 0.958),
    "ONH": ("Optic Nerve Head", 0.807, 0.874),
    "RPE": ("Retina Pigment Epithelium", 0.740, 0.871),
    "ADC": ("Adrenal Cortex", 0.787, 0.947),
    "ADM": ("Adrenal Medulla", 0.799, 0.899),
    "THR": ("Thyroid", 0.810, 0.723),
    "LUN": ("Lung", 0.803, 0.911),
    "AXL": ("Axilliary Lymphonodes", 0.724, 0.615),
    "MEL": ("Mesenteric Lymphonodes", 0.746, 0.833),
    "SPL": ("Spleen", 0.703, 0.741),
    "LIV": ("Liver", 0.730, 0.777),
    "PAN": ("Pancreas", 0.812, 0.960),
    "OMF": ("Omental Fat", 0.796, 0.892),
    "WAM": ("White Adipose Mesenteric", 0.814, 0.964),
    "WAP": ("White Adipose Mesenteric", 0.745, 0.609),   # ver nota de errata
    "WAR": ("White Adipose Perirenal", 0.813, 0.708),
    "WAS": ("White Adipose Subcutaneous", 0.811, 0.878),
    "WAT": ("White Adipose Tissue", 0.892, 0.961),
    "SKI": ("Skin", 0.788, 0.953),
    "BLA": ("Bladder", 0.745, 0.880),
    "KIC": ("Kidney Cortex", 0.811, 0.765),
    "KIM": ("Kidney Medulla", 0.764, 0.728),
    "AMY": ("Amygdala", 0.702, 0.686),
    "ARC": ("Arcuate Nucleus", 0.745, 0.680),
    "DMH": ("Dorsomedial Hypothalamus", 0.703, 0.728),
    "HIP": ("Hippocampus", 0.704, 0.734),
    "MMB": ("Mammilary Bodies", 0.735, 0.789),
    "PRA": ("Preoptic Area", 0.702, 0.741),
    "PUT": ("Putamen", 0.784, 0.813),
    "PVN": ("Paraventricular Nuclei", 0.728, 0.744),
    "SCN": ("Suprachiasmatic Nuclei", 0.711, 0.729),
    "VIC": ("Visual Cortex", 0.740, 0.706),
    "VMH": ("Ventromedial Hypothalamus", 0.776, 0.806),
    "AOR": ("Aorta", 0.742, 0.660),
    "HEA": ("Heart", 0.841, 0.935),
    "MUA": ("Muscle Abdominal", 0.799, 0.953),
    "MUG": ("Muscle Gastrocnemian", 0.701, 0.884),
    "ANT": ("Antrum", 0.835, 0.740),
    "ASC": ("Ascending Colon", 0.746, 0.959),
    "CEC": ("Cecum", 0.701, 0.734),
    "DEC": ("Descending Colon", 0.739, 0.846),
    "DUO": ("Duodenum", 0.816, 0.747),
    "ILE": ("Ileon", 0.709, 0.726),
    "OES": ("Oesophagus", 0.903, 0.973),
    "SMM": ("Smooth Muscle", 0.741, 0.674),
    "STF": ("Stomach Fundus", 0.824, 0.734),
}


def parsear_columnas(columnas: list[str]) -> dict[str, list[tuple]]:
    """Agrupa las columnas de muestra por tejido.

    Devuelve ``{codigo: [(hora_zt, es_aproximada, nombre_columna), ...]}``,
    ordenado por hora. Las columnas que no son muestra se ignoran.
    """
    por_tejido: dict[str, list[tuple]] = defaultdict(list)
    for col in columnas:
        m = PATRON.fullmatch(col)
        if m:
            por_tejido[m.group(1)].append(
                (int(m.group(2)), bool(m.group(3)), col)
            )
    return {k: sorted(v) for k, v in sorted(por_tejido.items())}


def preparar(origen: Path, destino: Path) -> pd.DataFrame:
    """Parte el atlas en un conjunto por tejido. Devuelve el indice."""
    print(f"Leyendo {origen} ...")
    bruto = pd.read_csv(origen)
    print(f"  {bruto.shape[0]} genes x {bruto.shape[1]} columnas")

    if "Symbol" not in bruto.columns:
        raise ValueError("El fichero no tiene columna 'Symbol'.")

    # El simbolo del gen pasa a ser el indice. EnsemblID y Symbol dejan de ser
    # columnas, de modo que no pueden colarse como muestras.
    simbolos = bruto["Symbol"].astype(str).str.strip()
    bruto = bruto.drop(columns=[c for c in ("EnsemblID", "Symbol")
                                if c in bruto.columns])
    bruto.index = simbolos
    bruto.index.name = "gene_symbol"

    por_tejido = parsear_columnas(list(bruto.columns))
    print(f"  {len(por_tejido)} tejidos reconocidos")

    destino.mkdir(parents=True, exist_ok=True)
    filas: list[dict] = []
    horas_vistas: set[int] = set()

    for codigo, muestras in por_tejido.items():
        horas = [h for h, _, _ in muestras]
        aprox = [c for _, a, c in muestras if a]
        cols = [c for _, _, c in muestras]

        sub = bruto[cols].copy()
        # Nombres de muestra homogeneos, con la hora explicita
        sub.columns = [f"ZT{h:02d}" for h in horas]

        carpeta = destino / codigo
        carpeta.mkdir(exist_ok=True)
        sub.to_csv(carpeta / "matriz.csv", index=True, index_label="gene_symbol")

        ref = TABLA_A.get(codigo)
        horas_vistas.update(horas)
        info = dict(
            codigo=codigo,
            nombre=ref[0] if ref else None,
            en_tabla_A=ref is not None,
            r2_ave_circust=ref[1] if ref else None,
            rho_circust=ref[2] if ref else None,
            origen=str(origen),
            referencia="Mure et al. (2018), GEO GSE98965",
            n_genes=int(sub.shape[0]),
            n_muestras=int(sub.shape[1]),
            horas_zt=horas,
            muestras_hora_aproximada=aprox,
            completo=(len(horas) == 12),
        )
        (carpeta / "info.json").write_text(json.dumps(info, indent=2,
                                                      ensure_ascii=False))
        filas.append(info)
        marca = "" if info["completo"] else f"  <- solo {len(horas)} muestras"
        aviso = f"  ({len(aprox)} con hora aproximada)" if aprox else ""
        print(f"  {codigo:<5} {sub.shape[1]:>3} muestras"
              f"{aviso}{marca}")

    pd.DataFrame({
        "sample_name": [f"ZT{h:02d}" for h in sorted(horas_vistas)],
        "true_time_h": [float(h % PERIOD_HOURS) for h in sorted(horas_vistas)],
    }).to_csv(destino / "tiempos.csv", index=False)

    indice = pd.DataFrame(filas)[
        ["codigo", "nombre", "n_genes", "n_muestras", "completo",
         "en_tabla_A", "r2_ave_circust", "rho_circust"]
    ]
    indice.to_csv(destino / "tejidos.csv", index=False)
    return indice


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Parte el atlas de babuino en un conjunto por tejido."
    )
    raiz = Path(__file__).resolve().parents[1]
    ap.add_argument("--origen", type=Path,
                    default=raiz / "data" / "GSE98965_baboon_tissue_expression_FPKM.csv")
    ap.add_argument("-o", "--destino", type=Path,
                    default=raiz / "data" / "babuinos")
    args = ap.parse_args()

    indice = preparar(args.origen, args.destino)
    print(f"\n{len(indice)} tejidos escritos en {args.destino}")
    print(f"  completos (12 muestras)      : {int(indice.completo.sum())}")
    print(f"  retenidos por CIRCUST        : {int(indice.en_tabla_A.sum())}")
    print(f"  descartados por su filtro    : {int((~indice.en_tabla_A).sum())}")


if __name__ == "__main__":
    main()
