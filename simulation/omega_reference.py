#!/usr/bin/env python3
"""
simulation/omega_reference.py
=============================
Mide la distribucion empirica del apuntamiento ``omega`` sobre atlas ya
calculados, y la resume como un reparto entre las clases de forma de onda que
usa el generador sintetico.

Para que sirve
--------------
Los valores de ``WAVEFORM_CLASSES`` y las proporciones de ``WAVEFORM_MIXES``
en ``simulation/generate.py`` no son elegidos a ojo, salen de esta medicion.
El script existe para que esa tabla sea reproducible y para poder rehacerla si
se anaden atlas nuevos.

Que compara
-----------
Tres fuentes, que responden a tres preguntas distintas:

1. **GTEx con 480 muestras y FMM base.** Es la referencia de verdad. Con ese
   tamano muestral el apuntamiento esta bien determinado, de modo que sirve
   para decidir como debe ser la mezcla realista del generador.
2. **Babuino con 12 muestras y FMM base.** Es el sintoma. Muestra a donde se va
   el estimador cuando las muestras escasean.
3. **Babuino con 12 muestras y FMM estabilizado.** Es la prueba de que lo
   anterior era en buena parte artefacto de estimacion y no biologia.

Advertencia metodologica
------------------------
La verdad del generador se calibra con la fuente 1 y nunca con las fuentes 2 o
3. Generar datos sinteticos con los omega que el babuino estima seria simular
un artefacto del estimador para despues medir como se corrige, que es un
argumento circular.

Uso
---
    python -m simulation.omega_reference
    python -m simulation.omega_reference --salida tabla_omega.csv
"""
from __future__ import annotations

import argparse
import glob
from pathlib import Path

import numpy as np
import pandas as pd

from simulation.generate import WAVEFORM_CLASSES

# Umbral de ajuste con el que se filtra que gen se considera bien descrito por
# el modelo. Es el mismo que usa CIRCUST para declarar ritmicidad.
R2_MIN: float = 0.5

# Fuentes por defecto, relativas a la raiz del repositorio.
FUENTES: dict[str, str] = {
    "GTEx 480 muestras, FMM base":
        "output/circust_atlas.csv",
    "Babuino 12 muestras, FMM base":
        "output/baboon/base_auto/*/run/circust_atlas.csv",
    "Babuino 12 muestras, FMM estabilizado":
        "output/baboon/base_stab/*/run/circust_atlas.csv",
}


def bordes_de_clase() -> tuple[list[str], np.ndarray]:
    """Nombres de clase de mayor a menor omega, y los cortes entre ellas.

    Cada gen se asigna a la clase cuyo omega tiene mas cerca, de modo que los
    cortes son los puntos medios entre valores consecutivos.
    """
    nombres = sorted(WAVEFORM_CLASSES, key=WAVEFORM_CLASSES.get)
    valores = np.array([WAVEFORM_CLASSES[c] for c in nombres])
    cortes = (valores[:-1] + valores[1:]) / 2.0
    return nombres, cortes


def omegas_de(patron: str, r2_min: float = R2_MIN) -> np.ndarray:
    """Todos los omega de los atlas que casan con ``patron``, ya filtrados."""
    rutas = sorted(glob.glob(patron))
    if not rutas:
        raise FileNotFoundError(f"ningun atlas casa con {patron!r}.")
    trozos = []
    for ruta in rutas:
        tabla = pd.read_csv(ruta)
        faltan = {"omega", "r2_fmm"} - set(tabla.columns)
        if faltan:
            raise ValueError(f"{ruta} no tiene las columnas {sorted(faltan)}.")
        trozos.append(tabla.loc[tabla["r2_fmm"] >= r2_min, "omega"].to_numpy())
    return np.concatenate(trozos)


def reparto(omegas: np.ndarray) -> dict[str, float]:
    """Porcentaje de genes que cae en cada clase de forma de onda."""
    nombres, cortes = bordes_de_clase()
    indice = np.digitize(omegas, cortes)
    return {c: 100.0 * float(np.mean(indice == i))
            for i, c in enumerate(nombres)}


def tabla(fuentes: dict[str, str], r2_min: float = R2_MIN) -> pd.DataFrame:
    """Una fila por clase, una columna por fuente, mas mediana y tamano."""
    nombres, _ = bordes_de_clase()
    columnas: dict[str, list[float]] = {}
    for etiqueta, patron in fuentes.items():
        om = omegas_de(patron, r2_min)
        pct = reparto(om)
        columnas[etiqueta] = ([pct[c] for c in nombres]
                              + [float(np.median(om)), float(om.size)])
    return pd.DataFrame(
        columnas,
        index=[f"{c} (omega {WAVEFORM_CLASSES[c]:.2f})" for c in nombres]
              + ["mediana de omega", "genes medidos"],
    )


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--raiz", type=Path, default=Path("."),
                    help="raiz del repositorio, por si se invoca desde fuera.")
    ap.add_argument("--r2-min", type=float, default=R2_MIN,
                    help="umbral de R2_FMM para considerar un gen ajustado.")
    ap.add_argument("--salida", type=Path, default=None,
                    help="si se indica, escribe la tabla como CSV.")
    args = ap.parse_args()

    fuentes = {k: str(args.raiz / v) for k, v in FUENTES.items()}
    resumen = tabla(fuentes, args.r2_min)

    print(f"Distribucion de omega, genes con R2_FMM >= {args.r2_min}")
    print(resumen.round(1).to_string())
    print("\nLas tres primeras filas son porcentajes de genes.")

    if args.salida is not None:
        resumen.to_csv(args.salida)
        print(f"\nEscrito en {args.salida}")


if __name__ == "__main__":
    main()
