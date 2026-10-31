"""
Precision de la implementacion en R frente a la de Python, sobre la misma
matriz sintetica y contra los tiempos verdaderos.

Por que hace falta
------------------
El bloque computacional mide que Python es unas 77 veces mas rapido que R.
Ese dato no vale de nada si la re-implementacion no es al menos igual de
precisa, asi que hay que comprobarlo. Como los datos son sinteticos tenemos
la verdad, de modo que no nos limitamos a ver si las dos implementaciones se
parecen: podemos decir cual acierta mas.

Que se compara con que
----------------------
El CIRCUST original resume las K repeticiones con una mediana de parametros
FMM y **no las colapsa en un orden de muestras**. La unica fase por muestra
que publica sale del orden preliminar, la Etapa 2. De ahi que haya tres
comparaciones y no una:

``preliminar R`` contra ``preliminar Python``
    Fidelidad de la re-implementacion. Deberian coincidir.
cada uno contra la verdad
    Precision. Es lo que acompana al argumento de velocidad.
``robusto Python`` contra ``preliminar R``
    Que anade la aportacion, porque en R ese orden no existe.

Uso
---
    python -m simulation.comparar_r --in output/simulacion
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from simulation.metrics import evaluate
from simulation.truth import TWO_PI

# Nombre legible de cada salida comparada.
SALIDAS = {
    "r_preliminar":  "R, orden preliminar",
    "py_preliminar": "Python, orden preliminar",
    "py_robusto":    "Python, orden robusto",
}


def _leer(celda: Path, fichero: str) -> pd.Series | None:
    """Fases indexadas por muestra, o ``None`` si esa salida no existe."""
    ruta = celda / fichero
    if not ruta.exists():
        return None
    d = pd.read_csv(ruta)
    return pd.Series(d["phase"].to_numpy(dtype=float) % TWO_PI,
                     index=d["sample"].astype(str))


def salidas_de_celda(celda: Path) -> dict[str, pd.Series]:
    """Las tres escalas estimadas de una celda, las que existan."""
    return {k: v for k, v in (
        ("r_preliminar",  _leer(celda, "r_fases.csv")),
        ("py_preliminar", _leer(celda, "py_fases_preliminar.csv")),
        ("py_robusto",    _leer(celda, "py_fases_robusto.csv")),
    ) if v is not None}


def _giro_a_la_verdad(theta_true: np.ndarray, theta_hat: np.ndarray) -> float:
    """
    Giro en horas que mejor alinea la salida con la verdad, con su signo.

    Es el diagnostico que distingue un fallo de anclaje de una diferencia de
    convenio. Si R sale desplazado la misma cantidad en todas las celdas, su
    origen es otro y basta con declararlo. Si el desplazamiento varia de celda
    a celda, entonces su anclaje falla de verdad.
    """
    d = (np.asarray(theta_hat, float) - np.asarray(theta_true, float))
    c, s = np.cos(d).mean(), np.sin(d).mean()
    # con signo y en (-12, 12]: un desfase de 22,7 h es en realidad de 1,3 h
    # hacia atras, y en el rango 0 a 24 se leeria como un desastre
    giro = np.arctan2(s, c) * 24.0 / TWO_PI
    return float(giro)


def contra_la_verdad(raiz: Path) -> pd.DataFrame:
    """
    Precision de cada salida frente a los tiempos verdaderos, celda a celda.

    Las metricas son las mismas del resto del estudio, calculadas sobre las
    muestras que esa salida conserva. R y Python pueden descartar muestras
    distintas como atipicas, asi que se informa cuantas se evaluan.
    """
    filas = []
    for celda in sorted(Path(raiz).glob("*/")):
        verdad = celda / "verdad.csv"
        if not verdad.exists():
            continue
        v = pd.read_csv(verdad)
        theta_true = pd.Series(
            v["t_true"].to_numpy(dtype=float) * TWO_PI / 24.0,
            index=v["sample"].astype(str))
        partes = celda.name.rstrip("/").split("_")
        for clave, fases in salidas_de_celda(celda).items():
            comunes = [s for s in fases.index if s in theta_true.index]
            if len(comunes) < 3:
                continue
            m = evaluate(theta_true.loc[comunes].to_numpy(),
                         fases.loc[comunes].to_numpy())
            filas.append(dict(
                celda=celda.name.rstrip("/"), salida=SALIDAS[clave],
                n_genes=int(partes[0][1:]), n_samples=int(partes[1][1:]),
                n_evaluadas=len(comunes),
                giro_h=_giro_a_la_verdad(theta_true.loc[comunes].to_numpy(),
                                         fases.loc[comunes].to_numpy()),
                **m))
    return pd.DataFrame(filas)


def entre_implementaciones(raiz: Path) -> pd.DataFrame:
    """
    Acuerdo entre dos salidas, sin usar la verdad.

    Responde a "dan lo mismo", que es distinto de "aciertan lo mismo". Se
    comparan sobre las muestras que ambas conservan.
    """
    pares = [("r_preliminar", "py_preliminar"),
             ("r_preliminar", "py_robusto")]
    filas = []
    for celda in sorted(Path(raiz).glob("*/")):
        s = salidas_de_celda(celda)
        for a, b in pares:
            if a not in s or b not in s:
                continue
            comunes = [x for x in s[a].index if x in s[b].index]
            if len(comunes) < 3:
                continue
            m = evaluate(s[a].loc[comunes].to_numpy(), s[b].loc[comunes].to_numpy())
            filas.append(dict(celda=celda.name.rstrip("/"),
                              a=SALIDAS[a], b=SALIDAS[b],
                              n_comunes=len(comunes), **m))
    return pd.DataFrame(filas)


def resumen(t: pd.DataFrame) -> pd.DataFrame:
    """Mediana por salida de las medidas que van al capitulo."""
    medidas = [c for c in ("error_fase_h", "error_orden_h", "rho_fl",
                           "p2", "giro_h") if c in t.columns]
    g = t.groupby("salida")
    salida = g[medidas].median()
    salida["n_celdas"] = g.size()
    return salida


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--in", dest="entrada", default="output/simulacion")
    args = ap.parse_args()

    entrada = Path(args.entrada)
    raiz = entrada / "rbench"
    tablas = entrada / "tablas"
    tablas.mkdir(parents=True, exist_ok=True)

    t = contra_la_verdad(raiz)
    if t.empty:
        print("No hay celdas con verdad y salidas volcadas. "
              "Corre antes: python -m simulation.rbench --reps 1 --run-r")
        return
    t.to_csv(tablas / "r_vs_python_precision.csv", index=False)
    print(f"  escrita r_vs_python_precision.csv ({len(t)} filas)\n")
    print(resumen(t).round(3).to_string())

    a = entre_implementaciones(raiz)
    if "giro_h" in t.columns:
        print("\n  giro respecto de la verdad, celda a celda (h):")
        print(t.pivot_table(index="celda", columns="salida",
                            values="giro_h").round(2).to_string())

    if not a.empty:
        a.to_csv(tablas / "r_vs_python_acuerdo.csv", index=False)
        print(f"\n  escrita r_vs_python_acuerdo.csv ({len(a)} filas)\n")
        print(a.groupby(["a", "b"])[["error_fase_h", "rho_fl"]]
               .median().round(3).to_string())


if __name__ == "__main__":
    main()
