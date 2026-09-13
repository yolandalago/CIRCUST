"""
Tablas y figuras del capitulo de simulaciones.

Lee el CSV largo que produce ``run.py`` y produce las tablas que van a la
memoria.

Tres decisiones de agregacion, todas tomadas midiendo y no por costumbre
-----------------------------------------------------------------------
1. **Las medidas continuas se agregan con MEDIANA, no con media.** ``rho_circ``
   cambia de signo en torno al 1-2 % de las replicas aunque la reconstruccion
   sea perfecta, porque su direccion media queda mal determinada cuando las
   muestras cubren todo el ciclo. Medido en el bloque sin ruido: 1 de 54
   reconstrucciones con orden exacto dio ``rho = -0,093``. Una sola replica
   asi arrastra la media de una celda de 20.

2. **La orientacion se agrega como TASA de acierto**, no como mediana, porque
   es una proporcion. Y va en columna propia porque ``medae_h`` mezcla dos
   cosas: si el sentido de giro sale invertido da 5 a 7 h aunque el orden sea
   perfecto, con lo que en las celdas donde la orientacion no es decidible esa
   metrica mide suerte y no reconstruccion.

3. **Las filas con fallo se excluyen de las medianas pero su proporcion se
   reporta** como columna. La tasa de fallo por celda es un resultado del
   estudio, es la frontera de aplicabilidad del metodo.

Uso
---
    python -m simulation.analyze --in output/simulacion --out output/simulacion
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

# Medidas continuas, se agregan con mediana.
MEDIDAS = ["medae_h", "mae_h", "mae_h_oracle", "tau_circ", "rho_circ", "tau",
           "p2", "r2_anchor", "evidencia", "desacuerdo_K", "n_top"]

# Proporciones, se agregan con media porque son tasas.
TASAS = ["orientacion_ok"]


def _con_nivel(tabla: pd.DataFrame, columnas: list[str],
               nivel: str) -> pd.DataFrame:
    """Convierte las columnas en un MultiIndex ``(medida, nivel)``."""
    salida = tabla.copy()
    salida.columns = pd.MultiIndex.from_product([columnas, [nivel]])
    return salida


def load(path: Path | str) -> pd.DataFrame:
    """Carga el CSV largo de resultados."""
    return pd.read_csv(path)


def _resumen(df: pd.DataFrame, por: list[str]) -> pd.DataFrame:
    """
    Mediana y dispersion de las medidas por grupo, tasas, y tasa de fallo.

    Las filas con fallo no entran en las medianas pero si en el denominador
    de la tasa de fallo.
    """
    fallo = df.groupby(por)["status"].apply(
        lambda s: float((s == "fallo").mean())
    )
    ok = df[df["status"] == "ok"]
    g  = ok.groupby(por)

    medidas = [m for m in MEDIDAS if m in ok.columns]
    tabla = _con_nivel(g[medidas].median(), medidas, "mediana")
    iqr   = _con_nivel(g[medidas].quantile(0.75) - g[medidas].quantile(0.25),
                       medidas, "iqr")
    tabla = pd.concat([tabla, iqr], axis=1)

    for t in (t for t in TASAS if t in ok.columns):
        tabla[(t, "tasa")] = g[t].mean()
    tabla[("tasa_fallo", "")] = fallo
    return tabla.sort_index(axis=1)


def tabla_sin_ruido(df: pd.DataFrame) -> pd.DataFrame:
    """
    Verificacion de que el pipeline funciona con datos perfectos.

    Cualquier error aqui no puede venir del ruido, asi que aisla los sesgos
    del propio algoritmo. Se agrupa por forma de onda y muestreo, que son los
    factores que discriminan.
    """
    s = df[df["bloque"] == "sin_ruido"]
    if s.empty:
        raise ValueError("no hay filas del bloque sin_ruido")
    por = [c for c in ("omega", "sampling", "n_samples") if c in s.columns]
    return _resumen(s, por[:2] if len(por) > 1 else por)


def tabla_principal(df: pd.DataFrame) -> pd.DataFrame:
    """Bloque P: reconstruccion por muestreo, morfologia y senal-ruido."""
    p = df[df["bloque"] == "principal"]
    if p.empty:
        raise ValueError("no hay filas del bloque principal")
    return _resumen(p, ["sampling", "omega", "snr"])


def tabla_regularizacion(df: pd.DataFrame) -> pd.DataFrame:
    """
    Bloque F: diferencia **pareada** entre stabilized y base.

    La diferencia se calcula replica a replica sobre la matriz identica y solo
    despues se resume. Promediar primero y restar despues perderia toda la
    potencia del diseno pareado, que es la razon de que los dos bloques
    compartan semilla de generacion.
    """
    clave = ["n_samples", "sampling", "omega", "snr", "rep"]
    base = df[(df["bloque"] == "principal") & (df["status"] == "ok")]
    stab = df[(df["bloque"] == "regularizacion") & (df["status"] == "ok")]
    if base.empty or stab.empty:
        raise ValueError("faltan filas de principal o de regularizacion")
    par = base.merge(stab, on=clave, suffixes=("_base", "_stab"))
    if par.empty:
        raise ValueError("el emparejamiento P contra F no cuadra")

    columnas = []
    for m in MEDIDAS + TASAS:
        a, b = f"{m}_base", f"{m}_stab"
        if a in par.columns and b in par.columns:
            par[f"delta_{m}"] = par[b] - par[a]
            columnas.append(f"delta_{m}")

    g = par.groupby(["sampling", "omega", "snr"])
    tabla = _con_nivel(g[columnas].median(), columnas, "mediana")
    iqr   = _con_nivel(g[columnas].quantile(0.75) - g[columnas].quantile(0.25),
                       columnas, "iqr")
    return pd.concat([tabla, iqr], axis=1).sort_index(axis=1)


def tabla_agregacion(df: pd.DataFrame) -> pd.DataFrame:
    """
    Bloque G: metodos de resumen de las K, separados por cantidad de senal.

    ``prop_rhythmic`` es el eje que gobierna el desacuerdo entre las K
    repeticiones. Sin variarlo los metodos empatan por construccion, porque
    donde las K coinciden cualquier resumen devuelve el mismo orden.
    """
    g = df[df["bloque"] == "agregacion"]
    if g.empty:
        raise ValueError("no hay filas del bloque agregacion")
    return _resumen(g, ["prop_rhythmic", "snr", "robust_method",
                        "selection_criterion", "aggregation"])


def tabla_gcore(df: pd.DataFrame) -> pd.DataFrame:
    """Bloque G_core: cuanto conocimiento a priori hace falta."""
    g = df[df["bloque"] == "gcore"]
    if g.empty:
        raise ValueError("no hay filas del bloque gcore")
    return _resumen(g, ["n_core", "snr"])


def tabla_invariancia(df: pd.DataFrame) -> pd.DataFrame:
    """
    Bloque de invariancia: comprobacion de que el pipeline solo ve A/sigma.

    Los dos miembros de cada par tienen el mismo cociente y distinta escala
    absoluta, asi que deben dar metricas identicas.
    """
    inv = df[(df["bloque"] == "invariancia") & (df["status"] == "ok")]
    if inv.empty:
        raise ValueError("no hay filas del bloque invariancia")
    return inv.groupby(["snr", "sigma"])[MEDIDAS].median()


def tabla_tiempos(df_tiempos: pd.DataFrame) -> pd.DataFrame:
    """Bloque C: tiempo por implementacion y tamano de matriz."""
    agg   = df_tiempos.groupby(["implementacion", "n_genes", "n_samples"])
    ancho = agg["t_total_s"].median().unstack("implementacion")
    if {"python", "r"} <= set(ancho.columns):
        ancho["factor_r_sobre_python"] = ancho["r"] / ancho["python"]
    return ancho


def calibracion_orientacion(df: pd.DataFrame, indicador: str = "evidencia",
                            n_tramos: int = 8) -> pd.DataFrame:
    """
    Relacion entre un indicador **observable sin la verdad** y el acierto de
    la orientacion.

    Es la base del aviso de fase absoluta: en un tejido nuevo no se puede
    medir el error, pero si se puede medir el indicador, y esta tabla dice
    que esperar para cada valor. El umbral no se fija aqui, se elige mirando
    la curva con un criterio explicito de precision y cobertura.
    """
    s = df[(df["status"] == "ok") & df[indicador].notna()].copy()
    if s.empty:
        raise ValueError(f"no hay filas con {indicador} definido")
    s["tramo"] = pd.qcut(s[indicador], n_tramos, duplicates="drop")
    g = s.groupby("tramo", observed=True)
    return pd.DataFrame({
        "n":                  g.size(),
        f"{indicador}_medio": g[indicador].median(),
        "acierto_orientacion": g["orientacion_ok"].mean(),
        "medae_h":            g["medae_h"].median(),
    })


def _figura(datos, x, y, hue, titulo, xlabel, ylabel, out: Path,
            logx: bool = True) -> None:
    """Grafico de lineas con una serie por nivel de ``hue``."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # Sin esto se escribian PNG vacios cuando el bloque todavia no se habia
    # ejecutado, y un fichero vacio parece una figura valida.
    if len(datos) == 0:
        raise ValueError(f"no hay datos para la figura {out.name}")

    fig, ax = plt.subplots(figsize=(7.0, 4.2))
    for nivel, sub in datos.groupby(hue):
        sub = sub.sort_values(x)
        ax.plot(sub[x], sub[y], marker="o", label=f"{hue} = {nivel}")
    if logx:
        ax.set_xscale("log", base=2)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(titulo)
    ax.legend(frameon=False)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=200)
    plt.close(fig)


def figura_snr(df: pd.DataFrame, out: Path) -> None:
    """Degradacion del error con la senal-ruido, una serie por muestreo."""
    p = df[(df["bloque"] == "principal") & (df["status"] == "ok")]
    datos = p.groupby(["sampling", "snr"])["medae_h"].median().reset_index()
    _figura(datos, "snr", "medae_h", "sampling",
            "Reconstruccion frente a la razon senal-ruido",
            "A / sigma", "Mediana del error absoluto (h)", out)


def figura_orientacion(df: pd.DataFrame, out: Path) -> None:
    """Tasa de acierto de orientacion, que va aparte del error."""
    p = df[(df["bloque"] == "principal") & (df["status"] == "ok")]
    datos = p.groupby(["omega", "snr"])["orientacion_ok"].mean().reset_index()
    _figura(datos, "snr", "orientacion_ok", "omega",
            "Acierto del sentido de giro",
            "A / sigma", "Proporcion de celdas bien orientadas", out)


def figura_regularizacion(df: pd.DataFrame, out: Path) -> None:
    """Ganancia pareada del FMM estabilizado, una serie por omega."""
    t = tabla_regularizacion(df)[("delta_medae_h", "mediana")].reset_index()
    t = t.rename(columns={t.columns[-1]: "delta_medae_h"})
    datos = t.groupby(["omega", "snr"])["delta_medae_h"].median().reset_index()
    _figura(datos, "snr", "delta_medae_h", "omega",
            "Ganancia del FMM estabilizado (negativo es mejor)",
            "A / sigma", "Diferencia pareada de mediana del error (h)", out)


def figura_agregacion(df: pd.DataFrame, out: Path) -> None:
    """Comparacion de metodos de resumen, una serie por metodo."""
    g = df[(df["bloque"] == "agregacion") & (df["status"] == "ok")].copy()
    g["metodo"] = np.where(
        g["robust_method"] == "best_k",
        "best_k:" + g["selection_criterion"].astype(str),
        "aggregate:" + g["aggregation"].astype(str),
    )
    datos = (g.groupby(["metodo", "prop_rhythmic"])["medae_h"]
               .median().reset_index()
               .rename(columns={"prop_rhythmic": "prop"}))
    _figura(datos, "prop", "medae_h", "metodo",
            "Metodos de resumen de las K repeticiones",
            "Proporcion de genes ritmicos", "Mediana del error absoluto (h)", out)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Tablas y figuras del capitulo de simulaciones.",
    )
    parser.add_argument("--in", dest="entrada", type=str,
                        default="output/simulacion")
    parser.add_argument("--out", dest="salida", type=str,
                        default="output/simulacion")
    parser.add_argument("--images", type=str, default="docs/images")
    args = parser.parse_args()

    entrada = Path(args.entrada)
    df = load(entrada / "resultados.csv")

    dir_tablas = Path(args.salida) / "tablas"
    dir_tablas.mkdir(parents=True, exist_ok=True)

    constructores = {
        "sin_ruido":      tabla_sin_ruido,
        "principal":      tabla_principal,
        "regularizacion": tabla_regularizacion,
        "agregacion":     tabla_agregacion,
        "gcore":          tabla_gcore,
        "invariancia":    tabla_invariancia,
    }
    for nombre, fn in constructores.items():
        try:
            fn(df).to_csv(dir_tablas / f"tabla_{nombre}.csv")
            print(f"  escrita tabla_{nombre}.csv")
        except (KeyError, ValueError) as exc:
            print(f"  tabla_{nombre} omitida: {exc}")

    for indicador in ("evidencia", "r2_anchor"):
        try:
            calibracion_orientacion(df, indicador).to_csv(
                dir_tablas / f"calibracion_{indicador}.csv")
            print(f"  escrita calibracion_{indicador}.csv")
        except (KeyError, ValueError) as exc:
            print(f"  calibracion_{indicador} omitida: {exc}")

    ruta_tiempos = entrada / "tiempos.csv"
    if ruta_tiempos.exists():
        tabla_tiempos(pd.read_csv(ruta_tiempos)).to_csv(
            dir_tablas / "tabla_tiempos.csv")
        print("  escrita tabla_tiempos.csv")

    imagenes = Path(args.images)
    for nombre, fn in (
        ("sim-snr",            figura_snr),
        ("sim-orientacion",    figura_orientacion),
        ("sim-regularizacion", figura_regularizacion),
        ("sim-agregacion",     figura_agregacion),
    ):
        try:
            fn(df, imagenes / f"{nombre}.png")
            print(f"  escrita {nombre}.png")
        except (KeyError, ValueError, IndexError) as exc:
            print(f"  {nombre}.png omitida: {exc}")


if __name__ == "__main__":
    main()
