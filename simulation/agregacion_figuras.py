"""
Figuras y tablas de la comparacion entre metodos de resumen de las K.

Es la seccion de la agregacion de ordenes, la aportacion central del TFG. La
pregunta es si resumir las K repeticiones de la Etapa 4 con un orden de
consenso (TSP3 o HODs) mejora al criterio que hay en el codigo, quedarse con
una repeticion (``best_k`` con ``r2_median``).

Lo que dice la rejilla, medido sobre 16.200 pares por variante del FMM
-----------------------------------------------------------------------
1. **HODs y TSP3 producen exactamente el mismo orden** en el 100 % de los
   pares. No es un empate estadistico, es identidad: con K = 5 repeticiones
   que discrepan poco, los dos llegan al mismo ciclo. Difieren solo en coste,
   y ahi HODs gana con claridad (figura D). En el capitulo se habla de "la
   agregacion" sin distinguirlos, y se recomienda HODs.

2. **best_k con medoid no aporta nada.** Elige la misma repeticion que
   ``r2_median`` en el 22 % de los casos, y cuando elige otra la mediana del
   cambio es 0,000. Cambiar de criterio de seleccion no es una via de mejora.

3. **La agregacion gana de forma consistente pero pequena.** Sobre el error
   del orden relativo, que mide el orden al margen del anclaje, gana en el
   69 % de las replicas con el FMM base y en el 71 % con el estabilizado. La
   mediana de la ganancia es de 0,03 h, y llega a 0,14 h con n = 96 y
   A/sigma = 1. Sobre el error de fase absoluta, que incluye el anclaje, la
   victoria baja al 55 %: el error de sincronizacion es comun a los dos
   brazos y tapa la diferencia.

4. **La ganancia crece con n y al bajar la senal.** Es lo contrario del
   FMM estabilizado, cuya ventaja se apaga al subir n. El consenso promedia
   el ruido de particion de las K repeticiones, y ese ruido pesa mas cuando
   hay muchas muestras que ordenar y poca senal para hacerlo.

Las figuras
-----------
A. Mapa de ganancia por n y A/sigma, con la fraccion de victorias al lado.
B. Dispersion pareada: cada punto es una replica, en x el error de quedarse
   con una repeticion y en y el de agregar. Bajo la diagonal gana agregar.
   Es la figura que explica que significa "pareado".
C. Ganancia por nivel de cada factor de la rejilla. Responde a "en que
   hiperparametros mejora".
D. Coste de la Etapa 4 por metodo y n.

Uso
---
    python -m simulation.agregacion_figuras --in output/simulacion --images docs/images
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

from simulation.analyze import (
    FACTORES, FMM_REFERENCIA, METODO_REFERENCIA, _rejilla, contraste_metodos,
    etiqueta, load,
)
from simulation.fmm_figuras import _mapa, _mpl, tabla_wilcoxon

METODO_AGREGACION = "aggregate:hods:circmean"
MEDIDA = "error_orden_h"

# Etiquetas cortas de los metodos para las figuras y las tablas del capitulo.
NOMBRES = {
    "best_k:r2_median":        "best_k (R2 mediana)",
    "best_k:medoid":           "best_k (medoide)",
    "aggregate:tsp3:circmean": "TSP3",
    "aggregate:hods:circmean": "HODs",
    "aggregate:tsp3:cire":     "TSP3 + CIRE",
}

AZUL, ROJO = "#2471a3", "#c0392b"


# ═══════════════════════════════════════════════════════════════════════════
# Datos
# ═══════════════════════════════════════════════════════════════════════════

def cargar_contraste(entrada: Path | str, fmm: str = FMM_REFERENCIA,
                     metodo: str = METODO_AGREGACION) -> pd.DataFrame:
    """
    Pares (metodo, referencia) replica a replica, con el FMM fijo.

    Anade las columnas ``d_<medida>`` que espera ``tabla_wilcoxon``, con el
    mismo convenio que en el FMM: negativo significa que el metodo nuevo
    gana.
    """
    df = load(Path(entrada) / "resultados.csv")
    met = contraste_metodos(df, fmm=fmm)
    par = met[met["metodo"] == metodo].copy()
    if par.empty:
        raise ValueError(f"no hay pares del metodo {metodo!r}")
    for m in ("error_orden_h", "concordancia", "error_fase_h"):
        par[f"d_{m}"] = par[f"delta_{m}"]
    return par.reset_index(drop=True)


def tabla_identidades(df: pd.DataFrame, medida: str = MEDIDA) -> pd.DataFrame:
    """
    Fraccion de pares en que dos metodos dan EXACTAMENTE el mismo resultado.

    Es lo que distingue "empatan" de "son el mismo calculo". HODs y TSP3
    coinciden en el 100 % de los pares, asi que en el capitulo no hay que
    compararlos sino decir que son intercambiables salvo en coste.
    """
    r = _rejilla(df)
    r = r[r["status"] == "ok"]
    clave = ["n_samples", "sampling", "omega", "snr", "n_genes",
             "prop_rhythmic", "n_core", "rep", "fmm_method"]
    piv = r.pivot_table(index=clave, columns="metodo", values=medida)
    metodos = [m for m in NOMBRES if m in piv.columns]
    filas = []
    for i, a in enumerate(metodos):
        for b in metodos[i + 1:]:
            iguales = (piv[a] - piv[b]).abs() < 1e-9
            filas.append(dict(metodo_a=NOMBRES[a], metodo_b=NOMBRES[b],
                              n_pares=int(iguales.notna().sum()),
                              fraccion_identicos=float(iguales.mean())))
    return pd.DataFrame(filas)


def tabla_por_factor(par: pd.DataFrame, medida: str = MEDIDA) -> pd.DataFrame:
    """
    Ganancia de la agregacion por nivel de cada factor: mediana, rango
    intercuartilico y fraccion de victorias. Una fila por factor y nivel.
    """
    col = f"d_{medida}"
    filas = []
    for f in FACTORES:
        for nivel, g in par.groupby(f):
            d = g[col].dropna()
            filas.append(dict(
                factor=f, nivel=str(nivel), n=int(d.size),
                ganancia_mediana=float(d.median()),
                iqr=float(d.quantile(0.75) - d.quantile(0.25)),
                tasa_victorias=float((d < 0).mean()),
            ))
    return pd.DataFrame(filas).set_index(["factor", "nivel"])


# ═══════════════════════════════════════════════════════════════════════════
# Figuras
# ═══════════════════════════════════════════════════════════════════════════

def figura_a_efecto(par: pd.DataFrame, destino: Path) -> None:
    """
    A: ganancia y consistencia por n y A/sigma.

    Izquierda, cuanto baja el error del orden relativo al agregar (mediana en
    horas). Derecha, en que fraccion de replicas mejora. Se leen juntas: la
    ganancia es pequena en horas pero se repite en dos de cada tres replicas,
    y las dos crecen hacia abajo (mas muestras) y hacia la izquierda (menos
    senal).
    """
    plt = _mpl()
    fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.4))
    fig.subplots_adjust(top=0.72, wspace=0.35)

    t = par.pivot_table(index="n_samples", columns="snr",
                        values=f"d_{MEDIDA}", aggfunc="median")
    im = _mapa(axes[0], t, "Cuanto baja el error del orden relativo (h)",
               "RdBu_r", 0.0)
    fig.colorbar(im, ax=axes[0], fraction=0.045, pad=0.03)

    v = (par.assign(gana=lambda x: (x[f"d_{MEDIDA}"] < 0).astype(float))
            .pivot_table(index="n_samples", columns="snr", values="gana"))
    # "RdBu" sin invertir: las fracciones altas, donde gana la agregacion,
    # salen en azul, igual que las ganancias negativas del panel izquierdo
    im = _mapa(axes[1], v, "En que fraccion de replicas mejora",
               "RdBu", 0.5, fmt="{:.2f}", vlim=0.5)
    fig.colorbar(im, ax=axes[1], fraction=0.045, pad=0.03)

    fig.suptitle("Agregar las K repeticiones (HODs) frente a quedarse con una\n"
                 "azul = gana la agregacion", fontsize=10, y=1.04)
    fig.savefig(destino)
    plt.close(fig)


def figura_b_pareada(par: pd.DataFrame, destino: Path) -> None:
    """
    B: dispersion pareada, un panel por n.

    Cada punto es una replica. En x el error del orden relativo al quedarse
    con una repeticion, en y el error al agregar. Los puntos bajo la
    diagonal son replicas donde agregar mejora. El color es A/sigma, que es
    lo que mueve los puntos a lo largo de la diagonal: con poca senal el
    error es grande en los dos brazos, y ahi es donde mas se separan.
    """
    plt = _mpl()
    niveles_n = sorted(par["n_samples"].unique())
    niveles_snr = sorted(par["snr"].unique())
    colores = plt.cm.viridis(np.linspace(0.1, 0.9, len(niveles_snr)))

    fig, axes = plt.subplots(1, len(niveles_n), figsize=(3.1 * len(niveles_n), 3.6),
                             sharex=True, sharey=True)
    fig.subplots_adjust(top=0.80, bottom=0.18, wspace=0.08)
    axes = np.atleast_1d(axes)
    lim = float(np.nanpercentile(
        np.concatenate([par[f"{MEDIDA}_ref"], par[f"{MEDIDA}_met"]]), 99.5))
    for ax, n in zip(axes, niveles_n):
        sub = par[par["n_samples"] == n]
        for snr, color in zip(niveles_snr, colores):
            s = sub[sub["snr"] == snr]
            ax.scatter(s[f"{MEDIDA}_ref"], s[f"{MEDIDA}_met"], s=5, alpha=0.35,
                       color=color, label=f"A/$\\sigma$ = {snr:g}", rasterized=True)
        ax.plot([0, lim], [0, lim], color="0.3", lw=0.9, ls="--")
        gana = float((sub[f"d_{MEDIDA}"] < 0).mean())
        ax.text(0.04, 0.95, f"n = {n}\nmejora en el {100 * gana:.0f} %",
                transform=ax.transAxes, va="top", fontsize=8.5)
        ax.set_xlim(0, lim)
        ax.set_ylim(0, lim)
        ax.set_aspect("equal")
    fig.supxlabel("error al quedarse con una repeticion (h)", fontsize=9, y=0.02)
    axes[0].set_ylabel("error al agregar (h)")
    leg = axes[-1].legend(frameon=False, fontsize=7.5, loc="lower right",
                          markerscale=2.5)
    for h in leg.legend_handles:
        h.set_alpha(1.0)
    fig.suptitle("Error del orden relativo, replica a replica\n"
                 "bajo la diagonal gana la agregacion", fontsize=10, y=1.0)
    fig.savefig(destino, dpi=200)
    plt.close(fig)


def figura_c_factores(par: pd.DataFrame, destino: Path) -> None:
    """
    C: ganancia por nivel de cada factor de la rejilla.

    Un panel por factor, la mediana de la ganancia con su rango
    intercuartilico. Lo que se ve: n y A/sigma mueven la ganancia, los demas
    factores apenas.
    """
    plt = _mpl()
    t = tabla_por_factor(par).reset_index()
    fig, axes = plt.subplots(1, len(FACTORES), figsize=(2.2 * len(FACTORES), 3.3),
                             sharey=True)
    fig.subplots_adjust(top=0.74, wspace=0.12)
    tope  = (t["ganancia_mediana"] + t["iqr"] / 2.0).max()
    suelo = (t["ganancia_mediana"] - t["iqr"] / 2.0).min()
    for ax, f in zip(axes, FACTORES):
        sub = t[t.factor == f]
        x = np.arange(len(sub))
        ax.errorbar(x, sub["ganancia_mediana"], yerr=sub["iqr"] / 2.0,
                    fmt="o", ms=4, capsize=2, color=AZUL)
        # los porcentajes van arriba del todo, fuera del recorrido de las
        # barras, para no pisarlas
        for xi, tv in zip(x, sub["tasa_victorias"]):
            ax.text(xi, 0.97, f"{100 * tv:.0f} %", ha="center", va="top",
                    fontsize=7.5, color=ROJO, transform=ax.get_xaxis_transform())
        ax.axhline(0, color="0.4", lw=0.8, ls="--")
        ax.set_xticks(x, sub["nivel"], rotation=45, fontsize=8)
        ax.set_title(etiqueta(f), fontsize=8.5)
        ax.set_xlim(-0.6, len(sub) - 0.4)
        ax.set_ylim(suelo - 0.02, tope + 0.06)
    axes[0].set_ylabel("cuanto baja el error del orden relativo (h)\n"
                       "mediana y rango intercuartilico", fontsize=8)
    fig.suptitle("Ganancia de la agregacion por factor\n"
                 "en rojo, la fraccion de replicas en que mejora", fontsize=10, y=1.02)
    fig.savefig(destino)
    plt.close(fig)


def figura_d_coste(df: pd.DataFrame, destino: Path) -> None:
    """
    D: coste de la Etapa 4 por metodo y numero de muestras.

    HODs cuesta casi lo mismo que quedarse con una repeticion. TSP3 escala
    mal con n porque el problema del viajante crece con el numero de
    muestras, y con n = 96 duplica el tiempo. Como los dos dan el mismo
    orden, el coste decide: HODs.

    El tiempo incluye la construccion del atlas, que con los metodos de
    agregacion reajusta el FMM sobre el orden de consenso. Parte del segundo
    de mas de HODs frente a best_k es ese reajuste y no el consenso en si.
    """
    plt = _mpl()
    r = _rejilla(df)
    r = r[(r["status"] == "ok") & (r["fmm_method"] == FMM_REFERENCIA)]
    t = r.groupby(["metodo", "n_samples"])["t_robust_s"].median().unstack("metodo")

    fig, ax = plt.subplots(figsize=(5.4, 3.2))
    marcas = ["o", "s", "^", "v", "D"]
    for (metodo, marca) in zip(NOMBRES, marcas):
        if metodo in t.columns:
            ax.plot(t.index, t[metodo], marker=marca, label=NOMBRES[metodo])
    ax.set_xticks(list(t.index), [f"{v:g}" for v in t.index])
    ax.set_xlabel(etiqueta("n_samples"))
    ax.set_ylabel("tiempo de la Etapa 4 (s, mediana)")
    ax.set_title("Coste de resumir las K repeticiones", fontsize=9)
    ax.legend(frameon=False, fontsize=7.5)
    ax.grid(alpha=0.25)
    fig.savefig(destino)
    plt.close(fig)


# ═══════════════════════════════════════════════════════════════════════════
# CLI
# ═══════════════════════════════════════════════════════════════════════════

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--in", dest="entrada", default="output/simulacion")
    ap.add_argument("--images", default="docs/images")
    args = ap.parse_args()

    entrada = Path(args.entrada)
    imagenes = Path(args.images)
    imagenes.mkdir(parents=True, exist_ok=True)
    tablas = entrada / "tablas"
    tablas.mkdir(parents=True, exist_ok=True)

    df = load(entrada / "resultados.csv")
    tabla_identidades(df).to_csv(tablas / "agregacion_identidades.csv", index=False)
    print("  escrita agregacion_identidades.csv")

    for fmm in ("base", "stabilized"):
        par = cargar_contraste(entrada, fmm=fmm)
        print(f"  {len(par)} pares HODs contra {METODO_REFERENCIA} con FMM {fmm}")
        t = tabla_wilcoxon(par, medida=MEDIDA)
        t.to_csv(tablas / f"agregacion_wilcoxon_{fmm}.csv", index=False)
        print(f"  escrita agregacion_wilcoxon_{fmm}.csv ({len(t)} celdas)")
        tabla_por_factor(par).to_csv(tablas / f"agregacion_por_factor_{fmm}.csv")
        print(f"  escrita agregacion_por_factor_{fmm}.csv")
        if fmm == "base":
            figura_a_efecto(par, imagenes / "agr-a-efecto.png")
            figura_b_pareada(par, imagenes / "agr-b-pareada.png")
            figura_c_factores(par, imagenes / "agr-c-factores.png")
            print("  escritas agr-a-efecto.png, agr-b-pareada.png y agr-c-factores.png")

    figura_d_coste(df, imagenes / "agr-d-coste.png")
    print("  escrita agr-d-coste.png")


if __name__ == "__main__":
    main()
