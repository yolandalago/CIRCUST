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
   y ahi HODs gana con claridad (ver figura C). En el capitulo se habla de
   "la agregacion" sin distinguirlos, y se recomienda HODs.

2. **best_k con medoid no aporta nada.** Elige la misma repeticion que
   ``r2_median`` en el 22 % de los casos, y cuando elige otra la mediana del
   delta es 0,000. Cambiar de criterio de seleccion no es una via de mejora.

3. **La agregacion gana de forma consistente pero pequena.** Sobre el error
   del oraculo, que mide el orden al margen del anclaje, gana en el 69 % de
   las replicas con el FMM base y en el 71 % con el estabilizado. La mediana
   de la ganancia es de 0,03 h, y llega a 0,14 h con n = 96 y A/sigma = 1.
   Sobre ``medae_h``, que incluye el anclaje, la victoria baja al 55 %: el
   error de sincronizacion es comun a los dos brazos y tapa la diferencia.

4. **La ganancia crece con n y al bajar la senal.** Es lo contrario del
   FMM estabilizado, cuya ventaja se apaga al subir n. El consenso promedia
   el ruido de particion de las K repeticiones, y ese ruido pesa mas cuando
   hay muchas muestras que ordenar y poca senal para hacerlo.

5. **No depende del desacuerdo entre las K tanto como se esperaba.** La
   hipotesis era que solo puede haber ganancia donde las K discrepan. La
   ganancia existe ya en el tramo de menor desacuerdo (73 % de victorias) y
   crece hasta la mitad de la escala, donde se estabiliza. La agregacion
   aprovecha discrepancias pequenas, del orden de 0,005 en la escala 0 a 1
   de ``desacuerdo_K``.

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
    FMM_REFERENCIA, METODO_REFERENCIA, _rejilla, contraste_metodos, load,
)
from simulation.fmm_figuras import _mapa, _mpl, tabla_wilcoxon

METODO_AGREGACION = "aggregate:hods:circmean"

# Etiquetas cortas para las figuras y las tablas del capitulo.
NOMBRES = {
    "best_k:r2_median":        "best_k (R2 mediana)",
    "best_k:medoid":           "best_k (medoide)",
    "aggregate:tsp3:circmean": "TSP3",
    "aggregate:hods:circmean": "HODs",
    "aggregate:tsp3:cire":     "TSP3 + CIRE",
}


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
    for m in ("mae_h_oracle", "abs_tau_circ", "medae_h"):
        par[f"d_{m}"] = par[f"delta_{m}"]
    return par.reset_index(drop=True)


def tabla_identidades(df: pd.DataFrame, medida: str = "mae_h_oracle") -> pd.DataFrame:
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


def tabla_por_desacuerdo(par: pd.DataFrame, n_tramos: int = 6) -> pd.DataFrame:
    """
    Ganancia de la agregacion por tramo de desacuerdo entre las K.

    ``desacuerdo_K`` es observable sin la verdad, asi que esta tabla dice
    tambien cuando merece la pena agregar en un tejido real.
    """
    p = par.copy()
    p["tramo"] = pd.qcut(p["desacuerdo_K_ref"], n_tramos, duplicates="drop")
    g = p.groupby("tramo", observed=True)
    return pd.DataFrame({
        "n":                  g.size(),
        "desacuerdo_mediana": g["desacuerdo_K_ref"].median(),
        "oracle_referencia":  g["mae_h_oracle_ref"].median(),
        "delta_mediana":      g["d_mae_h_oracle"].median(),
        "tasa_victorias":     g["d_mae_h_oracle"].apply(lambda x: float((x < 0).mean())),
    })


# ═══════════════════════════════════════════════════════════════════════════
# Figuras
# ═══════════════════════════════════════════════════════════════════════════

def figura_a_efecto(par: pd.DataFrame, destino: Path) -> None:
    """
    A: efecto pareado y consistencia, por n y A/sigma.

    Izquierda la mediana de la ganancia en horas del oraculo, derecha la
    fraccion de replicas en que gana la agregacion. Se leen juntas: la
    ganancia es pequena en horas pero se repite en dos de cada tres replicas,
    y las dos crecen hacia abajo (mas muestras) y hacia la izquierda (menos
    senal).
    """
    plt = _mpl()
    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.4))
    fig.subplots_adjust(top=0.72, wspace=0.35)

    t = par.pivot_table(index="n_samples", columns="snr",
                        values="d_mae_h_oracle", aggfunc="median")
    im = _mapa(axes[0], t, "Ganancia mediana (h del oraculo)", "RdBu_r", 0.0)
    fig.colorbar(im, ax=axes[0], fraction=0.045, pad=0.03)

    v = (par.assign(gana=lambda x: (x.d_mae_h_oracle < 0).astype(float))
            .pivot_table(index="n_samples", columns="snr", values="gana"))
    # "RdBu" sin invertir: las fracciones altas, donde gana la agregacion,
    # salen en azul, igual que las ganancias negativas del panel izquierdo
    im = _mapa(axes[1], v, "Fraccion de replicas en que gana", "RdBu", 0.5,
               fmt="{:.2f}", vlim=0.5)
    fig.colorbar(im, ax=axes[1], fraction=0.045, pad=0.03)

    fig.suptitle("Agregacion (HODs) frente a quedarse con una repeticion\n"
                 "azul = gana la agregacion", fontsize=10, y=1.04)
    fig.savefig(destino)
    plt.close(fig)


def figura_b_desacuerdo(par: pd.DataFrame, destino: Path) -> None:
    """
    B: la ganancia frente al desacuerdo entre las K repeticiones.

    Es la figura del mecanismo. La agregacion solo puede cambiar el orden
    donde las K discrepan, y ``desacuerdo_K`` mide cuanto discrepan. La
    ganancia aparece ya con desacuerdos minusculos y se estabiliza a mitad
    de escala: no hace falta que las K se contradigan mucho para que el
    consenso mejore a la mejor de ellas.
    """
    plt = _mpl()
    t = tabla_por_desacuerdo(par, n_tramos=8)
    x = t["desacuerdo_mediana"].to_numpy()

    fig, ax1 = plt.subplots(figsize=(5.6, 3.3))
    ax1.plot(x, t["delta_mediana"], marker="o", color="#2471a3",
             label="ganancia mediana (h)")
    ax1.axhline(0, color="0.4", lw=0.8, ls="--")
    ax1.set_xlabel("desacuerdo entre las K repeticiones (mediana del tramo)")
    ax1.set_ylabel("ganancia en el oraculo (h)\nnegativo = gana la agregacion",
                   color="#2471a3")
    ax1.set_xscale("log")

    ax2 = ax1.twinx()
    ax2.plot(x, t["tasa_victorias"], marker="s", color="#c0392b",
             label="fraccion de victorias")
    ax2.set_ylim(0.4, 1.0)
    ax2.set_ylabel("fraccion de replicas en que gana", color="#c0392b")

    ax1.set_title("La ganancia aparece con desacuerdos pequenos", fontsize=9)
    fig.savefig(destino)
    plt.close(fig)


def figura_c_coste(df: pd.DataFrame, destino: Path) -> None:
    """
    C: coste de la Etapa 4 por metodo y numero de muestras.

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
    ax.set_xlabel("n muestras")
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
        t = tabla_wilcoxon(par, medida="mae_h_oracle")
        t.to_csv(tablas / f"agregacion_wilcoxon_{fmm}.csv", index=False)
        print(f"  escrita agregacion_wilcoxon_{fmm}.csv ({len(t)} celdas)")
        tabla_por_desacuerdo(par).to_csv(tablas / f"agregacion_desacuerdo_{fmm}.csv")
        print(f"  escrita agregacion_desacuerdo_{fmm}.csv")
        if fmm == "base":
            figura_a_efecto(par, imagenes / "agr-a-efecto.png")
            figura_b_desacuerdo(par, imagenes / "agr-b-desacuerdo.png")
            print("  escritas agr-a-efecto.png y agr-b-desacuerdo.png")

    figura_c_coste(df, imagenes / "agr-c-coste.png")
    print("  escrita agr-c-coste.png")


if __name__ == "__main__":
    main()
