#!/usr/bin/env python3
"""
simulation/report.py
====================
Analisis de los tres bloques y figuras para la memoria.

Separa a proposito dos cosas que suelen ir mezcladas. Las **tablas** resumen lo
medido y llevan su contraste pareado, de modo que una diferencia pequena no se
presenta como un hallazgo. Las **figuras** solo existen donde una tabla no
transmite la forma del resultado, que aqui ocurre en tres sitios: el cruce del
sesgo de omega con el tamano muestral, el reparto de ganancias y perdidas de la
regularizacion, y la frontera de aplicabilidad.

Diseno pareado
--------------
Todas las comparaciones se hacen **dentro de la misma replica**, que comparte
datos, semilla y conjunto TOP. Comparar medianas entre grupos mezclaria la
diferencia entre metodos con la variabilidad entre celdas, que es mucho mayor.
El contraste es un Wilcoxon de rangos con signo, sin supuesto de normalidad
sobre unas diferencias que no tienen por que serlo.

Uso
---
    python -m simulation.report R
    python -m simulation.report A
    python -m simulation.report todo --figuras docs/images
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

from simulation.grid import RAIZ_POR_DEFECTO

# ===========================================================================
# Paleta
# ===========================================================================
#
# Validada con el comprobador de la guia de visualizacion: banda de luminosidad,
# suelo de croma, separacion bajo daltonismo (peor par delta E 24,7 en protan) y
# contraste sobre el fondo. Las figuras van a una memoria impresa, de modo que
# ademas de color cada serie lleva marcador propio y las barras llevan su valor
# escrito, para que nada dependa solo del color.
AZUL     = "#2a78d6"   # categorico 1, y polo "gana el estabilizado"
NARANJA  = "#eb6834"   # categorico 2
ROJO     = "#e34948"   # polo "gana el base"
GRIS     = "#c9c8c3"   # contexto, series no destacadas
FONDO    = "#fcfcfb"
TINTA    = "#0b0b0b"
TINTA_2  = "#52514e"
RAMPA_AZUL = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95"]

ORDEN_METODOS = ["median", "best_k:r2_median", "best_k:n_well_fitted",
                 "best_k:medoid", "aggregate:tsp3", "aggregate:hods"]

# Cortes del desacuerdo entre las K repeticiones. El de 0,03 no es arbitrario:
# es el que separa el regimen de babuino, donde las repeticiones coinciden y
# agregar no puede aportar nada, del de GSE71620.
CORTES_DESACUERDO = [-1, 0.01, 0.03, 9]
ETIQUETAS_DESACUERDO = ["< 0,01", "0,01 a 0,03", "> 0,03"]


def _estilo():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        "pdf.fonttype": 42, "ps.fonttype": 42,
        "figure.facecolor": FONDO, "axes.facecolor": FONDO,
        "savefig.facecolor": FONDO,
        "font.size": 9, "axes.labelsize": 9, "axes.titlesize": 10,
        "axes.edgecolor": "#d8d7d2", "axes.linewidth": 0.8,
        "axes.labelcolor": TINTA_2, "text.color": TINTA,
        "xtick.color": TINTA_2, "ytick.color": TINTA_2,
        "xtick.labelsize": 8, "ytick.labelsize": 8,
        "grid.color": "#e8e7e3", "grid.linewidth": 0.6,
        "legend.frameon": False, "legend.fontsize": 8,
        "axes.spines.top": False, "axes.spines.right": False,
    })
    return plt


# ===========================================================================
# Lectura y utilidades
# ===========================================================================

def leer(bloque: str, raiz: Path | str = RAIZ_POR_DEFECTO) -> pd.DataFrame:
    ruta = Path(raiz) / "resultados" / f"{bloque}.csv"
    if not ruta.exists():
        raise SystemExit(f"No existe {ruta}. Lanza antes ese bloque.")
    return pd.read_csv(ruta, keep_default_na=False, na_values=[""])


def validas(d: pd.DataFrame) -> pd.DataFrame:
    """Filas de replicas que llegaron a terminar."""
    return d[d["error"].isna() | (d["error"] == "")]


def _clave(d: pd.DataFrame) -> list[str]:
    return ["celda", "replica"]


def emparejar(d: pd.DataFrame, columna: str, a: str, b: str,
              metricas: list[str]) -> pd.DataFrame:
    """Empareja dos niveles de ``columna`` dentro de cada replica.

    Devuelve una fila por replica con las metricas de ambos lados y su
    diferencia ``b - a``, de modo que un valor negativo significa que ``b``
    tiene menos error, o sea que gana.
    """
    v = validas(d)
    izq = v[v[columna] == a].set_index(_clave(v))
    der = v[v[columna] == b].set_index(_clave(v))
    j = izq.join(der[metricas], rsuffix="__b", how="inner")
    for m in metricas:
        j[f"dif_{m}"] = j[f"{m}__b"] - j[m]
    return j.reset_index()


def contraste(dif: pd.Series) -> dict:
    """Wilcoxon de rangos con signo sobre las diferencias pareadas."""
    x = pd.to_numeric(dif, errors="coerce").dropna()
    x = x[x != 0.0]
    if len(x) < 6:
        return {"n": int(len(x)), "p": np.nan}
    return {"n": int(len(x)), "p": float(wilcoxon(x).pvalue)}


def _resumen(j: pd.DataFrame, por: str, metrica: str) -> pd.DataFrame:
    """Mediana de la diferencia pareada y su contraste, por nivel de ``por``."""
    filas = []
    for nivel, sub in j.groupby(por, sort=True):
        c = contraste(sub[f"dif_{metrica}"])
        filas.append({por: nivel,
                      "mediana_dif": float(sub[f"dif_{metrica}"].median()),
                      "gana": "estabilizado" if sub[f"dif_{metrica}"].median() < 0
                              else "base",
                      **c})
    return pd.DataFrame(filas)


# ===========================================================================
# Bloque R: regularizacion del FMM
# ===========================================================================

METRICAS_R = ["verdad_omega_medae", "verdad_omega_sesgo",
              "verdad_picos_medae_h"]


def tabla_sesgo_omega(d: pd.DataFrame) -> pd.DataFrame:
    """Sesgo del apuntamiento estimado, por tamano muestral y modelo.

    El signo es lo que importa. Negativo significa subestimar omega, es decir
    estimar formas mas apuntadas que las reales, que es la patologia observada
    en babuino al ajustar con doce muestras.
    """
    v = validas(d)
    return v.pivot_table(index="n_samples", columns="fmm_method",
                         values=["verdad_omega_sesgo", "verdad_omega_medae",
                                 "verdad_picos_medae_h"],
                         aggfunc="median")


def tabla_bloque_r(d: pd.DataFrame,
                   metrica: str = "verdad_omega_medae") -> pd.DataFrame:
    """Diferencia pareada estabilizado menos base, por cada eje del bloque."""
    j = emparejar(d, "fmm_method", "base", "stabilized", METRICAS_R)
    partes = []
    for eje in ("n_samples", "waveform_mix", "sampling"):
        t = _resumen(j, eje, metrica)
        t.insert(0, "eje", eje)
        t = t.rename(columns={eje: "nivel"})
        partes.append(t)
    return pd.concat(partes, ignore_index=True)


# ===========================================================================
# Bloque A: orden representativo
# ===========================================================================

def con_estrato(d: pd.DataFrame) -> pd.DataFrame:
    """Anade el estrato de desacuerdo entre las K repeticiones.

    Es la variable mediadora del bloque: si las repeticiones coinciden no hay
    desacuerdo que resolver y ningun metodo de consenso puede aportar nada, de
    modo que promediar sobre todas las celdas enmascararia el efecto.
    """
    v = validas(d).copy()
    v["estrato"] = pd.cut(v["sd_rho_reps"], CORTES_DESACUERDO,
                          labels=ETIQUETAS_DESACUERDO)
    return v


def tabla_bloque_a(d: pd.DataFrame, metrica: str = "rho") -> pd.DataFrame:
    """Mediana de la metrica por metodo y estrato de desacuerdo."""
    v = con_estrato(d)
    t = v.pivot_table(index="metodo", columns="estrato", values=metrica,
                      aggfunc="median", observed=True)
    return t.reindex([m for m in ORDEN_METODOS if m in t.index])


def tabla_agregacion(d: pd.DataFrame, referencia: str = "best_k:r2_median",
                     metrica: str = "rho") -> pd.DataFrame:
    """Diferencia pareada de cada metodo frente al de referencia, por estrato.

    Aqui un valor **positivo** significa mejor, porque la metrica es rho y en
    rho mas es mejor.
    """
    v = con_estrato(d)
    base = v[v["metodo"] == referencia].set_index(_clave(v))
    filas = []
    for metodo in ORDEN_METODOS:
        if metodo == referencia or metodo not in set(v["metodo"]):
            continue
        otro = v[v["metodo"] == metodo].set_index(_clave(v))
        j = base[[metrica, "estrato"]].join(
            otro[[metrica]], rsuffix="__b", how="inner").dropna()
        j["dif"] = j[f"{metrica}__b"] - j[metrica]
        for estrato, sub in j.groupby("estrato", observed=True):
            filas.append({"metodo": metodo, "estrato": estrato,
                          "mediana_dif": float(sub["dif"].median()),
                          **contraste(sub["dif"])})
    return pd.DataFrame(filas)


def tabla_seleccion(d: pd.DataFrame) -> pd.DataFrame:
    """Con que frecuencia dos criterios eligen la misma repeticion.

    Si un criterio elige lo mismo que otro casi siempre, no es un criterio
    distinto. Y si elige distinto pero no rinde mejor, no aporta.
    """
    v = validas(d)
    w = v.pivot_table(index=_clave(v), columns="metodo",
                      values="rep_elegida", aggfunc="first")
    filas = []
    criterios = [c for c in ORDEN_METODOS if c.startswith("best_k")]
    for i, a in enumerate(criterios):
        for b in criterios[i + 1:]:
            if a not in w or b not in w:
                continue
            par = w[[a, b]].dropna()
            filas.append({"criterio_a": a, "criterio_b": b, "n": len(par),
                          "coinciden_pct": 100.0 * float(
                              (par[a] == par[b]).mean())})
    return pd.DataFrame(filas)


# ===========================================================================
# Frontera de aplicabilidad
# ===========================================================================

def frontera(d: pd.DataFrame) -> pd.DataFrame:
    """Porcentaje de replicas que no llegan a terminar, por SNR y muestras.

    No es un fallo del estudio sino una medida: marca donde CIRCUST deja de ser
    aplicable. El mecanismo es que el conjunto TOP colapsa y la etapa robusta se
    queda sin genes que muestrear.
    """
    rep = d.groupby(_clave(d)).first().reset_index()
    rep["fallo"] = ~(rep["error"].isna() | (rep["error"] == ""))
    return rep.pivot_table(index="n_samples", columns="snr", values="fallo",
                           aggfunc="mean") * 100.0


# ===========================================================================
# Figuras
# ===========================================================================

def figura_sesgo_omega(d: pd.DataFrame, salida: Path) -> Path:
    """El sesgo del apuntamiento cruza el cero segun el tamano muestral.

    Es el resultado central del bloque R y una tabla no transmite bien que hay
    un cruce. Se dibuja como mancuerna, un par por tamano muestral, contra la
    linea de sesgo nulo.
    """
    plt = _estilo()
    v = validas(d)
    t = v.pivot_table(index="n_samples", columns="fmm_method",
                      values="verdad_omega_sesgo", aggfunc="median")
    q = v.groupby(["n_samples", "fmm_method"])["verdad_omega_sesgo"].quantile(
        [0.25, 0.75]).unstack()

    fig, ax = plt.subplots(figsize=(6.6, 3.1))
    ax.axvline(0, color=TINTA_2, lw=1.0, ls=(0, (4, 3)), zorder=1)
    ax.text(0, -0.62, "sin sesgo", ha="center", va="bottom", fontsize=7.5,
            color=TINTA_2)

    for y, n in enumerate(t.index):
        b, s = t.loc[n, "base"], t.loc[n, "stabilized"]
        ax.plot([b, s], [y, y], color=GRIS, lw=2.0, zorder=2,
                solid_capstyle="round")
        for x, color, marca in ((b, AZUL, "o"), (s, NARANJA, "D")):
            lo, hi = q.loc[(n, "base" if color == AZUL else "stabilized")]
            ax.plot([lo, hi], [y, y], color=color, lw=1.0, alpha=0.45, zorder=3)
            ax.scatter([x], [y], s=52, color=color, marker=marca, zorder=4,
                       edgecolor=FONDO, linewidth=1.2)
        ax.annotate(f"{b:+.3f}", (b, y), textcoords="offset points",
                    xytext=(0, 9), ha="center", fontsize=7.5, color=AZUL)
        ax.annotate(f"{s:+.3f}", (s, y), textcoords="offset points",
                    xytext=(0, -15), ha="center", fontsize=7.5, color=NARANJA)

    ax.set_yticks(range(len(t.index)))
    ax.set_yticklabels([f"{n} muestras" for n in t.index])
    ax.set_ylim(-0.95, 2.55)
    ax.invert_yaxis()
    ax.set_xlabel("sesgo de $\\omega$ frente al valor verdadero "
                  "(negativo = onda mas apuntada de lo real)")
    ax.grid(axis="x", zorder=0)
    ax.scatter([], [], s=52, color=AZUL, marker="o", label="FMM base")
    ax.scatter([], [], s=52, color=NARANJA, marker="D", label="FMM estabilizado")
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, -0.42), ncol=2)
    fig.tight_layout()
    return _guardar(fig, salida)


def figura_ganancia_r(d: pd.DataFrame, salida: Path) -> Path:
    """Donde gana y donde pierde la regularizacion, en sus tres ejes.

    Dos paneles con las mismas categorias, uno por metrica, porque el resultado
    interesante es que **no coinciden**: la regularizacion mejora la fase casi
    siempre y sin embargo empeora el apuntamiento en las ondas mas picudas.
    """
    plt = _estilo()
    j = emparejar(d, "fmm_method", "base", "stabilized", METRICAS_R)

    ejes = [("n_samples", ["12", "48", "192"], "muestras"),
            ("waveform_mix", ["simetrica", "gtex", "apuntada", "extrema"],
             "forma de onda"),
            ("sampling", ["equiespaciado", "aleatorio", "huecos",
                          "concentrado"], "muestreo")]
    etiquetas, valores = [], {m: [] for m in
                              ("verdad_omega_medae", "verdad_picos_medae_h")}
    for eje, niveles, titulo in ejes:
        etiquetas.append((f"— {titulo} —", True))
        for m in valores:
            valores[m].append(np.nan)
        for niv in niveles:
            sub = j[j[eje].astype(str) == niv]
            etiquetas.append((niv, False))
            for m in valores:
                valores[m].append(float(sub[f"dif_{m}"].median())
                                  if len(sub) else np.nan)

    y = np.arange(len(etiquetas))
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 4.4), sharey=True)
    titulos = [("verdad_omega_medae", "error de $\\omega$"),
               ("verdad_picos_medae_h", "error de la fase de pico (h)")]

    for ax, (metrica, titulo) in zip(axes, titulos):
        vals = np.array(valores[metrica], dtype=float)
        colores = [AZUL if v < 0 else ROJO for v in np.nan_to_num(vals)]
        ax.barh(y, vals, color=colores, height=0.62, zorder=3)
        ax.axvline(0, color=TINTA_2, lw=1.0, zorder=4)
        for yy, v in zip(y, vals):
            if np.isnan(v):
                continue
            ax.annotate(f"{v:+.3f}", (v, yy), textcoords="offset points",
                        xytext=(6 if v >= 0 else -6, 0),
                        ha="left" if v >= 0 else "right", va="center",
                        fontsize=7, color=TINTA_2)
        ax.set_title(titulo, pad=8)
        ax.set_xlabel("estabilizado menos base")
        ax.grid(axis="x", zorder=0)
        m = np.nanmax(np.abs(vals)) * 1.45
        ax.set_xlim(-m, m)

    axes[0].set_yticks(y)
    axes[0].set_yticklabels([t for t, _ in etiquetas])
    for tick, (_, cabecera) in zip(axes[0].get_yticklabels(), etiquetas):
        tick.set_color(TINTA if cabecera else TINTA_2)
        tick.set_fontsize(8 if cabecera else 8)
    axes[0].invert_yaxis()

    fig.suptitle("Azul, gana el estabilizado.   Rojo, gana el base.",
                 fontsize=8.5, color=TINTA_2, y=0.985)
    fig.tight_layout(rect=(0, 0, 1, 0.955))
    return _guardar(fig, salida)


def figura_bloque_a(d: pd.DataFrame, salida: Path) -> Path:
    """La agregacion no mejora, y el desacuerdo es un sintoma, no una ocasion.

    Panel izquierdo: rho por estrato de desacuerdo, con el mejor metodo
    destacado y el resto en gris, porque la lectura es que van juntos. Panel
    derecho: diferencia pareada de la agregacion frente a elegir la mejor
    repeticion, que es la comparacion que decide.
    """
    plt = _estilo()
    t = tabla_bloque_a(d)
    dif = tabla_agregacion(d)
    n_por_estrato = (con_estrato(d).query("metodo == 'best_k:r2_median'")
                     .groupby("estrato", observed=True).size())

    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.1))
    x = np.arange(len(t.columns))

    ax = axes[0]
    for metodo in t.index:
        if metodo == "median":
            continue
        destacado = metodo == "best_k:r2_median"
        ax.plot(x, t.loc[metodo], color=AZUL if destacado else GRIS,
                lw=2.0 if destacado else 1.4, zorder=4 if destacado else 2,
                marker="o" if destacado else None, ms=5,
                markeredgecolor=FONDO, markeredgewidth=1.0,
                label="best_k:r2_median" if destacado else None)
    ax.plot([], [], color=GRIS, lw=1.4, label="los otros cuatro metodos")
    ax.set_xticks(x)
    ax.set_xticklabels([f"{c}\n(n={n_por_estrato.get(c, 0)})"
                        for c in t.columns])
    ax.set_xlabel("desacuerdo entre las K repeticiones")
    ax.set_ylabel("$\\rho$ frente al orden verdadero")
    ax.grid(axis="y", zorder=0)
    ax.legend(loc="upper right")

    # La escala de este panel es deliberada. Las diferencias son del orden de
    # 0,003 en rho y con un eje ajustado a ellas pareceria haber un efecto
    # donde no lo hay. Se fija el rango al de la dispersion real de las
    # diferencias por replica, que es lo que permite ver que la mediana esta
    # sobre el cero y muy dentro del ruido.
    ax = axes[1]
    v = con_estrato(d)
    base = v[v["metodo"] == "best_k:r2_median"].set_index(_clave(v))
    otro = v[v["metodo"] == "aggregate:tsp3"].set_index(_clave(v))
    j = base[["rho", "estrato"]].join(otro[["rho"]], rsuffix="__b",
                                      how="inner").dropna()
    j["dif"] = j["rho__b"] - j["rho"]

    estratos = [c for c in t.columns if c in set(j["estrato"])]
    med = [float(j.loc[j.estrato == e, "dif"].median()) for e in estratos]
    q1 = [float(j.loc[j.estrato == e, "dif"].quantile(0.25)) for e in estratos]
    q3 = [float(j.loc[j.estrato == e, "dif"].quantile(0.75)) for e in estratos]
    x = np.arange(len(estratos))

    ax.axhline(0, color=TINTA_2, lw=1.0, zorder=4)
    ax.vlines(x, q1, q3, color=GRIS, lw=7, zorder=3)
    ax.scatter(x, med, s=46, color=AZUL, zorder=5, edgecolor=FONDO,
               linewidth=1.2)
    for xx, v_ in zip(x, med):
        ax.annotate(f"{v_:+.4f}", (xx, v_), textcoords="offset points",
                    xytext=(0, 11), ha="center", fontsize=7.5, color=TINTA_2)

    ax.set_xticks(x)
    ax.set_xticklabels(estratos)
    ax.set_xlim(-0.6, len(estratos) - 0.4)
    ax.set_xlabel("desacuerdo entre las K repeticiones")
    ax.set_ylabel("$\\rho$ de agregar menos $\\rho$ de elegir")
    ax.grid(axis="y", zorder=0)
    ax.set_title("Mediana y recorrido intercuartilico por replica",
                 fontsize=8.5, color=TINTA_2, pad=8)

    fig.tight_layout()
    return _guardar(fig, salida)


def figura_frontera(d: pd.DataFrame, salida: Path) -> Path:
    """Donde CIRCUST deja de ser aplicable.

    Rejilla de magnitud, luego una sola tinta de claro a oscuro. Cada celda
    lleva su valor escrito, de modo que la figura se lee tambien en blanco y
    negro.
    """
    from matplotlib.colors import LinearSegmentedColormap
    plt = _estilo()
    t = frontera(d)

    fig, ax = plt.subplots(figsize=(4.6, 2.7))
    cmap = LinearSegmentedColormap.from_list("azul", RAMPA_AZUL)
    im = ax.imshow(t.to_numpy(), cmap=cmap, vmin=0, vmax=100, aspect="auto")

    for i in range(t.shape[0]):
        for j_ in range(t.shape[1]):
            v = t.to_numpy()[i, j_]
            ax.text(j_, i, f"{v:.0f}", ha="center", va="center", fontsize=9,
                    color="#ffffff" if v > 55 else TINTA)

    ax.set_xticks(range(t.shape[1]))
    ax.set_xticklabels([f"{c:g}" for c in t.columns])
    ax.set_yticks(range(t.shape[0]))
    ax.set_yticklabels([f"{i}" for i in t.index])
    ax.set_xlabel("razon senal-ruido")
    ax.set_ylabel("numero de muestras")
    ax.set_title("Replicas que no llegan a terminar (%)", pad=8, fontsize=9)
    for lado in ("left", "bottom"):
        ax.spines[lado].set_visible(False)
    ax.tick_params(length=0)
    fig.colorbar(im, ax=ax, fraction=0.04, pad=0.03).outline.set_visible(False)
    fig.tight_layout()
    return _guardar(fig, salida)


def _guardar(fig, salida: Path) -> Path:
    salida = Path(salida)
    salida.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(salida, bbox_inches="tight")
    fig.savefig(salida.with_suffix(".png"), dpi=170, bbox_inches="tight")
    fig.clf()
    return salida


# ===========================================================================
# CLI
# ===========================================================================

def _informe_r(d: pd.DataFrame) -> None:
    print("=" * 74)
    print("BLOQUE R, regularizacion del FMM")
    print("=" * 74)
    v = validas(d)
    print(f"{len(v)} filas validas de {len(d)}, "
          f"{v.groupby(_clave(v)).ngroups} replicas\n")
    print("Sesgo y error de omega, y error de la fase de pico, por tamano")
    print(tabla_sesgo_omega(d).round(3).to_string())
    for metrica, titulo in (("verdad_omega_medae", "ERROR DE OMEGA"),
                            ("verdad_picos_medae_h", "ERROR DE LA FASE")):
        print(f"\n{titulo}, diferencia pareada estabilizado menos base")
        print(tabla_bloque_r(d, metrica).round(4).to_string(index=False))


def _informe_a(d: pd.DataFrame) -> None:
    print("=" * 74)
    print("BLOQUE A, orden representativo")
    print("=" * 74)
    v = validas(d)
    print(f"{len(v)} filas validas de {len(d)}, "
          f"{v.groupby(_clave(v)).ngroups} replicas\n")
    print("rho por metodo y estrato de desacuerdo")
    print(tabla_bloque_a(d).round(3).to_string())
    print("\nDiferencia pareada frente a best_k:r2_median (positivo = mejor)")
    print(tabla_agregacion(d).round(4).to_string(index=False))
    print("\nCon que frecuencia dos criterios eligen la misma repeticion")
    print(tabla_seleccion(d).round(1).to_string(index=False))
    print("\nFrontera de aplicabilidad, replicas que no terminan (%)")
    print(frontera(d).round(0).to_string())


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("bloque", choices=["R", "A", "todo"])
    ap.add_argument("--raiz", type=Path, default=RAIZ_POR_DEFECTO)
    ap.add_argument("--figuras", type=Path, default=None,
                    help="directorio donde escribir las figuras en PDF.")
    args = ap.parse_args()

    quiere = ["R", "A"] if args.bloque == "todo" else [args.bloque]
    for bloque in quiere:
        d = leer(bloque, args.raiz)
        (_informe_r if bloque == "R" else _informe_a)(d)
        print()
        if args.figuras is None:
            continue
        if bloque == "R":
            print("figuras:",
                  figura_sesgo_omega(d, args.figuras / "sim-sesgo-omega.pdf"),
                  figura_ganancia_r(d, args.figuras / "sim-ganancia-fmm.pdf"))
        else:
            print("figuras:",
                  figura_bloque_a(d, args.figuras / "sim-agregacion.pdf"),
                  figura_frontera(d, args.figuras / "sim-frontera.pdf"))


if __name__ == "__main__":
    main()
