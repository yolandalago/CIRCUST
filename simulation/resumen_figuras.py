"""
Figuras del resumen de resultados de una tanda de simulaciones.

Lee los CSV de ``run.py`` y ``rbench.py`` y dibuja las figuras del documento
de resultados. Todo sale de los CSV, ningun numero va escrito a mano, asi que
al regenerar la simulacion basta con volver a correrlo.

Uso
---
    python -m simulation.resumen_figuras --out output/simulacion_2026-10
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

ACENTO = "#2a6fdb"
GRIS = "#8f8f89"
GRIS_CLARO = "#c9c9c3"
TEXTO = "#222222"
TEXTO_2 = "#666666"
REJILLA = "#e8e8e5"
NARANJA = "#d9822b"

CLAVE = ["n_samples", "sampling", "omega", "snr", "n_genes", "prop_rhythmic",
         "n_core", "rep"]
REFERENCIA = "best_k:r2_median"
NOMBRES = {
    "best_k:r2_median":     "best_k, R² mediana",
    "best_k:n_well_fitted": "best_k, genes bien ajustados",
    "best_k:medoid":        "best_k, medoide",
    "circmean":             "media circular",
    "hod":                  "HODs",
    "tsp":                  "TSP3",
}
ORDEN = ["best_k:n_well_fitted", "best_k:medoid", "circmean", "hod", "tsp"]

plt.rcParams.update({
    "font.size": 9.5, "axes.edgecolor": GRIS_CLARO, "axes.labelcolor": TEXTO_2,
    "xtick.color": TEXTO_2, "ytick.color": TEXTO_2, "figure.dpi": 160,
    "axes.titlesize": 10, "axes.titlecolor": TEXTO,
})


def _coma(v, dec: int = 2) -> str:
    return f"{v:.{dec}f}".replace(".", ",").replace("-", "−")


def _nivel(v) -> str:
    """1.5 -> '1,5', 2.0 -> '2'."""
    return f"{v:g}".replace(".", ",")


def _coma_ejes(ax, ejes: str = "xy") -> None:
    """Coma decimal en las marcas de los ejes."""
    from matplotlib.ticker import FuncFormatter
    f = FuncFormatter(lambda v, _: f"{v:g}".replace(".", ",").replace("-", "−"))
    if "x" in ejes:
        ax.xaxis.set_major_formatter(f)
    if "y" in ejes:
        ax.yaxis.set_major_formatter(f)


def _limpiar(ax, eje: str = "y") -> None:
    for lado in ("top", "right"):
        ax.spines[lado].set_visible(False)
    ax.grid(axis=eje, color=REJILLA, linewidth=0.8)
    ax.set_axisbelow(True)


def _titulos(fig, titulo: str, subtitulo: str) -> None:
    fig.text(0.015, 0.985, titulo, fontsize=12.5, color=TEXTO, va="top")
    fig.text(0.015, 0.925, subtitulo, fontsize=9, color=TEXTO_2, va="top")


def cargar(raiz: Path) -> pd.DataFrame:
    d = pd.read_csv(raiz / "rejilla" / "resultados.csv")
    d["metodo"] = np.where(d.robust_method == "best_k",
                           "best_k:" + d.selection_criterion, d.aggregation)
    d["genes"] = ((d.n_genes * d.prop_rhythmic).round().astype(int).astype(str)
                  + " / " + d.n_genes.astype(str))
    return d


def _referencia(d: pd.DataFrame) -> pd.DataFrame:
    return d[(d.fmm_method == "base") & (d.metodo == REFERENCIA)]


def _pareado(d: pd.DataFrame, metodo: str, medida: str, fmm: str) -> pd.DataFrame:
    """Diferencia metodo - referencia sobre la misma replica y el mismo FMM."""
    b = d[d.fmm_method == fmm]
    ref = b[b.metodo == REFERENCIA].set_index(CLAVE)[medida]
    m = b[b.metodo == metodo].set_index(CLAVE)[medida]
    return (m - ref).rename("delta").reset_index()


# ── Dificultad ─────────────────────────────────────────────────────────────

def fig_mapas_dificultad(d: pd.DataFrame, destino: Path) -> None:
    r = _referencia(d)
    ns = sorted(r.n_samples.unique())
    vmax = r.groupby(["n_samples", "omega", "snr"]).error_orden_h.median().max()
    cmap = LinearSegmentedColormap.from_list("seq", ["#f3f6fc", ACENTO, "#123a7a"])
    fig, axes = plt.subplots(1, len(ns), figsize=(9, 3.6), sharey=True)
    fig.subplots_adjust(top=0.72, left=0.08, right=0.99, bottom=0.15, wspace=0.08)
    for ax, n in zip(axes, ns):
        t = r[r.n_samples == n].pivot_table(index="omega", columns="snr",
                                           values="error_orden_h", aggfunc="median")
        ax.imshow(t.to_numpy(), cmap=cmap, vmin=0, vmax=vmax, aspect="auto")
        for i in range(t.shape[0]):
            for j in range(t.shape[1]):
                v = t.to_numpy()[i, j]
                ax.text(j, i, _coma(v), ha="center", va="center", fontsize=9,
                        color="white" if v > 0.55 * vmax else TEXTO)
        ax.set_xticks(range(t.shape[1]), [_nivel(v) for v in t.columns])
        ax.set_yticks(range(t.shape[0]), [f"ω = {_nivel(v)}" for v in t.index])
        ax.set_xlabel("A/σ")
        ax.set_title(f"n = {n}")
        for lado in ax.spines.values():
            lado.set_visible(False)
        ax.tick_params(length=0)
    _titulos(fig, "El error del orden lo deciden ω y la señal. n apenas cambia el patrón",
             "Mediana del error del orden (h) por escenario. best_k por R² mediana con FMM base,\n"
             "sobre los tres muestreos, los cinco tamaños de matriz y los tres conjuntos core.")
    fig.savefig(destino / "fig1_dificultad_mapas.png")
    plt.close(fig)


def fig_cajas_factores(d: pd.DataFrame, destino: Path) -> None:
    r = _referencia(d)
    factores = [("n_samples", "muestras n"), ("sampling", "muestreo"),
                ("omega", "apuntamiento ω"), ("snr", "señal A/σ"),
                ("genes", "rítmicos / total"), ("n_core", "genes core declarados")]
    orden_genes = ["50 / 1000", "200 / 1000", "50 / 5000", "250 / 5000", "1000 / 5000"]
    fig, axes = plt.subplots(2, 3, figsize=(9, 5.6), sharey=True)
    fig.subplots_adjust(top=0.82, left=0.07, right=0.99, bottom=0.09, hspace=0.5, wspace=0.08)
    for ax, (col, nombre) in zip(axes.ravel(), factores):
        niveles = orden_genes if col == "genes" else sorted(r[col].unique())
        datos = [r.loc[r[col] == v, "error_orden_h"].dropna().to_numpy() for v in niveles]
        cajas = ax.boxplot(datos, widths=0.55, showfliers=False, whis=(5, 95),
                           patch_artist=True, medianprops=dict(color=ACENTO, lw=2))
        for caja in cajas["boxes"]:
            caja.set(facecolor="#efefec", edgecolor=GRIS)
        for parte in ("whiskers", "caps"):
            for linea in cajas[parte]:
                linea.set(color=GRIS)
        etiquetas = [v.replace(" / ", "\n/ ") if col == "genes"
                     else (v if isinstance(v, str) else _nivel(v)) for v in niveles]
        ax.set_xticks(range(1, len(niveles) + 1), etiquetas, fontsize=8.5)
        _coma_ejes(ax, "y")
        ax.set_title(nombre)
        _limpiar(ax)
    axes[0, 0].set_ylabel("error del orden (h)")
    axes[1, 0].set_ylabel("error del orden (h)")
    rango = max(r.groupby(c).error_orden_h.median().pipe(lambda s: s.max() - s.min())
                for c, _ in factores if c not in ("omega", "snr"))
    _titulos(fig, f"ω y A/σ dominan. El resto de factores mueve la mediana como mucho {_coma(rango, 1)} h",
             "Error del orden por nivel de cada factor, sobre el resto de factores. Caja = cuartiles, "
             "línea azul = mediana, bigotes = percentiles 5 y 95.\nbest_k por R² mediana con FMM base.")
    fig.savefig(destino / "fig2_dificultad_factores.png")
    plt.close(fig)


# ── Metodos de resumen ─────────────────────────────────────────────────────

def fig_bosque(d: pd.DataFrame, destino: Path) -> None:
    medidas = [("error_orden_h", "error del orden (min)", 60),
               ("cre", "CRE", 1),
               ("error_fase_h", "error de fase absoluta (min)", 60)]
    fig, axes = plt.subplots(1, 3, figsize=(9.5, 3.9), sharey=True)
    fig.subplots_adjust(top=0.76, left=0.24, right=0.99, bottom=0.14, wspace=0.12)
    y = np.arange(len(ORDEN))[::-1]
    for ax, (col, nombre, escala) in zip(axes, medidas):
        for k, (fmm, color, relleno) in enumerate([("base", ACENTO, ACENTO),
                                                   ("circafmm", GRIS, "white")]):
            for yi, met in zip(y, ORDEN):
                x = _pareado(d, met, col, fmm).delta.dropna() * escala
                q1, me, q3 = x.quantile([0.25, 0.5, 0.75])
                yy = yi + (0.14 if k == 0 else -0.14)
                ax.plot([q1, q3], [yy, yy], color=color, lw=2)
                ax.scatter(me, yy, s=34, color=relleno, edgecolor=color, zorder=3,
                           label=("FMM base" if fmm == "base" else "CircaFMM") if yi == y[0] else None)
        ax.axvline(0, color=GRIS, lw=1)
        ax.set_title(nombre)
        _limpiar(ax, "x")
        _coma_ejes(ax, "x")
    axes[0].set_yticks(y, [NOMBRES[m] for m in ORDEN])
    axes[0].tick_params(axis="y", length=0)
    axes[0].legend(frameon=False, loc="lower left", bbox_to_anchor=(-0.02, 1.08), ncol=2)
    _titulos(fig, "La media circular y HODs mejoran a la referencia en las tres medidas. TSP3 tiene colas largas",
             "Diferencia pareada con best_k por R² mediana en la misma réplica: punto = mediana, "
             "línea = cuartiles. Negativo = mejora.")
    fig.savefig(destino / "fig3_metodos_bosque.png")
    plt.close(fig)


def fig_violines(d: pd.DataFrame, destino: Path) -> None:
    lim = 60.0
    fig, ax = plt.subplots(figsize=(9, 4.4))
    fig.subplots_adjust(top=0.78, left=0.08, right=0.99, bottom=0.12)
    datos, fuera = [], []
    for met in ORDEN:
        x = _pareado(d, met, "error_orden_h", "base").delta.dropna().to_numpy() * 60
        fuera.append(np.mean(np.abs(x) > lim) * 100)
        datos.append(np.clip(x, -lim, lim))
    partes = ax.violinplot(datos, showextrema=False, widths=0.8, points=300)
    for cuerpo, met in zip(partes["bodies"], ORDEN):
        cuerpo.set(facecolor=ACENTO if met in ("circmean", "hod") else GRIS_CLARO,
                   edgecolor="none", alpha=0.75)
    for i, (x, f) in enumerate(zip(datos, fuera), start=1):
        ax.scatter(i, np.median(x), color=TEXTO, s=14, zorder=3)
        ax.text(i, lim + 4, f"{_coma(f, 1)} % fuera", ha="center", fontsize=8, color=TEXTO_2)
    ax.axhline(0, color=GRIS, lw=1)
    ax.set_xticks(range(1, len(ORDEN) + 1), [NOMBRES[m] for m in ORDEN])
    ax.set_ylabel("diferencia en el error del orden (min)")
    ax.set_ylim(-lim - 5, lim + 12)
    _limpiar(ax)
    _titulos(fig, "Casi todas las diferencias son de pocos minutos, pero TSP3 pierde a veces más de una hora",
             f"Distribución de la diferencia pareada con best_k por R² mediana, FMM base. Negativo = mejora.\n"
             f"Recortada a ±{lim:.0f} min, con el porcentaje de casos que queda fuera. Punto = mediana.")
    fig.savefig(destino / "fig4_metodos_violines.png")
    plt.close(fig)


# ── CircaFMM ───────────────────────────────────────────────────────────────

def _pares_fmm(d: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Los dos brazos alineados, solo pares con el sentido acertado en ambos."""
    clave = CLAVE + ["metodo"]
    b = d[d.fmm_method == "base"].set_index(clave)
    c = d[d.fmm_method == "circafmm"].set_index(clave).reindex(b.index)
    ok = (b.tau_circular > 0) & (c.tau_circular > 0)
    return b[ok], c[ok]


def fig_circafmm_victorias(d: pd.DataFrame, destino: Path) -> None:
    b, c = _pares_fmm(d)
    medidas = [("atlas_variabilidad_pico_h", "variabilidad del error en el pico"),
               ("atlas_error_onda", "error de la onda")]
    cmap = LinearSegmentedColormap.from_list("div", [NARANJA, "#f2f2ef", ACENTO])
    norm = TwoSlopeNorm(vmin=0, vcenter=50, vmax=100)
    fig, axes = plt.subplots(1, 2, figsize=(8.5, 3.8), sharey=True)
    fig.subplots_adjust(top=0.72, left=0.11, right=0.99, bottom=0.14, wspace=0.08)
    for ax, (col, nombre) in zip(axes, medidas):
        gana = ((c[col] - b[col]) < 0).rename("gana").reset_index()
        t = gana.pivot_table(index="snr", columns="n_samples", values="gana", aggfunc="mean") * 100
        ax.imshow(t.to_numpy(), cmap=cmap, norm=norm, aspect="auto")
        for i in range(t.shape[0]):
            for j in range(t.shape[1]):
                ax.text(j, i, f"{t.to_numpy()[i, j]:.0f} %", ha="center", va="center",
                        fontsize=9, color=TEXTO)
        ax.set_xticks(range(t.shape[1]), [f"n = {n}" for n in t.columns])
        ax.set_yticks(range(t.shape[0]), [f"A/σ = {_nivel(v)}" for v in t.index])
        ax.set_title(nombre)
        for lado in ax.spines.values():
            lado.set_visible(False)
        ax.tick_params(length=0)
    _titulos(fig, "CircaFMM mejora el atlas con poca información y lo empeora cuando sobra",
             "Porcentaje de pares en que CircaFMM reduce el error frente al FMM base: azul por encima del 50 %, "
             "naranja por debajo.\nSolo pares con el sentido de giro acertado en los dos modelos, todos los métodos de resumen.")
    fig.savefig(destino / "fig5_circafmm_victorias.png")
    plt.close(fig)


def fig_circafmm_omega(d: pd.DataFrame, destino: Path) -> None:
    b, c = _pares_fmm(d)
    omegas = sorted(b.reset_index().omega.unique())
    medidas = [("atlas_sesgo_omega", "sesgo de ω (ω estimado − ω real)", True),
               ("atlas_error_onda", "error de la onda", False)]
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.8))
    fig.subplots_adjust(top=0.76, left=0.07, right=0.99, bottom=0.15, wspace=0.22)
    x = np.arange(len(omegas))
    for ax, (col, nombre, cero) in zip(axes, medidas):
        for k, (brazo, etiqueta, color, relleno) in enumerate([(b, "FMM base", GRIS, "white"),
                                                               (c, "CircaFMM", ACENTO, ACENTO)]):
            v = brazo.reset_index()
            for xi, w in zip(x, omegas):
                s = v.loc[v.omega == w, col].dropna()
                q1, me, q3 = s.quantile([0.25, 0.5, 0.75])
                xx = xi + (-0.12 if k == 0 else 0.12)
                ax.plot([xx, xx], [q1, q3], color=color, lw=2)
                ax.scatter(xx, me, s=36, color=relleno, edgecolor=color, zorder=3,
                           label=etiqueta if xi == 0 else None)
                ax.annotate(_coma(me), (xx, me), xytext=(7 if k else -7, 0),
                            textcoords="offset points", va="center",
                            ha="left" if k else "right", fontsize=8, color=TEXTO)
        if cero:
            ax.axhline(0, color=GRIS, lw=1)
        ax.set_xticks(x, [f"ω real = {_nivel(w)}" for w in omegas])
        ax.set_xlim(-0.5, len(omegas) - 0.35)
        ax.set_title(nombre)
        _limpiar(ax)
        _coma_ejes(ax, "y")
    axes[0].legend(frameon=False, loc="lower left", bbox_to_anchor=(0, 1.08), ncol=2)
    _titulos(fig, "CircaFMM sobreestima ω en ondas apuntadas y solo mejora la onda cuando ω ≥ 0,5",
             "Punto = mediana, línea = cuartiles por conjunto de datos. Pares con el sentido acertado en los dos modelos.")
    fig.savefig(destino / "fig6_circafmm_omega.png")
    plt.close(fig)


# ── Tiempos ────────────────────────────────────────────────────────────────

def fig_tiempos(raiz: Path, destino: Path) -> None:
    t = pd.read_csv(raiz / "tiempos.csv")
    m = t.groupby(["n_genes", "n_samples", "implementacion"]).t_total_s.median().unstack()
    m = m.sort_index(ascending=[False, False])
    fig, ax = plt.subplots(figsize=(8, 4.6))
    fig.subplots_adjust(top=0.8, left=0.22, right=0.86, bottom=0.12)
    for i, (_, fila) in enumerate(m.iterrows()):
        py, rr = fila.get("python"), fila.get("r")
        if pd.notna(rr):
            ax.plot([py, rr], [i, i], color=GRIS_CLARO, lw=2, zorder=1)
            ax.scatter(rr, i, color=GRIS, s=40, zorder=2)
            ax.text(1.02, i, f"×{rr / py:.0f}", transform=ax.get_yaxis_transform(),
                    va="center", color=TEXTO)
        else:
            ax.text(1.02, i, "sin dato de R", transform=ax.get_yaxis_transform(),
                    va="center", color=TEXTO_2)
        ax.scatter(py, i, color=ACENTO, s=40, zorder=3)
    ax.set_yticks(range(len(m)), [f"{g} genes · n = {n}" for g, n in m.index])
    ax.set_xscale("log")
    from matplotlib.ticker import FixedLocator, NullFormatter
    ax.xaxis.set_major_locator(FixedLocator([1, 10, 100, 1000]))
    ax.set_xticklabels(["1 s", "10 s", "100 s", "1000 s"])
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_xlabel("tiempo hasta las K repeticiones (s, escala logarítmica)")
    ax.scatter([], [], color=ACENTO, label="Python")
    ax.scatter([], [], color=GRIS, label="R")
    ax.legend(frameon=False, loc="lower left", bbox_to_anchor=(0, 1.0), ncol=2)
    for lado in ("top", "right", "left"):
        ax.spines[lado].set_visible(False)
    ax.grid(axis="x", color=REJILLA)
    ax.tick_params(axis="y", length=0)
    fac = (m["r"] / m["python"]).dropna()
    _titulos(fig, f"Python tarda entre {fac.min():.0f} y {fac.max():.0f} veces menos que R",
             "Mediana del tiempo total por tamaño de matriz, un solo proceso en las dos implementaciones.")
    fig.savefig(destino / "fig7_tiempos_python_r.png")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument("--out", default="output/simulacion_2026-10")
    args = parser.parse_args()
    raiz = Path(args.out)
    destino = raiz / "figuras"
    destino.mkdir(parents=True, exist_ok=True)
    d = cargar(raiz)
    fig_mapas_dificultad(d, destino)
    fig_cajas_factores(d, destino)
    fig_bosque(d, destino)
    fig_violines(d, destino)
    fig_circafmm_victorias(d, destino)
    fig_circafmm_omega(d, destino)
    fig_tiempos(raiz, destino)
    print(f"Figuras en {destino}")


if __name__ == "__main__":
    main()
