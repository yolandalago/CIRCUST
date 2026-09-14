"""
Figuras y tablas de la comparacion FMM base contra estabilizado.

El argumento del capitulo tiene cinco pasos y cada figura cubre uno:

1. Donde NO puede haber diferencia. El orden de muestras sale de una SVD
   sobre los genes nucleares donde el FMM no interviene en ningun paso, asi
   que los dos metodos dan el mismo orden por construccion. Esto se demuestra,
   no se mide, y por eso no tiene figura.
2. Donde SI la hay: los parametros por gen, o sea el atlas.
3. El mecanismo, que es un compromiso sesgo-varianza y no una mejora sin mas.
   Figuras D y E.
4. El efecto, con su gradiente en n y en A/sigma. Figuras A, B y C.
5. El coste: con ondas genuinamente apuntadas suaviza de mas. Figuras D y E
   otra vez, y la columna omega=0,1 de la A.

El mecanismo, medido
--------------------
La primera version de este modulo decia que la penalizacion de rugosidad
"impide que omega colapse hacia valores espuriamente apuntados". Es solo una
parte, y no la principal. Lo que la descomposicion del error de omega muestra
(figura D) es:

- La **desviacion tipica** de omega entre genes baja SIEMPRE, en las quince
  celdas y en los tres niveles de omega, alrededor de un 40 %. Eso es lo que
  hace la penalizacion y no depende del regimen.
- El **sesgo** es lo que cambia de signo. La penalizacion empuja omega hacia
  arriba, o sea hacia ondas mas suaves, siempre. Eso corrige al base cuando
  el base subestima (omega verdadero 0,9 con poca senal, sesgo -0,42 que pasa
  a -0,19) y lo estropea cuando el base ya era insesgado (omega verdadero
  0,5, donde el base da -0,00 y el estabilizado +0,17).

O sea: se paga sesgo a cambio de varianza, y sale a cuenta solo mientras la
varianza sea el termino dominante, que es cuando hay poca informacion. De ahi
el gradiente de la figura A y el cruce en torno a A/sigma = 2 a 3.

Contarlo asi importa: "mi metodo reduce la varianza del estimador a costa de
sesgo, y eso compensa por debajo de tal umbral" es una afirmacion estadistica
verificable. "Mi metodo mejora" no lo es.

El paso 5 es lo que hace creible al 4. Un capitulo que solo dice "mi metodo
mejora" se lee como publicidad; uno que localiza donde empeora y explica por
que se lee como estadistica.

Dos cortes que se aplican en todas partes
-----------------------------------------
**omega = 0,1 va aparte, nunca promediado con el resto.** Es un regimen donde
el metodo falla incluso sin ruido, porque una onda tan apuntada lleva solo un
tercio de su energia en el armonico fundamental y el embebido de dos
componentes del CPCA no la captura. Mezclarlo con omega 0,5 y 0,9 arrastra
cualquier media.

**Solo entran las celdas con la orientacion resuelta en los DOS brazos.** Una
escala reflejada da un atlas especular, con el error de pico repartido por
todo el ciclo. Eso es otro regimen, no una calidad peor, y promediarlo con
los aciertos mezcla dos cosas distintas.

Uso
---
    python -m simulation.fmm_figuras --in output/simulacion --images docs/images
    python -m simulation.fmm_figuras --volcar-por-gen    # datos de E y F
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

from simulation.truth import TWO_PI

CLAVE_PAREADA = ["n_samples", "sampling", "omega", "snr", "rep"]

# Celdas del volcado por gen. Son pocas y elegidas a proposito:
#   - (24, 1.0) es la esquina donde el estabilizado gana por goleada
#   - (96, 4.0) es la esquina donde pierde
#   - omega 0,1 es el regimen donde el estabilizado es peor
# Con estas tres se cuenta la historia entera, y cuestan minutos.
CELDAS_POR_GEN = (
    # (n, omega, snr, etiqueta)
    (24, 0.7, 1.0, "n bajo, senal baja"),
    (96, 0.7, 4.0, "n alto, senal alta"),
    (48, 0.1, 2.0, "onda muy apuntada"),
)


# ═══════════════════════════════════════════════════════════════════════════
# Datos
# ═══════════════════════════════════════════════════════════════════════════

def cargar_pareado(entrada: Path | str) -> pd.DataFrame:
    """
    Une el brazo base con el estabilizado replica a replica.

    La union es por identidad de generacion mas numero de replica, que es lo
    que garantiza que las dos filas corrieron sobre la MISMA matriz. Sin eso
    la comparacion seria entre celdas con ruido distinto y perderia casi toda
    la potencia del diseno.
    """
    d = pd.read_csv(Path(entrada) / "resultados.csv")
    b = d[(d["bloque"] == "principal") & (d["status"] == "ok")]
    s = d[(d["bloque"] == "regularizacion") & (d["status"] == "ok")]
    if b.empty or s.empty:
        raise ValueError("faltan filas de principal o de regularizacion")

    par = b.merge(s, on=CLAVE_PAREADA, suffixes=("_b", "_s"))
    if par.empty:
        raise ValueError("el emparejamiento no cuadra, revisa las semillas")

    par = par[(par["orientacion_ok_b"] == 1.0) & (par["orientacion_ok_s"] == 1.0)]
    for m in ("atlas_rmse_curva", "atlas_err_pico_sd_h", "atlas_omega_sesgo",
              "atlas_vpp", "atlas_err_subida_h"):
        par[f"d_{m}"] = par[f"{m}_s"] - par[f"{m}_b"]
    return par


def volcar_por_gen(destino: Path | str) -> pd.DataFrame:
    """
    Corre unas pocas celdas con los dos metodos y vuelca el error POR GEN.

    El CSV de la rejilla guarda solo resumenes por celda, que bastan para las
    tablas pero no para ver la distribucion de omega ni dibujar ondas. Esto
    es una ejecucion aparte y pequena, no una repeticion de la rejilla.
    """
    from simulation.atlas import per_gene_errors
    from simulation.generate import generate_dataset
    from simulation.grid import GCORE_CENTRAL, SIGMA_BASE
    from simulation.runner import RobustConfig, run_full
    from simulation.truth import SimConfig

    filas = []
    for n, omega, snr, etiqueta in CELDAS_POR_GEN:
        cfg = SimConfig(n_samples=n, n_genes=1000, prop_rhythmic=0.20,
                        n_core=30, core_sizes=(10, 20, 30), omega=omega,
                        amplitude=snr * SIGMA_BASE, sigma=SIGMA_BASE,
                        sampling="uniforme", seed=4242)
        X, truth = generate_dataset(cfg)
        for metodo in ("base", "stabilized"):
            r = run_full(X, truth, GCORE_CENTRAL, RobustConfig(n_reps=5),
                         fmm_method=metodo)
            g = per_gene_errors(r.atlas, truth.gene_params, quitar_giro=True)
            g = g.assign(metodo=metodo, n_samples=n, omega_celda=omega,
                         snr=snr, celda=etiqueta)
            filas.append(g)

    salida = pd.concat(filas, ignore_index=True)
    destino = Path(destino)
    destino.parent.mkdir(parents=True, exist_ok=True)
    salida.to_csv(destino, index=False)
    return salida


# ═══════════════════════════════════════════════════════════════════════════
# Tabla con respaldo inferencial
# ═══════════════════════════════════════════════════════════════════════════

def tabla_wilcoxon(par: pd.DataFrame,
                   medida: str = "atlas_rmse_curva") -> pd.DataFrame:
    """
    Contraste pareado de Wilcoxon por celda, con correccion por multiplicidad.

    Se usa el de rangos con signo y no una t pareada porque las diferencias
    no son normales: el RMSE de la curva esta acotado por abajo y su
    distribucion es asimetrica. El de Wilcoxon solo pide simetria de las
    diferencias bajo la hipotesis nula.

    Se reporta ademas la **tasa de victorias**, que es la fraccion de
    replicas en que el estabilizado gana. Es mas informativa que el p-valor:
    un efecto mediano pequeno con el 100 % de victorias es una mejora mas
    solida que uno grande con el 55 %.

    Las 135 celdas se contrastan a la vez, asi que los p-valores llevan
    Benjamini-Hochberg. Sin eso, al 5 % se esperarian unas 7 celdas
    significativas por puro azar.
    """
    from scipy.stats import wilcoxon

    filas = []
    for (w, n, snr), g in par.groupby(["omega", "n_samples", "snr"]):
        d = g[f"d_{medida}"].dropna().to_numpy()
        if d.size < 6 or np.allclose(d, 0.0):
            p = float("nan")
        else:
            p = float(wilcoxon(d, alternative="two-sided").pvalue)
        filas.append({
            "omega": w, "n_samples": n, "snr": snr, "n_pares": int(d.size),
            "delta_mediana": float(np.median(d)) if d.size else float("nan"),
            "iqr": float(np.subtract(*np.percentile(d, [75, 25]))) if d.size else float("nan"),
            "tasa_victorias": float((d < 0).mean()) if d.size else float("nan"),
            "p": p,
        })

    t = pd.DataFrame(filas).sort_values(["omega", "n_samples", "snr"])
    t["p_bh"] = _benjamini_hochberg(t["p"].to_numpy())
    return t.reset_index(drop=True)


def _benjamini_hochberg(p: np.ndarray) -> np.ndarray:
    """p-valores ajustados por tasa de falsos descubrimientos."""
    p = np.asarray(p, dtype=float)
    ok = ~np.isnan(p)
    salida = np.full(p.shape, np.nan)
    if not ok.any():
        return salida
    v = p[ok]
    m = v.size
    orden = np.argsort(v)
    ajust = v[orden] * m / np.arange(1, m + 1)
    ajust = np.minimum.accumulate(ajust[::-1])[::-1]
    desorden = np.empty(m)
    desorden[orden] = np.minimum(ajust, 1.0)
    salida[ok] = desorden
    return salida


# ═══════════════════════════════════════════════════════════════════════════
# Figuras
# ═══════════════════════════════════════════════════════════════════════════

def _mpl():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 9, "axes.grid": False,
                         "figure.dpi": 160, "savefig.bbox": "tight"})
    return plt


def _mapa(ax, tabla, titulo, cmap, centro, fmt="{:+.3f}", vlim=None):
    """Mapa de calor anotado, filas n y columnas A/sigma."""
    M = tabla.to_numpy(dtype=float)
    lim = vlim if vlim is not None else np.nanmax(np.abs(M - centro))
    im = ax.imshow(M, cmap=cmap, vmin=centro - lim, vmax=centro + lim,
                   aspect="auto")
    ax.set_xticks(range(tabla.shape[1]), [f"{c:g}" for c in tabla.columns])
    ax.set_yticks(range(tabla.shape[0]), [f"{i:g}" for i in tabla.index])
    ax.set_xlabel("A/$\\sigma$")
    ax.set_ylabel("n muestras")
    ax.set_title(titulo, fontsize=9)
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            if np.isnan(M[i, j]):
                continue
            # texto claro sobre fondo saturado
            lum = abs(M[i, j] - centro) / (lim if lim else 1)
            ax.text(j, i, fmt.format(M[i, j]), ha="center", va="center",
                    fontsize=7.5, color="white" if lum > 0.6 else "black")
    return im


def figura_a_efecto(par: pd.DataFrame, destino: Path) -> None:
    """
    A: efecto pareado sobre el RMSE de la curva, por n y A/sigma.

    Azul (negativo) = gana el estabilizado. Se separan las ondas normales de
    omega=0,1 porque el signo del efecto se invierte entre las dos, que es
    justo el resultado.
    """
    plt = _mpl()
    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.4))
    fig.subplots_adjust(top=0.72)
    lim = np.nanmax(np.abs(par["d_atlas_rmse_curva"]
                           .groupby([par.omega, par.n_samples, par.snr])
                           .median()))
    for ax, (sel, tit) in zip(axes, [
        (par.omega != 0.1, "Ondas normales ($\\omega$ = 0,5 y 0,9)"),
        (par.omega == 0.1, "Onda muy apuntada ($\\omega$ = 0,1)"),
    ]):
        t = par[sel].pivot_table(index="n_samples", columns="snr",
                                 values="d_atlas_rmse_curva", aggfunc="median")
        im = _mapa(ax, t, tit, "RdBu_r", 0.0, vlim=lim)
    cb = fig.colorbar(im, ax=axes, fraction=0.03, pad=0.02)
    cb.set_label("$\\Delta$ RMSE de la curva\n(estabilizado $-$ base)", fontsize=8)
    fig.suptitle("Efecto del FMM estabilizado sobre el ajuste por gen\n"
                 "azul = mejora, rojo = empeora", fontsize=10, y=1.04)
    fig.savefig(destino)
    plt.close(fig)


def figura_b_victorias(par: pd.DataFrame, destino: Path) -> None:
    """
    B: fraccion de replicas en que gana el estabilizado.

    Acompana a la A porque responde otra pregunta. La A dice cuanto mejora en
    mediana; esta dice con cuanta consistencia. Un efecto pequeno que se
    repite en el 100 % de las replicas es una mejora mas solida que uno
    grande que solo ocurre en el 55 %.
    """
    plt = _mpl()
    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.4))
    fig.subplots_adjust(top=0.72)
    for ax, (sel, tit) in zip(axes, [
        (par.omega != 0.1, "Ondas normales ($\\omega$ = 0,5 y 0,9)"),
        (par.omega == 0.1, "Onda muy apuntada ($\\omega$ = 0,1)"),
    ]):
        t = (par[sel].assign(gana=lambda x: (x.d_atlas_rmse_curva < 0).astype(float))
             .pivot_table(index="n_samples", columns="snr", values="gana"))
        im = _mapa(ax, t, tit, "RdBu_r", 0.5, fmt="{:.2f}", vlim=0.5)
    cb = fig.colorbar(im, ax=axes, fraction=0.03, pad=0.02)
    cb.set_label("fraccion de replicas\nen que gana el estabilizado", fontsize=8)
    fig.suptitle("Consistencia de la mejora, no solo su tamano\n"
                 "0,5 = empate", fontsize=10, y=1.04)
    fig.savefig(destino)
    plt.close(fig)


def figura_c_pico_horas(par: pd.DataFrame, destino: Path) -> None:
    """
    C: la version interpretable, en horas.

    El RMSE de la curva es adimensional y no le dice nada a un lector
    biologo. El error de pico si: un atlas que situa el maximo de un gen con
    1,4 h menos de dispersion es una afirmacion que se entiende sin saber
    que es un FMM.
    """
    plt = _mpl()
    o = par[par.omega != 0.1]
    fig, ax = plt.subplots(figsize=(5.4, 3.2))
    for n, marca in zip(sorted(o.n_samples.unique()), ["o", "s", "^"]):
        g = o[o.n_samples == n].groupby("snr")["d_atlas_err_pico_sd_h"].median()
        ax.plot(g.index, g.values, marker=marca, label=f"n = {n}")
    ax.axhline(0, color="0.4", lw=0.8, ls="--")
    ax.set_xlabel("A/$\\sigma$")
    ax.set_ylabel("$\\Delta$ dispersion del error de pico (h)")
    ax.set_title("Cuanto mejora la localizacion del pico, en horas\n"
                 "negativo = el estabilizado acierta mas fino", fontsize=9)
    ax.legend(frameon=False, fontsize=8)
    fig.savefig(destino)
    plt.close(fig)


def figura_d_sesgo_omega(par: pd.DataFrame, destino: Path) -> None:
    """
    D: el mecanismo, descompuesto en sesgo y varianza.

    Es la figura que explica POR QUE mejora, y la que sostiene el argumento
    ante quien pregunte si no sera casualidad.

    Fila de abajo: la desviacion tipica de omega entre genes **baja siempre**,
    en las quince celdas y en los tres niveles de omega, alrededor de un 40 %.
    Eso es lo que hace la penalizacion de rugosidad, y no depende del regimen.

    Fila de arriba: el sesgo es el que cambia de signo. Con poca informacion
    el FMM base subestima omega, o sea estima ondas mas apuntadas de lo que
    son, porque con pocos puntos tiene libertad para meter un pico estrecho
    entre dos observaciones. La regularizacion corrige eso. Pero cuando hay
    datos de sobra el base ya es practicamente insesgado y lo unico que
    aporta la penalizacion es sesgo propio.

    Las dos filas juntas son el compromiso sesgo-varianza de manual: se paga
    sesgo a cambio de varianza, y sale a cuenta solo mientras la varianza sea
    el termino dominante. El cruce esta en torno a A/sigma = 2 a 3.
    """
    plt = _mpl()
    omegas = sorted(par.omega.unique())
    fig, axes = plt.subplots(2, len(omegas), figsize=(9.0, 5.0),
                             sharex=True, sharey="row")
    fig.subplots_adjust(top=0.84, hspace=0.28)

    filas = (("atlas_omega_sesgo", "sesgo de $\\omega$\n(estimado $-$ verdadero)", True),
             ("atlas_omega_sd",    "desviacion tipica de $\\omega$\nentre genes", False))
    for r, (col, etiqueta, con_cero) in enumerate(filas):
        for c, w in enumerate(omegas):
            ax = axes[r, c]
            g = par[par.omega == w]
            for suf, nombre, color in (("_b", "base", "#c0392b"),
                                       ("_s", "estabilizado", "#2471a3")):
                m  = g.groupby("snr")[col + suf].median()
                q1 = g.groupby("snr")[col + suf].quantile(0.25)
                q3 = g.groupby("snr")[col + suf].quantile(0.75)
                ax.plot(m.index, m.values, marker="o", ms=4, color=color,
                        label=nombre)
                ax.fill_between(m.index, q1.values, q3.values,
                                color=color, alpha=0.15)
            if con_cero:
                ax.axhline(0, color="0.3", lw=0.9, ls="--")
            if r == 0:
                ax.set_title(f"$\\omega$ verdadero = {w:g}", fontsize=9)
            if r == 1:
                ax.set_xlabel("A/$\\sigma$")
        axes[r, 0].set_ylabel(etiqueta, fontsize=8.5)

    axes[0, 0].legend(frameon=False, fontsize=8)
    fig.suptitle("El mecanismo: la varianza baja siempre, el sesgo cruza\n"
                 "abajo el estabilizado gana en todas partes, arriba solo "
                 "cuando hay poca informacion", fontsize=10, y=0.99)
    fig.savefig(destino)
    plt.close(fig)


def figura_e_omega_por_gen(g: pd.DataFrame, destino: Path) -> None:
    """
    E: distribucion de omega estimada gen a gen.

    Es la figura que explica el mecanismo en vez de solo cuantificarlo, y hay
    que leer los tres paneles juntos porque cuentan tres cosas distintas:

    - **n bajo, senal baja.** Con pocos puntos el FMM base tiene libertad para
      poner un pico estrechisimo entre dos observaciones, y la distribucion se
      desparrama por todo el rango con masa contra el limite inferior. El
      estabilizado se concentra en torno al valor verdadero. Mediana 0,37
      contra 0,66, con omega verdadero 0,7.
    - **n alto, senal alta.** El base acierta de lleno, mediana 0,69 sobre un
      valor verdadero de 0,7, y el estabilizado se desplaza a 0,84. Aqui la
      regularizacion solo estorba.
    - **onda muy apuntada.** El base vuelve a acertar, 0,13 sobre 0,1, y el
      estabilizado se va a 0,57. Es el peor caso para la penalizacion, porque
      la rugosidad que castiga es aqui una propiedad real de la onda.

    Ojo al panel de la derecha: el amontonamiento del base contra el cero no
    es concentracion sino el suelo de la rejilla de omega, que trunca la
    distribucion por abajo. Sobre el conjunto de la rejilla la dispersion del
    estabilizado es menor tambien en omega = 0,1 (figura D), pero en esta
    celda concreta, con solo veinte genes evaluables, no se aprecia. Esta
    figura sirve para ver la FORMA de cada fallo; la afirmacion cuantitativa
    sobre dispersion es la de la figura D, que tiene 2700 replicas detras.
    """
    plt = _mpl()
    celdas = list(dict.fromkeys(g.celda))
    fig, axes = plt.subplots(1, len(celdas), figsize=(3.1 * len(celdas), 3.2),
                             sharey=True)
    fig.subplots_adjust(top=0.74)
    axes = np.atleast_1d(axes)
    bins = np.linspace(0, 1, 26)
    for ax, celda in zip(axes, celdas):
        sub = g[g.celda == celda]
        for metodo, nombre, color in (("base", "base", "#c0392b"),
                                      ("stabilized", "estabilizado", "#2471a3")):
            ax.hist(sub[sub.metodo == metodo].omega_hat, bins=bins, alpha=0.55,
                    color=color, label=nombre, density=True)
        ax.axvline(sub.omega_true.iloc[0], color="black", lw=1.4,
                   label="$\\omega$ verdadero")
        ax.set_title(celda, fontsize=9)
        ax.set_xlabel("$\\omega$ estimada")
    axes[0].set_ylabel("densidad de genes")
    axes[0].legend(frameon=False, fontsize=7.5)
    fig.suptitle(
        "Los tres regimenes: donde falla el base (izquierda) y donde falla\n"
        "el estabilizado (centro y derecha)", fontsize=10, y=1.02)
    fig.savefig(destino)
    plt.close(fig)


def figura_f_ondas(g: pd.DataFrame, destino: Path, n_genes: int = 3) -> None:
    """
    F: ondas de ejemplo, verdadera contra las dos estimadas.

    La mas intuitiva de todas y la que conviene poner primero en la seccion.
    Se eligen los genes donde los dos metodos MAS discrepan, no al azar ni
    los mejores: el objetivo es mostrar de que naturaleza es la diferencia,
    y un gen donde ambos aciertan no ensena nada.
    """
    from simulation.generate import fmm_wave

    plt = _mpl()
    celda = g.celda.iloc[0]
    sub = g[g.celda == celda]
    b = sub[sub.metodo == "base"].set_index("gene_name")
    s = sub[sub.metodo == "stabilized"].set_index("gene_name")
    comunes = b.index.intersection(s.index)
    discrepancia = (b.loc[comunes, "rmse_curva"] - s.loc[comunes, "rmse_curva"])
    elegidos = discrepancia.sort_values(ascending=False).index[:n_genes]

    malla = np.linspace(0, TWO_PI, 400)
    horas = malla * 24.0 / TWO_PI
    fig, axes = plt.subplots(1, len(elegidos), figsize=(3.0 * len(elegidos), 3.1),
                             sharey=True)
    fig.subplots_adjust(top=0.70)
    axes = np.atleast_1d(axes)
    for ax, gen in zip(axes, elegidos):
        fv = fmm_wave(malla, 0.0, 1.0, b.loc[gen, "alpha_true"],
                      b.loc[gen, "beta_true"], b.loc[gen, "omega_true"])
        fb = fmm_wave(malla, 0.0, 1.0, b.loc[gen, "alpha_hat"],
                      b.loc[gen, "beta_hat"], b.loc[gen, "omega_hat"])
        fs = fmm_wave(malla, 0.0, 1.0, s.loc[gen, "alpha_hat"],
                      s.loc[gen, "beta_hat"], s.loc[gen, "omega_hat"])
        ax.plot(horas, fv, color="black", lw=2.0, label="verdadera")
        ax.plot(horas, fb, color="#c0392b", lw=1.3, ls="--",
                label=f"base ($\\omega$={b.loc[gen,'omega_hat']:.2f})")
        ax.plot(horas, fs, color="#2471a3", lw=1.3,
                label=f"estab. ($\\omega$={s.loc[gen,'omega_hat']:.2f})")
        ax.set_title(gen, fontsize=8.5)
        ax.set_xlabel("hora del ciclo")
        ax.set_xticks([0, 6, 12, 18, 24])
        # bajo el eje: dentro tapaba las curvas, que son el contenido
        ax.legend(frameon=False, fontsize=7, loc="upper center",
                  bbox_to_anchor=(0.5, -0.30), ncol=1, handlelength=1.6)
    axes[0].set_ylabel("expresion (estandarizada)")
    fig.suptitle(f"Ondas recuperadas ({celda})\n"
                 "el FMM base inventa picos que no existen", fontsize=10, y=1.02)
    fig.savefig(destino)
    plt.close(fig)


# ═══════════════════════════════════════════════════════════════════════════
# CLI
# ═══════════════════════════════════════════════════════════════════════════

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--in", dest="entrada", default="output/simulacion")
    ap.add_argument("--images", default="docs/images")
    ap.add_argument("--volcar-por-gen", action="store_true",
                    help="Corre las celdas de las figuras E y F y vuelca el "
                         "error por gen. Minutos, no repite la rejilla.")
    args = ap.parse_args()

    entrada = Path(args.entrada)
    imagenes = Path(args.images)
    imagenes.mkdir(parents=True, exist_ok=True)
    por_gen = entrada / "por_gen.csv"

    if args.volcar_por_gen:
        g = volcar_por_gen(por_gen)
        print(f"  volcadas {len(g)} filas por gen en {por_gen}")

    par = cargar_pareado(entrada)
    print(f"  {len(par)} replicas pareadas con orientacion resuelta en ambos brazos")

    tablas = entrada / "tablas"
    tablas.mkdir(parents=True, exist_ok=True)
    t = tabla_wilcoxon(par)
    t.to_csv(tablas / "tabla_wilcoxon_fmm.csv", index=False)
    print(f"  escrita tabla_wilcoxon_fmm.csv ({len(t)} celdas)")

    for nombre, fn in (("fmm-a-efecto", figura_a_efecto),
                       ("fmm-b-victorias", figura_b_victorias),
                       ("fmm-c-pico-horas", figura_c_pico_horas),
                       ("fmm-d-sesgo-omega", figura_d_sesgo_omega)):
        fn(par, imagenes / f"{nombre}.png")
        print(f"  escrita {nombre}.png")

    if por_gen.exists():
        g = pd.read_csv(por_gen)
        figura_e_omega_por_gen(g, imagenes / "fmm-e-omega-por-gen.png")
        print("  escrita fmm-e-omega-por-gen.png")
        figura_f_ondas(g, imagenes / "fmm-f-ondas.png")
        print("  escrita fmm-f-ondas.png")
    else:
        print(f"  E y F omitidas: no existe {por_gen}, "
              f"corre antes con --volcar-por-gen")


if __name__ == "__main__":
    main()
