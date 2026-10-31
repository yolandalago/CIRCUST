"""
Ejemplos de datos simulados con el ajuste del FMM base y de CircaFMM.

Toma matrices de la rejilla (las mismas que corrio la simulacion, por su
semilla), elige unos genes ritmicos y dibuja sus puntos frente al tiempo
verdadero junto a la onda verdadera y los dos ajustes. Los ajustes se hacen
sobre los tiempos verdaderos, para ver solo como ajusta cada modelo, sin el
error del orden.

Los datos se muestran normalizados a [-1, 1] gen a gen, que es lo que recibe
el FMM dentro del pipeline. La onda verdadera pasa por la misma transformacion
afin, asi que es comparable con los puntos y con los ajustes.

Uso
---
    python -m simulation.ejemplos_ajuste --omega 0.1 --snr 2 --n 48 --genes 3
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from circust.fitting.circafmm import CircaFMMModel, prepare_dataset
from circust.fitting.fmm import FMMModel
from simulation.generate import fmm_wave, generate_dataset
from simulation.grid import build_grid
from simulation.truth import TWO_PI

MUESTREOS = ("uniforme", "agrupado", "hueco")
NOMBRE_MUESTREO = {"uniforme": "Muestreo uniforme",
                   "agrupado": "Muestreo agrupado (4, 12 y 20 h)",
                   "hueco": "Muestreo con hueco (10-16 h)"}
ACENTO, NARANJA, GRIS, TEXTO = "#2a6fdb", "#d9822b", "#9a9a94", "#222222"


def _coma(v, dec=2) -> str:
    return f"{v:.{dec}f}".replace(".", ",")


def _celda(muestreo: str, omega: float, snr: float, n: int, rep: int):
    """La configuracion de la rejilla, con su semilla, para estos factores."""
    for c in build_grid("rejilla", rep + 1):
        k = c.cfg
        if (c.rep == rep and k.sampling == muestreo and k.omega == omega
                and k.snr == snr and k.n_samples == n
                and k.n_genes == 1000 and k.prop_rhythmic == 0.20):
            return k
    raise ValueError("No hay esa celda en la rejilla")


def _normalizar(x: np.ndarray, lo: float, hi: float) -> np.ndarray:
    return 2.0 * (x - lo) / (hi - lo) - 1.0


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    p.add_argument("--omega", type=float, default=0.1)
    p.add_argument("--snr", type=float, default=2.0)
    p.add_argument("--n", type=int, default=48)
    p.add_argument("--genes", type=int, default=3)
    p.add_argument("--rep", type=int, default=0)
    p.add_argument("--semilla", type=int, default=7)
    p.add_argument("--out", default="output/simulacion_2026-10/figuras")
    a = p.parse_args()

    malla = np.linspace(0.0, TWO_PI, 400)
    horas = malla * 24.0 / TWO_PI
    fig, axes = plt.subplots(len(MUESTREOS), a.genes, figsize=(3.4 * a.genes, 2.9 * len(MUESTREOS)),
                             sharex=True, sharey=True)
    fig.subplots_adjust(top=0.84, left=0.06, right=0.99, bottom=0.07, hspace=0.42, wspace=0.08)
    rng = np.random.default_rng(a.semilla)

    for fila, muestreo in zip(axes, MUESTREOS):
        cfg = _celda(muestreo, a.omega, a.snr, a.n, a.rep)
        matriz, verdad = generate_dataset(cfg)
        theta, t = verdad.theta_true, verdad.t_true
        base = FMMModel()
        circa = CircaFMMModel(prepared=prepare_dataset(theta))
        params = verdad.gene_params
        ritmicos = params.index[params["kind"] == "ritmico"].to_numpy()
        for ax, gen in zip(fila, rng.choice(ritmicos, a.genes, replace=False)):
            x = matriz.loc[gen].to_numpy(float)
            lo, hi = x.min(), x.max()
            y = _normalizar(x, lo, hi)
            g = params.loc[gen]
            onda = _normalizar(fmm_wave(malla, 0.0, cfg.amplitude, g.alpha, g.beta, g.omega), lo, hi)
            fb, fc = base.fit(y, theta), circa.fit(y, theta)
            curva = lambda f: fmm_wave(malla, f.params["M"], f.params["A"], f.params["alpha"],
                                       f.params["beta"], f.params["omega"])
            ax.scatter(t, y, s=14, color=GRIS, alpha=0.85, zorder=2, label="datos")
            ax.plot(horas, onda, color=TEXTO, lw=1.6, ls="--", zorder=3, label="onda verdadera")
            ax.plot(horas, curva(fb), color=NARANJA, lw=1.8, zorder=4, label="FMM base")
            ax.plot(horas, curva(fc), color=ACENTO, lw=1.8, zorder=5, label="CircaFMM")
            ax.text(0.02, 0.03,
                    f"ω base = {_coma(fb.params['omega'])} (R² {_coma(fb.r2)})\n"
                    f"ω CircaFMM = {_coma(fc.params['omega'])} (R² {_coma(fc.r2)})",
                    transform=ax.transAxes, fontsize=7.5, color=TEXTO, va="bottom",
                    bbox=dict(facecolor="white", edgecolor="none", alpha=0.8, pad=1.5))
            ax.set_title(gen, fontsize=9, color=TEXTO)
            ax.set_xlim(0, 24)
            ax.set_xticks([0, 6, 12, 18, 24])
            ax.set_ylim(-2.1, 2.1)
            for lado in ("top", "right"):
                ax.spines[lado].set_visible(False)
            ax.grid(color="#ececea", lw=0.7)
            ax.set_axisbelow(True)
        fila[0].set_ylabel("expresión normalizada")
        fila[0].annotate(NOMBRE_MUESTREO[muestreo], (0, 1.18), xycoords="axes fraction",
                         fontsize=10.5, color=TEXTO, fontweight="bold")
    for ax in axes[-1]:
        ax.set_xlabel("tiempo verdadero (h)")
    manejadores, etiquetas = axes[0, 0].get_legend_handles_labels()
    fig.legend(manejadores, etiquetas, frameon=False, fontsize=8.5, ncol=4,
               loc="upper left", bbox_to_anchor=(0.055, 0.94))
    fig.text(0.06, 0.985,
             f"Genes simulados con ω = {_coma(a.omega, 1)}, A/σ = {a.snr:g} y n = {a.n}: "
             "datos, onda verdadera y ajustes", fontsize=12, color=TEXTO, va="top")
    fig.text(0.06, 0.955,
             "Matrices de la rejilla (réplica 0, 1000 genes, 20 % rítmicos). Ajustes sobre los "
             "tiempos verdaderos y los datos normalizados a [−1, 1].",
             fontsize=8.5, color="#666666", va="top")
    destino = Path(a.out)
    destino.mkdir(parents=True, exist_ok=True)
    nombre = destino / f"ejemplos_omega{str(a.omega).replace('.', '')}_snr{a.snr:g}_n{a.n}.png"
    fig.savefig(nombre, dpi=160)
    print(nombre)


if __name__ == "__main__":
    main()
