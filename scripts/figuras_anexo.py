#!/usr/bin/env python3
"""
figuras_anexo.py
================
Figuras explicativas de las secciones A.8 y A.9 del Anexo de la memoria:

  - anexo-spline.pdf : base de B-splines periodicos, ajuste a un gen core y
                       residuos estandarizados (deteccion de muestras atipicas).
  - anexo-seno.pdf   : por que la comparacion lineal de instantes falla sobre
                       la circunferencia y como lo resuelve el seno de la diferencia.

Usa el propio codigo del detector (``circust.fitting.cyclic_spline``), de modo
que la figura refleja exactamente lo que hace la implementacion.

    PYTHONPATH=src python3 scripts/figuras_anexo.py [--salida docs/images]
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from circust.fitting.cyclic_spline import (
    cyclic_spline_basis, fit_cyclic_spline_residuals,
)

# Paleta categorica de referencia (modo claro), en orden fijo.
AZUL, NARANJA, AGUA, AMARILLO = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
TINTA, TINTA_2, GRIS = "#0b0b0b", "#52514e", "#c9c8c3"

plt.rcParams.update({
    "font.size": 7.5, "axes.titlesize": 8, "axes.labelsize": 7.5,
    "xtick.labelsize": 7, "ytick.labelsize": 7,
    "font.family": "serif", "mathtext.fontset": "cm",
    "axes.edgecolor": TINTA_2, "axes.labelcolor": TINTA, "text.color": TINTA,
    "xtick.color": TINTA_2, "ytick.color": TINTA_2,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": 0.6, "figure.facecolor": "white", "axes.facecolor": "white",
})

DOS_PI = 2 * np.pi
PI_TICKS = [0, np.pi, DOS_PI]
PI_LABELS = ["0", r"$\pi$", r"$2\pi$"]


def figura_spline(n_knots: int = 5) -> plt.Figure:
    rng = np.random.default_rng(7)
    fig, axs = plt.subplots(1, 3, figsize=(6.5, 2.2), constrained_layout=True)

    # (a) Base de B-splines periodicos con nodos uniformes en todo el ciclo
    ax = axs[0]
    t = np.linspace(0.0, DOS_PI, 600)
    B = cyclic_spline_basis(t, n_knots)
    knots = np.linspace(0.0, DOS_PI, n_knots)
    for k in knots:
        ax.axvline(k, color=GRIS, lw=0.6, ls=":", zorder=0)
    colores = [AZUL, NARANJA, AGUA, AMARILLO]
    for j in range(B.shape[1]):
        ax.plot(t, B[:, j], color=colores[j], lw=1.5)
        jmax = int(np.argmax(B[:, j]))
        ax.text(t[jmax], B[jmax, j] + 0.04, rf"$B_{{{j + 1}}}$", ha="center",
                va="bottom", fontsize=7, color=TINTA)
    ax.plot(t, B.sum(axis=1), color=TINTA_2, lw=0.8, ls="--")
    ax.text(DOS_PI * 0.30, 1.03, r"$\sum_j B_j=1$", fontsize=6.5, color=TINTA_2,
            va="bottom", ha="left")
    ax.set_ylim(-0.03, 1.2)
    ax.set_xticks(knots)
    ax.set_xticklabels(["0", r"$\pi/2$", r"$\pi$", r"$3\pi/2$", r"$2\pi$"])
    ax.set_xlabel(r"Fase $\varphi$ (rad)")
    ax.set_title(r"(a) Base de B-splines ($m=5$ nodos)", loc="left")

    # Gen sintetico con perfil asimetrico y una muestra atipica
    n = 36
    phi = np.sort(rng.uniform(0.0, DOS_PI, n))
    omega, alpha, beta = 0.55, 1.2, 0.0
    señal = 0.7 * np.cos(beta + 2 * np.arctan(omega * np.tan((phi - alpha) / 2)))
    y = señal + rng.normal(0.0, 0.10, n)
    i_out = int(np.argmin(np.abs(phi - 4.2)))
    y[i_out] += 1.25

    basis = cyclic_spline_basis(phi, n_knots)
    coef, *_ = np.linalg.lstsq(basis, y, rcond=None)
    ajuste = cyclic_spline_basis(t, n_knots) @ coef
    r = fit_cyclic_spline_residuals(y, phi, n_knots)
    assert abs(r[i_out]) > 3.5 and np.all(np.abs(np.delete(r, i_out)) < 3.5), r

    # (b) Ajuste del spline a un gen core
    ax = axs[1]
    mask = np.ones(n, dtype=bool)
    mask[i_out] = False
    ax.plot(t, ajuste, color=AZUL, lw=1.6, zorder=2)
    ax.scatter(phi[mask], y[mask], s=9, color=TINTA_2, zorder=3)
    ax.scatter(phi[i_out], y[i_out], s=26, color=NARANJA, edgecolor="white",
               linewidth=0.8, zorder=4)
    ax.annotate("muestra atípica", (phi[i_out], y[i_out]),
                xytext=(phi[i_out] + 0.25, y[i_out] + 0.32), fontsize=6.5,
                color=TINTA, arrowprops=dict(arrowstyle="-", color=TINTA_2, lw=0.6))
    ax.text(0.15, ajuste.min() - 0.05, "ajuste spline", fontsize=6.5, color=AZUL,
            va="top")
    ax.set_xticks(PI_TICKS)
    ax.set_xticklabels(PI_LABELS)
    ax.set_xlabel(r"Fase $\varphi$ (rad)")
    ax.set_ylabel("Expresión normalizada")
    ax.set_title("(b) Ajuste a un gen core", loc="left")

    # (c) Residuos estandarizados y umbral
    ax = axs[2]
    for u in (3.5, -3.5):
        ax.axhline(u, color=TINTA_2, lw=0.7, ls="--")
    ax.axhline(0.0, color=GRIS, lw=0.6)
    ax.scatter(phi[mask], r[mask], s=9, color=TINTA_2, zorder=3)
    ax.scatter(phi[i_out], r[i_out], s=26, color=NARANJA, edgecolor="white",
               linewidth=0.8, zorder=4)
    ax.text(0.1, 3.6, r"$|r|=3{,}5$", fontsize=6.5, color=TINTA_2, va="bottom")
    ax.set_ylim(-5.5, 6.2)
    ax.set_xticks(PI_TICKS)
    ax.set_xticklabels(PI_LABELS)
    ax.set_xlabel(r"Fase $\varphi$ (rad)")
    ax.set_ylabel(r"Residuo estandarizado $r$")
    ax.set_title("(c) Candidata si $|r|>3{,}5$", loc="left")
    return fig


def figura_seno() -> plt.Figure:
    fig, axs = plt.subplots(1, 2, figsize=(6.5, 2.7), constrained_layout=True,
                            gridspec_kw={"width_ratios": [1.0, 1.35]})

    # (a) Reloj de 24 h: el corte de la escala separa 23 h y 1 h
    ax = axs[0]
    ax.set_aspect("equal")
    ax.axis("off")
    th = np.linspace(0, DOS_PI, 400)
    ax.plot(np.sin(th), np.cos(th), color=GRIS, lw=1.0)

    def xy(h, r=1.0):
        a = DOS_PI * h / 24.0
        return r * np.sin(a), r * np.cos(a)

    for h, lab in [(0, "0 h"), (6, "6 h"), (12, "12 h"), (18, "18 h")]:
        x0, y0 = xy(h, 1.0)
        x1, y1 = xy(h, 1.07)
        ax.plot([x0, x1], [y0, y1], color=TINTA_2, lw=0.7)
        xt, yt = xy(h, 1.2)
        ax.text(xt, yt, lab, ha="center", va="center", fontsize=6.5, color=TINTA_2)
    ax.plot([0, 0], [0.0, 1.0], color=TINTA_2, lw=0.8, ls="--")
    ax.text(0.04, 0.08, "corte de\nla escala", fontsize=6, color=TINTA_2,
            ha="left", va="bottom")
    # arco de g1 (23 h) a g2 (1 h) a traves del corte
    arco = np.linspace(23.0, 25.0, 60)
    ax.plot(*zip(*[xy(h, 0.8) for h in arco]), color=TINTA, lw=1.0)
    ax.annotate("", xy=xy(25.0, 0.8), xytext=xy(24.85, 0.8),
                arrowprops=dict(arrowstyle="->", color=TINTA, lw=1.0))
    g1, g2 = xy(23.0), xy(1.0)
    ax.scatter(*g1, s=34, color=AZUL, zorder=3, edgecolor="white", linewidth=0.8)
    ax.scatter(*g2, s=34, color=NARANJA, zorder=3, edgecolor="white", linewidth=0.8)
    ax.text(g1[0] - 0.14, g1[1] + 0.02, r"$g_1$ (23 h)", ha="right", fontsize=6.8)
    ax.text(g2[0] + 0.14, g2[1] + 0.02, r"$g_2$ (1 h)", ha="left", fontsize=6.8)
    ax.text(0.0, 0.52, r"$\delta=2$ h", ha="center", fontsize=7)
    ax.text(0.0, -0.22, r"$23<1$ es falso" + "\n" + r"pero $\sin\delta>0$",
            ha="center", va="top", fontsize=7)
    ax.set_xlim(-1.45, 1.45)
    ax.set_ylim(-1.3, 1.35)
    ax.set_title("(a) «Después» sobre la circunferencia", loc="left", y=0.99)

    # (b) El seno de la diferencia
    ax = axs[1]
    d = np.linspace(0.0, 24.0, 600)
    s = np.sin(DOS_PI * d / 24.0)
    ax.fill_between(d, 0, s, where=s >= 0, color=AZUL, alpha=0.14, lw=0)
    ax.fill_between(d, 0, s, where=s <= 0, color=NARANJA, alpha=0.14, lw=0)
    ax.axhline(0, color=TINTA_2, lw=0.6)
    ax.plot(d, -s, color=TINTA_2, lw=1.0, ls="--")
    ax.plot(d, s, color=AZUL, lw=1.8)
    ax.text(0.5, 1.0, r"$\sin\delta$", color=AZUL, fontsize=8, va="bottom",
            ha="left")
    ax.text(0.5, -1.02, r"tras invertir: $\sin(-\delta)$", color=TINTA_2,
            fontsize=6.5, va="top")
    ax.text(6, 0.45, r"$g_2$ por delante de $g_1$", ha="center", fontsize=6.8,
            color=TINTA)
    ax.text(18, -0.45, r"$g_2$ por detrás de $g_1$", ha="center", fontsize=6.8,
            color=TINTA)
    ax.scatter([6], [1.0], s=22, color=AZUL, zorder=4)
    ax.annotate("máxima información\n(un cuarto de ciclo)", (6, 1.0),
                xytext=(9.2, 0.98), fontsize=6.3, color=TINTA, va="center",
                arrowprops=dict(arrowstyle="-", color=TINTA_2, lw=0.6))
    ax.scatter([12], [0.0], s=22, color=NARANJA, zorder=4)
    ax.annotate("antifase: sin información", (12, 0.0), xytext=(12.8, 0.28),
                fontsize=6.3, color=TINTA,
                arrowprops=dict(arrowstyle="-", color=TINTA_2, lw=0.6))
    ax.scatter([2], [np.sin(DOS_PI * 2 / 24)], s=22, color=NARANJA, zorder=4,
               edgecolor="white", linewidth=0.6)
    ax.annotate(r"ejemplo de (a): $\delta=2$ h", (2, np.sin(DOS_PI * 2 / 24)),
                xytext=(0.4, 0.1), fontsize=6.3, color=TINTA,
                arrowprops=dict(arrowstyle="-", color=TINTA_2, lw=0.6))
    ax.set_xlim(0, 24)
    ax.set_ylim(-1.25, 1.3)
    ax.set_xticks([0, 6, 12, 18, 24])
    ax.set_xlabel(r"Diferencia $\delta=\hat t_U^{(g_2)}-\hat t_U^{(g_1)}$ (h)")
    ax.set_ylabel(r"$\sin\delta$")
    ax.set_title("(b) Aportación de un par a la puntuación", loc="left")
    return fig


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--salida", default="docs/images")
    args = ap.parse_args()
    out = Path(args.salida)
    out.mkdir(parents=True, exist_ok=True)
    for nombre, fig in [("anexo-spline", figura_spline()),
                        ("anexo-seno", figura_seno())]:
        fig.savefig(out / f"{nombre}.pdf")
        fig.savefig(out / f"{nombre}.png", dpi=160)
        plt.close(fig)
        print("guardada:", out / f"{nombre}.pdf")


if __name__ == "__main__":
    main()
