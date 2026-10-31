"""
Regla de seleccion en la rejilla de CircaFMM: argmax de J (articulo)
frente a la meseta del codigo publicado de CircaFMM.

Por que hace falta
------------------
La Seccion S3.4.3 del suplemento elige ``(alpha, omega)`` como el maximo de
``J = P - lambda g R``, y asi lo hace ``fit_gene``. El codigo publicado
(``fmm_fit_power_rough`` en R y su version C++, que es la que usan
``runCircaFMM.R`` y la calibracion de lambda) hace otra cosa:

1. toma la meseta ``{g : J_g >= max J - 0.01}``;
2. dentro de ella elige el punto que minimiza ``WRSS + lambda g R``;
3. si hay empate (tolerancia relativa 1e-10) prefiere el omega mayor.

Comparado directamente con el codigo R sobre 300 genes, el argmax coincide
siempre, pero la seleccion final cambia en el 93 %. Aqui se mide cual de las
dos reglas estima mejor. Cada regla calibra su propio ``lambda*`` con el mismo
Monte Carlo del articulo, porque en R la calibracion tambien usa la meseta.

Los genes se escalan a [-1, 1] como en ``preprocessing.py``. La tolerancia de
0,01 es absoluta sobre J, asi que su efecto depende de esa escala.

Uso
---
    python -m simulation.seleccion_meseta --out output/seleccion_meseta
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMBA_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

from simulation.rugosidad_dense_vs_obs import (
    DISENOS, HORAS, NS, OMEGAS, PROJECT_ROOT, SNRS, T_DENSE, TWO_PI,
    alpha_para_pico, diff_circular, fs, muestrear_tiempos, onda_fmm,
)
from simulation.escala_lambda import escalar

REGLAS = ("argmax", "meseta")
OMEGA0, POT = 0.10, 3.0
TOL_SCORE, TOL_OBJ = 1e-2, 1e-10
LAMBDAS = 2.0 ** np.linspace(-1, 4, 24)


# ── Ajuste vectorizado de B series con cualquiera de las dos reglas ─────────

def _potencia_y_rss(X, w, grid, cache, chunk=250):
    """``P`` (B, nG) y ``WRSS`` (B, nG) de la WLS en cada punto de la rejilla.

    Ninguna de las dos depende de lambda, asi que se calculan una vez.
    """
    B = X.shape[0]
    sw = w.sum()
    P = np.empty((B, grid.C.shape[0]))
    RSS = np.empty_like(P)
    for i0 in range(0, B, chunk):
        Xb = X[i0:i0 + chunk]
        mu = (Xb @ w) / sw
        wx0 = w[None, :] * (Xb - mu[:, None])
        U, V = wx0 @ grid.C.T, wx0 @ grid.S.T
        P[i0:i0 + chunk] = (U * U + V * V) / cache.den[None, :]
        wx = w[None, :] * Xb
        b = np.stack([wx.sum(1)[:, None].repeat(grid.C.shape[0], 1),
                      wx @ grid.C.T, wx @ grid.S.T], axis=2)          # (b, nG, 3)
        coef = np.einsum("gij,bgj->bgi", cache.inv_XWX, b)
        fit = (coef[:, :, :1] + coef[:, :, 1:2] * grid.C[None]
               + coef[:, :, 2:3] * grid.S[None])                      # (b, nG, n)
        RSS[i0:i0 + chunk] = (w[None, None, :] * (Xb[:, None, :] - fit) ** 2).sum(2)
    return P, RSS


def seleccionar(P, RSS, pen, omega_vec, regla):
    """Indice de rejilla elegido para cada fila."""
    J = P - pen[None, :]
    if regla == "argmax":
        return np.argmax(J, axis=1)
    meseta = J >= J.max(axis=1, keepdims=True) - TOL_SCORE
    obj = np.where(meseta, RSS + pen[None, :], np.inf)
    best = obj.min(axis=1, keepdims=True)
    empate = obj <= best + TOL_OBJ * np.maximum(np.abs(best), 1.0)
    return np.argmax(np.where(empate, omega_vec[None, :], -np.inf), axis=1)


def calibrar(t, w, grid, cache, regla, seed0, B=1000, R=5, q=0.975):
    g_om = (OMEGA0 / (grid.omega_vec + OMEGA0)) ** POT * grid.roughness
    Q = np.empty((R, LAMBDAS.shape[0]))
    for r in range(R):
        rng = np.random.default_rng(seed0 + r)
        X = rng.standard_normal((B, t.shape[0]))
        X = X - ((X @ w) / w.sum())[:, None]
        P, RSS = _potencia_y_rss(X, w, grid, cache)
        for j, lam in enumerate(LAMBDAS):
            idx = seleccionar(P, RSS, lam * g_om, grid.omega_vec, regla)
            C, S = grid.C[idx], grid.S[idx]
            b = np.stack([(X * w).sum(1), (X * w * C).sum(1), (X * w * S).sum(1)], 1)
            coef = np.einsum("bij,bj->bi", cache.inv_XWX[idx], b)
            fo = (coef[:, :1] + coef[:, 1:2] * C + coef[:, 2:3] * S)[:, grid.order]
            Q[r, j] = np.quantile(np.abs(np.diff(fo, axis=1, append=fo[:, :1])).max(1), q)
    return fs._pick_lambda_elbow(LAMBDAS, np.median(Q, axis=0), smooth_k=3)


def params_en(y, idx, grid, w, t):
    fit = fs._wls_fit_at(y, float(grid.alpha_vec[idx]), float(grid.omega_vec[idx]), w, t)
    pico = (fit["alpha"] + 2.0 * np.arctan2(np.sin(-fit["beta"] / 2.0) / fit["omega"],
                                            np.cos(-fit["beta"] / 2.0))) % TWO_PI
    ss = ((y - y.mean()) ** 2).sum()
    r2 = 1.0 - ((y - fit["fitted"]) ** 2).sum() / ss
    return fit, pico, r2


# ── Una extracción de tiempos ───────────────────────────────────────────────

def correr_extraccion(diseno, n, draw, n_rep, n_null, seed):
    rng = np.random.default_rng([seed, DISENOS.index(diseno), n, draw])
    t = muestrear_tiempos(diseno, n, rng)
    w = fs._time_density_weights(t)
    grid = fs.precompute_grid(t)
    cache = fs.prepare_cache(w, grid)
    g_om = (OMEGA0 / (grid.omega_vec + OMEGA0)) ** POT * grid.roughness
    lam = {r: calibrar(t, w, grid, cache, r, seed + 1000 * draw) for r in REGLAS}
    base = dict(diseno=diseno, n=n, draw=draw,
                lambda_argmax=lam["argmax"], lambda_meseta=lam["meseta"])

    # Señales
    meta, Y, verdad = [], [], []
    for omega in OMEGAS:
        for snr in SNRS:
            for r in range(n_rep):
                pico = rng.uniform(0.0, TWO_PI)
                beta = rng.uniform(0.0, TWO_PI)
                alpha = alpha_para_pico(pico, beta, omega)
                y0 = onda_fmm(t, 0.0, 1.0, alpha, beta, omega) + rng.standard_normal(n) / snr
                y, a, b = escalar(y0, "minmax")
                Y.append(y)
                meta.append(dict(rep=r, omega_true=omega, snr=snr, pico=pico,
                                 alpha=alpha, beta=beta, a=a, b=b))
    Y = np.array(Y)
    P, RSS = _potencia_y_rss(Y, w, grid, cache)

    filas = []
    for regla in REGLAS:
        pen = lam[regla] * g_om
        idx = seleccionar(P, RSS, pen, grid.omega_vec, regla)
        n_meseta = (P - pen >= (P - pen).max(1, keepdims=True) - TOL_SCORE).sum(1)
        for k, m in enumerate(meta):
            fit, pk, r2 = params_en(Y[k], idx[k], grid, w, t)
            verdad_den = m["a"] * onda_fmm(T_DENSE, 0.0, 1.0, m["alpha"], m["beta"],
                                           m["omega_true"]) + m["b"]
            curva = onda_fmm(T_DENSE, fit["M"], fit["A"], fit["alpha"], fit["beta"], fit["omega"])
            filas.append(dict(
                base, tipo="senal", regla=regla, rep=m["rep"],
                omega_true=m["omega_true"], snr=m["snr"], n_meseta=int(n_meseta[k]),
                err_pico_h=diff_circular(pk, m["pico"]) * HORAS,
                err_alpha=diff_circular(fit["alpha"], m["alpha"]),
                err_beta=diff_circular(fit["beta"], m["beta"]),
                omega_hat=fit["omega"], err_omega=fit["omega"] - m["omega_true"],
                err_A_rel=fit["A"] / m["a"] - 1.0,
                err_M=(fit["M"] - m["b"]) / m["a"],
                rmse_curva=float(np.sqrt(np.mean((curva - verdad_den) ** 2)) / m["a"]),
                r2=r2,
            ))

    # Ruido puro
    Z = np.array([escalar(rng.standard_normal(n), "minmax")[0] for _ in range(n_null)])
    P0, RSS0 = _potencia_y_rss(Z, w, grid, cache)
    for regla in REGLAS:
        idx = seleccionar(P0, RSS0, lam[regla] * g_om, grid.omega_vec, regla)
        for k in range(n_null):
            fit, _, r2 = params_en(Z[k], idx[k], grid, w, t)
            fo = fit["fitted"][grid.order]
            filas.append(dict(base, tipo="nulo", regla=regla, rep=k,
                              omega_hat=fit["omega"], r2=r2,
                              salto=float(np.abs(np.diff(fo, append=fo[:1])).max())))
    return filas


# ── Resúmenes ───────────────────────────────────────────────────────────────

def tabla(g):
    return pd.Series({
        "pico_h":      np.median(g.err_pico_h),
        "pico<2h":     np.mean(g.err_pico_h < 2.0),
        "|d_alpha|":   np.median(g.err_alpha),
        "|d_beta|":    np.median(g.err_beta),
        "|d_omega|":   np.median(np.abs(g.err_omega)),
        "sesgo_omega": np.median(g.err_omega),
        "|dA|/A":      np.median(np.abs(g.err_A_rel)),
        "sesgo_A":     np.median(g.err_A_rel),
        "|dM|/A":      np.median(np.abs(g.err_M)),
        "rmse/A":      np.median(g.rmse_curva),
        "R2":          np.median(g.r2),
    })


def resumir(df):
    pd.set_option("display.width", 220)
    pd.set_option("display.max_columns", 30)
    s = df[df.tipo == "senal"]
    print("\n=== lambda* (mediana por diseno y n) ===")
    print(df.drop_duplicates(["diseno", "n", "draw"]).groupby(["diseno", "n"])[
        ["lambda_argmax", "lambda_meseta"]].median().round(3))
    print("\n=== Senal: todo agregado ===")
    print(s.groupby("regla").apply(tabla, include_groups=False).round(3))
    print("\n=== Senal: por omega verdadera ===")
    print(s.groupby(["omega_true", "regla"]).apply(tabla, include_groups=False).round(3))
    print("\n=== Senal: por diseno y n ===")
    print(s.groupby(["diseno", "n", "regla"]).apply(tabla, include_groups=False).round(3))
    print("\n=== Senal: por SNR ===")
    print(s.groupby(["snr", "regla"]).apply(tabla, include_groups=False).round(3))

    clave = ["diseno", "n", "draw", "omega_true", "snr", "rep"]
    w = s.pivot_table(index=clave, columns="regla", values=["rmse_curva", "err_pico_h", "omega_hat"])
    e = 1e-9
    print("\n=== Pareado: fraccion de victorias ===")
    print(pd.Series({
        "misma omega":         np.mean(np.isclose(w.omega_hat.argmax, w.omega_hat.meseta)),
        "rmse: gana meseta":   np.mean(w.rmse_curva.meseta < w.rmse_curva.argmax - e),
        "rmse: gana argmax":   np.mean(w.rmse_curva.argmax < w.rmse_curva.meseta - e),
        "pico: gana meseta":   np.mean(w.err_pico_h.meseta < w.err_pico_h.argmax - e),
        "pico: gana argmax":   np.mean(w.err_pico_h.argmax < w.err_pico_h.meseta - e),
    }).round(3))
    print("\nTamano de la meseta (senal): mediana",
          s[s.regla == "meseta"].n_meseta.median(), " p90",
          s[s.regla == "meseta"].n_meseta.quantile(0.9))

    z = df[df.tipo == "nulo"]
    print("\n=== H0 (ruido puro, escala [-1, 1]) ===")
    print(z.groupby("regla").apply(lambda g: pd.Series({
        "omega_hat_med":   np.median(g.omega_hat),
        "frac_omega<0.1":  np.mean(g.omega_hat < 0.1),
        "frac_omega<0.25": np.mean(g.omega_hat < 0.25),
        "salto_med":       np.median(g.salto),
        "salto_p975":      np.quantile(g.salto, 0.975),
        "R2_med":          np.median(g.r2),
    }), include_groups=False).round(3))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--out", default="output/seleccion_meseta")
    ap.add_argument("--draws", type=int, default=3)
    ap.add_argument("--reps", type=int, default=100)
    ap.add_argument("--nulls", type=int, default=200)
    ap.add_argument("--jobs", type=int, default=12)
    ap.add_argument("--seed", type=int, default=20260925)
    ap.add_argument("--solo-resumen", action="store_true")
    args = ap.parse_args()
    out = PROJECT_ROOT / args.out
    out.mkdir(parents=True, exist_ok=True)
    csv = out / "resultados.csv"
    if not args.solo_resumen:
        tareas = [(d, n, k) for d in DISENOS for n in NS for k in range(args.draws)]
        res = Parallel(n_jobs=args.jobs, verbose=5)(
            delayed(correr_extraccion)(d, n, k, args.reps, args.nulls, args.seed)
            for d, n, k in tareas)
        pd.DataFrame([f for filas in res for f in filas]).to_csv(csv, index=False)
    resumir(pd.read_csv(csv))


if __name__ == "__main__":
    main()
