"""
Efecto de la escala del gen sobre la penalizacion de rugosidad de CircaFMM, y de corregirla dividiendo P por la varianza del gen.

Por que hace falta
------------------
CircaFMM calibra ``lambda*`` con ruido nulo de varianza 1 (Algoritmo S3,
"sigma = 1 without loss of generality"), pero luego ajusta los genes sin
estandarizar. La potencia ``P`` crece con la varianza del gen y la
penalizacion ``lambda g(omega) R`` no depende de los datos, asi que la fuerza
efectiva de la penalizacion depende de la escala del gen. En CIRCUST los genes
estan escalados a [-1, 1], con varianzas muy por debajo de 1.

Variantes
---------
``original``
    Como el articulo y ``circafmm``: ``J = P - lambda g R``.
``corregido``
    ``J = P / s2_w - lambda g R``, con ``s2_w`` la varianza ponderada del gen.
    Equivale a estandarizar el gen solo para elegir (alpha, omega). La WLS
    final se hace sobre los datos tal cual. La calibracion de ``lambda`` usa
    el mismo criterio corregido, para que calibracion y ajuste coincidan.

Como ``argmax(P / s2 - lambda g R) = argmax(P - s2 lambda g R)``, el ajuste
corregido se obtiene con ``fit_gene`` pasando la penalizacion multiplicada
por ``s2_w``: no hay que tocar el modulo.

Escalas del gen
---------------
``minmax``  la del pipeline: cada gen se lleva a [-1, 1] como en
            ``preprocessing.py``.
``crudo``   A = 1, sigma = 1 / SNR.
``x10``     lo mismo multiplicado por 10.

Con la correccion el ajuste deberia ser identico en las tres escalas. Los
errores de A, M y la curva se miden en unidades de la A verdadera en esa
escala.

Uso
---
    python -m simulation.escala_lambda --out output/escala_lambda
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMBA_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse
from dataclasses import replace
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

from simulation.rugosidad_dense_vs_obs import (
    DISENOS, HORAS, NS, OMEGAS, PROJECT_ROOT, SNRS, T_DENSE, TWO_PI,
    alpha_para_pico, diff_circular, fs, muestrear_tiempos, onda_fmm,
)

ESCALAS    = ("minmax", "crudo", "x10")
VARIANTES  = ("original", "corregido")
OMEGA0, POT = 0.10, 3.0


def var_ponderada(y: np.ndarray, w: np.ndarray) -> float:
    mu = (w * y).sum() / w.sum()
    return float((w * (y - mu) ** 2).sum() / w.sum())


def escalar(y: np.ndarray, escala: str) -> tuple[np.ndarray, float, float]:
    """Devuelve ``(y', a, b)`` con ``y' = a*y + b``."""
    if escala == "minmax":
        lo, hi = y.min(), y.max()
        a = 2.0 / (hi - lo)
        return a * y - a * lo - 1.0, a, -a * lo - 1.0
    a = 1.0 if escala == "crudo" else 10.0
    return a * y, a, 0.0


# ── Calibración de lambda, original o corregida ─────────────────────────────

def _saltos(Xnull, w, grid, cache, lam, corregido):
    """Copia de ``fs._simulate_spikes_vectorized`` con el criterio corregido opcional."""
    sw = float(w.sum())
    mu_b = (Xnull @ w) / sw
    Xc = Xnull - mu_b[:, None]
    wx0 = w[None, :] * Xc
    U = wx0 @ grid.C.T
    V = wx0 @ grid.S.T
    pow_ = (U * U + V * V) / cache.den[None, :]
    if corregido:
        s2 = (w[None, :] * Xc * Xc).sum(axis=1) / sw
        pow_ = pow_ / s2[:, None]
    g_om = (OMEGA0 / (grid.omega_vec + OMEGA0)) ** POT
    score = pow_ - lam * (g_om * grid.roughness)[None, :]
    idx = np.argmax(score, axis=1)
    C_at, S_at = grid.C[idx], grid.S[idx]
    b = np.stack([(Xnull * w).sum(1), (Xnull * w * C_at).sum(1),
                  (Xnull * w * S_at).sum(1)], axis=1)
    coef = np.einsum("bij,bj->bi", cache.inv_XWX[idx], b)
    fit = coef[:, :1] + coef[:, 1:2] * C_at + coef[:, 2:3] * S_at
    fo = fit[:, grid.order]
    return np.abs(np.diff(fo, axis=1, append=fo[:, :1])).max(axis=1)


def calibrar(t, w, grid, cache, corregido, seed0, B=1000, R=5, q=0.975):
    lam_grid = 2.0 ** np.linspace(-1, 4, 24)
    Q = []
    for r in range(R):
        rng = np.random.default_rng(seed0 + r)
        X = rng.standard_normal((B, t.shape[0]))
        X = X - ((X * w).sum(1) / w.sum())[:, None]
        Q.append([np.quantile(_saltos(X, w, grid, cache, lam, corregido), q)
                  for lam in lam_grid])
    return fs._pick_lambda_elbow(lam_grid, np.median(np.array(Q), axis=0), smooth_k=3)


# ── Una extracción de tiempos ───────────────────────────────────────────────

def ajustar(y, prep, corregido):
    if corregido:
        prep = replace(prep, penalty=prep.penalty * var_ponderada(y, prep.weights))
    return fs.fit_gene(y, prep)


def correr_extraccion(diseno, n, draw, n_rep, n_null, seed):
    rng = np.random.default_rng([seed, DISENOS.index(diseno), n, draw])
    t = muestrear_tiempos(diseno, n, rng)
    w = fs._time_density_weights(t)
    grid = fs.precompute_grid(t)
    cache = fs.prepare_cache(w, grid)

    preps = {}
    for var in VARIANTES:
        lam = calibrar(t, w, grid, cache, var == "corregido", seed + 1000 * draw)
        preps[var] = fs.FMMPreparedDataset(
            times=t, weights=w, grid=grid, cache=cache, lambda_star=float(lam),
            omega0=OMEGA0, penalty_power=POT,
            penalty=fs._penalty_vector(grid, float(lam), OMEGA0, POT))

    base = dict(diseno=diseno, n=n, draw=draw,
                lambda_original=preps["original"].lambda_star,
                lambda_corregido=preps["corregido"].lambda_star)
    filas = []

    for omega in OMEGAS:
        for snr in SNRS:
            for r in range(n_rep):
                pico = rng.uniform(0.0, TWO_PI)
                beta = rng.uniform(0.0, TWO_PI)
                alpha = alpha_para_pico(pico, beta, omega)
                y0 = onda_fmm(t, 0.0, 1.0, alpha, beta, omega) \
                     + rng.standard_normal(n) / snr
                verdad_den = onda_fmm(T_DENSE, 0.0, 1.0, alpha, beta, omega)
                for esc in ESCALAS:
                    y, a, b = escalar(y0, esc)
                    for var, prep in preps.items():
                        fit = ajustar(y, prep, var == "corregido")
                        p = fit.params
                        curva = onda_fmm(T_DENSE, p["M"], p["A"], p["alpha"],
                                         p["beta"], p["omega"])
                        err = (curva - (a * verdad_den + b)) / a
                        filas.append(dict(
                            base, tipo="senal", escala=esc, variante=var, rep=r,
                            omega_true=omega, snr=snr,
                            err_pico_h=diff_circular(fit.peak_time, pico) * HORAS,
                            err_alpha=diff_circular(p["alpha"], alpha),
                            err_beta=diff_circular(p["beta"], beta),
                            omega_hat=p["omega"],
                            err_omega=p["omega"] - omega,
                            err_A_rel=p["A"] / a - 1.0,
                            err_M=(p["M"] - b) / a,
                            rmse_curva=float(np.sqrt(np.mean(err ** 2))),
                            r2=fit.r2,
                        ))

    for r in range(n_null):
        y0 = rng.standard_normal(n)
        for esc in ESCALAS:
            y, a, _ = escalar(y0, esc)
            for var, prep in preps.items():
                fit = ajustar(y, prep, var == "corregido")
                fo = fit.fitted[grid.order]
                salto = np.abs(np.diff(fo, append=fo[:1])).max()
                filas.append(dict(
                    base, tipo="nulo", escala=esc, variante=var, rep=r,
                    omega_hat=fit.params["omega"],
                    salto_rel=float(salto / np.sqrt(var_ponderada(y, w))),
                    r2=fit.r2,
                ))
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


def resumir(df: pd.DataFrame) -> None:
    pd.set_option("display.width", 220)
    pd.set_option("display.max_columns", 30)
    s = df[df.tipo == "senal"]

    print("\n=== lambda* (mediana por diseno y n) ===")
    print(df.drop_duplicates(["diseno", "n", "draw"]).groupby(["diseno", "n"])[
        ["lambda_original", "lambda_corregido"]].median().round(3))

    print("\n=== Senal: por escala y variante (todo agregado) ===")
    print(s.groupby(["escala", "variante"]).apply(tabla, include_groups=False).round(3))

    print("\n=== Senal, escala minmax (la de CIRCUST): por omega verdadera ===")
    sm = s[s.escala == "minmax"]
    print(sm.groupby(["omega_true", "variante"]).apply(tabla, include_groups=False).round(3))

    print("\n=== Senal, escala minmax: por diseno y n ===")
    print(sm.groupby(["diseno", "n", "variante"]).apply(tabla, include_groups=False).round(3))

    print("\n=== Senal, escala minmax: por SNR ===")
    print(sm.groupby(["snr", "variante"]).apply(tabla, include_groups=False).round(3))

    clave = ["diseno", "n", "draw", "omega_true", "snr", "rep"]
    w = sm.pivot_table(index=clave, columns="variante",
                       values=["rmse_curva", "err_pico_h"])
    eps = 1e-9
    print("\n=== Pareado en minmax: fraccion de victorias ===")
    print(pd.Series({
        "rmse: gana corregido": np.mean(w.rmse_curva.corregido < w.rmse_curva.original - eps),
        "rmse: gana original":  np.mean(w.rmse_curva.original < w.rmse_curva.corregido - eps),
        "pico: gana corregido": np.mean(w.err_pico_h.corregido < w.err_pico_h.original - eps),
        "pico: gana original":  np.mean(w.err_pico_h.original < w.err_pico_h.corregido - eps),
        "misma omega":          np.mean(np.isclose(
            sm[sm.variante == "original"].omega_hat.values,
            sm[sm.variante == "corregido"].omega_hat.values)),
    }).round(3))

    print("\n=== Invariancia: omega_hat identica entre escalas (misma senal) ===")
    for var in VARIANTES:
        piv = s[s.variante == var].pivot_table(index=clave, columns="escala",
                                               values="omega_hat")
        print(f"{var:10s} minmax=crudo: {np.mean(np.isclose(piv.minmax, piv.crudo)):.3f}"
              f"   crudo=x10: {np.mean(np.isclose(piv.crudo, piv.x10)):.3f}")

    z = df[df.tipo == "nulo"]
    print("\n=== H0 (ruido puro): por escala y variante ===")
    print(z.groupby(["escala", "variante"]).apply(lambda g: pd.Series({
        "omega_hat_med":   np.median(g.omega_hat),
        "frac_omega<0.1":  np.mean(g.omega_hat < 0.1),
        "frac_omega<0.25": np.mean(g.omega_hat < 0.25),
        "salto/sd_med":    np.median(g.salto_rel),
        "salto/sd_p975":   np.quantile(g.salto_rel, 0.975),
        "R2_med":          np.median(g.r2),
    }), include_groups=False).round(3))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--out", default="output/escala_lambda")
    ap.add_argument("--draws", type=int, default=4)
    ap.add_argument("--reps", type=int, default=100)
    ap.add_argument("--nulls", type=int, default=200)
    ap.add_argument("--jobs", type=int, default=18)
    ap.add_argument("--seed", type=int, default=20260925)
    ap.add_argument("--solo-resumen", action="store_true",
                    help="no simular, resumir el CSV existente")
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
