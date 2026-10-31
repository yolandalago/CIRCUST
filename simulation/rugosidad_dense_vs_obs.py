"""
Rugosidad densa frente a rugosidad sobre las observaciones en CircaFMM.

Por que hace falta
------------------
``circafmm`` penaliza la rugosidad ``R_dense``, evaluada sobre una
rejilla uniforme de 512 puntos, y no ``R_obs``, evaluada sobre los instantes
observados y sus huecos circulares. El argumento es que ``R_obs`` no ve un
pico que cae entero dentro de un hueco, y que se dispara con instantes casi
duplicados. Aqui se comprueba con datos sinteticos si esa eleccion da ajustes
mejores, iguales o peores.

Diseno
------
Para cada diseno de muestreo (``regular``, ``hueco`` de 120 grados,
``agrupado`` en tres cumulos von Mises), cada ``n`` y cada extraccion de
tiempos se prepara el dataset dos veces, identico salvo por la rugosidad:
mismos pesos de densidad, misma rejilla (alpha, omega), mismo cache WLS.
``lambda*`` se calibra por separado para cada rugosidad con el Monte Carlo del
articulo (B = 1000, R = 5, q = 0,975), porque las dos escalas no son
comparables. Sobre esos tiempos se ajustan las mismas senales FMM y los mismos
ruidos con las dos variantes, de modo que la comparacion es pareada.

``R_obs`` se discretiza igual que ``R_dense`` (misma ``_curvature_roughness``),
cambiando solo la rejilla: los instantes observados ordenados y sus huecos
circulares.

Metricas
--------
Parametros: error circular de la hora del pico (horas) y fraccion con error
menor de 2 h, errores circulares de alpha y beta (rad), error absoluto,
relativo y sesgo de omega, error relativo y sesgo de A, y error de M. Con
omega cerca de 1 alpha y beta dejan de ser identificables por separado (beta
absorbe alpha), asi que ahi la hora del pico es la medida de fase fiable. Calidad del ajuste: RMSE entre la curva ajustada y
la verdadera sobre 512 puntos uniformes (en unidades de A), el mismo RMSE
restringido a los huecos grandes del muestreo, y el R2 sobre las
observaciones. Bajo H0 (ruido puro): el pico oculto, cuanto sube la curva
ajustada entre muestras por encima de su maximo en las muestras (en unidades
de sigma), y la fraccion de ajustes con omega < 0,1.

Uso
---
    python -m simulation.rugosidad_dense_vs_obs --out output/rugosidad
"""
from __future__ import annotations

import os

# Un hilo BLAS/numba por proceso: el paralelismo va por extracciones de tiempos.
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMBA_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse
import sys
from dataclasses import replace
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT / "src"))

from circust.fitting import circafmm as fs

TWO_PI = 2.0 * np.pi
HORAS = 24.0 / TWO_PI

DISENOS   = ("regular", "hueco", "agrupado")
NS        = (12, 24, 48)
OMEGAS    = (0.1, 0.25, 0.5, 0.9)
SNRS      = (2.0, 4.0)          # A / sigma
N_DENSE   = 512
T_DENSE   = np.linspace(0.0, TWO_PI, N_DENSE, endpoint=False)


# ── Diseños de muestreo ─────────────────────────────────────────────────────

def muestrear_tiempos(diseno: str, n: int, rng: np.random.Generator) -> np.ndarray:
    if diseno == "regular":
        t = np.linspace(0.0, TWO_PI, n, endpoint=False)
        t = t + rng.uniform(-0.25, 0.25, n) * (TWO_PI / n)
    elif diseno == "hueco":
        # Uniforme en la circunferencia salvo un hueco de 120 grados
        inicio = rng.uniform(0.0, TWO_PI)
        t = inicio + TWO_PI / 3.0 + rng.uniform(0.0, 2.0 * TWO_PI / 3.0, n)
    elif diseno == "agrupado":
        centros = rng.uniform(0.0, TWO_PI, 3)
        t = rng.vonmises(centros[rng.integers(0, 3, n)], 8.0)
    else:
        raise ValueError(diseno)
    return np.sort(t % TWO_PI)


def huecos_circulares(t_ord: np.ndarray) -> np.ndarray:
    """``dt[k] = t[k+1] - t[k]`` con el ultimo salto cerrando la circunferencia."""
    return np.diff(np.append(t_ord, t_ord[0] + TWO_PI))


def mascara_huecos_grandes(t_ord: np.ndarray, puntos: np.ndarray) -> np.ndarray:
    """Que ``puntos`` caen en un hueco mayor que el doble del espaciado medio."""
    dt = huecos_circulares(t_ord)
    umbral = 2.0 * TWO_PI / t_ord.shape[0]
    en_hueco = np.zeros(puntos.shape[0], dtype=bool)
    for k in np.flatnonzero(dt > umbral):
        en_hueco |= ((puntos - t_ord[k]) % TWO_PI) < dt[k]
    return en_hueco


# ── Modelo FMM ──────────────────────────────────────────────────────────────

def onda_fmm(t, M, A, alpha, beta, omega):
    return M + A * np.cos(beta + 2.0 * np.arctan(omega * np.tan((t - alpha) / 2.0)))


def alpha_para_pico(pico: float, beta: float, omega: float) -> float:
    """Inversa de la formula de ``peak_time`` de ``fit_gene``."""
    return (pico - 2.0 * np.arctan2(np.sin(-beta / 2.0) / omega,
                                    np.cos(-beta / 2.0))) % TWO_PI


def diff_circular(a, b):
    return np.abs(np.arctan2(np.sin(a - b), np.cos(a - b)))


# ── Preparación con una u otra rugosidad ────────────────────────────────────

def preparar_ambas(t: np.ndarray, seed: int) -> dict[str, fs.FMMPreparedDataset]:
    weights = fs._time_density_weights(t)
    grid_dense = fs.precompute_grid(t, n_dense=N_DENSE)
    R_obs = fs._curvature_roughness(grid_dense.C, grid_dense.S,
                                    huecos_circulares(t[grid_dense.order]),
                                    grid_dense.order)
    grid_obs = replace(grid_dense, roughness=R_obs)
    cache = fs.prepare_cache(weights, grid_dense)

    out = {}
    for nombre, grid in (("dense", grid_dense), ("obs", grid_obs)):
        lam = fs.calibrate_lambda_repeat(
            n_repeats=5, seed0=seed, agg="median",
            times=t, weights=weights, grid=grid, cache=cache,
            n_simulations=1000, quantile=0.975, sigma=1.0,
            omega0=0.10, penalty_power=3.0,
        )
        out[nombre] = fs.FMMPreparedDataset(
            times=t, weights=weights, grid=grid, cache=cache,
            lambda_star=float(lam), omega0=0.10, penalty_power=3.0,
            penalty=fs._penalty_vector(grid, float(lam), 0.10, 3.0),
        )
    return out


# ── Una extracción de tiempos: calibra y ajusta todas las señales ──────────

def correr_extraccion(diseno: str, n: int, draw: int, n_rep: int, n_null: int,
                      seed: int) -> list[dict]:
    rng = np.random.default_rng([seed, DISENOS.index(diseno), n, draw])
    t = muestrear_tiempos(diseno, n, rng)
    preps = preparar_ambas(t, seed=seed + 1000 * draw)
    hueco_dense = mascara_huecos_grandes(t, T_DENSE)

    filas = []
    base = dict(diseno=diseno, n=n, draw=draw,
                lambda_dense=preps["dense"].lambda_star,
                lambda_obs=preps["obs"].lambda_star)

    # Señales FMM
    for omega in OMEGAS:
        for snr in SNRS:
            for r in range(n_rep):
                pico = rng.uniform(0.0, TWO_PI)
                beta = rng.uniform(0.0, TWO_PI)
                alpha = alpha_para_pico(pico, beta, omega)
                sigma = 1.0 / snr
                verdad_obs = onda_fmm(t, 0.0, 1.0, alpha, beta, omega)
                verdad_den = onda_fmm(T_DENSE, 0.0, 1.0, alpha, beta, omega)
                y = verdad_obs + sigma * rng.standard_normal(n)
                pico_en_hueco = bool(mascara_huecos_grandes(t, np.array([pico]))[0])

                for rug, prep in preps.items():
                    fit = fs.fit_gene(y, prep)
                    p = fit.params
                    curva = onda_fmm(T_DENSE, p["M"], p["A"], p["alpha"],
                                     p["beta"], p["omega"])
                    err = curva - verdad_den
                    filas.append(dict(
                        base, tipo="senal", rugosidad=rug, rep=r,
                        omega_true=omega, snr=snr, pico_en_hueco=pico_en_hueco,
                        err_pico_h=diff_circular(fit.peak_time, pico) * HORAS,
                        err_alpha=diff_circular(p["alpha"], alpha),
                        err_beta=diff_circular(p["beta"], beta),
                        err_omega=p["omega"] - omega,
                        err_omega_rel=(p["omega"] - omega) / omega,
                        err_A_rel=p["A"] - 1.0,
                        err_M=p["M"],
                        rmse_curva=float(np.sqrt(np.mean(err ** 2))),
                        rmse_curva_hueco=(float(np.sqrt(np.mean(err[hueco_dense] ** 2)))
                                          if hueco_dense.any() else np.nan),
                        r2=fit.r2,
                    ))

    # Ruido puro (H0)
    for r in range(n_null):
        y = rng.standard_normal(n)
        for rug, prep in preps.items():
            fit = fs.fit_gene(y, prep)
            p = fit.params
            curva = onda_fmm(T_DENSE, p["M"], p["A"], p["alpha"], p["beta"], p["omega"])
            obs = onda_fmm(t, p["M"], p["A"], p["alpha"], p["beta"], p["omega"])
            filas.append(dict(
                base, tipo="nulo", rugosidad=rug, rep=r,
                omega_hat=p["omega"],
                pico_oculto=float(np.abs(curva - p["M"]).max()
                                  - np.abs(obs - p["M"]).max()),
                A_hat=p["A"],
            ))
    return filas


# ── Resúmenes ───────────────────────────────────────────────────────────────

def resumir(df: pd.DataFrame) -> None:
    pd.set_option("display.width", 200)
    pd.set_option("display.max_columns", 30)
    s = df[df.tipo == "senal"]

    print("\n=== lambda* calibrado (mediana por diseno y n) ===")
    lam = df.drop_duplicates(["diseno", "n", "draw"]).groupby(["diseno", "n"])[
        ["lambda_dense", "lambda_obs"]].median()
    print(lam.round(3))

    def tabla(g):
        return pd.Series({
            "pico_h":      np.median(g.err_pico_h),
            "pico<2h":     np.mean(g.err_pico_h < 2.0),
            "|d_alpha|":   np.median(g.err_alpha),
            "|d_beta|":    np.median(g.err_beta),
            "|d_omega|":   np.median(np.abs(g.err_omega)),
            "|d_om|/om":   np.median(np.abs(g.err_omega_rel)),
            "sesgo_omega": np.median(g.err_omega),
            "|dA|/A":      np.median(np.abs(g.err_A_rel)),
            "sesgo_A":     np.median(g.err_A_rel),
            "|dM|":        np.median(np.abs(g.err_M)),
            "rmse":        np.median(g.rmse_curva),
            "rmse_hueco":  np.nanmedian(g.rmse_curva_hueco),
            "R2":          np.median(g.r2),
        })

    print("\n=== Senal: medianas por diseno, n y rugosidad (todas las omega y SNR) ===")
    print(s.groupby(["diseno", "n", "rugosidad"]).apply(tabla).round(3))

    print("\n=== Senal: medianas por omega verdadera y rugosidad (todos los disenos) ===")
    print(s.groupby(["omega_true", "rugosidad"]).apply(tabla).round(3))

    print("\n=== Senal, disenos con huecos: pico dentro o fuera de un hueco grande ===")
    sh = s[s.diseno != "regular"]
    print(sh.groupby(["pico_en_hueco", "rugosidad"]).apply(tabla).round(3))

    # Comparación pareada
    clave = ["diseno", "n", "draw", "omega_true", "snr", "rep"]
    w = s.pivot_table(index=clave, columns="rugosidad",
                      values=["rmse_curva", "err_pico_h"])
    par = pd.DataFrame({
        "gana_dense_rmse": (w["rmse_curva"]["dense"] < w["rmse_curva"]["obs"] - 1e-9),
        "gana_obs_rmse":   (w["rmse_curva"]["obs"] < w["rmse_curva"]["dense"] - 1e-9),
        "gana_dense_pico": (w["err_pico_h"]["dense"] < w["err_pico_h"]["obs"] - 1e-9),
        "gana_obs_pico":   (w["err_pico_h"]["obs"] < w["err_pico_h"]["dense"] - 1e-9),
        "d_rmse":          w["rmse_curva"]["obs"] - w["rmse_curva"]["dense"],
    }).reset_index()
    print("\n=== Pareado (misma senal y ruido): fraccion de victorias, d_rmse = obs - dense ===")
    print(par.groupby(["diseno", "n"])[["gana_dense_rmse", "gana_obs_rmse",
                                        "gana_dense_pico", "gana_obs_pico",
                                        "d_rmse"]].mean().round(3))

    z = df[df.tipo == "nulo"]
    print("\n=== H0 (ruido puro): pico oculto entre muestras (en sigma) y omega pequena ===")
    print(z.groupby(["diseno", "n", "rugosidad"]).apply(lambda g: pd.Series({
        "pico_oculto_med": np.median(g.pico_oculto),
        "pico_oculto_p90": np.quantile(g.pico_oculto, 0.9),
        "frac_omega<0.1":  np.mean(g.omega_hat < 0.1),
        "A_hat_med":       np.median(g.A_hat),
    })).round(3))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--out", default="output/rugosidad")
    ap.add_argument("--draws", type=int, default=4,
                    help="extracciones de tiempos por diseno y n")
    ap.add_argument("--reps", type=int, default=100,
                    help="senales por (extraccion, omega, SNR)")
    ap.add_argument("--nulls", type=int, default=200,
                    help="ajustes de ruido puro por extraccion")
    ap.add_argument("--jobs", type=int, default=12)
    ap.add_argument("--seed", type=int, default=20260925)
    args = ap.parse_args()

    out = PROJECT_ROOT / args.out
    out.mkdir(parents=True, exist_ok=True)

    tareas = [(d, n, k) for d in DISENOS for n in NS for k in range(args.draws)]
    res = Parallel(n_jobs=args.jobs, verbose=5)(
        delayed(correr_extraccion)(d, n, k, args.reps, args.nulls, args.seed)
        for d, n, k in tareas
    )
    df = pd.DataFrame([f for filas in res for f in filas])
    df.to_csv(out / "resultados.csv", index=False)
    resumir(df)


if __name__ == "__main__":
    main()
