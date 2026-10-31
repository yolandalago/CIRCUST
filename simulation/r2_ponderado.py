"""
R2 ponderado frente a R2 sin pesos en CircaFMM.

Por que hace falta
------------------
``fit_gene`` estima M, A y beta por minimos cuadrados PONDERADOS con los pesos
de densidad temporal, pero hasta el 28/09/2026 devolvia un R2 sin pesos,
``1 - sum(e^2) / sum((y - mean(y))^2)``. Ahora devuelve el coherente con el
ajuste,

    R2_w = 1 - sum(w e^2) / sum(w (y - mean_w(y))^2),

y este script compara las dos definiciones.

El R2 de CircaFMM decide cosas en el pipeline: el filtro y el orden
dentro de cada sector del conjunto TOP (``R2_FMM >= 0.5``), el umbral
``selection_r2_threshold = 0.65`` de la Etapa 4 y la puerta de R2 de los pares
de consistencia en la sincronizacion. Antes de cambiarlo, y de tener que
repetir la rejilla unificada, se mide cuanto cambia.

Partes
------
``genes``     Gen a gen: diferencia entre los dos R2, correlacion de rangos y
              cuantos genes cruzan los umbrales 0,5 y 0,65.
``pipeline``  Celdas de la rejilla con ``fmm_method="circafmm"``, corridas
              dos veces sobre la misma matriz, una con cada R2. Se compara el
              conjunto TOP, el orden de muestras y el atlas. La variante sin
              pesos se obtiene sustituyendo ``circafmm.fit_gene`` dentro
              de cada tarea.

Uso
---
    python -m simulation.r2_ponderado --parte genes
    python -m simulation.r2_ponderado --parte pipeline --jobs 12
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMBA_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse
from dataclasses import replace

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from scipy import stats

from simulation.rugosidad_dense_vs_obs import (
    DISENOS, OMEGAS, PROJECT_ROOT, SNRS, TWO_PI, alpha_para_pico, fs,
    muestrear_tiempos, onda_fmm,
)
from simulation.escala_lambda import escalar

UMBRALES = (0.5, 0.65)
VARIANTES = ("sin_pesos", "ponderado")


def r2_ponderado(y: np.ndarray, fitted: np.ndarray, w: np.ndarray) -> float:
    mu = (w * y).sum() / w.sum()
    ss = (w * (y - mu) ** 2).sum()
    return float(1.0 - (w * (y - fitted) ** 2).sum() / ss) if ss > 0 else 0.0


def r2_sin_pesos(y: np.ndarray, fitted: np.ndarray) -> float:
    ss = ((y - y.mean()) ** 2).sum()
    return float(1.0 - ((y - fitted) ** 2).sum() / ss) if ss > 0 else 0.0


# ── Parte 1: gen a gen ──────────────────────────────────────────────────────

def genes_extraccion(diseno, n, draw, n_rep, n_null, seed):
    rng = np.random.default_rng([seed, DISENOS.index(diseno), n, draw])
    t = muestrear_tiempos(diseno, n, rng)
    prep = fs.prepare_dataset(t, seed=seed + draw)
    w = prep.cache.weights
    filas = []

    def uno(y, **extra):
        f = _FIT_GENE(y, prep)
        filas.append(dict(diseno=diseno, n=n, draw=draw, r2=r2_sin_pesos(y, f.fitted),
                          r2w=r2_ponderado(y, f.fitted, w),
                          w_max=float(w.max()), w_cv=float(w.std() / w.mean()), **extra))

    for omega in OMEGAS:
        for snr in SNRS:
            for _ in range(n_rep):
                pk, be = rng.uniform(0.0, TWO_PI, 2)
                y0 = onda_fmm(t, 0.0, 1.0, alpha_para_pico(pk, be, omega), be, omega) \
                     + rng.standard_normal(n) / snr
                uno(escalar(y0, "minmax")[0], tipo="senal", omega_true=omega, snr=snr)
    for _ in range(n_null):
        uno(escalar(rng.standard_normal(n), "minmax")[0], tipo="nulo")
    return filas


def resumir_genes(df):
    pd.set_option("display.width", 220)
    df["d"] = df.r2w - df.r2
    print("\n=== R2 ponderado - R2 sin pesos, por diseno y n ===")
    def t(g):
        out = {"d_media": g.d.mean(), "|d|_mediana": g.d.abs().median(),
               "|d|_p90": g.d.abs().quantile(0.9),
               "spearman": stats.spearmanr(g.r2, g.r2w).statistic,
               "cv_pesos": g.w_cv.median()}
        for u in UMBRALES:
            a, b = g.r2 >= u, g.r2w >= u
            out[f"cruza_{u}"] = np.mean(a != b)
        return pd.Series(out)
    print(df.groupby(["diseno", "n"]).apply(t, include_groups=False).round(3))
    print("\n=== Por tipo de gen (todos los disenos) ===")
    print(df.groupby("tipo").apply(t, include_groups=False).round(3))


# ── Parte 2: pipeline completo ──────────────────────────────────────────────

_FIT_GENE = fs.fit_gene      # el del modulo, que ya devuelve el R2 ponderado


def _fit_gene_sin_pesos(data, prepared, *a, **k):
    """``fit_gene`` con el R2 sin pesos que tenia antes del cambio."""
    res = _FIT_GENE(data, prepared, *a, **k)
    y = np.asarray(data, dtype=float)
    return replace(res, r2=r2_sin_pesos(y, res.fitted))


def _fijar_variante(variante):
    """Fija ``fs.fit_gene`` para esta tarea.

    Se fija SIEMPRE, en las dos variantes: joblib reutiliza los procesos entre
    tareas, y una sustitucion hecha en una tarea seguiria activa en la
    siguiente que cayera en el mismo proceso.
    """
    fs.fit_gene = _fit_gene_sin_pesos if variante == "sin_pesos" else _FIT_GENE


def pipeline_celda(cfg, n_core, variante):
    from simulation.generate import generate_dataset
    from simulation.grid import METODO_REFERENCIA
    from simulation.metrics import evaluate
    from simulation.atlas import atlas_metrics
    from simulation.runner import run_stages_0_3, run_stage_4

    _fijar_variante(variante)
    fila = dict(n=cfg.n_samples, muestreo=cfg.sampling, omega=cfg.omega,
                snr=cfg.amplitude / cfg.sigma, n_genes=cfg.n_genes,
                prop=cfg.prop_rhythmic, variante=variante, status="ok")
    try:
        matriz, truth = generate_dataset(cfg)
        estado = run_stages_0_3(matriz, truth, n_core, "circafmm")
        e4 = run_stage_4(estado, truth, METODO_REFERENCIA, "circafmm")
        fila.update(evaluate(truth.theta_for(list(e4.theta_hat.index)), e4.theta_hat.values))
        fila.update(atlas_metrics(e4.atlas, truth.gene_params, quitar_giro=True))
        fila.update(n_top=estado.n_top, top="|".join(sorted(estado.top_matrix.index)),
                    r2_anchor=estado.r2_anchor, pares_ok=estado.pares_ok,
                    theta=";".join(f"{s}:{v:.6f}" for s, v in e4.theta_hat.items()))
    except Exception as exc:                                  # noqa: BLE001
        fila.update(status="fallo", error=f"{type(exc).__name__}: {exc}")
    return fila


def celdas_pipeline(max_celdas):
    from simulation.grid import _con_semilla, escenarios_rejilla, GCORE_CENTRAL
    escen = [c for c in escenarios_rejilla()
             if (c.n_genes, c.prop_rhythmic) == (1000, 0.20)
             and c.amplitude / c.sigma in (1.5, 4.0)]
    return [_con_semilla("rejilla", c, 0) for c in escen][:max_celdas], GCORE_CENTRAL


def resumir_pipeline(df):
    pd.set_option("display.width", 220)
    pd.set_option("display.max_columns", 30)
    ok = df[df.status == "ok"]
    print("\nfallos por variante:", df[df.status != "ok"].groupby("variante").size().to_dict())
    clave = ["n", "muestreo", "omega", "snr"]
    a = ok[ok.variante == "sin_pesos"].set_index(clave)
    b = ok[ok.variante == "ponderado"].set_index(clave)
    comun = a.index.intersection(b.index)
    a, b = a.loc[comun], b.loc[comun]

    def jacc(x, y):
        x, y = set(x.split("|")), set(y.split("|"))
        return len(x & y) / len(x | y) if x | y else 1.0

    def dif_orden(x, y):
        px = dict(s.split(":") for s in x.split(";"))
        py = dict(s.split(":") for s in y.split(";"))
        k = sorted(set(px) & set(py))
        d = np.angle(np.exp(1j * (np.array([float(px[s]) for s in k])
                                  - np.array([float(py[s]) for s in k]))))
        return np.median(np.abs(d)) * 24 / TWO_PI

    r = pd.DataFrame({
        "jaccard_TOP":   [jacc(x, y) for x, y in zip(a.top, b.top)],
        "d_n_top":       b.n_top.values - a.n_top.values,
        "mismo_orden":   [x == y for x, y in zip(a.theta, b.theta)],
        "dif_fase_h":    [dif_orden(x, y) for x, y in zip(a.theta, b.theta)],
    }, index=comun)
    metr = ["error_orden_h", "error_fase_medio_h", "tau_circular",
            "atlas_error_onda", "atlas_err_pico_h", "atlas_sesgo_omega", "atlas_r2_mediana"]
    for m in metr:
        if m in a.columns:
            r[f"d_{m}"] = b[m].values - a[m].values

    print(f"\n=== Celdas pareadas: {len(r)} ===")
    print("TOP identico en", round(np.mean(r.jaccard_TOP == 1), 3),
          "| Jaccard mediano", round(r.jaccard_TOP.median(), 3), " minimo", round(r.jaccard_TOP.min(), 3))
    print("orden de muestras identico en", round(r.mismo_orden.mean(), 3),
          "| diferencia de fase mediana entre variantes (h):", round(r.dif_fase_h.median(), 4),
          " max", round(r.dif_fase_h.max(), 3))
    print("\n=== Diferencias (ponderado - sin pesos) en las metricas del estudio ===")
    cols = [c for c in r.columns if c.startswith("d_")]
    print(r[cols].describe(percentiles=[0.1, 0.5, 0.9]).T.round(4))
    print("\n=== Por muestreo ===")
    print(r.groupby(level="muestreo")[["jaccard_TOP", "dif_fase_h", "d_error_orden_h",
                                       "d_atlas_err_pico_h", "d_atlas_error_onda"]].mean().round(4))
    for m in ("error_orden_h", "atlas_err_pico_h", "atlas_error_onda"):
        x = r[f"d_{m}"].dropna()
        x = x[x != 0]
        if len(x) > 5:
            print(f"Wilcoxon {m} (celdas con diferencia, n={len(x)}): p={stats.wilcoxon(x).pvalue:.3g}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--parte", choices=("genes", "pipeline"), required=True)
    ap.add_argument("--out", default="output/r2_ponderado")
    ap.add_argument("--draws", type=int, default=3)
    ap.add_argument("--reps", type=int, default=100)
    ap.add_argument("--nulls", type=int, default=200)
    ap.add_argument("--max-celdas", type=int, default=1000)
    ap.add_argument("--jobs", type=int, default=12)
    ap.add_argument("--seed", type=int, default=20260928)
    ap.add_argument("--solo-resumen", action="store_true")
    args = ap.parse_args()
    out = PROJECT_ROOT / args.out
    out.mkdir(parents=True, exist_ok=True)
    csv = out / f"{args.parte}.csv"

    if not args.solo_resumen:
        if args.parte == "genes":
            tareas = [(d, n, k) for d in DISENOS for n in (12, 24, 48, 96)
                      for k in range(args.draws)]
            res = Parallel(n_jobs=args.jobs, verbose=5)(
                delayed(genes_extraccion)(d, n, k, args.reps, args.nulls, args.seed)
                for d, n, k in tareas)
            pd.DataFrame([f for filas in res for f in filas]).to_csv(csv, index=False)
        else:
            celdas, gcore = celdas_pipeline(args.max_celdas)
            res = Parallel(n_jobs=args.jobs, verbose=5)(
                delayed(pipeline_celda)(c, gcore, v) for c in celdas for v in VARIANTES)
            pd.DataFrame(res).to_csv(csv, index=False)

    df = pd.read_csv(csv)
    (resumir_genes if args.parte == "genes" else resumir_pipeline)(df)


if __name__ == "__main__":
    main()
