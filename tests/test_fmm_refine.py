"""
Refinamiento local opcional de (alpha, omega) de CircaFMM, Seccion
S3.4.4 y Algoritmo S2 del suplemento de CircaFMM.

``_refine_jit`` se compara con una transcripcion literal del Algoritmo S2 en
Python puro, que calcula la rugosidad exacta sobre la rejilla densa en lugar
de leerla de la tabla en omega.
"""
import math

import numpy as np
import pytest

from circust.fitting import circafmm as fs

TWO_PI = 2.0 * math.pi
OMEGA0, POT = 0.10, 3.0
LAMBDA = 1.5


def _r_dense_exacta(alpha, omega, n_dense=512):
    td = np.linspace(0.0, TWO_PI, n_dense, endpoint=False)
    C, S = fs._mobius_cos_sin(td, alpha, omega)
    dt = np.full(n_dense, TWO_PI / n_dense)
    return float(fs._curvature_roughness(C[None], S[None], dt, np.arange(n_dense))[0])


def _J(x, t, w, alpha, omega):
    """Criterio J = P - lambda g(omega) R_dense, Secciones S3.3 y S3.4.3."""
    C, S = fs._mobius_cos_sin(t, alpha, omega)
    x0 = x - (w * x).sum() / w.sum()
    P = ((w * x0 * C).sum() ** 2 + (w * x0 * S).sum() ** 2) / (w * (C * C + S * S)).sum()
    return P - LAMBDA * (OMEGA0 / (omega + OMEGA0)) ** POT * _r_dense_exacta(alpha, omega)


def _secundario(x, t, w, alpha, omega):
    """WRSS(alpha, omega) + lambda g(omega) R(alpha, omega)."""
    fit = fs._wls_fit_at(x, alpha, omega, w, t)
    return fit["rss"] + LAMBDA * (OMEGA0 / (omega + OMEGA0)) ** POT * _r_dense_exacta(alpha, omega)


def _algoritmo_s2_literal(x, t, w, a0, o0, S=3, d_alpha=TWO_PI / 48, eta=1.25,
                          omega_min=1e-4, omega_max=1.0, tau=1e-4):
    a, o = a0, o0
    j_cur = _J(x, t, w, a, o)
    for _ in range(S):
        vecinos = [(a + d_alpha, o), (a - d_alpha, o), (a, o * eta), (a, o / eta),
                   (a + d_alpha, o * eta), (a - d_alpha, o * eta),
                   (a + d_alpha, o / eta), (a - d_alpha, o / eta)]
        vecinos = [(ai % TWO_PI, min(max(oi, omega_min), omega_max)) for ai, oi in vecinos]
        js = [_J(x, t, w, ai, oi) for ai, oi in vecinos]
        jmax = max(js)
        retenidos = [k for k, j in enumerate(js) if j >= jmax - tau * max(1.0, abs(jmax))]
        # entre los retenidos, minimo del objetivo secundario y empate a omega mayor
        best = min(retenidos, key=lambda k: (_secundario(x, t, w, *vecinos[k]), -vecinos[k][1]))
        if js[best] > j_cur:
            a, o = vecinos[best]
            j_cur = js[best]
        else:
            d_alpha /= 2.0
            eta = math.sqrt(eta)
    return a, o


@pytest.fixture(scope="module")
def dataset():
    rng = np.random.default_rng(7)
    t = np.sort(rng.vonmises(rng.uniform(0, TWO_PI, 3)[rng.integers(0, 3, 30)], 6.0) % TWO_PI)
    prep = fs.prepare_dataset(t, lambda_star=LAMBDA)
    genes = []
    for om in (0.05, 0.15, 0.4, 0.8):
        for _ in range(8):
            a, b = rng.uniform(0, TWO_PI, 2)
            y = np.cos(b + 2 * np.arctan(om * np.tan((t - a) / 2))) + 0.3 * rng.standard_normal(30)
            genes.append(y)
    return t, prep, genes


def test_refinado_igual_que_algoritmo_s2_literal(dataset):
    t, prep, genes = dataset
    w = prep.cache.weights
    for y in genes:
        grueso = fs.fit_gene(y, prep, refine=False).params
        esperado = _algoritmo_s2_literal(y, t, w, grueso["alpha"], grueso["omega"])
        fino = fs.fit_gene(y, prep, refine=True).params
        assert fino["omega"] == pytest.approx(esperado[1], rel=1e-9)
        assert abs(math.remainder(fino["alpha"] - esperado[0], TWO_PI)) < 1e-9


def test_refinado_nunca_empeora_J(dataset):
    t, prep, genes = dataset
    w = prep.cache.weights
    for y in genes:
        g = fs.fit_gene(y, prep, refine=False).params
        f = fs.fit_gene(y, prep, refine=True).params
        assert _J(y, t, w, f["alpha"], f["omega"]) >= _J(y, t, w, g["alpha"], g["omega"]) - 1e-9
        assert 1e-4 <= f["omega"] <= 1.0


def test_refinado_desactivado_por_defecto():
    from circust.cpca import CPCA
    from circust.robust_order import RobustOrderEstimator
    from circust.top_genes import TopGeneSelector

    assert fs.CircaFMMModel().refine is False
    assert CPCA(fmm_method="circafmm", verbose=False).fmm_refine is False
    assert TopGeneSelector(fmm_method="circafmm", verbose=False).fmm_refine is False
    assert RobustOrderEstimator(fmm_method="circafmm", verbose=False).fmm_refine is False


def test_cpca_pasa_el_refinado_al_modelo():
    from circust.cpca import CPCA

    t = np.linspace(0, TWO_PI, 24, endpoint=False)
    circa = {"lambda_star": LAMBDA}
    assert CPCA(fmm_method="circafmm", circafmm_kwargs=circa, fmm_refine=True,
                verbose=False)._make_fmm_model(t).refine is True
    assert CPCA(fmm_method="circafmm", circafmm_kwargs=circa,
                verbose=False)._make_fmm_model(t).refine is False
