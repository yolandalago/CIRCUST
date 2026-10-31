"""
Tests de la eleccion del recorrido en el TSP.

Como en ``isocir``, se preseleccionan los recorridos mas cortos y se desempata
por MSCE. El orden de la media circular no entra entre los candidatos: es el
metodo ``circmean``, aparte, y meterlo en el TSP mezclaria los dos metodos.
"""
import numpy as np

from circust.circular_aggregation import (
    _build_distance_matrix_tsp, _canonical_cycle, _msce_cire,
    _solve_tsp_multi_start, _tour_length, aggregate_circular_orders,
)

TWO_PI = 2.0 * np.pi


def _ruidoso(n=100, sigma=0.3, seed=0):
    rng = np.random.default_rng(seed)
    phi = np.sort(rng.uniform(0, TWO_PI, n))
    Theta = (phi + rng.normal(0, sigma, (5, n))) % TWO_PI
    return Theta, _build_distance_matrix_tsp(Theta)


def test_solo_se_evaluan_los_cinco_mas_cortos():
    Theta, D = _ruidoso()
    evaluados = []
    def score(t):
        evaluados.append(_canonical_cycle(t))
        return _msce_cire(Theta, t)
    _solve_tsp_multi_start(D, score_fn=score)
    assert 1 <= len(evaluados) <= 5
    assert len(set(evaluados)) == len(evaluados)


def test_el_criterio_elige_entre_los_cinco_y_no_fuera():
    Theta, D = _ruidoso()
    evaluados = []
    def score(t):
        evaluados.append(_canonical_cycle(t))
        return -len(evaluados)                 # gana el ultimo evaluado
    elegido, _ = _solve_tsp_multi_start(D, score_fn=score)
    assert _canonical_cycle(elegido) == evaluados[-1]


def test_sin_criterio_devuelve_el_mas_corto():
    Theta, D = _ruidoso(seed=1)
    tour, L = _solve_tsp_multi_start(D)
    assert L == _tour_length(D, tour)
    elegido, L_msce = _solve_tsp_multi_start(D, score_fn=lambda t: _msce_cire(Theta, t))
    assert L <= L_msce + 1e-12


def test_el_tsp_no_usa_la_media_circular():
    """Con ruido alto el orden de la media circular es mejor que cualquier
    recorrido del TSP. Si el TSP lo devolviera, se estaria colando ``circmean``
    dentro de ``tsp``."""
    Theta, _ = _ruidoso()
    cm = aggregate_circular_orders(Theta, method="circmean")["order"]
    tsp = aggregate_circular_orders(Theta, method="tsp")["order"]
    assert _canonical_cycle(tsp) != _canonical_cycle(cm)


# ── Parametro alpha de la ecuacion (7) ──────────────────────────────────────

def test_alpha_por_defecto_es_tsp3():
    Theta, _ = _ruidoso(n=30)
    from circust.circular_aggregation import _distance_matrix_tsp_single
    E3 = _distance_matrix_tsp_single(Theta[0], 3.0)
    assert np.array_equal(_distance_matrix_tsp_single(Theta[0]), E3)


def test_alpha_uno_es_simetrica_e_infinito_es_d_r():
    from circust.circular_aggregation import _distance_matrix_tsp_single
    rng = np.random.default_rng(4)
    th = rng.uniform(0, TWO_PI, 20)
    E1 = _distance_matrix_tsp_single(th, 1.0)
    assert np.allclose(E1, E1.T)
    dif = (th[None, :] - th[:, None]) % TWO_PI
    d_R = np.where(dif <= np.pi, 1 - np.cos(dif), 3 - np.cos(dif - np.pi))
    np.fill_diagonal(d_R, 0.0)
    assert np.allclose(_distance_matrix_tsp_single(th, np.inf), d_R)


def test_alpha_menor_que_uno_se_rechaza():
    import pytest
    from circust.robust_order import RobustOrderEstimator
    Theta, _ = _ruidoso(n=20)
    with pytest.raises(ValueError):
        aggregate_circular_orders(Theta, method="tsp", tsp_alpha=0.5)
    with pytest.raises(ValueError):
        RobustOrderEstimator(tsp_alpha=0.5, verbose=False)


def test_tsp3_sigue_aceptandose_como_sinonimo():
    from circust.robust_order import RobustOrderEstimator
    Theta, _ = _ruidoso(n=20)
    a = aggregate_circular_orders(Theta, method="tsp3")
    b = aggregate_circular_orders(Theta, method="tsp")
    assert a["method"] == b["method"] == "tsp"
    assert np.array_equal(a["order"], b["order"])
    assert RobustOrderEstimator(aggregation_method="tsp3", verbose=False).aggregation_method == "tsp"
