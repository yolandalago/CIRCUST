"""
Tests de los tres metodos de consenso de ``aggregate_circular_orders``.

``circmean`` da la fase de cada item como media circular ponderada de sus K
angulos y el orden que inducen. ``tsp`` y ``hod`` dan el orden de Barragan et
al. y reconstruyen las fases con el CIRE de esa media, ponderado por la
longitud resultante de cada item, de modo que las fases verifican el orden.
"""
import numpy as np
import pytest

from circust.circular_aggregation import (
    _canonical_cycle, aggregate_circular_orders, circular_isotonic_regression,
)

TWO_PI = 2.0 * np.pi


def _theta(n=40, sigma=0.3, K=5, seed=0):
    rng = np.random.default_rng(seed)
    phi = rng.uniform(0, TWO_PI, n)
    return (phi + rng.normal(0, sigma, (K, n))) % TWO_PI, phi


def _descensos(fases, orden):
    s = fases[orden]
    k = int(np.argmin(s))
    return int(np.sum(np.diff(np.roll(s, -k)) < -1e-12))


def test_circmean_es_el_metodo_por_defecto():
    Theta, _ = _theta()
    assert aggregate_circular_orders(Theta)["method"] == "circmean"


def test_circmean_da_la_media_ponderada_y_su_orden():
    Theta, _ = _theta(seed=1)
    out = aggregate_circular_orders(Theta, method="circmean")
    w = out["weights"]
    media = np.arctan2((w[:, None] * np.sin(Theta)).sum(0),
                       (w[:, None] * np.cos(Theta)).sum(0)) % TWO_PI
    assert out["phases"] == pytest.approx(media)
    assert np.array_equal(out["order"], np.argsort(media, kind="stable"))
    assert out["tour_length"] is None


@pytest.mark.parametrize("metodo", ["tsp", "hod"])
def test_las_fases_de_tsp_y_hod_verifican_su_orden(metodo):
    Theta, _ = _theta(n=60, sigma=0.6, seed=2)
    out = aggregate_circular_orders(Theta, method=metodo)
    assert _descensos(out["phases"], out["order"]) == 0


def test_las_fases_de_tsp_son_el_cire_ponderado_de_la_media():
    Theta, _ = _theta(seed=3)
    out = aggregate_circular_orders(Theta, method="tsp")
    w = out["weights"]
    S = (w[:, None] * np.sin(Theta)).sum(0)
    C = (w[:, None] * np.cos(Theta)).sum(0)
    esperado, _ = circular_isotonic_regression(np.arctan2(S, C) % TWO_PI, out["order"],
                                               weights=np.hypot(S, C))
    assert out["phases"] == pytest.approx(esperado)


def test_sin_ruido_los_tres_metodos_recuperan_la_verdad():
    rng = np.random.default_rng(4)
    phi = rng.uniform(0, TWO_PI, 30)
    Theta = np.tile(phi, (5, 1))
    for metodo in ("circmean", "tsp", "hod"):
        out = aggregate_circular_orders(Theta, method=metodo)
        assert _canonical_cycle(out["order"]) == _canonical_cycle(np.argsort(phi))
        assert np.allclose(np.angle(np.exp(1j * (out["phases"] - phi))), 0.0, atol=1e-9)


def test_cire_ponderado_con_pesos_uno_es_el_cire():
    rng = np.random.default_rng(5)
    for _ in range(30):
        n = int(rng.integers(3, 25))
        x = rng.uniform(0, TWO_PI, n)
        o = rng.permutation(n)
        assert np.array_equal(circular_isotonic_regression(x, o, np.ones(n))[0],
                              circular_isotonic_regression(x, o)[0])


def test_metodo_auto_ya_no_existe():
    Theta, _ = _theta()
    with pytest.raises(ValueError):
        aggregate_circular_orders(Theta, method="auto")
