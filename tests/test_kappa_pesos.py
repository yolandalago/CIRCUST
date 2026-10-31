"""
Tests de los pesos por repeticion de la agregacion circular.

Barragan et al. (2015) ponderan cada experimento con
``w_j = kappa_j / sum(kappa)``, con ``theta_ij ~ vM(phi_ij, kappa_j)``. Los
kappa salen de Fernandez et al. (2012, Supp. 1.3), cuya maxima verosimilitud
no esta acotada y diverge al iterarla, asi que aqui se estiman por pares: la
diferencia entre dos repeticiones sobre el mismo item cancela ``phi_i`` y su
longitud resultante esperada es ``A(kappa_j) A(kappa_l)``.
"""
import itertools

import numpy as np
import pytest

from circust.circular_aggregation import (
    KAPPA_MAX, _a1, _a1_inv, aggregate_circular_orders,
    estimate_kappa_pairwise, kappa_weights,
)

TWO_PI = 2.0 * np.pi


def _theta(kappas, n, seed=0):
    """K repeticiones de las mismas n fases con ruido vM(0, kappa_j)."""
    rng = np.random.default_rng(seed)
    phi = rng.uniform(0, TWO_PI, n)
    return np.array([(phi + rng.vonmises(0, k, n)) % TWO_PI for k in kappas]), phi


def test_a1_es_estable_y_su_inversa_la_deshace():
    for k in (0.3, 2.0, 29.08, 700.0, 5000.0):
        assert np.isfinite(_a1(k))
        assert _a1_inv(_a1(k)) == pytest.approx(k, rel=1e-6)
    assert _a1_inv(0.0) == 0.0
    assert _a1_inv(1.0) == KAPPA_MAX


def test_recupera_kappas_conocidos():
    kappas = np.array([200.0, 50.0, 10.0, 2.0])
    Theta, _ = _theta(kappas, n=5000)
    est = estimate_kappa_pairwise(Theta)
    assert est == pytest.approx(kappas, rel=0.15)


def test_un_giro_por_repeticion_no_cambia_nada():
    sin_giro, _ = _theta([100.0, 20.0, 5.0], n=500, seed=3)
    con_giro = (sin_giro + np.array([[1.3], [4.0], [5.9]])) % TWO_PI
    assert estimate_kappa_pairwise(con_giro) == pytest.approx(
        estimate_kappa_pairwise(sin_giro), rel=1e-9)


def test_con_tres_repeticiones_coincide_con_la_formula_cerrada():
    Theta, _ = _theta([80.0, 15.0, 4.0], n=400, seed=5)
    rho = {}
    for j, l in itertools.combinations(range(3), 2):
        d = Theta[j] - Theta[l]
        rho[j, l] = np.hypot(np.cos(d).mean(), np.sin(d).mean())
    a1 = np.sqrt(rho[0, 1] * rho[0, 2] / rho[1, 2])
    assert estimate_kappa_pairwise(Theta)[0] == pytest.approx(_a1_inv(a1), rel=1e-9)


def test_kappas_iguales_dan_pesos_casi_iguales():
    Theta, _ = _theta([50.0] * 5, n=3000, seed=1)
    w = kappa_weights(Theta)
    assert w.sum() == pytest.approx(1.0)
    assert w == pytest.approx(np.full(5, 0.2), abs=0.02)


def test_la_repeticion_mala_casi_no_pesa():
    Theta, _ = _theta([300.0, 300.0, 300.0, 300.0, 1.0], n=100, seed=2)
    w = kappa_weights(Theta)
    assert w[4] < 0.01


def test_necesita_tres_repeticiones():
    Theta, _ = _theta([50.0, 10.0], n=100)
    with pytest.raises(ValueError):
        estimate_kappa_pairwise(Theta)


def test_aggregate_usa_kappa_por_defecto_y_respeta_pesos_dados():
    Theta, _ = _theta([300.0, 300.0, 300.0, 300.0, 1.0], n=40, seed=4)
    por_defecto = aggregate_circular_orders(Theta, method="hods")
    assert por_defecto["weights"] == pytest.approx(kappa_weights(Theta))

    dados = aggregate_circular_orders(Theta, weights=np.array([2, 2, 2, 2, 2.]),
                                      method="hods")
    assert dados["weights"] == pytest.approx(np.full(5, 0.2))

    with pytest.raises(ValueError):
        aggregate_circular_orders(Theta, weights=np.array([1.0, 1.0]), method="hods")


def test_con_dos_repeticiones_aggregate_cae_a_pesos_iguales():
    Theta, _ = _theta([50.0, 10.0], n=30)
    out = aggregate_circular_orders(Theta, method="hods")
    assert out["weights"] == pytest.approx([0.5, 0.5])
