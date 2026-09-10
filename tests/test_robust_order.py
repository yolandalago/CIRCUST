"""
Tests de la correccion de orientaciones de la Etapa 4.

``_fix_orientations`` detecta las repeticiones cuyo sentido de recorrido
salio invertido y las corrige antes de que se usen para nada. Decidia por el
**signo** de la correlacion de Jammalamadaka, y ese signo no es fiable cuando
las muestras cubren todo el ciclo: sobre repeticiones genuinamente invertidas
del pipeline se midieron valores positivos de +0,11 a +0,36, con lo que la
condicion ``corr < 0`` no se cumplia nunca y no se corregia ninguna.
"""
import numpy as np
import pytest

from circust.robust_order import (
    RobustOrderEstimator, _build_theta, _tau_circular,
)

TWO_PI = 2.0 * np.pi


def _rep(theta_por_muestra):
    """Convierte fases por muestra en (orden, escala) como los produce la etapa."""
    orden = np.argsort(theta_por_muestra % TWO_PI)
    return orden, (theta_por_muestra % TWO_PI)[orden]


def _escenario(n=48, seed=0):
    """Preliminar mas dos reps: la 0 bien orientada y la 1 reflejada."""
    rng = np.random.default_rng(seed)
    theta_pre = np.sort(rng.uniform(0, TWO_PI, n))
    bien = (theta_pre + rng.normal(0, 0.05, n)) % TWO_PI
    mal  = (TWO_PI - bien) % TWO_PI
    o0, e0 = _rep(bien)
    o1, e1 = _rep(mal)
    return theta_pre, np.vstack([o0, o1]), np.vstack([e0, e1])


def test_tau_circular_identifica_la_rep_invertida():
    theta_pre, ordenes, escalas = _escenario()
    Theta = _build_theta(ordenes, escalas, theta_pre.size)
    assert _tau_circular(Theta[0], theta_pre) > 0.8
    assert _tau_circular(Theta[1], theta_pre) < -0.8


def test_fix_orientations_corrige_la_rep_invertida():
    theta_pre, ordenes, escalas = _escenario()
    est = RobustOrderEstimator(n_reps=2, verbose=False)
    o_fix, e_fix = est._fix_orientations(ordenes, escalas, theta_pre,
                                         theta_pre.size)
    Theta = _build_theta(o_fix, e_fix, theta_pre.size)
    assert _tau_circular(Theta[0], theta_pre) > 0.8
    assert _tau_circular(Theta[1], theta_pre) > 0.8, (
        "la repeticion reflejada deberia haber quedado corregida"
    )


def test_fix_orientations_no_toca_lo_que_ya_esta_bien():
    rng = np.random.default_rng(1)
    n = 48
    theta_pre = np.sort(rng.uniform(0, TWO_PI, n))
    reps = [_rep((theta_pre + rng.normal(0, 0.05, n)) % TWO_PI) for _ in range(3)]
    ordenes = np.vstack([o for o, _ in reps])
    escalas = np.vstack([e for _, e in reps])

    est = RobustOrderEstimator(n_reps=3, verbose=False)
    o_fix, e_fix = est._fix_orientations(ordenes, escalas, theta_pre, n)
    assert np.array_equal(o_fix, ordenes)
    assert np.allclose(e_fix, escalas)


def test_tau_circular_es_invariante_a_rotacion():
    """El origen de cada repeticion es arbitrario y no debe influir."""
    rng = np.random.default_rng(2)
    theta = np.sort(rng.uniform(0, TWO_PI, 40))
    for delta in (0.0, 1.0, 3.0, 5.5):
        assert _tau_circular(theta, (theta + delta) % TWO_PI) == pytest.approx(1.0)


def test_tau_circular_con_menos_de_tres_muestras():
    assert np.isnan(_tau_circular(np.array([0.0, 1.0]), np.array([0.0, 1.0])))
