"""Tests de simulation/metrics.py."""
import numpy as np
import pytest

from simulation.metrics import (
    circular_correlation, circular_distance, evaluate,
    kendall_tau_fixed_cut, oracle_mae_h, tau_circular, to_hours,
)
from simulation.truth import TWO_PI


def _ang(horas):
    return np.asarray(horas, dtype=float) * TWO_PI / 24.0


def test_ejemplo_de_la_spec_23h_contra_1h_da_2h():
    d = circular_distance(_ang([23.0]), _ang([1.0]))
    assert to_hours(d)[0] == pytest.approx(2.0)


def test_distancia_circular_es_simetrica_y_acotada():
    rng = np.random.default_rng(0)
    a, b = rng.uniform(0, TWO_PI, 50), rng.uniform(0, TWO_PI, 50)
    assert np.allclose(circular_distance(a, b), circular_distance(b, a))
    assert circular_distance(a, b).max() <= np.pi + 1e-12


def test_correlacion_circular_vale_uno_consigo_misma():
    theta = np.random.default_rng(0).uniform(0, TWO_PI, 40)
    assert circular_correlation(theta, theta) == pytest.approx(1.0)


def test_correlacion_circular_invariante_a_rotacion():
    theta = np.random.default_rng(0).uniform(0, TWO_PI, 40)
    rotado = (theta + 1.234) % TWO_PI
    assert circular_correlation(theta, rotado) == pytest.approx(1.0)


def test_correlacion_circular_no_esta_definida_con_rejilla_exacta():
    """
    Con muestras exactamente equiespaciadas la resultante es del orden de
    1e-17, la direccion media es arbitraria y deja de ser equivariante bajo
    rotacion. Devolver nan en vez de un numero sin sentido, que ahi daba 0,26
    al comparar una escala consigo misma rotada.

    Solo afecta al muestreo equiespaciado, que usa el bloque computacional,
    donde la medida de interes es el tiempo y no la calidad.
    """
    rejilla = _ang(np.linspace(0, 24, 40, endpoint=False))
    assert np.isnan(circular_correlation(rejilla, rejilla))
    # una muestra uniforme aleatoria del mismo tamano si esta definida
    aleatoria = np.random.default_rng(1).uniform(0, TWO_PI, 40)
    assert not np.isnan(circular_correlation(aleatoria, aleatoria))


def test_kendall_perfecto_y_perfectamente_invertido():
    theta = _ang(np.linspace(0, 24, 24, endpoint=False))
    assert kendall_tau_fixed_cut(theta, theta) == pytest.approx(1.0)
    assert kendall_tau_fixed_cut(theta, (-theta) % TWO_PI) < 0.0


def test_tau_circular_es_invariante_a_cualquier_rotacion():
    """
    La propiedad que justifica los tripletes. Con corte fijo, un desfase de
    3 rad deja tau en 0,015, que se lee como si no hubiera relacion ninguna
    cuando la reconstruccion es perfecta.
    """
    theta = np.sort(np.random.default_rng(3).uniform(0, TWO_PI, 40))
    for delta in (0.0, 0.4, 1.5, 3.0, 5.7):
        rotado = (theta + delta) % TWO_PI
        assert tau_circular(theta, rotado) == pytest.approx(1.0)
    # con corte fijo la de 3 rad se lee como ausencia de relacion
    assert abs(kendall_tau_fixed_cut(theta, (theta + 3.0) % TWO_PI)) < 0.1


def test_tau_circular_da_menos_uno_al_reflejar_aunque_este_rotado():
    """
    El signo separa las dos patologias: +1 con rotacion es que fallo el
    origen, -1 es que fallo la orientacion. Con corte fijo, reflejado y
    rotado a la vez da -0,015, indistinguible del azar.
    """
    theta = np.sort(np.random.default_rng(3).uniform(0, TWO_PI, 40))
    for delta in (0.0, 2.0, 4.5):
        reflejado = (-theta + delta) % TWO_PI
        assert tau_circular(theta, reflejado) == pytest.approx(-1.0)
    assert abs(kendall_tau_fixed_cut(theta, (-theta + 2.0) % TWO_PI)) < 0.1


def test_tau_circular_vale_cero_sin_relacion_y_degrada_con_el_ruido():
    rng = np.random.default_rng(11)
    theta = np.sort(rng.uniform(0, TWO_PI, 40))
    assert abs(tau_circular(theta, rng.uniform(0, TWO_PI, 40))) < 0.2
    suave = tau_circular(theta, (theta + rng.normal(0, 0.2, 40)) % TWO_PI)
    fuerte = tau_circular(theta, (theta + rng.normal(0, 0.8, 40)) % TWO_PI)
    assert 1.0 > suave > fuerte > 0.0


def test_oracle_encuentra_la_rotacion_y_la_reflexion():
    theta = _ang(np.linspace(0, 24, 48, endpoint=False))
    estimado = (-theta + 2.0) % TWO_PI     # reflejado y rotado
    assert oracle_mae_h(theta, estimado) < 0.1
    assert oracle_mae_h(theta, theta) == pytest.approx(0.0, abs=1e-9)


def test_evaluate_devuelve_todas_las_claves():
    theta = np.random.default_rng(2).uniform(0, TWO_PI, 24)
    m = evaluate(theta, theta)
    esperadas = {"mae_h", "medae_h", "rho_circ", "tau",
                 "tau_circ", "p1", "p2", "p4", "mae_h_oracle"}
    assert set(m) == esperadas
    assert m["mae_h"] == pytest.approx(0.0)
    assert m["p1"] == pytest.approx(1.0)
    assert m["rho_circ"] == pytest.approx(1.0)


def test_proporciones_cuentan_bien():
    theta    = _ang([0.0, 0.0, 0.0, 0.0])
    estimado = _ang([0.5, 1.5, 3.0, 6.0])
    m = evaluate(theta, estimado)
    assert m["p1"] == pytest.approx(0.25)   # solo 0.5 h
    assert m["p2"] == pytest.approx(0.50)   # 0.5 y 1.5 h
    assert m["p4"] == pytest.approx(0.75)   # 0.5, 1.5 y 3.0 h


def test_longitudes_distintas_fallan():
    with pytest.raises(ValueError, match="misma longitud"):
        evaluate(_ang([0.0, 1.0]), _ang([0.0]))


def test_rho_es_solido_en_mediana_pero_tiene_cola_de_signo_contrario():
    """
    Fragilidad medida de la correlacion circular, no un bug.

    Cuando las muestras cubren todo el ciclo la resultante es pequena
    (R ~ 0,13 de media con n=48). En las semillas donde sale especialmente
    baja, un ruido pequeno manda la direccion media de la escala estimada
    casi al lado opuesto de la verdadera y rho cambia de signo. Con R=0,008
    y ruido de 0,2 rad se midio rho = -0,86.

    Ocurre en torno al 1 % de las semillas. La mediana aguanta y la media no,
    asi que analyze.py agrega rho con mediana.
    """
    rhos = []
    for s in range(200):
        rng = np.random.default_rng(s)
        theta = rng.uniform(0, TWO_PI, 48)
        est   = (theta + rng.normal(0, 0.2, 48)) % TWO_PI
        rhos.append(circular_correlation(theta, est))
    rhos = np.asarray(rhos)

    assert np.median(rhos) > 0.90, "la mediana de rho deberia aguantar"
    assert (rhos < 0).mean() < 0.05, "la cola no deberia pasar del 5 %"
    # y la cola existe de verdad, no es una precaucion teorica
    assert rhos.min() < 0.0
