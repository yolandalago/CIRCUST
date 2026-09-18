"""Tests de simulation/metrics.py."""
import numpy as np
import pytest

from simulation.metrics import (
    circular_correlation, circular_distance, evaluate,
    fisher_lee_correlation, kendall_tau_fixed_cut, error_orden_h,
    tau_circular, to_hours,
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
    assert error_orden_h(theta, estimado) < 0.1
    assert error_orden_h(theta, theta) == pytest.approx(0.0, abs=1e-9)


def test_evaluate_devuelve_todas_las_claves():
    theta = np.random.default_rng(2).uniform(0, TWO_PI, 24)
    m = evaluate(theta, theta)
    esperadas = {"error_fase_medio_h", "error_fase_h", "rho_circ", "rho_fl", "tau",
                 "tau_circular", "acierto_sentido", "p1", "p2", "p4",
                 "error_orden_h"}
    assert set(m) == esperadas
    assert m["error_fase_medio_h"] == pytest.approx(0.0)
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


def test_acierto_sentido_marca_el_sentido_de_giro():
    """
    Columna separada del error porque son cosas distintas: el orden puede
    estar bien recuperado y el sentido al reves. Promediada sobre replicas da
    directamente la tasa de acierto de orientacion, que es lo que el capitulo
    necesita para no mezclar reconstruccion con suerte.
    """
    rng = np.random.default_rng(5)
    theta = np.sort(rng.uniform(0, TWO_PI, 40))
    ruido = (theta + rng.normal(0, 0.2, 40)) % TWO_PI

    assert evaluate(theta, ruido)["acierto_sentido"] == 1.0
    assert evaluate(theta, (-ruido) % TWO_PI)["acierto_sentido"] == 0.0
    # rotar no cambia el sentido de giro
    assert evaluate(theta, (ruido + 2.5) % TWO_PI)["acierto_sentido"] == 1.0


# ---------------------------------------------------------------------------
# Correlacion circular de Fisher y Lee
# ---------------------------------------------------------------------------

def test_fisher_lee_funciona_donde_jammalamadaka_no_esta_definida():
    """
    La razon de adoptarla. Jammalamadaka resta a cada serie su direccion
    media, y con muestreo exactamente equiespaciado esa direccion no existe.
    Fisher-Lee se construye sobre diferencias por pares y no necesita
    ninguna media.
    """
    rejilla = _ang(np.linspace(0, 24, 48, endpoint=False))
    est = (rejilla + np.random.default_rng(0).normal(0, 0.2, 48)) % TWO_PI
    assert np.isnan(circular_correlation(rejilla, est))
    assert fisher_lee_correlation(rejilla, est) > 0.9


def test_fisher_lee_conserva_las_propiedades_necesarias():
    theta = np.sort(np.random.default_rng(0).uniform(0, TWO_PI, 40))
    assert fisher_lee_correlation(theta, theta) == pytest.approx(1.0)
    # invariante a rotar cualquiera de las dos escalas
    for d in (0.4, 1.5, 3.0, 5.7):
        assert fisher_lee_correlation(theta, (theta + d) % TWO_PI) == pytest.approx(1.0)
    # antisimetrica bajo reflexion, incluso combinada con rotacion
    for d in (0.0, 2.0, 4.5):
        assert fisher_lee_correlation(theta, (-theta + d) % TWO_PI) == pytest.approx(-1.0)
    # sin relacion
    aleatoria = np.random.default_rng(1).uniform(0, TWO_PI, 40)
    assert abs(fisher_lee_correlation(theta, aleatoria)) < 0.2


def test_fisher_lee_no_cambia_de_signo_donde_jammalamadaka_si():
    """
    Medido: sobre 400 reconstrucciones perfectas con muestreo uniforme,
    Jammalamadaka cambia de signo en el 1,2 % de los casos y Fisher-Lee en
    ninguno.
    """
    jam, fl = [], []
    for s in range(300):
        r = np.random.default_rng(s)
        th = r.uniform(0, TWO_PI, 48)
        est = (th + r.normal(0, 0.2, 48)) % TWO_PI
        jam.append(circular_correlation(th, est))
        fl.append(fisher_lee_correlation(th, est))
    jam, fl = np.asarray(jam), np.asarray(fl)
    assert (fl < 0).sum() == 0, "Fisher-Lee no deberia equivocarse de signo"
    assert np.median(fl) > 0.9
    # y la fragilidad de la otra existe de verdad, no es una precaucion
    assert (jam < 0).sum() > 0


def test_fisher_lee_con_menos_de_dos_muestras():
    assert np.isnan(fisher_lee_correlation(np.array([0.0]), np.array([0.0])))


def _fisher_lee_doble_suma(a, b):
    """
    Definicion literal de Fisher y Lee, en O(n^2). Vive aqui como referencia
    contra la que verificar el algebra de la forma computacional que usa el
    modulo, que es O(n).
    """
    a = np.asarray(a, float); b = np.asarray(b, float)
    iu = np.triu_indices(a.size, 1)
    da = np.sin(a[:, None] - a[None, :])[iu]
    db = np.sin(b[:, None] - b[None, :])[iu]
    den = np.sqrt(np.sum(da**2) * np.sum(db**2))
    return float(np.sum(da*db)/den) if den else float("nan")


def test_la_forma_computacional_es_la_definicion_literal():
    """
    El modulo usa la descomposicion en sumas simples de Fisher (1993):

        rho = 4(AB - CD) / sqrt[(n^2 - E^2 - F^2)(n^2 - G^2 - H^2)]

    que es algebraicamente identica a la doble suma pero cuesta O(n) en vez
    de O(n^2). Con n=480, que es el tamano del tejido real de GTEx, la
    diferencia es de 108 us contra 11,6 ms.
    """
    rng = np.random.default_rng(0)
    for _ in range(300):
        n = int(rng.integers(5, 150))
        a, b = rng.uniform(0, TWO_PI, n), rng.uniform(0, TWO_PI, n)
        assert fisher_lee_correlation(a, b) == pytest.approx(
            _fisher_lee_doble_suma(a, b), abs=1e-12)
