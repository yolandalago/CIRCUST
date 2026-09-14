"""
Tests de la correccion de orientaciones de la Etapa 4.

``_fix_orientations`` detecta las repeticiones cuyo sentido de recorrido
salio invertido y las corrige antes de que se usen para nada. Decidia por el
**signo** de la correlacion de Jammalamadaka, y ese signo no es fiable cuando
las muestras cubren todo el ciclo: sobre repeticiones genuinamente invertidas
del pipeline se midieron valores positivos de +0,11 a +0,36, con lo que la
condicion ``corr < 0`` no se cumplia nunca y no se corregia ninguna.

Se usa la de Fisher y Lee, que se construye sobre diferencias por pares y no
necesita ninguna direccion media. Medido sobre 720 escenarios de inversion
conocida: Jammalamadaka acierta el signo el 95,3 % de las veces y Fisher-Lee
el 100 %. Frente a la alternativa por tripletes, que tambien acierta el
100 %, es 49 veces mas barata con n=96 (0,22 ms contra 10,63 ms).
"""
import numpy as np
import pytest

from circust.robust_order import (
    RobustOrderEstimator, _build_theta, _fisher_lee,
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


def test_fisher_lee_identifica_la_rep_invertida():
    theta_pre, ordenes, escalas = _escenario()
    Theta = _build_theta(ordenes, escalas, theta_pre.size)
    assert _fisher_lee(Theta[0], theta_pre) > 0.8
    assert _fisher_lee(Theta[1], theta_pre) < -0.8


def test_fix_orientations_corrige_la_rep_invertida():
    theta_pre, ordenes, escalas = _escenario()
    est = RobustOrderEstimator(n_reps=2, verbose=False)
    o_fix, e_fix = est._fix_orientations(ordenes, escalas, theta_pre,
                                         theta_pre.size)
    Theta = _build_theta(o_fix, e_fix, theta_pre.size)
    assert _fisher_lee(Theta[0], theta_pre) > 0.8
    assert _fisher_lee(Theta[1], theta_pre) > 0.8, (
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


def test_fisher_lee_es_invariante_a_rotacion():
    """El origen de cada repeticion es arbitrario y no debe influir."""
    rng = np.random.default_rng(2)
    theta = np.sort(rng.uniform(0, TWO_PI, 40))
    for delta in (0.0, 1.0, 3.0, 5.5):
        assert _fisher_lee(theta, (theta + delta) % TWO_PI) == pytest.approx(1.0)


def test_fisher_lee_con_una_sola_muestra():
    assert np.isnan(_fisher_lee(np.array([0.0]), np.array([0.0])))


# ---------------------------------------------------------------------------
# Calibracion de lambda del FMM estabilizado
# ---------------------------------------------------------------------------

def test_lambda_no_se_recalibra_por_repeticion_por_defecto():
    """
    Recalibrar lambda en cada una de las K repeticiones era el valor por
    defecto y costaba el 99 % del sobrecoste del FMM estabilizado: 27 s por
    calibracion, seis por celda, y el bloque de regularizacion pasaba de 4 a
    20 horas.

    Y no aportaba nada, medido:

    1. La calibracion no determina lambda de forma estable. Sobre la MISMA
       escala, cambiando solo la semilla del Monte Carlo, devuelve 1,23,
       1,23, 4,12, 4,12 y 8,76. Esa dispersion (193 %) es MAYOR que la que
       hay entre repeticiones distintas (150 %), asi que lo que varia entre
       repeticiones es ruido del estimador, no informacion sobre los datos.
    2. El valor tampoco cambia el resultado donde el metodo funciona: con
       omega=0,5 y A/sigma=2, pasar de lambda 1,23 a 8,76 da el mismo orden
       bit a bit, 0,000 h de diferencia de fase.
    """
    est = RobustOrderEstimator(verbose=False)
    assert est.recalibrate_lambda_per_rep is False


def test_run_pipeline_no_recalibra_lambda():
    """El estudio y la herramienta tienen que ejecutar lo mismo."""
    import inspect
    import sys
    from pathlib import Path
    raiz = Path(__file__).resolve().parent.parent
    if str(raiz / "scripts") not in sys.path:
        sys.path.insert(0, str(raiz / "scripts"))
    import run_pipeline
    fuente = inspect.getsource(run_pipeline)
    assert "recalibrate_lambda_per_rep" in fuente, (
        "run_pipeline deberia fijarlo explicitamente, no depender del default"
    )
