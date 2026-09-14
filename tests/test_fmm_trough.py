"""
Tests del tiempo de valle del modelo FMM.

El pico ``t_U`` es donde el argumento del coseno vale 0, y el valle ``t_L``
donde vale pi. El codigo calculaba el valle como ``t_U + pi``, o sea a media
vuelta del pico, lo cual solo es cierto para una sinusoide: la onda FMM
compone el coseno con una deformacion del eje temporal, y al deshacerla esa
media vuelta deja de ser media vuelta en tiempo de reloj.
"""
import numpy as np
import pytest

from circust.fitting.fmm import FMMModel

TWO_PI = 2.0 * np.pi
REJILLA = np.linspace(0.0, TWO_PI, 200_000, endpoint=False)


def _onda(theta, M, A, alpha, beta, omega):
    mob = 2.0 * np.arctan(omega * np.tan((np.asarray(theta) - alpha) / 2.0))
    return M + A * np.cos(beta + mob)


def _dcirc(a, b):
    d = abs(a - b) % TWO_PI
    return min(d, TWO_PI - d)


def test_el_valle_cae_en_el_minimo_de_la_onda():
    """Comprobacion directa contra el argmin sobre una rejilla fina."""
    rng = np.random.default_rng(0)
    for _ in range(200):
        a, b = rng.uniform(0, TWO_PI, 2)
        w = float(rng.uniform(0.1, 1.0))
        real = REJILLA[int(np.argmin(_onda(REJILLA, 0.0, 1.0, a, b, w)))]
        assert _dcirc(FMMModel.trough_time(a, b, w), real) < 1e-4


def test_en_el_valle_la_onda_vale_M_menos_A_y_en_el_pico_M_mas_A():
    """El argumento del coseno vale pi en el valle y 0 en el pico."""
    rng = np.random.default_rng(1)
    for _ in range(200):
        a, b = rng.uniform(0, TWO_PI, 2)
        w = float(rng.uniform(0.1, 1.0))
        M, A = rng.normal(0, 2), float(rng.uniform(0.5, 3))
        tL = FMMModel.trough_time(a, b, w)
        tU = FMMModel.peak_time(a, b, w)
        assert _onda(tL, M, A, a, b, w) == pytest.approx(M - A, abs=1e-9)
        assert _onda(tU, M, A, a, b, w) == pytest.approx(M + A, abs=1e-9)


def test_solo_con_omega_uno_el_valle_esta_a_media_vuelta_del_pico():
    """
    La razon del bug. Con omega=1 no hay deformacion y ``t_U + pi`` acierta.
    Cuanto mas apuntada la onda, mas se equivoca: hasta 7,6 h con omega=0,2.
    """
    rng = np.random.default_rng(2)
    errores = {}
    for w in (0.2, 0.5, 1.0):
        e = []
        for _ in range(80):
            a, b = rng.uniform(0, TWO_PI, 2)
            tU = FMMModel.peak_time(a, b, w)
            e.append(_dcirc((tU + np.pi) % TWO_PI,
                            FMMModel.trough_time(a, b, w)) * 24 / TWO_PI)
        errores[w] = float(np.median(e))

    assert errores[1.0] < 1e-6, "con omega=1 la regla antigua es exacta"
    assert errores[0.5] > 2.0, "con omega=0,5 se equivoca en horas"
    assert errores[0.2] > 6.0, "con omega=0,2 el error supera las 6 h"


def test_la_rejilla_de_la_etapa_4_usa_el_valle_correcto():
    """
    robust_order calculaba el valle como pico mas media vuelta al rellenar la
    columna 6 del stats_table, que alimenta ``trough_phase_rad``.
    """
    import inspect
    import circust.robust_order as ro
    fuente = inspect.getsource(ro)
    assert "pkL = (pkU + pi)" not in fuente, (
        "el valle se sigue calculando como pico mas media vuelta"
    )
    assert "trough_time" in fuente
