"""
Desempate determinista de la rejilla del FMM base (``_best_step1``).

Con omega = 1 el modelo FMM es un coseno y todos los alpha de la rejilla dan la
misma curva, de modo que sus RSS son iguales en aritmetica exacta y solo
difieren en el ultimo bit. Ese bit depende del camino de calculo (Numba o
NumPy), y antes decidia el punto de partida de Nelder-Mead: sobre 2000 genes
reales, el pico de un 1-2 % de los genes del TOP se movia hasta 0,76 h entre
los dos caminos.
"""
import numpy as np
import pytest

import circust.fitting.fmm as fmm


def test_entre_candidatos_empatados_gana_el_primero_de_la_rejilla():
    # Mismo candidato dos veces, con el segundo 1e-16 relativo por debajo:
    # sin redondeo ganaria el segundo; con el desempate, el primero.
    fila = np.array([0.0, 1.0, 0.5, 2.0, 1.0, 0.25])
    rejilla = np.vstack([fila, fila, fila])
    rejilla[1, 2] = 1.5
    rejilla[1, 5] = 0.25 * (1 - 1e-16)
    rejilla[2, 2] = 2.5
    rejilla[2, 5] = 0.40                       # claramente peor
    datos = np.linspace(-1, 1, 30)
    elegido = fmm._best_step1(datos, rejilla)
    assert elegido[2] == 0.5


@pytest.mark.skipif(not fmm._HAS_NUMBA, reason="hace falta Numba para comparar")
def test_numba_y_numpy_dan_el_mismo_pico_en_ondas_casi_cosenoidales(monkeypatch):
    rng = np.random.default_rng(3)
    n = 60
    t = np.sort(rng.uniform(0, 2 * np.pi, n))
    genes = []
    for _ in range(40):
        a, b = rng.uniform(0, 2 * np.pi, 2)
        w = rng.choice([1.0, 0.95, 0.9])
        y = np.cos(b + 2 * np.arctan(w * np.tan((t - a) / 2)))
        genes.append(y + rng.normal(0, 0.3, n))

    con_numba = [fmm.FMMModel().fit(g, t).peak_time for g in genes]
    monkeypatch.setattr(fmm, "_HAS_NUMBA", False)    # rejilla y objetivo en NumPy
    con_numpy = [fmm.FMMModel().fit(g, t).peak_time for g in genes]

    dif_h = np.abs((np.array(con_numba) - np.array(con_numpy) + np.pi)
                   % (2 * np.pi) - np.pi) * 24 / (2 * np.pi)
    assert dif_h.max() < 0.01
