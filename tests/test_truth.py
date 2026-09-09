"""Tests de simulation/truth.py."""
import numpy as np
import pytest

from simulation.truth import (
    ANCHOR_NAME, OPPOSITE_NAME, TWO_PI,
    core_layout, nested_core_sets, consistency_chain,
    consistency_pairs_from_chain, _spread_order,
)


def test_anclas_en_pi_y_cero():
    names, peaks = core_layout(30, 0.60)
    peak_of = dict(zip(names, peaks))
    assert peak_of[ANCHOR_NAME] == pytest.approx(np.pi)
    assert peak_of[OPPOSITE_NAME] == pytest.approx(0.0)


def test_reparto_dia_noche_sigue_day_fraction():
    names, peaks = core_layout(30, 0.60)
    dia   = [p for p in peaks if 0.0 <= p < np.pi]
    noche = [p for p in peaks if np.pi <= p < TWO_PI]
    assert len(dia) == 18      # CORE_B mas 17
    assert len(noche) == 12    # CORE_A mas 11


def test_picos_dentro_de_su_semicirculo_y_distintos():
    names, peaks = core_layout(30, 0.60)
    assert len(set(np.round(peaks, 12))) == 30
    for nm, p in zip(names, peaks):
        if nm.startswith("CORE_D"):
            assert 0.0 < p < np.pi
        if nm.startswith("CORE_N"):
            assert np.pi < p < TWO_PI


def test_spread_order_es_permutacion():
    for total in (1, 2, 5, 11, 17):
        order = _spread_order(total)
        assert sorted(order) == list(range(total))


def test_conjuntos_core_anidados_y_con_anclas():
    names, peaks = core_layout(30, 0.60)
    sets = nested_core_sets(names, (10, 20, 30), 0.60)
    assert [len(sets[k]) for k in (10, 20, 30)] == [10, 20, 30]
    assert set(sets[10]) <= set(sets[20]) <= set(sets[30])
    for k in (10, 20, 30):
        assert ANCHOR_NAME in sets[k] and OPPOSITE_NAME in sets[k]


def test_conjuntos_core_conservan_la_proporcion():
    names, peaks = core_layout(30, 0.60)
    peak_of = dict(zip(names, peaks))
    sets = nested_core_sets(names, (10, 20, 30), 0.60)
    esperado = {10: (6, 4), 20: (12, 8), 30: (18, 12)}
    for k, (n_dia, n_noche) in esperado.items():
        dia = sum(1 for g in sets[k] if peak_of[g] < np.pi)
        assert (dia, len(sets[k]) - dia) == (n_dia, n_noche)


def test_cadena_de_consistencia_va_de_cero_a_pi():
    names, peaks = core_layout(30, 0.60)
    peak_of = dict(zip(names, peaks))
    sets = nested_core_sets(names, (10, 20, 30), 0.60)
    chain = consistency_chain(names, peaks, sets, n_chain=3)
    assert chain[0] == OPPOSITE_NAME
    assert chain[-1] == ANCHOR_NAME
    assert len(chain) == 3
    picos = [peak_of[g] for g in chain]
    assert picos == sorted(picos)
    # los intermedios estan en el conjunto mas pequeno, para estar
    # disponibles en los tres niveles de G_core
    assert all(g in sets[10] for g in chain)


def test_pares_son_consecutivos_de_la_cadena():
    pares = consistency_pairs_from_chain(["A", "B", "C", "D"])
    assert pares == [("A", "B"), ("B", "C"), ("C", "D")]
