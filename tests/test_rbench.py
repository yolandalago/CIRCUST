"""Tests de simulation/rbench.py."""
import numpy as np
import pandas as pd
import pytest

from simulation.generate import generate_dataset
from simulation.rbench import (
    R_CORE_NAMES, consolidar_r, export_for_r, rename_map_for_r,
)
from simulation.truth import ANCHOR_NAME, OPPOSITE_NAME, SimConfig


def _celda():
    cfg = SimConfig(n_samples=24, n_genes=400, prop_rhythmic=0.20, n_core=30,
                    seed=2)
    return generate_dataset(cfg)


def test_export_escribe_los_cuatro_ficheros(tmp_path):
    X, truth = _celda()
    rutas = export_for_r(X, truth, tmp_path)

    assert set(rutas) == {"matriz", "core_genes", "pares", "matriz_r"}
    for ruta in rutas.values():
        assert ruta.exists() and ruta.stat().st_size > 0

    leida = pd.read_csv(rutas["matriz"], index_col=0)
    pd.testing.assert_frame_equal(leida, X, check_exact=False, rtol=1e-9)

    core = rutas["core_genes"].read_text().split()
    assert ANCHOR_NAME in core and OPPOSITE_NAME in core

    pares = pd.read_csv(rutas["pares"])
    assert list(pares.columns) == ["antes", "despues"]
    assert len(pares) == len(truth.consistency_pairs())


def test_la_copia_para_r_solo_cambia_las_etiquetas(tmp_path):
    X, truth = _celda()
    rutas = export_for_r(X, truth, tmp_path)

    para_r = pd.read_csv(rutas["matriz_r"], index_col=0)
    assert para_r.shape == X.shape
    assert np.allclose(para_r.values, X.values)
    for nombre in R_CORE_NAMES:
        assert nombre in para_r.index, nombre


def test_el_renombrado_preserva_dbp_menor_que_cry1_menor_que_arntl():
    _, truth = _celda()
    mapa    = rename_map_for_r(truth)
    peak_of = dict(zip(truth.core_names, truth.core_peaks))

    assert mapa[ANCHOR_NAME] == "ARNTL"
    assert mapa[OPPOSITE_NAME] == "DBP"
    intermedio = truth.chain[1]
    assert mapa[intermedio] == "CRY1"
    assert peak_of[OPPOSITE_NAME] < peak_of[intermedio] < peak_of[ANCHOR_NAME]
    assert len(set(mapa.values())) == len(mapa)


def test_el_renombrado_exige_al_menos_doce_core():
    cfg = SimConfig(n_samples=24, n_genes=400, prop_rhythmic=0.20, n_core=10,
                    core_sizes=(10,), seed=2)
    _, truth = generate_dataset(cfg)
    with pytest.raises(ValueError, match="al menos 12"):
        rename_map_for_r(truth)


def test_consolidar_r_pasa_a_formato_ancho(tmp_path):
    celda = tmp_path / "g1000_n24_r2"
    celda.mkdir()
    pd.DataFrame({
        "rep": [1] * 5,
        "stage": ["preprocess", "preliminary_order", "top_genes",
                  "robust_order", "elapsed_total"],
        "seconds": [1.0, 2.0, 3.0, 4.0, 10.0],
        "implementacion": ["r"] * 5,
        "celda": ["g1000_n24_r2"] * 5,
    }).to_csv(celda / "timing_r.csv", index=False)

    t = consolidar_r(tmp_path)
    assert len(t) == 1
    fila = t.iloc[0]
    assert fila["implementacion"] == "r"
    assert (fila["n_genes"], fila["n_samples"], fila["rep"]) == (1000, 24, 2)
    assert fila["t_total_s"] == 10.0
    assert fila["t_robust_order_s"] == 4.0


def test_consolidar_r_sin_ficheros_devuelve_vacio(tmp_path):
    assert consolidar_r(tmp_path).empty
