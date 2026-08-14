"""Tests del volcado de celdas a disco.

La pieza que traduce un ``SyntheticDataset`` en memoria a la forma de unidad
que el evaluador de ``validation/`` espera encontrar en disco.
"""
import json

import numpy as np
import pandas as pd
import pytest

from circust.preprocessing import load_expression_matrix
from simulation.cells import (
    GENE_COLUMN, drop_matrix, matrix_on_disk, read_truth, regenerate,
    write_cell,
)
from simulation.generate import generate


@pytest.fixture
def ds():
    return generate("celda", n_samples=8, n_genes=600, prop_rhythmic=0.05,
                    sampling="aleatorio", waveform_mix="gtex", seed=4)


# ===========================================================================
# Escritura
# ===========================================================================

def test_write_cell_deja_los_cuatro_ficheros(ds, tmp_path):
    ruta = write_cell(ds, tmp_path / "rep00")
    assert ruta.name == "matriz.parquet"
    assert ruta.exists()
    for nombre in ("tiempos.csv", "params.csv", "meta.json"):
        assert (tmp_path / "rep00" / nombre).exists(), nombre


def test_write_cell_crea_el_directorio_si_no_existe(ds, tmp_path):
    ruta = write_cell(ds, tmp_path / "a" / "b" / "rep00")
    assert ruta.exists()


@pytest.mark.parametrize("formato", ["parquet", "csv"])
def test_la_matriz_se_relee_identica_con_el_cargador_del_pipeline(
        ds, tmp_path, formato):
    ruta = write_cell(ds, tmp_path / "rep00", formato=formato)
    leida = load_expression_matrix(str(ruta), gene_column=GENE_COLUMN)
    assert list(leida.index) == list(ds.matrix.index)
    assert list(leida.columns) == list(ds.matrix.columns)
    assert np.allclose(leida.to_numpy(), ds.matrix.to_numpy())


def test_tambien_csv_escribe_los_dos_formatos(ds, tmp_path):
    ruta = write_cell(ds, tmp_path / "rep00", tambien_csv=True)
    assert ruta.suffix == ".parquet"
    assert (tmp_path / "rep00" / "matriz.csv").exists()


def test_los_tiempos_siguen_el_formato_de_los_datos_reales(ds, tmp_path):
    write_cell(ds, tmp_path / "rep00")
    t = pd.read_csv(tmp_path / "rep00" / "tiempos.csv")
    assert list(t.columns) == ["sample_name", "true_time_h"]
    assert list(t["sample_name"]) == list(ds.matrix.columns)
    assert np.allclose(t["true_time_h"], ds.true_times_h)


def test_params_guarda_solo_los_ritmicos(ds, tmp_path):
    write_cell(ds, tmp_path / "rep00")
    p = pd.read_csv(tmp_path / "rep00" / "params.csv")
    esperados = ds.gene_params[ds.gene_params["ritmico"]]
    assert len(p) == len(esperados)
    assert p["ritmico"].all()
    for col in ("gene", "clase", "alpha", "beta", "omega", "t_U_h"):
        assert col in p.columns, col


def test_el_meta_permite_saber_cuantos_genes_de_fondo_habia(ds, tmp_path):
    write_cell(ds, tmp_path / "rep00")
    meta = json.loads((tmp_path / "rep00" / "meta.json").read_text())
    assert meta["n_genes"] == 600
    assert meta["n_rhythmic"] == 30
    assert meta["name"] == "celda"
    assert meta["sampling"] == "aleatorio"
    assert meta["waveform_mix"] == "gtex"


# ===========================================================================
# Lectura del ground truth
# ===========================================================================

def test_read_truth_devuelve_las_dos_verdades_y_el_manifiesto(ds, tmp_path):
    write_cell(ds, tmp_path / "rep00")
    verdad = read_truth(tmp_path / "rep00")
    assert list(verdad.tiempos.index) == list(ds.matrix.columns)
    assert np.allclose(verdad.tiempos.to_numpy(), ds.true_times_h)
    assert verdad.params["ritmico"].all()
    assert verdad.meta["seed"] == 4


def test_read_truth_funciona_aunque_la_matriz_ya_no_este(ds, tmp_path):
    # Es el caso normal en la rejilla: la matriz se borra tras evaluar y el
    # ground truth se conserva para el analisis posterior.
    write_cell(ds, tmp_path / "rep00")
    drop_matrix(tmp_path / "rep00")
    verdad = read_truth(tmp_path / "rep00")
    assert len(verdad.params) == 30


# ===========================================================================
# Borrado y ciclo efimero
# ===========================================================================

def test_drop_matrix_borra_la_matriz_y_conserva_el_resto(ds, tmp_path):
    write_cell(ds, tmp_path / "rep00", tambien_csv=True)
    drop_matrix(tmp_path / "rep00")
    quedan = sorted(p.name for p in (tmp_path / "rep00").iterdir())
    assert quedan == ["meta.json", "params.csv", "tiempos.csv"]


def test_drop_matrix_no_falla_si_ya_no_hay_matriz(ds, tmp_path):
    write_cell(ds, tmp_path / "rep00")
    drop_matrix(tmp_path / "rep00")
    drop_matrix(tmp_path / "rep00")


def test_el_ciclo_efimero_borra_la_matriz_al_salir(ds, tmp_path):
    with matrix_on_disk(ds, tmp_path / "rep00") as ruta:
        assert ruta.exists()
    assert not ruta.exists()
    assert (tmp_path / "rep00" / "params.csv").exists()


def test_el_ciclo_efimero_conserva_la_matriz_si_se_pide(ds, tmp_path):
    with matrix_on_disk(ds, tmp_path / "rep00", conservar=True) as ruta:
        pass
    assert ruta.exists()


def test_el_ciclo_efimero_borra_la_matriz_aunque_falle_la_evaluacion(
        ds, tmp_path):
    with pytest.raises(RuntimeError):
        with matrix_on_disk(ds, tmp_path / "rep00") as ruta:
            raise RuntimeError("el evaluador se cayo")
    assert not ruta.exists()


# ===========================================================================
# Reproducibilidad
# ===========================================================================

def test_regenerate_reconstruye_la_celda_bit_a_bit(ds, tmp_path):
    write_cell(ds, tmp_path / "rep00")
    otra = regenerate(tmp_path / "rep00")
    assert np.array_equal(otra.matrix.to_numpy(), ds.matrix.to_numpy())
    assert np.array_equal(otra.true_times_h, ds.true_times_h)


def test_regenerate_respeta_los_parametros_del_esquema_de_muestreo(tmp_path):
    # Los sampling_kwargs no son argumentos con nombre de generate, viajan en
    # **kwargs. Si no quedasen anotados en el manifiesto, la celda regenerada
    # no seria la misma.
    ds = generate("hueco", n_samples=12, n_genes=600, sampling="huecos",
                  gap=(3.0, 11.0), seed=1)
    write_cell(ds, tmp_path / "rep00")
    otra = regenerate(tmp_path / "rep00")
    assert np.array_equal(otra.true_times_h, ds.true_times_h)


def test_el_manifiesto_lleva_una_huella_de_la_matriz(ds, tmp_path):
    write_cell(ds, tmp_path / "rep00")
    meta = json.loads((tmp_path / "rep00" / "meta.json").read_text())
    assert len(meta["matriz_sha1"]) == 40
    otro = write_cell(ds, tmp_path / "rep01")
    meta_b = json.loads((tmp_path / "rep01" / "meta.json").read_text())
    assert meta["matriz_sha1"] == meta_b["matriz_sha1"]
    assert otro.exists()
