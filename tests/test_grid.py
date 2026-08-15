"""Tests de la definicion de la rejilla.

Es codigo de datos, sin pipeline, de modo que los tests son rapidos y pueden
ser exhaustivos. Vigilan sobre todo tres cosas que si se rompen no dan error
sino resultados silenciosamente mal: que las semillas no se repitan, que las
rutas no colisionen, y que los valores anclados en datos reales sigan siendo
los que se midieron.
"""
from collections import Counter

import pytest

from simulation.generate import SAMPLING_SCHEMES, WAVEFORM_MIXES
from simulation.grid import (
    GRID_SNR, REFERENCIA, REPLICAS, SEED_BASE, bloque_a, bloque_fase0,
    bloque_p, bloque_r, bloque_tau, materializar, ruta_de, todas_las_celdas,
    vecinos_snr,
)


# ===========================================================================
# Celda de referencia
# ===========================================================================

def test_la_referencia_usa_los_tamanos_anclados_en_datos_reales():
    # 20.000 genes es el orden de los expresados y evita que las K
    # repeticiones degeneren, que es lo que dejaria al bloque A sin nada que
    # agregar.
    assert REFERENCIA["n_genes"] == 20_000
    assert REFERENCIA["prop_rhythmic"] == 0.05
    assert REFERENCIA["waveform_mix"] == "gtex"


def test_el_muestreo_de_referencia_no_es_degenerado():
    assert REFERENCIA["sampling"] in SAMPLING_SCHEMES
    assert REFERENCIA["sampling"] != "equiespaciado"


# ===========================================================================
# Estructura de los bloques
# ===========================================================================

def test_la_fase_cero_cruza_los_cuatro_factores():
    celdas = bloque_fase0()
    assert len(celdas) == 36
    assert {c.generate["snr"] for c in celdas} == {0.5, 1.0, 2.0}
    assert {c.generate["n_samples"] for c in celdas} == {12, 48}


def test_la_fase_cero_mueve_los_factores_que_gobiernan_el_desacuerdo():
    """Regresion del error de diseno del 14/08/2026.

    La primera fase 0 solo cruzaba SNR por numero de muestras y dejaba
    congelados la proporcion de ritmicos y el esquema de muestreo, que son
    justamente los que mueven el desacuerdo entre las K repeticiones. Ninguna
    celda cumplio el criterio porque la rejilla no contenia los mandos
    relevantes.
    """
    celdas = bloque_fase0()
    assert len({c.generate["prop_rhythmic"] for c in celdas}) > 1
    assert len({c.generate["sampling"] for c in celdas}) > 1


def test_la_fase_cero_conserva_un_ancla_del_regimen_saturado():
    # Con SNR 2 el desacuerdo medido baja a 0,002. Se deja en la rejilla para
    # que el criterio tenga que rechazarlo de forma visible.
    assert 2.0 in {c.generate["snr"] for c in bloque_fase0()}


def test_el_bloque_p_llega_a_la_escala_de_gtex():
    celdas = bloque_p()
    assert len(celdas) == 12
    assert (56_200, 480) in {(c.generate["n_genes"], c.generate["n_samples"])
                             for c in celdas}
    assert 12 in {c.generate["n_samples"] for c in celdas}


def test_el_bloque_p_conserva_las_matrices_para_R():
    for c in bloque_p():
        assert c.conservar_matriz
        assert c.tambien_csv


def test_los_demas_bloques_no_conservan_las_matrices():
    for c in bloque_fase0() + bloque_r() + bloque_a():
        assert not c.conservar_matriz


def test_el_bloque_r_cruza_muestreo_por_mezcla_por_tamano_muestral():
    celdas = bloque_r()
    assert len(celdas) == 48
    assert {c.generate["sampling"] for c in celdas} == set(SAMPLING_SCHEMES)
    assert {c.generate["waveform_mix"] for c in celdas} == set(WAVEFORM_MIXES)
    assert {c.generate["n_samples"] for c in celdas} == {12, 48, 192}


def test_el_bloque_r_usa_su_propia_ancla_y_no_la_calibrada():
    """El ancla de la fase 0 se calibro para el bloque A, no para este.

    Con SNR 1 y 1 por ciento de ritmicos las celdas de 192 muestras no corren,
    porque el conjunto TOP colapsa, de modo que heredarla dejaria en blanco un
    tercio del bloque. Sondeado el 15/08/2026.
    """
    calibrada = {**REFERENCIA, "snr": 1.0, "prop_rhythmic": 0.01}
    for c in bloque_r(calibrada):
        assert c.generate["snr"] == 2.0
        assert c.generate["prop_rhythmic"] == 0.05


def test_el_bloque_r_congela_el_resumen_robusto_y_mueve_el_fmm():
    for c in bloque_r():
        assert c.fmm_methods == ("base", "stabilized")
        assert c.metodos == ("best_k:r2_median",)


def test_el_bloque_a_congela_el_fmm_y_mueve_el_resumen():
    celdas = bloque_a()
    assert len(celdas) == 27
    for c in celdas:
        assert c.fmm_methods == ("base",)
        assert len(c.metodos) == 6


def test_el_bloque_a_usa_el_snr_calibrado_y_sus_vecinos():
    celdas = bloque_a({**REFERENCIA, "snr": 2.0})
    assert {c.generate["snr"] for c in celdas} == {1.0, 2.0, 3.0}


def test_los_vecinos_de_snr_no_se_salen_de_la_rejilla():
    assert vecinos_snr(0.25) == (0.25, 0.5, 1.0)
    assert vecinos_snr(3.0) == (1.0, 2.0, 3.0)
    assert vecinos_snr(1.0) == (0.5, 1.0, 2.0)


def test_un_snr_fuera_de_la_rejilla_falla():
    with pytest.raises(ValueError, match="0.7"):
        vecinos_snr(0.7)


# ===========================================================================
# Barrido de tau
# ===========================================================================

def test_el_barrido_de_tau_solo_toca_el_criterio_afectado():
    celdas = bloque_tau()
    etiquetas = {e for c in celdas for e, _ in c.metodos}
    assert all("n_well_fitted" in e for e in etiquetas)
    assert len(etiquetas) == 5


def test_el_barrido_de_tau_lleva_el_umbral_en_los_kwargs():
    umbrales = sorted({kw["selection_r2_threshold"]
                       for c in bloque_tau() for _, kw in c.metodos})
    assert umbrales == [0.60, 0.65, 0.70, 0.75, 0.80]


def test_el_barrido_de_tau_usa_semillas_de_calibracion():
    # Calibrar y evaluar con las mismas semillas seria sobreajuste.
    calibracion = {c.seed_base for c in bloque_fase0() + bloque_tau()}
    evaluacion  = {c.seed_base for c in bloque_r() + bloque_a()}
    assert not (calibracion & evaluacion)


# ===========================================================================
# Semillas y rutas
# ===========================================================================

def test_ninguna_semilla_se_repite_en_toda_la_rejilla():
    semillas = [c.seed_base + r
                for c in todas_las_celdas()
                for r in range(REPLICAS[c.bloque])]
    repetidas = [s for s, n in Counter(semillas).items() if n > 1]
    assert repetidas == []


def test_ninguna_ruta_se_repite_en_toda_la_rejilla():
    rutas = [ruta_de(c, r) for c in todas_las_celdas()
             for r in range(REPLICAS[c.bloque])]
    repetidas = [p for p, n in Counter(rutas).items() if n > 1]
    assert repetidas == []


def test_la_ruta_separa_bloque_celda_y_replica():
    c = bloque_fase0()[0]
    partes = ruta_de(c, 3).parts
    assert partes[-3] == c.bloque
    assert partes[-2] == c.nombre
    assert partes[-1] == "rep03"


def test_cada_bloque_tiene_su_propia_base_de_semillas():
    assert len(set(SEED_BASE.values())) == len(SEED_BASE)
    for c in todas_las_celdas():
        assert c.seed_base >= SEED_BASE[c.bloque]


# ===========================================================================
# Materializacion
# ===========================================================================

def test_materializar_respeta_los_parametros_de_la_celda():
    c = [x for x in bloque_r()
         if x.generate["sampling"] == "huecos"
         and x.generate["waveform_mix"] == "extrema"][0]
    ds = materializar(c, 0, n_genes=800)
    assert ds.meta["sampling"] == "huecos"
    assert ds.meta["waveform_mix"] == "extrema"
    assert ds.meta["n_genes"] == 800


def test_materializar_es_reproducible_y_distinto_entre_replicas():
    import numpy as np
    # 2.000 genes al 1 por ciento dan 20 ritmicos, suficientes para los 12
    # core que el CPCA necesita.
    c = bloque_fase0()[0]
    a = materializar(c, 0, n_genes=2000)
    b = materializar(c, 0, n_genes=2000)
    d = materializar(c, 1, n_genes=2000)
    assert np.array_equal(a.matrix.to_numpy(), b.matrix.to_numpy())
    assert not np.array_equal(a.matrix.to_numpy(), d.matrix.to_numpy())


def test_el_nombre_de_la_celda_es_legible_y_valido_como_directorio():
    for c in todas_las_celdas():
        assert c.nombre
        assert not set(c.nombre) & set("/\\ :")


# ===========================================================================
# Coste
# ===========================================================================

def test_el_numero_de_replicas_es_el_acordado():
    assert REPLICAS["fase0"] == 5
    assert REPLICAS["P"] == 3
    assert REPLICAS["R"] == 20
    assert REPLICAS["A"] == 20
