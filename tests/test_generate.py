"""Tests del generador sintetico.

Cubren las dos piezas que el estudio de simulacion necesita del generador y
que no existian antes: la mezcla de clases de forma de onda y el esquema de
muestreo aleatorio.
"""
from collections import Counter

import numpy as np
import pytest

from simulation.generate import (
    SAMPLING_SCHEMES, WAVEFORM_CLASSES, WAVEFORM_MIXES,
    assign_waveform_classes, generate, sampling_times,
)


# ===========================================================================
# Clases de forma de onda
# ===========================================================================

def test_las_clases_tienen_los_omega_anclados_en_los_datos_reales():
    # Los cuatro valores son los representantes de los tramos en los que se
    # midio la distribucion empirica de omega sobre GTEx y babuino. Si alguien
    # los cambia sin rehacer esa medicion, este test lo caza.
    assert WAVEFORM_CLASSES == {
        "sinusoidal": 0.85, "moderado": 0.65,
        "apuntado": 0.45, "muy_apuntado": 0.20,
    }


def test_todas_las_mezclas_suman_uno_y_usan_clases_conocidas():
    for nombre, mezcla in WAVEFORM_MIXES.items():
        assert set(mezcla) <= set(WAVEFORM_CLASSES), nombre
        assert abs(sum(mezcla.values()) - 1.0) < 1e-9, nombre


def test_la_mezcla_respeta_las_proporciones_exactamente():
    rng = np.random.default_rng(0)
    clases, omegas = assign_waveform_classes(
        100, {"sinusoidal": 0.2, "moderado": 0.6, "apuntado": 0.2},
        rng, jitter=0.0,
    )
    assert Counter(clases.tolist()) == {
        "sinusoidal": 20, "moderado": 60, "apuntado": 20}
    assert set(np.round(omegas[clases == "sinusoidal"], 9)) == {0.85}
    assert set(np.round(omegas[clases == "moderado"], 9)) == {0.65}


def test_el_reparto_por_resto_mayor_no_pierde_ni_inventa_genes():
    # 7 genes sobre proporciones 0,2 / 0,6 / 0,2 no da enteros. El reparto por
    # resto mayor tiene que devolver exactamente 7, ni 6 ni 8.
    rng = np.random.default_rng(1)
    clases, omegas = assign_waveform_classes(7, "gtex", rng, jitter=0.0)
    assert clases.size == 7
    assert omegas.size == 7
    assert set(clases.tolist()) <= set(WAVEFORM_CLASSES)


def test_la_mezcla_de_una_sola_clase_es_homogenea():
    rng = np.random.default_rng(2)
    clases, omegas = assign_waveform_classes(50, "simetrica", rng, jitter=0.0)
    assert set(clases.tolist()) == {"sinusoidal"}
    assert np.allclose(omegas, 0.85)


def test_el_jitter_no_saca_omega_del_rango_valido():
    rng = np.random.default_rng(3)
    _, omegas = assign_waveform_classes(500, "extrema", rng, jitter=0.5)
    assert np.all(omegas > 0.0)
    assert np.all(omegas <= 1.0)


def test_el_jitter_dispersa_pero_no_desplaza_la_clase():
    rng = np.random.default_rng(4)
    clases, omegas = assign_waveform_classes(5000, "simetrica", rng, jitter=0.05)
    assert omegas.std() > 0.01
    assert abs(np.median(omegas) - 0.85) < 0.02


def test_una_mezcla_desconocida_falla_con_mensaje_claro():
    rng = np.random.default_rng(5)
    with pytest.raises(ValueError, match="no_existe"):
        assign_waveform_classes(10, "no_existe", rng)


def test_una_clase_desconocida_dentro_de_una_mezcla_tambien_falla():
    rng = np.random.default_rng(6)
    with pytest.raises(ValueError, match="cuadrado"):
        assign_waveform_classes(10, {"cuadrado": 1.0}, rng)


def test_una_mezcla_que_no_suma_uno_falla():
    rng = np.random.default_rng(7)
    with pytest.raises(ValueError, match="suma"):
        assign_waveform_classes(10, {"sinusoidal": 0.3, "moderado": 0.3}, rng)


# ===========================================================================
# Muestreo
# ===========================================================================

def test_el_esquema_aleatorio_esta_disponible():
    assert "aleatorio" in SAMPLING_SCHEMES


def test_el_muestreo_aleatorio_cae_en_el_ciclo_y_va_ordenado():
    rng = np.random.default_rng(8)
    t = sampling_times(48, "aleatorio", rng)
    assert t.size == 48
    assert np.all((t >= 0.0) & (t < 24.0))
    assert np.all(np.diff(t) >= 0.0)


def test_el_muestreo_aleatorio_no_es_equiespaciado():
    # Es la condicion no equiespaciada del suplemento de CIRCUST: uniforme al
    # azar. Si saliese regular no estariamos probando nada distinto.
    rng = np.random.default_rng(9)
    t = sampling_times(96, "aleatorio", rng)
    d = np.diff(t)
    assert d.std() > 0.1 * d.mean()


def test_todos_los_esquemas_producen_el_numero_pedido_de_muestras():
    for esquema in SAMPLING_SCHEMES:
        rng = np.random.default_rng(10)
        t = sampling_times(36, esquema, rng)
        assert t.size == 36, esquema
        assert np.all((t >= 0.0) & (t < 24.0)), esquema


# ===========================================================================
# Integracion con generate
# ===========================================================================

def test_generate_propaga_la_clase_a_gene_params():
    ds = generate("t", n_samples=24, n_genes=2000, prop_rhythmic=0.05,
                  waveform_mix="gtex", seed=7)
    ritmicos = ds.gene_params[ds.gene_params["ritmico"]]
    assert set(ritmicos["clase"]) == set(WAVEFORM_MIXES["gtex"])
    fondo = ds.gene_params[~ds.gene_params["ritmico"]]
    assert (fondo["clase"] == "-").all()


def test_generate_deja_constancia_de_la_mezcla_en_el_manifiesto():
    ds = generate("t", n_samples=24, n_genes=2000, waveform_mix="apuntada",
                  seed=7)
    assert ds.meta["waveform_mix"] == "apuntada"


def test_generate_con_omega_escalar_conserva_el_comportamiento_anterior():
    ds = generate("t", n_samples=24, n_genes=2000, prop_rhythmic=0.05,
                  omega=0.9, omega_jitter=0.0, seed=7)
    ritmicos = ds.gene_params[ds.gene_params["ritmico"]]
    assert np.allclose(ritmicos["omega"], 0.9)
    assert set(ritmicos["clase"]) == {"manual"}


def test_generate_rechaza_omega_y_mezcla_a_la_vez():
    with pytest.raises(ValueError, match="waveform_mix"):
        generate("t", n_samples=24, n_genes=2000, omega=0.9,
                 waveform_mix="gtex", seed=7)


def test_la_mezcla_apuntada_produce_omegas_menores_que_la_simetrica():
    ds_sim = generate("s", n_samples=24, n_genes=4000, waveform_mix="simetrica",
                      seed=11)
    ds_apu = generate("a", n_samples=24, n_genes=4000, waveform_mix="extrema",
                      seed=11)
    om_sim = ds_sim.gene_params.loc[ds_sim.gene_params["ritmico"], "omega"]
    om_apu = ds_apu.gene_params.loc[ds_apu.gene_params["ritmico"], "omega"]
    assert np.median(om_apu) < np.median(om_sim)
