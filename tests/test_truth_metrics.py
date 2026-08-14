"""Tests de las metricas contra el ground truth.

Se construyen atlas artificiales a partir de la propia verdad, aplicandoles
transformaciones conocidas, de modo que se sabe de antemano que tiene que
salir. Asi los tests comprueban la metrica y no el algoritmo.
"""
import numpy as np
import pandas as pd
import pytest

from simulation.truth_metrics import (
    deteccion_top, error_omega, error_omega_por_clase, error_picos,
    evaluar_contra_verdad,
)

DOS_PI = 2.0 * np.pi


@pytest.fixture
def params():
    """Veinte genes ritmicos con fases repartidas y omegas de tres clases.

    El primero es ARNTL con pico en la hora cero, como en el generador, porque
    es el gen ancla y sin el no puede medirse el error de sincronizacion.
    """
    n = 20
    rng = np.random.default_rng(0)
    t_u = np.linspace(0.0, DOS_PI, n, endpoint=False)
    clases = np.array(["sinusoidal"] * 6 + ["moderado"] * 8 + ["apuntado"] * 6)
    omega = np.where(clases == "sinusoidal", 0.85,
             np.where(clases == "moderado", 0.65, 0.45))
    return pd.DataFrame({
        "gene":    ["ARNTL"] + [f"RIT{i:03d}" for i in range(1, n)],
        "ritmico": True,
        "clase":   clases,
        "omega":   omega + rng.normal(0.0, 1e-9, n),
        "t_U_rad": t_u,
        "t_U_h":   t_u * 24.0 / DOS_PI,
    })


def atlas_anclado(params, *, error_sinc_h=0.0, reflejar=False,
                  sesgo_omega=0.0):
    """Atlas artificial tal y como lo produciria el pipeline.

    Reproduce el convenio de CIRCUST: el pico del gen ancla se situa en pi, o
    sea en la hora 12 de la escala reconstruida, con independencia de donde
    pique de verdad. ``error_sinc_h`` desplaza el origen respecto de ese
    convenio, que es lo que la metrica tiene que recuperar.
    """
    t_a = float(params.loc[params["gene"] == "ARNTL", "t_U_h"].iloc[0])
    signo = -1.0 if reflejar else 1.0
    est_h = (signo * params["t_U_h"].to_numpy() - signo * t_a + 12.0
             + error_sinc_h) % 24.0
    return pd.DataFrame({
        "gene_name":      params["gene"],
        "peak_phase_rad": est_h * DOS_PI / 24.0,
        "omega":          params["omega"].to_numpy() + sesgo_omega,
    })


# ===========================================================================
# Fases de pico y sincronizacion
# ===========================================================================

def test_un_atlas_perfecto_no_tiene_error_ni_de_picos_ni_de_sincronizacion(
        params):
    r = error_picos(atlas_anclado(params), params)
    assert r["verdad_picos_medae_h"] < 1e-6
    assert r["verdad_sinc_error_h"] < 0.1
    assert r["verdad_sinc_reflejado"] == 0
    assert r["verdad_picos_n"] == 20


def test_el_anclaje_por_convenio_no_se_confunde_con_un_error(params):
    # ARNTL pica en la hora 0 pero el pipeline lo situa en la 12, asi que el
    # desplazamiento crudo de la alineacion vale -12 aunque no haya ningun
    # error. Es justo la confusion que la metrica tiene que evitar.
    r = error_picos(atlas_anclado(params), params)
    assert abs(abs(r["verdad_sinc_delta_h"]) - 12.0) < 0.1
    assert r["verdad_sinc_error_h"] < 0.1


def test_un_error_de_sincronizacion_se_absorbe_y_queda_medido(params):
    r = error_picos(atlas_anclado(params, error_sinc_h=5.0), params)
    # Girar el origen es una indeterminacion del metodo, no un error de
    # estimacion de las fases, asi que el error de picos sigue siendo nulo...
    assert r["verdad_picos_medae_h"] < 0.1
    # ...pero queda registrado como lo que es, el error de la etapa de
    # sincronizacion, que ninguna validacion previa de CIRCUST mide.
    assert abs(r["verdad_sinc_error_h"] - 5.0) < 0.1


def test_el_error_de_sincronizacion_no_depende_del_sentido_del_giro(params):
    izq = error_picos(atlas_anclado(params, error_sinc_h=3.0), params)
    der = error_picos(atlas_anclado(params, error_sinc_h=-3.0), params)
    assert abs(izq["verdad_sinc_error_h"] - 3.0) < 0.1
    assert abs(der["verdad_sinc_error_h"] - 3.0) < 0.1


def test_un_ciclo_recorrido_al_reves_se_detecta(params):
    r = error_picos(atlas_anclado(params, reflejar=True), params)
    assert r["verdad_picos_medae_h"] < 0.1
    assert r["verdad_sinc_reflejado"] == 1
    # El anclaje siguio siendo correcto, lo que fallo es la direccion.
    assert r["verdad_sinc_error_h"] < 0.1


def test_sin_gen_ancla_el_error_de_sincronizacion_es_nan(params):
    sin_ancla = params[params["gene"] != "ARNTL"]
    r = error_picos(atlas_anclado(params), sin_ancla)
    assert np.isnan(r["verdad_sinc_error_h"])
    assert r["verdad_picos_medae_h"] < 1e-6


def test_el_ruido_en_las_fases_se_traduce_en_error(params):
    rng = np.random.default_rng(1)
    atlas = atlas_anclado(params)
    atlas["peak_phase_rad"] = (atlas["peak_phase_rad"]
                               + rng.normal(0.0, 0.3, len(atlas))) % DOS_PI
    r = error_picos(atlas, params)
    assert 0.3 < r["verdad_picos_medae_h"] < 3.0
    assert r["verdad_picos_p75_h"] >= r["verdad_picos_medae_h"]


def test_sin_genes_en_comun_devuelve_nan_y_no_revienta(params):
    atlas = atlas_anclado(params)
    atlas["gene_name"] = [f"OTRO{i}" for i in range(len(atlas))]
    r = error_picos(atlas, params)
    assert np.isnan(r["verdad_picos_medae_h"])
    assert r["verdad_picos_n"] == 0


# ===========================================================================
# Omega
# ===========================================================================

def test_un_atlas_perfecto_no_tiene_error_de_omega(params):
    r = error_omega(atlas_anclado(params), params)
    assert r["verdad_omega_medae"] < 1e-6
    assert r["verdad_omega_n"] == 20


def test_el_sesgo_de_omega_conserva_el_signo(params):
    # Signo negativo significa subestimar omega, o sea estimar formas mas
    # apuntadas que las reales. Es exactamente la patologia medida en babuino
    # con el FMM base y pocas muestras, asi que el signo tiene que verse.
    r = error_omega(atlas_anclado(params, sesgo_omega=-0.2), params)
    assert abs(r["verdad_omega_sesgo"] + 0.2) < 1e-6
    assert abs(r["verdad_omega_medae"] - 0.2) < 1e-6


def test_omega_devuelve_nan_si_el_atlas_no_la_trae(params):
    atlas = atlas_anclado(params).drop(columns=["omega"])
    r = error_omega(atlas, params)
    assert np.isnan(r["verdad_omega_medae"])
    assert r["verdad_omega_n"] == 0


def test_el_error_de_omega_se_puede_desglosar_por_clase(params):
    tabla = error_omega_por_clase(atlas_anclado(params, sesgo_omega=-0.2), params)
    assert set(tabla["clase"]) == {"sinusoidal", "moderado", "apuntado"}
    assert np.allclose(tabla["verdad_omega_sesgo"], -0.2)
    assert tabla["verdad_omega_n"].sum() == 20


# ===========================================================================
# Deteccion del conjunto TOP
# ===========================================================================

def test_un_top_todo_ritmico_tiene_precision_uno(params):
    r = deteccion_top(list(params["gene"][:10]), params)
    assert r["verdad_top_n"] == 10
    assert r["verdad_top_precision"] == 1.0
    assert r["verdad_top_n_ritmicos"] == 10


def test_el_top_contaminado_con_genes_de_fondo_baja_la_precision(params):
    genes = list(params["gene"][:8]) + ["BG00001", "BG00002"]
    r = deteccion_top(genes, params)
    assert r["verdad_top_n"] == 10
    assert abs(r["verdad_top_precision"] - 0.8) < 1e-9


def test_el_hueco_maximo_detecta_un_top_mal_repartido(params):
    # Diez genes consecutivos cubren menos de medio ciclo, luego dejan un
    # hueco grande. Es lo que el cupo por sector del conjunto TOP evita.
    apinado = deteccion_top(list(params["gene"][:10]), params)
    repartido = deteccion_top(list(params["gene"][::2]), params)
    assert apinado["verdad_top_hueco_max_h"] > repartido["verdad_top_hueco_max_h"]
    assert repartido["verdad_top_hueco_max_h"] < 4.0


def test_un_top_vacio_no_revienta(params):
    r = deteccion_top([], params)
    assert r["verdad_top_n"] == 0
    assert np.isnan(r["verdad_top_precision"])


# ===========================================================================
# Agregado
# ===========================================================================

def test_evaluar_contra_verdad_junta_las_tres_familias(params):
    fila = evaluar_contra_verdad(atlas_anclado(params, error_sinc_h=3.0),
                                 params, genes_top=list(params["gene"]))
    for clave in ("verdad_picos_medae_h", "verdad_sinc_error_h",
                  "verdad_omega_medae", "verdad_top_precision"):
        assert clave in fila, clave
    assert fila["verdad_picos_medae_h"] < 0.1


def test_evaluar_contra_verdad_funciona_sin_lista_de_top(params):
    fila = evaluar_contra_verdad(atlas_anclado(params), params)
    assert np.isnan(fila["verdad_top_precision"])
    assert fila["verdad_picos_medae_h"] < 1e-6
