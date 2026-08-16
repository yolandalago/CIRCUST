"""Tests del analisis de los bloques y de las figuras.

Las funciones de tabla se prueban sobre datos construidos a mano, donde se
sabe el resultado. De las figuras se comprueba que se generan y que producen
los dos formatos, no su apariencia, que se revisa mirandolas.
"""
import numpy as np
import pandas as pd
import pytest

from simulation.report import (
    con_estrato, contraste, emparejar, figura_bloque_a, figura_frontera,
    figura_ganancia_r, figura_sesgo_omega, frontera, tabla_agregacion,
    tabla_bloque_a, tabla_bloque_r, tabla_seleccion, validas,
)


def filas_r(celda, replica, *, n=48, mezcla="gtex", muestreo="huecos",
            om_base=0.20, om_stab=0.15, pic_base=1.0, pic_stab=0.8):
    comun = {"celda": celda, "replica": replica, "n_samples": n,
             "waveform_mix": mezcla, "sampling": muestreo, "error": ""}
    return [
        {**comun, "fmm_method": "base", "verdad_omega_medae": om_base,
         "verdad_omega_sesgo": -0.1, "verdad_picos_medae_h": pic_base},
        {**comun, "fmm_method": "stabilized", "verdad_omega_medae": om_stab,
         "verdad_omega_sesgo": 0.02, "verdad_picos_medae_h": pic_stab},
    ]


@pytest.fixture
def bloque_r():
    d = []
    for rep in range(12):
        d += filas_r("n12", rep, n=12, om_base=0.20, om_stab=0.15)
        d += filas_r("n192", rep, n=192, om_base=0.10, om_stab=0.17)
    return pd.DataFrame(d)


# ===========================================================================
# Emparejado y contraste
# ===========================================================================

def test_emparejar_compara_dentro_de_la_misma_replica(bloque_r):
    j = emparejar(bloque_r, "fmm_method", "base", "stabilized",
                  ["verdad_omega_medae"])
    assert len(j) == 24
    n12 = j[j["n_samples"] == 12]["dif_verdad_omega_medae"]
    assert np.allclose(n12, -0.05)


def test_la_diferencia_negativa_significa_que_gana_el_segundo(bloque_r):
    t = tabla_bloque_r(bloque_r).set_index(["eje", "nivel"])
    assert t.loc[("n_samples", 12), "gana"] == "estabilizado"
    assert t.loc[("n_samples", 192), "gana"] == "base"


def test_el_contraste_no_se_pronuncia_con_pocos_pares():
    assert np.isnan(contraste(pd.Series([0.1, -0.2, 0.3]))["p"])


def test_el_contraste_ignora_las_diferencias_nulas():
    c = contraste(pd.Series([0.0] * 20 + [0.1] * 8))
    assert c["n"] == 8


def test_una_diferencia_sistematica_sale_significativa():
    c = contraste(pd.Series(np.linspace(0.05, 0.30, 20)))
    assert c["p"] < 0.01


# ===========================================================================
# Estratos de desacuerdo
# ===========================================================================

def test_los_estratos_cortan_donde_esta_documentado():
    d = pd.DataFrame({"celda": "c", "replica": range(4),
                      "sd_rho_reps": [0.005, 0.02, 0.05, 0.20],
                      "metodo": "best_k:r2_median", "rho": 0.8, "error": ""})
    assert list(con_estrato(d)["estrato"]) == [
        "< 0,01", "0,01 a 0,03", "> 0,03", "> 0,03"]


def test_el_estrato_de_babuino_queda_por_debajo_del_umbral():
    # 0,012 es el valor medido en los 47 tejidos de babuino, donde ya se sabe
    # que agregar no aporta. Tiene que caer en el estrato bajo.
    d = pd.DataFrame({"celda": ["c"], "replica": [0], "sd_rho_reps": [0.012],
                      "metodo": ["best_k:r2_median"], "rho": [0.8],
                      "error": [""]})
    assert con_estrato(d)["estrato"].iloc[0] == "0,01 a 0,03"


# ===========================================================================
# Frontera de aplicabilidad
# ===========================================================================

def test_la_frontera_cuenta_sobre_todas_las_replicas_no_solo_las_validas():
    # Una celda que falla en tres de cinco replicas tiene que dar 60, no
    # desaparecer del recuento por no tener filas validas.
    d = pd.DataFrame([
        {"celda": "c", "replica": r, "n_samples": 192, "snr": 1.0,
         "error": "ValueError: colapso" if r < 3 else "",
         "metodo": "best_k:r2_median", "rho": np.nan if r < 3 else 0.8}
        for r in range(5)])
    assert frontera(d).loc[192, 1.0] == 60.0


def test_una_celda_que_falla_siempre_sale_al_cien_por_cien():
    d = pd.DataFrame([
        {"celda": "c", "replica": r, "n_samples": 192, "snr": 0.5,
         "error": "ValueError", "metodo": "", "rho": np.nan}
        for r in range(5)])
    assert frontera(d).loc[192, 0.5] == 100.0


def test_validas_descarta_las_filas_con_error():
    d = pd.DataFrame({"celda": "c", "replica": [0, 1],
                      "error": ["", "ValueError: x"]})
    assert len(validas(d)) == 1


# ===========================================================================
# Bloque A
# ===========================================================================

@pytest.fixture
def bloque_a():
    d = []
    for rep in range(30):
        sd = 0.005 if rep < 15 else 0.05
        for metodo, rho, elegida in [("best_k:r2_median", 0.80, 1),
                                     ("best_k:medoid", 0.79, 2),
                                     ("aggregate:tsp3", 0.78, np.nan)]:
            d.append({"celda": "c", "replica": rep, "metodo": metodo,
                      "rho": rho, "sd_rho_reps": sd, "rep_elegida": elegida,
                      "n_samples": 48, "snr": 1.0, "error": ""})
    return pd.DataFrame(d)


def test_la_tabla_del_bloque_a_respeta_el_orden_canonico(bloque_a):
    t = tabla_bloque_a(bloque_a)
    assert list(t.index) == ["best_k:r2_median", "best_k:medoid",
                             "aggregate:tsp3"]


def test_la_diferencia_frente_a_la_referencia_lleva_su_signo(bloque_a):
    t = tabla_agregacion(bloque_a)
    peor = t[t["metodo"] == "aggregate:tsp3"]
    assert (peor["mediana_dif"] < 0).all()


def test_la_tabla_de_seleccion_mide_cuanto_coinciden_los_criterios(bloque_a):
    t = tabla_seleccion(bloque_a)
    fila = t.iloc[0]
    assert fila["coinciden_pct"] == 0.0
    assert fila["n"] == 30


# ===========================================================================
# Figuras
# ===========================================================================

def test_las_figuras_del_bloque_r_se_generan(bloque_r, tmp_path):
    for f in (figura_sesgo_omega, figura_ganancia_r):
        salida = f(bloque_r, tmp_path / f"{f.__name__}.pdf")
        assert salida.exists()
        assert salida.with_suffix(".png").exists()


def test_las_figuras_del_bloque_a_se_generan(bloque_a, tmp_path):
    for f in (figura_bloque_a, figura_frontera):
        salida = f(bloque_a, tmp_path / f"{f.__name__}.pdf")
        assert salida.exists()
        assert salida.with_suffix(".png").exists()
