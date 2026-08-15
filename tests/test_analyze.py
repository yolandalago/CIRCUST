"""Tests del analisis de la fase de calibracion.

Trabajan sobre tablas de resultados construidas a mano, de modo que se sabe de
antemano que celda tiene que ganar. Asi se comprueba el criterio y no el
pipeline.
"""
import json

import numpy as np
import pandas as pd
import pytest

from simulation.analyze import (
    CONFIG_CRITERIO, SD_RHO_MINIMA, VENTANA_MEDAE_H, VENTANA_RHO,
    criterio_fase0, elegir_referencia, escribir_referencia, resumen_tau,
)
from simulation.grid import REFERENCIA


def filas(celda, snr, n, *, rho, mad, sd, reps=5, fallos=0, n_top=100):
    """Replicas de una celda con la configuracion del criterio."""
    fmm, metodo = CONFIG_CRITERIO
    out = []
    for r in range(reps):
        if r < fallos:
            out.append({"celda": celda, "replica": r, "snr": snr,
                        "n_samples": n, "fmm_method": "", "metodo": "",
                        "rho": np.nan, "mad_h": np.nan, "sd_rho_reps": np.nan,
                        "n_top": np.nan, "error": "ValueError: colapso"})
            continue
        out.append({"celda": celda, "replica": r, "snr": snr, "n_samples": n,
                    "fmm_method": fmm, "metodo": metodo, "rho": rho,
                    "mad_h": mad, "sd_rho_reps": sd, "n_top": n_top,
                    "error": ""})
        # Una segunda configuracion, que el criterio tiene que ignorar.
        out.append({"celda": celda, "replica": r, "snr": snr, "n_samples": n,
                    "fmm_method": fmm, "metodo": "aggregate:tsp3",
                    "rho": 0.99, "mad_h": 0.01, "sd_rho_reps": 0.9,
                    "n_top": n_top, "error": ""})
    return out


@pytest.fixture
def fase0():
    d = []
    d += filas("saturada",  3.0, 48, rho=0.995, mad=0.25, sd=0.001)
    d += filas("buena",     1.0, 48, rho=0.806, mad=1.15, sd=0.038)
    d += filas("degenerada", 1.0, 12, rho=0.610, mad=1.67, sd=0.012)
    d += filas("dificil",   0.5, 48, rho=0.348, mad=2.73, sd=0.035)
    d += filas("rota",      0.5, 192, rho=0.80, mad=1.5, sd=0.05, fallos=3)
    return pd.DataFrame(d)


# ===========================================================================
# Criterio
# ===========================================================================

def test_el_criterio_resume_una_fila_por_celda(fase0):
    t = criterio_fase0(fase0)
    assert len(t) == 5
    assert set(t["celda"]) == {"saturada", "buena", "degenerada", "dificil",
                               "rota"}


def test_el_criterio_solo_mira_la_configuracion_de_referencia(fase0):
    # La segunda configuracion de cada celda tiene rho 0,99 y sd 0,9. Si el
    # criterio la mirase, la celda saturada pasaria.
    t = criterio_fase0(fase0).set_index("celda")
    assert abs(t.loc["buena", "rho"] - 0.806) < 1e-9
    assert abs(t.loc["saturada", "sd_rho_reps"] - 0.001) < 1e-9


def test_la_celda_saturada_no_pasa_por_falta_de_desacuerdo(fase0):
    t = criterio_fase0(fase0).set_index("celda")
    assert not t.loc["saturada", "cumple_sd"]
    assert not t.loc["saturada", "cumple"]


def test_la_celda_degenerada_no_pasa_aunque_su_rho_este_bien(fase0):
    # sd 0,012 es el regimen de babuino, donde ya se sabe que las K
    # repeticiones coinciden y agregar no puede aportar nada.
    t = criterio_fase0(fase0).set_index("celda")
    assert t.loc["degenerada", "sd_rho_reps"] < SD_RHO_MINIMA
    assert not t.loc["degenerada", "cumple"]


def test_una_celda_que_no_corre_siempre_queda_descartada(fase0):
    # No es un criterio mas, es una condicion previa: una celda que falla en
    # parte de sus replicas no puede ser la referencia.
    t = criterio_fase0(fase0).set_index("celda")
    assert t.loc["rota", "tasa_exito"] < 1.0
    assert not t.loc["rota", "cumple"]


def test_la_celda_buena_cumple_los_tres_criterios(fase0):
    t = criterio_fase0(fase0).set_index("celda")
    fila = t.loc["buena"]
    assert fila["cumple_rho"] and fila["cumple_medae"] and fila["cumple_sd"]
    assert fila["cumple"]


def test_las_ventanas_son_las_fijadas_de_antemano():
    assert VENTANA_RHO == (0.70, 0.90)
    assert VENTANA_MEDAE_H == (1.0, 2.0)
    assert SD_RHO_MINIMA == 0.03


# ===========================================================================
# Eleccion de la referencia
# ===========================================================================

def test_se_elige_la_celda_que_cumple(fase0):
    ref = elegir_referencia(fase0)
    assert ref["snr"] == 1.0
    assert ref["n_samples"] == 48


def test_la_referencia_conserva_los_demas_parametros(fase0):
    ref = elegir_referencia(fase0)
    for clave in ("n_genes", "prop_rhythmic", "sampling", "waveform_mix"):
        assert ref[clave] == REFERENCIA[clave]


def test_con_varias_candidatas_gana_la_de_rho_mas_cercano_al_objetivo():
    d = pd.DataFrame(filas("a", 1.0, 48, rho=0.72, mad=1.5, sd=0.05)
                     + filas("b", 2.0, 48, rho=0.79, mad=1.5, sd=0.05))
    assert elegir_referencia(d)["snr"] == 2.0


def test_si_ninguna_celda_cumple_se_dice_en_vez_de_inventar():
    d = pd.DataFrame(filas("saturada", 3.0, 48, rho=0.99, mad=0.2, sd=0.001))
    with pytest.raises(ValueError, match="ninguna celda"):
        elegir_referencia(d)


def test_escribir_referencia_deja_el_fichero_que_lee_el_lanzador(tmp_path):
    ruta = escribir_referencia({**REFERENCIA, "snr": 1.0}, tmp_path)
    assert ruta.name == "referencia.json"
    assert json.loads(ruta.read_text())["snr"] == 1.0


# ===========================================================================
# Barrido de tau
# ===========================================================================

def test_el_resumen_de_tau_ordena_por_umbral():
    d = pd.DataFrame([
        {"celda": "snr1", "replica": r, "metodo": f"best_k:n_well_fitted@{t}",
         "rho": 0.8, "mad_h": 1.2, "rep_elegida": k, "score_seleccion": 10,
         "error": ""}
        for t in (0.6, 0.7, 0.8) for r, k in enumerate([0, 1, 0, 2, 1])
    ])
    t = resumen_tau(d)
    assert list(t["tau"]) == [0.6, 0.7, 0.8]


def test_el_resumen_de_tau_avisa_cuando_el_criterio_degenera():
    # Con tau alto ninguna repeticion supera el umbral, todas empatan a cero y
    # el criterio desempata por mediana de R2, o sea colapsa en r2_median.
    d = pd.DataFrame(
        [{"celda": "c", "replica": r, "metodo": "best_k:n_well_fitted@0.6",
          "rho": 0.8, "mad_h": 1.2, "rep_elegida": r % 3,
          "score_seleccion": 12, "error": ""} for r in range(5)]
        + [{"celda": "c", "replica": r, "metodo": "best_k:n_well_fitted@0.9",
            "rho": 0.8, "mad_h": 1.2, "rep_elegida": 0,
            "score_seleccion": 0, "error": ""} for r in range(5)]
    )
    t = resumen_tau(d).set_index("tau")
    assert not t.loc[0.6, "degenerado"]
    assert t.loc[0.9, "degenerado"]
