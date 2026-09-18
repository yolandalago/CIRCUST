"""
Tests de las tablas de la seccion de agregacion.

Se fijan las decisiones de analisis, no el aspecto de los dibujos.
"""
import numpy as np
import pandas as pd
import pytest

from simulation.agregacion_figuras import (
    cargar_contraste, tabla_identidades, tabla_por_desacuerdo,
)

METODOS = (
    ("best_k", "r2_median", "auto", "circmean"),
    ("best_k", "medoid", "auto", "circmean"),
    ("aggregate", "r2_median", "tsp3", "circmean"),
    ("aggregate", "r2_median", "hods", "circmean"),
)


def _falso(reps: int = 8) -> pd.DataFrame:
    """
    Rejilla falsa donde HODs y TSP3 son identicos, medoid coincide con la
    referencia en la mitad de las replicas, y la agregacion mejora 0,05 h
    con desacuerdo alto y nada con desacuerdo bajo.
    """
    filas = []
    for snr in (1.0, 4.0):
        for rep in range(reps):
            desacuerdo = 0.02 if rep % 2 == 0 else 0.001
            for rm, sel, agg, esc in METODOS:
                oracle = 1.0
                if agg in ("tsp3", "hods") and desacuerdo > 0.01:
                    oracle -= 0.05
                if sel == "medoid" and rep % 2 == 1:
                    oracle += 0.01          # elige otra repeticion
                filas.append(dict(
                    bloque="rejilla", fmm_method="base", sampling="uniforme",
                    omega=0.5, snr=snr, n_samples=48, rep=rep, n_core=20,
                    n_genes=1000, prop_rhythmic=0.2, status="ok",
                    robust_method=rm, selection_criterion=sel,
                    aggregation=agg, consensus_scale=esc,
                    mae_h_oracle=oracle, tau_circ=0.9, medae_h=oracle + 0.5,
                    orientacion_ok=1.0, desacuerdo_K=desacuerdo, t_robust_s=1.0,
                ))
    return pd.DataFrame(filas)


def _csv(tmp_path, df):
    df.to_csv(tmp_path / "resultados.csv", index=False)
    return tmp_path


def test_identidades_distinguen_iguales_de_empates():
    t = tabla_identidades(_falso()).set_index(["metodo_a", "metodo_b"])
    assert t.loc[("TSP3", "HODs"), "fraccion_identicos"] == pytest.approx(1.0)
    assert t.loc[("best_k (R2 mediana)", "best_k (medoide)"),
                 "fraccion_identicos"] == pytest.approx(0.5)
    assert t.loc[("best_k (R2 mediana)", "HODs"),
                 "fraccion_identicos"] == pytest.approx(0.5)


def test_el_contraste_es_pareado_y_negativo_cuando_gana(tmp_path):
    par = cargar_contraste(_csv(tmp_path, _falso()))
    assert len(par) == 2 * 8
    ganan = par[par.desacuerdo_K_ref > 0.01]
    empatan = par[par.desacuerdo_K_ref <= 0.01]
    assert np.allclose(ganan["d_mae_h_oracle"], -0.05)
    assert np.allclose(empatan["d_mae_h_oracle"], 0.0)


def test_tabla_por_desacuerdo_separa_los_tramos(tmp_path):
    par = cargar_contraste(_csv(tmp_path, _falso()))
    t = tabla_por_desacuerdo(par, n_tramos=2)
    assert len(t) == 2
    bajo, alto = t.iloc[0], t.iloc[1]
    assert bajo["delta_mediana"] == pytest.approx(0.0)
    assert alto["delta_mediana"] == pytest.approx(-0.05)
    assert alto["tasa_victorias"] == pytest.approx(1.0)


def test_falla_si_no_hay_pares_del_metodo(tmp_path):
    df = _falso()
    df = df[df.aggregation != "hods"]
    with pytest.raises(ValueError, match="no hay pares"):
        cargar_contraste(_csv(tmp_path, df))
