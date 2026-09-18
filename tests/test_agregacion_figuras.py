"""
Tests de las tablas de la seccion de agregacion.

Se fijan las decisiones de analisis, no el aspecto de los dibujos.
"""
import numpy as np
import pandas as pd
import pytest

from simulation.agregacion_figuras import (
    cargar_contraste, tabla_identidades, tabla_por_factor,
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
    con A/sigma 1 y nada con A/sigma 4.
    """
    filas = []
    for snr in (1.0, 4.0):
        for rep in range(reps):
            for rm, sel, agg, esc in METODOS:
                error = 1.0
                if agg in ("tsp3", "hods") and snr == 1.0:
                    error -= 0.05
                if sel == "medoid" and rep % 2 == 1:
                    error += 0.01          # elige otra repeticion
                filas.append(dict(
                    bloque="rejilla", fmm_method="base", sampling="uniforme",
                    omega=0.5, snr=snr, n_samples=48, rep=rep, n_core=20,
                    n_genes=1000, prop_rhythmic=0.2, status="ok",
                    robust_method=rm, selection_criterion=sel,
                    aggregation=agg, consensus_scale=esc,
                    error_orden_h=error, tau_circular=0.9,
                    error_fase_h=error + 0.5, acierto_sentido=1.0,
                    t_robust_s=1.0,
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
    ganan = par[par.snr == 1.0]
    empatan = par[par.snr == 4.0]
    assert np.allclose(ganan["d_error_orden_h"], -0.05)
    assert np.allclose(empatan["d_error_orden_h"], 0.0)


def test_tabla_por_factor_recupera_el_efecto_por_nivel(tmp_path):
    par = cargar_contraste(_csv(tmp_path, _falso()))
    t = tabla_por_factor(par)
    assert t.loc[("snr", "1.0"), "ganancia_mediana"] == pytest.approx(-0.05)
    assert t.loc[("snr", "1.0"), "tasa_victorias"] == pytest.approx(1.0)
    assert t.loc[("snr", "4.0"), "ganancia_mediana"] == pytest.approx(0.0)
    assert t.loc[("snr", "4.0"), "tasa_victorias"] == pytest.approx(0.0)
    # los factores sin variacion en el falso salen con un solo nivel
    assert len(t.loc["n_samples"]) == 1


def test_falla_si_no_hay_pares_del_metodo(tmp_path):
    df = _falso()
    df = df[df.aggregation != "hods"]
    with pytest.raises(ValueError, match="no hay pares"):
        cargar_contraste(_csv(tmp_path, df))
