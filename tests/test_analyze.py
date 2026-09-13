"""Tests de simulation/analyze.py."""
import numpy as np
import pandas as pd
import pytest

from simulation.analyze import (
    MEDIDAS, tabla_agregacion, tabla_gcore, tabla_principal,
    tabla_regularizacion, tabla_sin_ruido,
)


def _fila(**kw):
    base = dict(
        bloque="principal", fmm_method="base", sampling="uniforme", omega=0.5,
        snr=2.0, n_samples=48, rep=0, n_core=20, prop_rhythmic=0.20,
        robust_method="best_k", selection_criterion="r2_median",
        aggregation="auto", consensus_scale="circmean", status="ok",
        medae_h=1.0, mae_h=1.2, rho_circ=0.9, tau=0.8, tau_circ=0.9,
        orientacion_ok=1.0, p1=0.4, p2=0.7, p4=0.9, mae_h_oracle=0.5,
        r2_anchor=0.8, pares_ok=1.0, evidencia=1.7, desacuerdo_K=0.01,
        n_top=100, t_total_s=5.0,
    )
    base.update(kw)
    return base


def _falso():
    filas = []
    for bloque, fmm in (("principal", "base"), ("regularizacion", "stabilized")):
        for mu in ("uniforme", "hueco"):
            for w in (0.1, 0.9):
                for snr in (1.0, 3.0):
                    for rep in range(3):
                        filas.append(_fila(
                            bloque=bloque, fmm_method=fmm, sampling=mu,
                            omega=w, snr=snr, rep=rep,
                            # stabilized es 0,5 h mejor en el falso
                            medae_h=1.5 if fmm == "base" else 1.0))
    for prop in (0.01, 0.20):
        for metodo in ("best_k", "aggregate"):
            for rep in range(3):
                filas.append(_fila(bloque="agregacion", prop_rhythmic=prop,
                                   robust_method=metodo, aggregation="tsp3",
                                   rep=rep))
    for gcore in (10, 20, 30):
        for rep in range(3):
            filas.append(_fila(bloque="gcore", n_core=gcore, rep=rep))
    for w in (0.1, 0.9):
        for rep in range(3):
            filas.append(_fila(bloque="sin_ruido", omega=w, snr=1e6,
                               tau_circ=0.2 if w == 0.1 else 1.0,
                               orientacion_ok=0.0 if w == 0.1 else 1.0))
    return pd.DataFrame(filas)


def test_tabla_principal_resume_por_los_factores():
    t = tabla_principal(_falso())
    assert list(t.index.names) == ["sampling", "omega", "snr"]
    # una fila por combinacion de los tres factores del falso
    assert len(t) == 2 * 2 * 2
    assert {"medae_h", "tau_circ", "orientacion_ok"} <= set(
        t.columns.get_level_values(0))


def test_rho_se_agrega_con_mediana_y_no_con_media():
    """
    Decision tomada tras medir su fragilidad: rho cambia de signo en torno al
    1-2 % de las replicas aunque la reconstruccion sea perfecta, y una sola
    replica mala arrastra la media de la celda. La mediana aguanta.
    """
    df = _falso()
    p = df[df.bloque == "principal"].copy()
    # una replica de cada celda con rho catastrofico
    idx = p.groupby(["sampling", "omega", "snr"]).head(1).index
    p.loc[idx, "rho_circ"] = -0.9
    t = tabla_principal(pd.concat([df[df.bloque != "principal"], p]))
    # la mediana de (0.9, 0.9, -0.9) es 0.9; la media seria 0.3
    assert (t[("rho_circ", "mediana")] > 0.8).all()


def test_la_orientacion_se_agrega_como_TASA_no_como_mediana():
    """Es una proporcion de aciertos, no una medida continua."""
    df = _falso()
    p = df[df.bloque == "principal"].copy()
    idx = p.groupby(["sampling", "omega", "snr"]).head(1).index
    p.loc[idx, "orientacion_ok"] = 0.0        # 1 de cada 3 replicas falla
    t = tabla_principal(pd.concat([df[df.bloque != "principal"], p]))
    assert np.allclose(t[("orientacion_ok", "tasa")].values, 2 / 3)


def test_tabla_regularizacion_es_pareada():
    """La diferencia se calcula por replica, no entre medias de celda."""
    t = tabla_regularizacion(_falso())
    assert "delta_medae_h" in t.columns.get_level_values(0)
    # en el falso, stabilized es 0,5 h mejor en todas las celdas
    assert np.allclose(t[("delta_medae_h", "mediana")].values, -0.5)


def test_tabla_agregacion_separa_por_prop_rhythmic():
    t = tabla_agregacion(_falso()).reset_index()
    assert set(t["prop_rhythmic"]) == {0.01, 0.20}
    assert set(t["robust_method"]) == {"best_k", "aggregate"}


def test_tabla_gcore_tiene_los_tres_niveles():
    t = tabla_gcore(_falso()).reset_index()
    assert set(t["n_core"]) == {10, 20, 30}


def test_tabla_sin_ruido_separa_orden_de_orientacion():
    """
    Es la verificacion que se le ensena a la tutora: con datos perfectos y
    omega moderada el orden tiene que salir exacto, y con omega 0,1 se rompe.
    """
    t = tabla_sin_ruido(_falso()).reset_index()
    assert set(t["omega"]) == {0.1, 0.9}
    buena = t[t.omega == 0.9].iloc[0]
    mala  = t[t.omega == 0.1].iloc[0]
    assert buena[("tau_circ", "mediana")] == pytest.approx(1.0)
    assert buena[("orientacion_ok", "tasa")] == pytest.approx(1.0)
    assert mala[("orientacion_ok", "tasa")] == pytest.approx(0.0)


def test_las_filas_con_fallo_no_entran_en_las_medianas_pero_si_se_cuentan():
    df = _falso()
    malas = df[df.bloque == "principal"].head(4).copy()
    malas["status"] = "fallo"
    for m in MEDIDAS:
        malas[m] = np.nan
    t = tabla_principal(pd.concat([df, malas], ignore_index=True))
    assert not t[("medae_h", "mediana")].isna().any()
    assert "tasa_fallo" in t.columns.get_level_values(0)
    assert (t[("tasa_fallo", "")] > 0).any()


def test_no_se_escriben_figuras_vacias(tmp_path):
    """
    Antes se escribian PNG vacios cuando el bloque no se habia ejecutado, y un
    fichero vacio parece una figura valida.
    """
    from simulation.analyze import figura_snr
    solo_sin_ruido = _falso().query("bloque == 'sin_ruido'")
    with pytest.raises(ValueError, match="no hay datos"):
        figura_snr(solo_sin_ruido, tmp_path / "vacia.png")
    assert not (tmp_path / "vacia.png").exists()
