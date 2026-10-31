"""Tests de simulation/analyze.py."""
import numpy as np
import pandas as pd
import pytest

from simulation.analyze import (
    FACTORES, MEDIDAS, contraste_fmm, contraste_metodos, efectos_principales,
    interaccion, modelo_lineal, tabla_atlas, tabla_por_factor,
    tabla_reconstruccion, tabla_sin_ruido, tabla_tiempos_python,
)

METODOS = (
    ("best_k", "r2_median", "auto", "circmean"),
    ("best_k", "medoid", "auto", "circmean"),
    ("aggregate", "r2_median", "tsp3", "circmean"),
    ("aggregate", "r2_median", "hods", "circmean"),
)


def _fila(**kw):
    base = dict(
        bloque="rejilla", fmm_method="base", sampling="uniforme", omega=0.5,
        snr=2.0, n_samples=48, rep=0, n_core=20, n_genes=1000,
        prop_rhythmic=0.20,
        robust_method="best_k", selection_criterion="r2_median",
        aggregation="auto", consensus_scale="circmean", status="ok",
        error_fase_h=1.0, error_fase_medio_h=1.2, rho_circ=0.9, rho_fl=0.9, tau=0.8, tau_circular=0.9,
        acierto_sentido=1.0, p1=0.4, p2=0.7, p4=0.9, error_orden_h=0.5,
        r2_anchor=0.8, pares_ok=1.0, evidencia=1.7, desacuerdo_K=0.01,
        n_top=100, t_total_s=5.0, t_robust_s=2.0,
        atlas_error_onda=0.3, atlas_variabilidad_pico_h=1.0,
        atlas_sesgo_omega=-0.2, atlas_vpp=0.95,
    )
    base.update(kw)
    return base


def _falso(reps: int = 3) -> pd.DataFrame:
    """
    Rejilla falsa: 2 muestreos x 2 omegas x 2 SNR x 2 FMM x 4 metodos x reps.

    Efectos plantados, para que los tests sepan que esperar:
    - CircaFMM tiene 0,1 menos de RMSE de curva, y la ventaja es el
      doble con snr 1 que con snr 3
    - la agregacion (tsp3 y hods) tiene 0,1 h menos de error del orden relativo que la
      referencia, medoid 0,05 h menos
    """
    filas = []
    for mu in ("uniforme", "hueco"):
        for w in (0.1, 0.9):
            for snr in (1.0, 3.0):
                for rep in range(reps):
                    for fmm in ("base", "circafmm"):
                        for rm, sel, agg, esc in METODOS:
                            ventaja_fmm = (0.2 if snr == 1.0 else 0.1) if fmm == "circafmm" else 0.0
                            ventaja_met = {"tsp3": 0.1, "hods": 0.1}.get(agg, 0.0)
                            if sel == "medoid":
                                ventaja_met = 0.05
                            filas.append(_fila(
                                fmm_method=fmm, sampling=mu, omega=w, snr=snr,
                                rep=rep, robust_method=rm, selection_criterion=sel,
                                aggregation=agg, consensus_scale=esc,
                                atlas_error_onda=0.4 - ventaja_fmm,
                                error_orden_h=0.6 - ventaja_met,
                                error_fase_h=1.5 - ventaja_met,
                            ))
    for w in (0.1, 0.9):
        for rep in range(reps):
            filas.append(_fila(bloque="sin_ruido", omega=w, snr=1e6,
                               tau_circular=0.2 if w == 0.1 else 1.0,
                               acierto_sentido=0.0 if w == 0.1 else 1.0))
    return pd.DataFrame(filas)


# ── reconstruccion ─────────────────────────────────────────────────────────

def test_tabla_reconstruccion_usa_solo_la_referencia():
    t = tabla_reconstruccion(_falso())
    assert list(t.index.names) == ["sampling", "omega", "snr"]
    assert len(t) == 2 * 2 * 2
    # 3 replicas por celda: solo base y best_k:r2_median entran
    assert (t[("n", "")] == 3).all()
    assert {"error_fase_h", "concordancia", "acierto_sentido"} <= set(
        t.columns.get_level_values(0))


def test_tabla_por_factor_tiene_una_fila_por_nivel():
    t = tabla_por_factor(_falso()).reset_index()
    assert set(t["factor"]) == set(FACTORES)
    assert set(t[t.factor == "sampling"]["nivel"]) == {"uniforme", "hueco"}


def test_rho_se_agrega_con_mediana_y_no_con_media():
    """
    rho cambia de signo en torno al 1-2 % de las replicas aunque la
    reconstruccion sea perfecta. Una sola replica mala arrastra la media de
    la celda, la mediana aguanta.
    """
    df = _falso()
    ref = (df.bloque == "rejilla") & (df.fmm_method == "base") \
        & (df.selection_criterion == "r2_median") & (df.robust_method == "best_k")
    idx = df[ref].groupby(["sampling", "omega", "snr"]).head(1).index
    df.loc[idx, "rho_circ"] = -0.9
    t = tabla_reconstruccion(df)
    assert (t[("rho_circ", "mediana")] > 0.8).all()


def test_la_orientacion_se_agrega_como_TASA_no_como_mediana():
    df = _falso()
    ref = (df.bloque == "rejilla") & (df.fmm_method == "base") \
        & (df.selection_criterion == "r2_median") & (df.robust_method == "best_k")
    idx = df[ref].groupby(["sampling", "omega", "snr"]).head(1).index
    df.loc[idx, "acierto_sentido"] = 0.0
    t = tabla_reconstruccion(df)
    assert np.allclose(t[("acierto_sentido", "tasa")].values, 2 / 3)


def test_las_filas_con_fallo_no_entran_en_las_medianas_pero_si_se_cuentan():
    df = _falso()
    ref = (df.bloque == "rejilla") & (df.fmm_method == "base") \
        & (df.selection_criterion == "r2_median") & (df.robust_method == "best_k")
    malas = df[ref].head(4).copy()
    malas["rep"] += 100
    malas["status"] = "fallo"
    for m in MEDIDAS:
        if m in malas.columns:
            malas[m] = np.nan
    t = tabla_reconstruccion(pd.concat([df, malas], ignore_index=True))
    assert not t[("error_fase_h", "mediana")].isna().any()
    assert (t[("tasa_fallo", "")] > 0).any()


def test_la_calidad_del_orden_se_resume_en_valor_absoluto():
    """
    tau_circular cambia de signo cuando el ciclo se recorre al reves. Cuatro
    replicas +0,8 +0,8 -0,8 -0,8 resumidas con signo dan 0,0, que se lee
    como "no se recupero nada" cuando la calidad real es 0,8 en las cuatro.
    """
    filas = [_fila(rep=i, tau_circular=s * 0.8, acierto_sentido=float(s > 0))
             for i, s in enumerate([1, 1, -1, -1])]
    t = tabla_reconstruccion(pd.DataFrame(filas))
    assert t[("concordancia", "mediana")].iloc[0] == pytest.approx(0.8)
    assert t[("acierto_sentido", "tasa")].iloc[0] == pytest.approx(0.5)
    assert ("tau_circular", "mediana") not in t.columns


def test_tabla_sin_ruido_separa_orden_de_orientacion():
    t = tabla_sin_ruido(_falso()).reset_index()
    assert set(t["omega"]) == {0.1, 0.9}
    buena = t[t.omega == 0.9].iloc[0]
    mala  = t[t.omega == 0.1].iloc[0]
    assert buena[("concordancia", "mediana")] == pytest.approx(1.0)
    assert buena[("acierto_sentido", "tasa")] == pytest.approx(1.0)
    assert mala[("acierto_sentido", "tasa")] == pytest.approx(0.0)


# ── contraste FMM ──────────────────────────────────────────────────────────

def test_contraste_fmm_es_pareado_y_con_el_metodo_fijo():
    par = contraste_fmm(_falso())
    # 8 celdas x 3 replicas, un par por replica, solo la referencia
    assert len(par) == 8 * 3
    esperado = np.where(par["snr"] == 1.0, -0.2, -0.1)
    assert np.allclose(par["delta_atlas_error_onda"], esperado)


def test_contraste_fmm_descarta_las_replicas_con_orientacion_invertida():
    df = _falso()
    malas = (df.bloque == "rejilla") & (df.fmm_method == "circafmm") & (df.rep == 0)
    df.loc[malas, "acierto_sentido"] = 0.0
    par = contraste_fmm(df)
    assert 0 not in set(par["rep"])
    assert len(par) == 8 * 2


def test_contraste_fmm_falla_si_las_replicas_no_casan():
    df = _falso()
    df.loc[df.fmm_method == "circafmm", "rep"] += 50
    with pytest.raises(ValueError, match="no cuadra"):
        contraste_fmm(df)


def test_efectos_principales_detectan_la_interaccion_plantada():
    par = contraste_fmm(_falso())
    t = efectos_principales(par, "atlas_error_onda").reset_index()
    snr = t[t.factor == "snr"].set_index("nivel")
    assert snr.loc["1.0", "mediana"] == pytest.approx(-0.2)
    assert snr.loc["3.0", "mediana"] == pytest.approx(-0.1)
    assert (t["tasa_victorias"] == 1.0).all()
    assert set(t["factor"]) == set(FACTORES)


def test_interaccion_devuelve_el_mapa_de_dos_factores():
    par = contraste_fmm(_falso())
    m = interaccion(par, "atlas_error_onda", "omega", "snr")
    assert m.shape == (2, 2)
    assert m.loc[0.1, 1.0] == pytest.approx(-0.2)
    v = interaccion(par, "atlas_error_onda", "omega", "snr", "victorias")
    assert (v.values == 1.0).all()


def test_modelo_lineal_recupera_el_efecto_plantado():
    par = contraste_fmm(_falso())
    t = modelo_lineal(par, "atlas_error_onda", factores=["snr", "omega", "sampling"])
    # pasar de snr 1 a snr 3 reduce la ventaja en 0,1
    assert t.loc["C(snr)[T.3.0]", "coef"] == pytest.approx(0.1)
    assert abs(t.loc["C(omega)[T.0.9]", "coef"]) < 1e-9
    assert t.attrs["r2"] == pytest.approx(1.0)


def test_tabla_atlas_separa_los_dos_metodos_en_niveles():
    t = tabla_atlas(_falso())
    assert "fmm_method" in t.index.names
    col = ("atlas_error_onda", "mediana")
    base = t.xs("base", level="fmm_method")[col]
    circa = t.xs("circafmm", level="fmm_method")[col]
    assert (circa.values < base.values).all()


# ── contraste de metodos ───────────────────────────────────────────────────

def test_contraste_metodos_es_pareado_contra_la_referencia():
    met = contraste_metodos(_falso())
    assert set(met["metodo"]) == {"best_k:medoid", "aggregate:tsp3:circmean",
                                  "aggregate:hods:circmean"}
    tsp = met[met.metodo == "aggregate:tsp3:circmean"]
    assert len(tsp) == 8 * 3
    assert np.allclose(tsp["delta_error_orden_h"], -0.1)
    med = met[met.metodo == "best_k:medoid"]
    assert np.allclose(med["delta_error_orden_h"], -0.05)


def test_contraste_metodos_solo_usa_el_fmm_pedido():
    a = contraste_metodos(_falso(), fmm="base")
    b = contraste_metodos(_falso(), fmm="circafmm")
    assert len(a) == len(b)


def test_efectos_principales_por_metodo():
    met = contraste_metodos(_falso())
    t = efectos_principales(met, "error_orden_h", por=["metodo"]).reset_index()
    assert set(t["metodo"]) == set(met["metodo"])
    hods = t[(t.metodo == "aggregate:hods:circmean") & (t.factor == "snr")]
    assert np.allclose(hods["mediana"], -0.1)


def test_tabla_tiempos_python_separa_fmm_y_metodo():
    t = tabla_tiempos_python(_falso())
    assert list(t.index.names) == ["fmm_method", "metodo", "n_genes", "n_samples"]
    assert "t_total_s" in t.columns


def test_sin_rejilla_falla_con_mensaje():
    solo = _falso().query("bloque == 'sin_ruido'")
    with pytest.raises(ValueError, match="rejilla"):
        contraste_fmm(solo)
