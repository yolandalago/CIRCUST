"""
Tests de las metricas de calidad del atlas.

Comparan los parametros FMM estimados por gen contra los verdaderos. Son las
metricas del bloque de regularizacion, porque el FMM estabilizado no cambia
el orden de muestras (ese sale de una SVD donde el FMM no interviene) sino
la estimacion de los parametros por gen, que es lo que el atlas publica.

Los errores se guardan CON SIGNO para que el analisis pueda derivar el sesgo,
la dispersion y el error de subida sin ocupar columnas extra:

    error de subida = |error del pico - error del valle|
"""
import numpy as np
import pandas as pd
import pytest

from simulation.atlas import atlas_metrics
from simulation.truth import TWO_PI


def _ajuste(genes, M=0.0, A=1.0, alpha=1.0, beta=2.0, omega=0.7):
    return pd.DataFrame({
        "gene_name": genes,
        "r2_fmm":    [0.9] * len(genes),
        "M":         [M] * len(genes),
        "amplitude": [A] * len(genes),
        "alpha":     [alpha] * len(genes),
        "beta":      [beta] * len(genes),
        "omega":     [omega] * len(genes),
    })


def _verdad(genes, kinds, M=0.0, A=1.0, alpha=1.0, beta=2.0, omega=0.7):
    return pd.DataFrame({
        "gene": genes, "kind": kinds,
        "M": M, "A": A, "alpha": alpha, "beta": beta, "omega": omega,
    }).set_index("gene")


def test_ajuste_perfecto_da_error_cero():
    genes = ["RIT0001", "RIT0002", "CORE_A"]
    m = atlas_metrics(_ajuste(genes), _verdad(genes, ["ritmico"]*2 + ["core"]))
    assert m["atlas_rmse_curva"] == pytest.approx(0.0, abs=1e-9)
    assert m["atlas_err_pico_h"] == pytest.approx(0.0, abs=1e-9)
    assert m["atlas_err_valle_h"] == pytest.approx(0.0, abs=1e-9)
    assert m["atlas_omega_sesgo"] == pytest.approx(0.0, abs=1e-9)
    assert m["atlas_vpp"] == pytest.approx(1.0)


def test_los_errores_van_con_signo():
    """
    Guardarlos con signo permite derivar sesgo, dispersion y error de subida
    en el analisis, sin gastar columnas.
    """
    genes = ["RIT0001"]
    m = atlas_metrics(_ajuste(genes, omega=0.5), _verdad(genes, ["ritmico"], omega=0.7))
    assert m["atlas_omega_sesgo"] == pytest.approx(-0.2, abs=1e-9)
    m = atlas_metrics(_ajuste(genes, omega=0.9), _verdad(genes, ["ritmico"], omega=0.7))
    assert m["atlas_omega_sesgo"] == pytest.approx(+0.2, abs=1e-9)


def test_el_vpp_cuenta_los_genes_de_fondo_colados_en_el_atlas():
    """
    La pregunta practica del atlas: de los genes que contiene, cuantos son
    ritmicos de verdad. Un gen de fondo colado es un resultado falso
    publicado.
    """
    genes = [f"RIT{i:04d}" for i in range(6)] + ["GEN0001", "GEN0002"]
    kinds = ["ritmico"] * 6 + ["fondo"] * 2
    m = atlas_metrics(_ajuste(genes), _verdad(genes, kinds))
    assert m["atlas_vpp"] == pytest.approx(6 / 8)
    assert m["atlas_n_evaluados"] == 6      # el fondo no tiene onda que comparar


def test_los_genes_de_fondo_no_entran_en_los_errores():
    """No tienen onda verdadera, asi que compararlos no significa nada."""
    solo_rit = ["RIT0001", "RIT0002"]
    con_fondo = solo_rit + ["GEN0001"]
    a = atlas_metrics(_ajuste(solo_rit), _verdad(solo_rit, ["ritmico"]*2))
    b = atlas_metrics(_ajuste(con_fondo),
                      _verdad(con_fondo, ["ritmico"]*2 + ["fondo"]))
    assert a["atlas_rmse_curva"] == pytest.approx(b["atlas_rmse_curva"])
    assert a["atlas_vpp"] == 1.0 and b["atlas_vpp"] == pytest.approx(2/3)


def test_el_rmse_de_la_curva_no_depende_de_la_parametrizacion():
    """
    Dos parametrizaciones que dan la MISMA onda tienen que dar error cero.
    Es la propiedad que motiva usar la curva como medida principal: alpha,
    beta y omega interactuan y compararlos uno a uno penaliza diferencias
    que no existen en la onda.
    """
    genes = ["RIT0001"]
    # alpha y alpha + 2*pi describen exactamente la misma onda
    a = _ajuste(genes, alpha=1.0)
    b = _ajuste(genes, alpha=1.0 + TWO_PI)
    v = _verdad(genes, ["ritmico"], alpha=1.0)
    assert atlas_metrics(a, v)["atlas_rmse_curva"] == pytest.approx(
        atlas_metrics(b, v)["atlas_rmse_curva"], abs=1e-9)


def test_un_error_de_pico_se_mide_en_horas():
    genes = ["RIT0001"]
    from circust.fitting.fmm import FMMModel
    # se desplaza alpha, que desplaza el pico la misma cantidad
    dh = 2.0
    d_rad = dh * TWO_PI / 24.0
    m = atlas_metrics(_ajuste(genes, alpha=1.0 + d_rad),
                      _verdad(genes, ["ritmico"], alpha=1.0))
    assert abs(m["atlas_err_pico_h"]) == pytest.approx(dh, abs=1e-6)


def test_quitar_giro_descuenta_el_desplazamiento_comun_y_lo_registra():
    """
    Si la escala estimada esta girada, todos los picos se desplazan lo mismo.
    Ese es un error de anclaje, ya medido por mae_h sobre las muestras, y
    dejarlo dentro taparia la calidad del ajuste que es lo que aqui se mide.
    """
    genes = [f"RIT{i:04d}" for i in range(8)]
    d_rad = 3.0 * TWO_PI / 24.0
    a = _ajuste(genes, alpha=1.0 + d_rad)
    v = _verdad(genes, ["ritmico"] * 8, alpha=1.0)

    crudo = atlas_metrics(a, v)
    assert abs(crudo["atlas_err_pico_h"]) == pytest.approx(3.0, abs=1e-6)
    assert crudo["atlas_giro_h"] == pytest.approx(0.0)

    limpio = atlas_metrics(a, v, quitar_giro=True)
    assert limpio["atlas_err_pico_h"] == pytest.approx(0.0, abs=1e-9)
    assert limpio["atlas_rmse_curva"] == pytest.approx(0.0, abs=1e-9)
    assert abs(limpio["atlas_giro_h"]) == pytest.approx(3.0, abs=1e-6)


def test_quitar_giro_no_borra_el_error_que_difiere_entre_genes():
    """Quitar el giro comun no puede convertir un mal ajuste en uno bueno."""
    genes = ["RIT0001", "RIT0002"]
    a = _ajuste(genes)
    a.loc[0, "alpha"] = 1.0 + 2.0 * TWO_PI / 24.0    # +2 h
    a.loc[1, "alpha"] = 1.0 - 2.0 * TWO_PI / 24.0    # -2 h
    v = _verdad(genes, ["ritmico"] * 2, alpha=1.0)

    m = atlas_metrics(a, v, quitar_giro=True)
    assert m["atlas_err_pico_h"] == pytest.approx(0.0, abs=1e-6)  # se compensan
    assert m["atlas_err_pico_sd_h"] == pytest.approx(2.0, abs=1e-6)
    assert m["atlas_rmse_curva"] > 0.0


def test_sin_genes_evaluables_devuelve_nan_y_no_revienta():
    genes = ["GEN0001", "GEN0002"]
    m = atlas_metrics(_ajuste(genes), _verdad(genes, ["fondo"]*2))
    assert np.isnan(m["atlas_rmse_curva"])
    assert m["atlas_vpp"] == pytest.approx(0.0)
    assert m["atlas_n_evaluados"] == 0


def test_la_amplitud_no_se_compara_porque_el_min_max_la_borra():
    """
    El preprocesado normaliza min-max cada gen a [-1, 1]. Despues de eso M y
    A absolutos no son recuperables ni en principio, asi que compararlos
    daria un error que crece con el ruido aunque la onda este perfecta. La
    curva se compara estandarizada: mismo alpha, beta y omega tiene que dar
    error cero por mucho que la amplitud declarada sea otra.
    """
    genes = ["RIT0001"]
    v = _verdad(genes, ["ritmico"], M=3.0, A=4.0)
    m = atlas_metrics(_ajuste(genes, M=0.0, A=1.0), v)
    assert m["atlas_rmse_curva"] == pytest.approx(0.0, abs=1e-9)
    assert "atlas_err_A" not in m


def test_se_registra_el_r2_que_el_atlas_declara_de_si_mismo():
    """Es el indicador de calidad que si se calcula sin conocer la verdad."""
    genes = ["RIT0001", "RIT0002", "RIT0003"]
    a = _ajuste(genes)
    a["r2_fmm"] = [0.4, 0.8, 0.9]
    m = atlas_metrics(a, _verdad(genes, ["ritmico"] * 3))
    assert m["atlas_r2_mediana"] == pytest.approx(0.8)


# ---------------------------------------------------------------------------
# Tabla por gen, la que alimenta las figuras del capitulo
# ---------------------------------------------------------------------------

def test_per_gene_errors_da_una_fila_por_gen_evaluable():
    from simulation.atlas import per_gene_errors

    genes = ["RIT0001", "RIT0002", "CORE_A", "GEN0001"]
    kinds = ["ritmico", "ritmico", "core", "fondo"]
    g = per_gene_errors(_ajuste(genes), _verdad(genes, kinds))

    assert len(g) == 3                     # el de fondo queda fuera
    assert "GEN0001" not in set(g.gene_name)
    assert {"omega_true", "omega_hat", "err_omega", "rmse_curva",
            "pico_true_h", "pico_hat_h", "err_pico_h"} <= set(g.columns)


def test_per_gene_errors_y_atlas_metrics_no_pueden_discrepar():
    """
    Los resumenes salen de esta tabla, asi que tienen que reproducirse desde
    ella. Es lo que impide que las figuras por gen y las tablas del capitulo
    cuenten cosas distintas.
    """
    from simulation.atlas import atlas_metrics, per_gene_errors

    genes = [f"RIT{i:04d}" for i in range(5)]
    a = _ajuste(genes)
    a["omega"] = [0.5, 0.6, 0.7, 0.8, 0.9]
    a["alpha"] = [1.0, 1.1, 0.9, 1.2, 0.8]
    v = _verdad(genes, ["ritmico"] * 5)

    g = per_gene_errors(a, v, quitar_giro=True)
    m = atlas_metrics(a, v, quitar_giro=True)

    assert m["atlas_omega_sesgo"] == pytest.approx(g.err_omega.mean())
    assert m["atlas_rmse_curva"] == pytest.approx(g.rmse_curva.median())
    assert m["atlas_err_pico_h"] == pytest.approx(g.err_pico_h.mean())
    assert m["atlas_err_pico_sd_h"] == pytest.approx(g.err_pico_h.std(ddof=0))


def test_el_pico_por_gen_sale_en_horas_del_ciclo():
    from simulation.atlas import per_gene_errors

    genes = ["RIT0001"]
    g = per_gene_errors(_ajuste(genes), _verdad(genes, ["ritmico"]))
    assert 0.0 <= g.pico_true_h.iloc[0] < 24.0
    assert 0.0 <= g.pico_hat_h.iloc[0] < 24.0
