"""Tests de simulation/runner.py."""
import inspect

import numpy as np
import pandas as pd
import pytest

from simulation.generate import generate_dataset
from simulation.runner import (
    RobustConfig, run_full, run_stage_4, run_stages_0_3,
)
from simulation.truth import SimConfig


@pytest.fixture(scope="module")
def celda_pequena():
    """Celda minima con senal fuerte, para que el pipeline no falle."""
    cfg = SimConfig(n_samples=48, n_genes=400, prop_rhythmic=0.20,
                    n_core=20, core_sizes=(10, 20), omega=0.5,
                    amplitude=4.0, sigma=1.0, sampling="uniforme", seed=101)
    return generate_dataset(cfg)


def test_run_full_devuelve_fases_por_nombre(celda_pequena):
    X, truth = celda_pequena
    res = run_full(X, truth, n_core_declared=20, robust_cfg=RobustConfig(n_reps=3))
    assert isinstance(res.theta_hat, pd.Series)
    assert set(res.theta_hat.index) <= set(truth.sample_names)
    assert len(res.theta_hat) >= 0.9 * len(truth.sample_names)
    assert ((res.theta_hat >= 0.0) & (res.theta_hat < 2 * np.pi)).all()
    assert res.n_top > 0


def test_partir_el_pipeline_da_el_mismo_resultado(celda_pequena):
    """
    Condicion para que el atajo del bloque de agregacion sea legitimo:
    correr las etapas 0 a 3 y luego la 4 tiene que dar exactamente lo mismo
    que correrlo entero de una vez.
    """
    X, truth = celda_pequena
    rc = RobustConfig(n_reps=3, seed=7)
    entero  = run_full(X, truth, 20, rc)
    estado  = run_stages_0_3(X, truth, 20)
    partido = run_stage_4(estado, truth, rc)
    pd.testing.assert_series_equal(entero.theta_hat, partido.theta_hat)


def test_el_estado_intermedio_se_reutiliza_entre_metodos(celda_pequena):
    X, truth = celda_pequena
    estado = run_stages_0_3(X, truth, 20)
    a = run_stage_4(estado, truth, RobustConfig(robust_method="best_k", n_reps=3, seed=7))
    b = run_stage_4(estado, truth, RobustConfig(robust_method="aggregate",
                                                aggregation="tsp3", n_reps=3, seed=7))
    assert len(a.theta_hat) == len(b.theta_hat)
    assert set(a.theta_hat.index) == set(b.theta_hat.index)


def test_runner_no_redeclara_la_configuracion_del_algoritmo():
    """Decision 1 de la spec: los parametros vienen de run_pipeline."""
    import simulation.runner as runner
    fuente = inspect.getsource(runner)
    assert "from run_pipeline import DEFAULT_CONFIG" in fuente
    for clave in ("sparse_threshold", "tight_radius", "loose_radius",
                  "fmm_alpha_grid", "fmm_omega_grid", "fmm_reps"):
        # solo puede aparecer leyendose del diccionario, nunca asignandose
        assert f'"{clave}":' not in fuente, f"{clave} redeclarado en runner.py"


def test_los_pares_de_consistencia_sinteticos_llegan_a_las_dos_etapas():
    import simulation.runner as runner
    fuente = inspect.getsource(runner)
    # una vez para CircularSynchronizer, otra para RobustOrderEstimator
    assert fuente.count("consistency_pairs = ") >= 2


# ---------------------------------------------------------------------------
# Indicadores observables sin la verdad
# ---------------------------------------------------------------------------

def test_los_observables_salen_en_rango(celda_pequena):
    X, truth = celda_pequena
    res = run_full(X, truth, 20, RobustConfig(n_reps=3))
    assert 0.0 <= res.r2_anchor <= 1.0
    assert 0.0 <= res.pares_ok <= 1.0
    assert 0.0 <= res.desacuerdo_K <= 1.0
    assert 0 <= res.direction_flips <= 3
    assert res.n_top > 0
    assert res.n_outliers >= 0


def test_con_senal_fuerte_el_ancla_se_ajusta_bien(celda_pequena):
    """
    Con A/sigma = 4 el gen ancla tiene que quedar bien ajustado. Si r2_anchor
    saliera bajo aqui, el indicador no serviria para distinguir regimenes.
    """
    X, truth = celda_pequena
    estado = run_stages_0_3(X, truth, 20)
    assert estado.r2_anchor > 0.8


def test_el_ancla_se_ajusta_peor_cuando_baja_la_senal():
    """
    Lo que hace util a r2_anchor: tiene que separar regimenes sin usar la
    verdad. Se compara senal fuerte contra senal debil sobre la misma semilla.
    """
    def r2_de(snr):
        cfg = SimConfig(n_samples=48, n_genes=400, prop_rhythmic=0.20,
                        n_core=20, core_sizes=(10, 20), omega=0.5,
                        amplitude=snr, sigma=1.0, seed=202)
        X, truth = generate_dataset(cfg)
        return run_stages_0_3(X, truth, 20).r2_anchor

    assert r2_de(4.0) > r2_de(0.25)


def test_la_evidencia_a_priori_se_registra_y_esta_acotada(celda_pequena):
    """
    ``evidencia`` es la suma de |sin(delta)| sobre los pares evaluables, o sea
    cuanta informacion sobre la orientacion habia realmente disponible. Su
    maximo es el numero de pares, que se alcanza cuando todos estan separados
    6 h. Es la variable candidata a predecir si se puede confiar en el sentido
    de giro sin conocer la verdad.
    """
    X, truth = celda_pequena
    n_pares = len(truth.consistency_pairs())
    res = run_full(X, truth, 20, RobustConfig(n_reps=3))
    assert 0.0 <= res.evidencia <= n_pares


def test_la_evidencia_cae_cuando_los_core_se_ajustan_mal():
    """
    Lo que la hace util: con senal baja los genes core que portan el
    conocimiento a priori estan mal ajustados, no superan el filtro de R2 y
    la evidencia se desploma. Es la senal de que la orientacion no es
    decidible, y se calcula sin usar la verdad.
    """
    def evidencia_de(snr):
        cfg = SimConfig(n_samples=48, n_genes=400, prop_rhythmic=0.20,
                        n_core=20, core_sizes=(10, 20), omega=0.5,
                        amplitude=snr, sigma=1.0, seed=303)
        X, truth = generate_dataset(cfg)
        return run_stages_0_3(X, truth, 20).evidencia

    fuerte, debil = evidencia_de(4.0), evidencia_de(0.25)
    assert fuerte > debil or np.isnan(debil)


# ---------------------------------------------------------------------------
# Atlas: calidad del ajuste por gen
# ---------------------------------------------------------------------------

def test_run_full_publica_el_atlas(celda_pequena):
    """
    El atlas es lo que se compara entre FMM base y estabilizado, porque el
    orden de muestras sale de una SVD donde el FMM no interviene.
    """
    X, truth = celda_pequena
    res = run_full(X, truth, 20, RobustConfig(n_reps=3))
    assert res.atlas is not None
    assert len(res.atlas) == res.n_top
    for c in ("gene_name", "M", "amplitude", "alpha", "beta", "omega"):
        assert c in res.atlas.columns
    # todos los genes del atlas existen en la verdad, que es lo que permite
    # compararlos gen a gen
    assert set(res.atlas["gene_name"]) <= set(truth.gene_params.index)


def test_las_metricas_del_atlas_salen_en_rango(celda_pequena):
    from simulation.atlas import atlas_metrics

    X, truth = celda_pequena
    res = run_full(X, truth, 20, RobustConfig(n_reps=3))
    m = atlas_metrics(res.atlas, truth.gene_params, quitar_giro=True)

    assert 0.0 <= m["atlas_vpp"] <= 1.0
    assert m["atlas_n_evaluados"] > 0
    # con A/sigma = 4 el ajuste tiene que ser bueno, si no la medida no
    # discriminaria nada mas abajo
    assert m["atlas_rmse_curva"] < 0.5
    assert abs(m["atlas_err_pico_h"]) < 3.0


def test_se_puede_pedir_el_pipeline_sin_atlas(celda_pequena):
    """Valvula de escape: con aggregate el atlas se re-ajusta y cuesta."""
    X, truth = celda_pequena
    res = run_full(X, truth, 20, RobustConfig(n_reps=3), con_atlas=False)
    assert res.atlas is None
