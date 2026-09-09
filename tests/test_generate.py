"""Tests de simulation/generate.py."""
import numpy as np
import pandas as pd
import pytest

from simulation.generate import (
    alpha_for_peak, assign_omega, fmm_wave, generate_dataset, sample_times,
)
from simulation.truth import ANCHOR_NAME, OPPOSITE_NAME, SimConfig, TWO_PI


def _cfg(**kw):
    base = dict(n_samples=200, seed=7)
    base.update(kw)
    return SimConfig(**base)


def test_uniforme_cubre_el_ciclo():
    t = sample_times(_cfg(sampling="uniforme"), np.random.default_rng(1))
    assert t.shape == (200,)
    assert t.min() >= 0.0 and t.max() < 24.0
    # cada cuarto del dia recibe muestras
    for lo in (0, 6, 12, 18):
        assert ((t >= lo) & (t < lo + 6)).sum() > 0


def test_equiespaciado_es_determinista_y_regular():
    cfg = _cfg(n_samples=24, sampling="equiespaciado")
    t = sample_times(cfg, np.random.default_rng(1))
    assert np.allclose(t, np.arange(24, dtype=float))
    # no consume aleatoriedad, dos rng distintos dan lo mismo
    assert np.allclose(t, sample_times(cfg, np.random.default_rng(999)))


def test_hueco_no_produce_muestras_en_la_ventana():
    t = sample_times(_cfg(sampling="hueco", gap=(10.0, 16.0)),
                     np.random.default_rng(3))
    assert not ((t >= 10.0) & (t < 16.0)).any()
    assert t.min() >= 0.0 and t.max() < 24.0


def _mayor_hueco(t):
    """Mayor hueco circular en horas, que es lo que ve el algoritmo."""
    s = np.sort(np.asarray(t) % 24.0)
    return float(max(np.diff(np.r_[s, s[0] + 24.0])))


def test_los_tres_escenarios_dejan_huecos_crecientes():
    """
    El eje de muestreo tiene que tener tres niveles **reales**, no solo tres
    nombres. Con cluster_sd = 2 h los tres grupos teselaban el circulo y el
    escenario agrupado salia indistinguible del uniforme (2,46 h contra
    2,11 h de mayor hueco). Este test es el que impide que eso vuelva.
    """
    def mediana(sampling):
        return float(np.median([
            _mayor_hueco(sample_times(_cfg(n_samples=48, sampling=sampling),
                                      np.random.default_rng(s)))
            for s in range(50)
        ]))

    uniforme  = mediana("uniforme")
    agrupado  = mediana("agrupado")
    con_hueco = mediana("hueco")

    assert uniforme < agrupado < con_hueco
    assert agrupado >= 1.5 * uniforme, (
        f"el escenario agrupado no se distingue del uniforme: "
        f"{agrupado:.2f} h contra {uniforme:.2f} h"
    )


def test_agrupado_concentra_alrededor_de_los_centros():
    cfg = _cfg(n_samples=3000, sampling="agrupado",
               cluster_centers=(4.0, 12.0, 20.0), cluster_sd=2.0)
    t = sample_times(cfg, np.random.default_rng(5))
    ang = t * TWO_PI / 24.0
    centros = np.asarray(cfg.cluster_centers) * TWO_PI / 24.0
    # distancia circular al centro mas cercano
    d = np.abs(ang[:, None] - centros[None, :]) % TWO_PI
    d = np.minimum(d, TWO_PI - d).min(axis=1)
    horas = d * 24.0 / TWO_PI
    # con sigma_t = 2 h, la mayoria cae a menos de 4 h de un centro
    assert (horas <= 4.0).mean() > 0.90
    assert t.min() >= 0.0 and t.max() < 24.0


def test_misma_semilla_mismos_tiempos():
    cfg = _cfg(sampling="agrupado")
    a = sample_times(cfg, np.random.default_rng(11))
    b = sample_times(cfg, np.random.default_rng(11))
    assert np.array_equal(a, b)


def test_muestreo_desconocido_falla():
    with pytest.raises(ValueError, match="Muestreo desconocido"):
        sample_times(_cfg(sampling="inventado"), np.random.default_rng(1))


# ---------------------------------------------------------------------------
# Onda FMM, inversion del pico y ensamblado de la matriz
# ---------------------------------------------------------------------------

def test_inversion_del_pico_ida_y_vuelta_es_exacta():
    """alpha_for_peak invierte FMMModel.peak_time analiticamente."""
    from circust.fitting.fmm import FMMModel
    rng = np.random.default_rng(0)
    for _ in range(200):
        theta_peak = rng.uniform(0.0, TWO_PI)
        beta       = rng.uniform(0.0, TWO_PI)
        omega      = float(rng.choice([0.1, 0.3, 0.5, 0.7, 0.9]))
        alpha      = alpha_for_peak(theta_peak, beta, omega)
        recuperado = FMMModel.peak_time(alpha, beta, omega)
        d = abs(recuperado - theta_peak) % TWO_PI
        assert min(d, TWO_PI - d) < 1e-9


def test_el_maximo_de_la_onda_cae_en_el_pico_pedido():
    """Comprobacion independiente: el argmax sobre una rejilla fina."""
    rejilla = np.linspace(0.0, TWO_PI, 200_000, endpoint=False)
    rng = np.random.default_rng(1)
    for _ in range(20):
        theta_peak = rng.uniform(0.0, TWO_PI)
        beta       = rng.uniform(0.0, TWO_PI)
        omega      = float(rng.choice([0.1, 0.5, 0.9]))
        alpha      = alpha_for_peak(theta_peak, beta, omega)
        y          = fmm_wave(rejilla, 0.0, 1.0, alpha, beta, omega)
        obtenido   = rejilla[int(np.argmax(y))]
        d = abs(obtenido - theta_peak) % TWO_PI
        assert min(d, TWO_PI - d) < 1e-3


def test_assign_omega_escalar_y_mezcla():
    assert np.allclose(assign_omega(0.5, 7), 0.5)
    mezcla = assign_omega((0.1, 0.9), 10)
    assert sorted(mezcla.tolist()) == [0.1] * 5 + [0.9] * 5
    # reparto por resto mayor, exacto tambien con n no divisible
    assert len(assign_omega((0.1, 0.5, 0.9), 10)) == 10


def test_generate_dataset_forma_y_composicion():
    cfg = SimConfig(n_samples=48, n_genes=500, prop_rhythmic=0.20, n_core=30,
                    omega=0.5, amplitude=1.0, sigma=1.0, seed=42)
    X, truth = generate_dataset(cfg)
    assert X.shape == (500, 48)
    assert len(truth.theta_true) == 48
    kinds = truth.gene_params["kind"].value_counts()
    assert kinds["core"] == 30
    assert kinds["ritmico"] == 70      # 100 ritmicos menos 30 core
    assert kinds["fondo"] == 400
    assert ANCHOR_NAME in X.index and OPPOSITE_NAME in X.index


def test_generate_dataset_es_reproducible():
    cfg = SimConfig(n_samples=24, n_genes=200, seed=3)
    a, ta = generate_dataset(cfg)
    b, tb = generate_dataset(cfg)
    pd.testing.assert_frame_equal(a, b)
    assert np.array_equal(ta.theta_true, tb.theta_true)


def test_las_filas_estan_barajadas():
    """Ningun indice de fila debe filtrar el tipo de gen al pipeline."""
    cfg = SimConfig(n_samples=24, n_genes=300, seed=9)
    X, _ = generate_dataset(cfg)
    primeras = list(X.index[:30])
    assert not all(n.startswith("CORE_") for n in primeras)


def test_sin_ruido_los_core_ajustan_su_pico():
    """Con sigma=0 el FMM recupera el pico que se fijo por construccion."""
    from circust.fitting.fmm import FMMModel
    cfg = SimConfig(n_samples=96, n_genes=200, prop_rhythmic=0.30, n_core=10,
                    core_sizes=(10,), omega=0.5, amplitude=1.0, sigma=0.0,
                    sampling="equiespaciado", seed=4)
    X, truth = generate_dataset(cfg)
    theta = truth.theta_true
    orden = np.argsort(theta)
    peak_of = dict(zip(truth.core_names, truth.core_peaks))
    fit = FMMModel().fit(X.loc[ANCHOR_NAME].values[orden], theta[orden])
    d = abs(fit.peak_time - peak_of[ANCHOR_NAME]) % TWO_PI
    assert min(d, TWO_PI - d) < 0.05


def test_prop_rhythmic_menor_que_n_core_falla():
    cfg = SimConfig(n_samples=24, n_genes=100, prop_rhythmic=0.05, n_core=30)
    with pytest.raises(ValueError, match="genes ritmicos"):
        generate_dataset(cfg)
