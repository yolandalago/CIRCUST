"""Tests de simulation/generate.py."""
import numpy as np
import pytest

from simulation.generate import sample_times
from simulation.truth import SimConfig, TWO_PI


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
