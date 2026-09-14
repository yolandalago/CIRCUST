"""
Tests de las figuras y la tabla de la comparacion FMM base vs estabilizado.

Lo que se fija aqui no es el aspecto de los dibujos sino las decisiones de
analisis que hay detras, que son las que pueden hacer que el capitulo afirme
algo falso.
"""
import numpy as np
import pandas as pd
import pytest

from simulation.fmm_figuras import (
    _benjamini_hochberg, cargar_pareado, tabla_wilcoxon,
)


def _fila(bloque, fmm, rep, snr=2.0, omega=0.5, n=48, rmse=0.3, **kw):
    base = dict(
        bloque=bloque, fmm_method=fmm, rep=rep, snr=snr, omega=omega,
        n_samples=n, sampling="uniforme", status="ok", orientacion_ok=1.0,
        atlas_rmse_curva=rmse, atlas_err_pico_sd_h=1.0,
        atlas_omega_sesgo=-0.1, atlas_vpp=0.95, atlas_err_subida_h=1.0,
    )
    base.update(kw)
    return base


def _csv(tmp_path, filas):
    d = tmp_path / "resultados.csv"
    pd.DataFrame(filas).to_csv(d, index=False)
    return tmp_path


def test_el_emparejamiento_une_la_misma_replica_de_los_dos_brazos(tmp_path):
    filas = []
    for rep in range(4):
        filas.append(_fila("principal", "base", rep, rmse=0.40))
        filas.append(_fila("regularizacion", "stabilized", rep, rmse=0.25))
    par = cargar_pareado(_csv(tmp_path, filas))

    assert len(par) == 4
    # negativo = gana el estabilizado
    assert par["d_atlas_rmse_curva"].tolist() == pytest.approx([-0.15] * 4)


def test_no_se_emparejan_replicas_distintas(tmp_path):
    """
    Toda la potencia del diseno viene de que las dos filas corrieron sobre la
    MISMA matriz. Si las replicas no casan, no hay nada que emparejar y hay
    que fallar en vez de comparar cosas distintas.
    """
    filas = [_fila("principal", "base", 0),
             _fila("regularizacion", "stabilized", 99)]
    with pytest.raises(ValueError, match="no cuadra"):
        cargar_pareado(_csv(tmp_path, filas))


def test_se_descartan_las_replicas_con_la_orientacion_invertida(tmp_path):
    """
    Una escala reflejada da un atlas especular, que es otro regimen y no una
    calidad peor. Basta con que falle UN brazo para que el par no sirva.
    """
    filas = []
    for rep in range(4):
        # en la replica 0 falla el brazo estabilizado
        filas.append(_fila("principal", "base", rep))
        filas.append(_fila("regularizacion", "stabilized", rep,
                           orientacion_ok=0.0 if rep == 0 else 1.0))
    par = cargar_pareado(_csv(tmp_path, filas))
    assert len(par) == 3
    assert 0 not in set(par["rep"])


def test_la_tasa_de_victorias_cuenta_replicas_no_magnitud():
    """
    Es la medida que distingue una mejora consistente de una que depende de
    unas pocas replicas con efecto grande. Un metodo que gana poquito pero
    siempre es preferible a uno que gana mucho la mitad de las veces.
    """
    filas = []
    # nueve victorias pequenas y una derrota enorme: la mediana favorece al
    # estabilizado y la tasa de victorias tambien, pero la MEDIA no
    deltas = [-0.01] * 9 + [+5.0]
    for rep, d in enumerate(deltas):
        filas.append(_fila("principal", "base", rep, rmse=0.5))
        filas.append(_fila("regularizacion", "stabilized", rep, rmse=0.5 + d))
    par = pd.DataFrame(filas)
    par = cargar_pareado_desde(par)

    t = tabla_wilcoxon(par)
    assert len(t) == 1
    assert t["tasa_victorias"].iloc[0] == pytest.approx(0.9)
    assert t["delta_mediana"].iloc[0] < 0


def cargar_pareado_desde(df: pd.DataFrame) -> pd.DataFrame:
    """Version en memoria de ``cargar_pareado``, para no tocar disco."""
    import tempfile
    from pathlib import Path
    with tempfile.TemporaryDirectory() as tmp:
        df.to_csv(Path(tmp) / "resultados.csv", index=False)
        return cargar_pareado(tmp)


def test_wilcoxon_no_inventa_p_valores_con_pocos_pares():
    """Con menos de seis pares el contraste no tiene potencia ninguna."""
    filas = []
    for rep in range(3):
        filas.append(_fila("principal", "base", rep, rmse=0.5))
        filas.append(_fila("regularizacion", "stabilized", rep, rmse=0.3))
    t = tabla_wilcoxon(cargar_pareado_desde(pd.DataFrame(filas)))
    assert np.isnan(t["p"].iloc[0])
    assert np.isnan(t["p_bh"].iloc[0])
    # pero la descriptiva si se reporta
    assert t["delta_mediana"].iloc[0] == pytest.approx(-0.2)


def test_benjamini_hochberg_es_monotono_y_no_baja_ningun_p():
    p = np.array([0.001, 0.01, 0.02, 0.04, 0.5])
    q = _benjamini_hochberg(p)
    assert np.all(q >= p - 1e-12), "un ajuste nunca puede hacer mas significativo"
    assert np.all(np.diff(q) >= -1e-12), "tiene que respetar el orden"
    assert np.all(q <= 1.0)


def test_benjamini_hochberg_ignora_los_nan():
    q = _benjamini_hochberg(np.array([0.01, np.nan, 0.5]))
    assert np.isnan(q[1])
    assert not np.isnan(q[0]) and not np.isnan(q[2])


def test_las_celdas_se_contrastan_por_separado():
    """
    Las 135 celdas son regimenes distintos y el efecto cambia de signo entre
    ellos. Agregarlas antes de contrastar borraria justo el resultado.
    """
    filas = []
    for rep in range(10):
        # celda A: gana el estabilizado. celda B: gana el base
        filas.append(_fila("principal", "base", rep, snr=1.0, rmse=0.5))
        filas.append(_fila("regularizacion", "stabilized", rep, snr=1.0, rmse=0.3))
        filas.append(_fila("principal", "base", rep, snr=4.0, rmse=0.2))
        filas.append(_fila("regularizacion", "stabilized", rep, snr=4.0, rmse=0.3))
    t = tabla_wilcoxon(cargar_pareado_desde(pd.DataFrame(filas)))

    assert len(t) == 2
    baja = t[t.snr == 1.0].iloc[0]
    alta = t[t.snr == 4.0].iloc[0]
    assert baja["delta_mediana"] < 0 and baja["tasa_victorias"] == 1.0
    assert alta["delta_mediana"] > 0 and alta["tasa_victorias"] == 0.0
