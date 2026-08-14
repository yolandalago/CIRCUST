"""Tests del lanzador de la rejilla.

Las piezas rapidas (reanudacion, resolucion de bloques, referencia calibrada)
se prueban sin tocar el pipeline. La integracion corre de verdad, sobre una
celda diminuta construida a mano, porque lo que hay que comprobar es que el
ciclo completo encaja.
"""
import json

import pandas as pd
import pytest

from simulation.grid import REFERENCIA, Celda, bloque_r
from simulation.run import (
    cargar_referencia, celdas_de, columnas_de_factores, lanzar,
    ruta_resultados, ya_hechas,
)


def celda_minima(**cambios) -> Celda:
    """Una celda de juguete, suficiente para que el pipeline arranque."""
    generate = {**REFERENCIA, "n_genes": 800, "n_samples": 24,
                "prop_rhythmic": 0.05, "snr": 2.0}
    generate.update(cambios.pop("generate", {}))
    base = dict(bloque="R", nombre="prueba", generate=generate,
                fmm_methods=("base",), metodos=("best_k:r2_median",),
                seed_base=999_000)
    base.update(cambios)
    return Celda(**base)


# ===========================================================================
# Resolucion de bloques y referencia calibrada
# ===========================================================================

def test_celdas_de_devuelve_las_celdas_del_bloque_pedido():
    assert len(celdas_de(["R"])) == len(bloque_r())


def test_celdas_de_admite_varios_bloques():
    assert len(celdas_de(["R", "A"])) == len(bloque_r()) + 27


def test_un_bloque_desconocido_falla_con_mensaje_claro():
    with pytest.raises(ValueError, match="inventado"):
        celdas_de(["inventado"])


def test_sin_referencia_calibrada_se_usa_la_provisional(tmp_path):
    assert cargar_referencia(tmp_path) == REFERENCIA


def test_la_referencia_calibrada_manda_sobre_la_provisional(tmp_path):
    (tmp_path / "referencia.json").write_text(
        json.dumps({**REFERENCIA, "snr": 0.5, "n_samples": 192}))
    ref = cargar_referencia(tmp_path)
    assert ref["snr"] == 0.5
    assert ref["n_samples"] == 192


def test_los_bloques_de_evaluacion_heredan_la_referencia_calibrada(tmp_path):
    (tmp_path / "referencia.json").write_text(
        json.dumps({**REFERENCIA, "snr": 3.0}))
    celdas = celdas_de(["A"], cargar_referencia(tmp_path))
    assert {c.generate["snr"] for c in celdas} == {1.0, 2.0, 3.0}


# ===========================================================================
# Reanudacion
# ===========================================================================

def test_sin_fichero_previo_no_hay_nada_hecho(tmp_path):
    assert ya_hechas(ruta_resultados("R", tmp_path)) == set()


def test_ya_hechas_lee_los_pares_celda_replica_del_fichero(tmp_path):
    ruta = ruta_resultados("R", tmp_path)
    ruta.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"celda": ["a", "a", "b"], "replica": [0, 1, 0],
                  "rho": [0.1, 0.2, 0.3]}).to_csv(ruta, index=False)
    assert ya_hechas(ruta) == {("a", 0), ("a", 1), ("b", 0)}


# ===========================================================================
# Integracion
# ===========================================================================

def test_lanzar_evalua_escribe_y_borra_la_matriz(tmp_path):
    d = lanzar([celda_minima()], raiz=tmp_path, replicas=2, verbose=False)
    assert len(d) == 2
    assert set(d["replica"]) == {0, 1}
    assert ruta_resultados("R", tmp_path).exists()
    # El ground truth se conserva, la matriz no.
    celda_dir = tmp_path / "R" / "prueba" / "rep00"
    assert (celda_dir / "params.csv").exists()
    assert not (celda_dir / "matriz.parquet").exists()


def test_las_filas_llevan_los_factores_y_las_metricas_contra_la_verdad(tmp_path):
    d = lanzar([celda_minima()], raiz=tmp_path, replicas=1, verbose=False)
    for col in columnas_de_factores():
        assert col in d.columns, col
    for col in ("verdad_picos_medae_h", "verdad_omega_sesgo",
                "verdad_top_precision", "verdad_sinc_reflejado",
                "sd_rho_reps", "rho"):
        assert col in d.columns, col
    assert d["waveform_mix"].iloc[0] == "gtex"


def test_relanzar_no_repite_lo_ya_hecho(tmp_path):
    celda = celda_minima()
    lanzar([celda], raiz=tmp_path, replicas=2, verbose=False)
    antes = pd.read_csv(ruta_resultados("R", tmp_path))
    nuevas = lanzar([celda], raiz=tmp_path, replicas=2, verbose=False)
    despues = pd.read_csv(ruta_resultados("R", tmp_path))
    assert nuevas.empty
    assert len(despues) == len(antes)


def test_relanzar_completa_las_replicas_que_faltaban(tmp_path):
    celda = celda_minima()
    lanzar([celda], raiz=tmp_path, replicas=1, verbose=False)
    lanzar([celda], raiz=tmp_path, replicas=3, verbose=False)
    d = pd.read_csv(ruta_resultados("R", tmp_path))
    assert sorted(d["replica"].unique()) == [0, 1, 2]


def test_forzar_rehace_lo_ya_hecho(tmp_path):
    celda = celda_minima()
    lanzar([celda], raiz=tmp_path, replicas=1, verbose=False)
    nuevas = lanzar([celda], raiz=tmp_path, replicas=1, forzar=True,
                    verbose=False)
    assert len(nuevas) == 1


def test_las_filas_de_error_no_desalinean_las_columnas(tmp_path):
    """Una fila corta anadida a un CSV existente no puede correr las columnas.

    Las filas de error traen muchas menos claves que las de exito. Escribirlas
    en modo anadir sin respetar la cabecera meteria el mensaje de error en la
    columna equivocada, y el fichero quedaria corrupto sin que nadie lo note
    hasta el analisis.
    """
    ok = celda_minima()
    rota = celda_minima(nombre="rota", generate={"prop_rhythmic": 0.001})
    lanzar([ok], raiz=tmp_path, replicas=1, verbose=False)
    lanzar([rota], raiz=tmp_path, replicas=1, verbose=False)

    d = pd.read_csv(ruta_resultados("R", tmp_path), keep_default_na=False)
    fallo = d[d["celda"] == "rota"].iloc[0]
    assert fallo["error"].startswith("ValueError")
    assert fallo["n_samples"] == 24
    # Las columnas que la fila de error no rellena quedan vacias, no ocupadas
    # por valores de otra columna.
    assert fallo["metodo"] == ""


def test_una_columna_nueva_reescribe_el_fichero_en_vez_de_encajarla(tmp_path):
    ruta = ruta_resultados("R", tmp_path)
    ruta.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([{"celda": "a", "replica": 0, "rho": 0.5}]).to_csv(
        ruta, index=False)
    from simulation.run import _anadir
    _anadir(ruta, [{"celda": "b", "replica": 0, "rho": 0.6, "nueva": 7}])

    d = pd.read_csv(ruta)
    assert list(d.columns) == ["celda", "replica", "rho", "nueva"]
    assert len(d) == 2
    assert d.loc[d["celda"] == "b", "nueva"].iloc[0] == 7
    assert pd.isna(d.loc[d["celda"] == "a", "nueva"].iloc[0])


def test_una_celda_rota_no_tumba_la_ejecucion(tmp_path):
    # prop_rhythmic tan baja que no salen ni los doce genes core, luego el
    # generador se niega. En una rejilla de horas eso no puede abortarlo todo.
    rota = celda_minima(nombre="rota", generate={"prop_rhythmic": 0.001})
    d = lanzar([rota, celda_minima()], raiz=tmp_path, replicas=1,
               verbose=False)
    assert set(d["celda"]) == {"rota", "prueba"}
    fallo = d[d["celda"] == "rota"].iloc[0]
    assert isinstance(fallo["error"], str) and fallo["error"]
    assert pd.isna(fallo["rho"])
