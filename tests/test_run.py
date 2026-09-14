"""Tests de simulation/run.py."""
import numpy as np
import pandas as pd
import pytest

from simulation.grid import build_grid
from simulation.run import COLUMNAS, existing_ids, run_cell_group


def test_el_esquema_del_csv_tiene_lo_que_pide_la_spec():
    necesarias = {
        "bloque", "row_id", "rep", "seed",
        "n_samples", "sampling", "omega", "A", "sigma", "snr",
        "n_genes", "prop_rhythmic", "n_core",
        "fmm_method", "robust_method", "selection_criterion",
        "aggregation", "consensus_scale",
        "status", "error",
        # metricas
        "mae_h", "medae_h", "rho_circ", "tau", "tau_circ", "orientacion_ok",
        "p1", "p2", "p4", "mae_h_oracle",
        # observables sin la verdad
        "r2_anchor", "pares_ok", "evidencia", "desacuerdo_K",
        "n_top", "n_outliers", "n_evaluadas", "direction_flips",
        "t_total_s",
    }
    assert necesarias <= set(COLUMNAS)


def test_existing_ids_lee_un_csv_previo(tmp_path):
    destino = tmp_path / "resultados.csv"
    assert existing_ids(destino) == set()
    pd.DataFrame({"row_id": ["a|r0", "b|r1"]}).to_csv(destino, index=False)
    assert existing_ids(destino) == {"a|r0", "b|r1"}


def test_una_celda_sin_ruido_produce_una_fila_completa():
    """Celda facil: omega 0.5, que sabemos que reconstruye exacto."""
    celda = next(c for c in build_grid("sin_ruido", reps=1)
                 if c.cfg.omega == 0.5 and c.cfg.n_samples == 48
                 and c.cfg.sampling == "uniforme")
    filas = run_cell_group(celda.cfg, celda.n_core_declared, celda.fmm_method,
                           [celda.robust_cfg], bloque=celda.bloque, rep=celda.rep)
    assert len(filas) == 1
    fila = filas[0]
    assert set(fila) == set(COLUMNAS)
    assert fila["status"] == "ok"
    assert fila["n_evaluadas"] > 0
    # sin ruido y con omega moderada la reconstruccion tiene que ser buena
    assert fila["tau_circ"] > 0.9
    assert fila["orientacion_ok"] == 1.0


def test_un_fallo_de_las_etapas_0_a_3_se_registra_y_no_revienta(monkeypatch):
    """La tasa de fallo es un resultado, no un bug que esconder."""
    import simulation.run as modulo
    from simulation.runner import RobustConfig
    from simulation.truth import SimConfig

    def revienta(*a, **kw):
        raise ValueError("Quedan 1 genes tras filtrar el conjunto TOP")

    monkeypatch.setattr(modulo, "run_stages_0_3", revienta)
    cfg = SimConfig(n_samples=24, n_genes=300, prop_rhythmic=0.20, n_core=10,
                    core_sizes=(10,), seed=1)
    filas = run_cell_group(cfg, 10, "base", [RobustConfig(n_reps=3)],
                           bloque="prueba", rep=0)
    assert len(filas) == 1
    assert filas[0]["status"] == "fallo"
    assert "conjunto TOP" in filas[0]["error"]
    assert np.isnan(filas[0]["mae_h"])
    assert not np.isnan(filas[0]["t_total_s"])


def test_un_fallo_de_la_etapa_4_solo_afecta_a_su_metodo(monkeypatch):
    """Si un metodo de resumen falla, los demas de la misma celda siguen."""
    import simulation.run as modulo
    from simulation.runner import RobustConfig
    from simulation.truth import SimConfig

    real = modulo.run_stage_4
    llamadas = {"n": 0}

    def a_veces(estado, truth, robust_cfg, fmm_method="base"):
        llamadas["n"] += 1
        if llamadas["n"] == 1:
            raise RuntimeError("fallo simulado del primer metodo")
        return real(estado, truth, robust_cfg, fmm_method)

    monkeypatch.setattr(modulo, "run_stage_4", a_veces)
    cfg = SimConfig(n_samples=48, n_genes=400, prop_rhythmic=0.20, n_core=20,
                    core_sizes=(10, 20), omega=0.5, amplitude=4.0, sigma=1.0,
                    seed=77)
    filas = run_cell_group(cfg, 20, "base",
                           [RobustConfig(n_reps=3), RobustConfig(n_reps=3, seed=8)],
                           bloque="prueba", rep=0)
    assert [f["status"] for f in filas] == ["fallo", "ok"]


def test_los_metodos_comparten_las_etapas_0_a_3():
    """Un solo grupo produce una fila por metodo, con el mismo n_top."""
    from simulation.grid import METODOS_RESUMEN
    from simulation.truth import SimConfig
    cfg = SimConfig(n_samples=48, n_genes=400, prop_rhythmic=0.20, n_core=20,
                    core_sizes=(10, 20), omega=0.5, amplitude=4.0, sigma=1.0,
                    seed=55)
    filas = run_cell_group(cfg, 20, "base", METODOS_RESUMEN[:3],
                           bloque="prueba", rep=0)
    assert len(filas) == 3
    assert len({f["n_top"] for f in filas}) == 1
    # y la configuracion del metodo queda registrada en su fila
    assert {f["selection_criterion"] for f in filas} >= {"r2_median", "medoid"}


def test_el_csv_se_escribe_por_lotes_y_no_al_final(tmp_path, monkeypatch):
    """
    Un bloque tarda horas. Antes el CSV no se tocaba hasta que terminaba todo,
    de modo que una interrupcion tiraba el trabajo entero. Se comprueba que
    hay filas en disco antes de que acabe la ejecucion.
    """
    import simulation.run as run_mod

    destino = tmp_path / "resultados.csv"
    vistos: list[int] = []

    # el Parallel falso no mira el contenido de los grupos, solo cuenta
    grupos = [("sin_ruido", None, 10, "base", i, []) for i in range(6)]

    def falso_parallel(*a, **kw):
        assert kw.get("return_as") == "generator_unordered", (
            "sin generador las filas no pueden escribirse segun llegan"
        )
        def corre(tareas):
            for _ in tareas:
                # cuantas filas hay en disco en este momento
                vistos.append(
                    sum(1 for _ in destino.open()) - 1 if destino.exists() else 0
                )
                yield [{c: None for c in run_mod.COLUMNAS}
                       | {"status": "ok", "t_total_s": 1.0, "row_id": str(len(vistos))}]
        return corre

    monkeypatch.setattr(run_mod, "LOTE_ESCRITURA", 2)
    monkeypatch.setattr(run_mod, "Parallel", falso_parallel)
    monkeypatch.setattr(run_mod, "build_grid", lambda b, r: [])
    monkeypatch.setattr(run_mod, "_agrupar", lambda c: grupos)
    monkeypatch.setattr(
        "sys.argv",
        ["run.py", "--bloque", "sin_ruido", "--out", str(tmp_path)],
    )
    run_mod.main()

    assert destino.exists()
    assert sum(1 for _ in destino.open()) - 1 == 6
    assert max(vistos) > 0, (
        "el CSV seguia vacio al llegar al ultimo grupo, no se escribio por lotes"
    )
