"""Tests de simulation/grid.py."""
import pytest

from simulation.grid import (
    BLOQUES, METODOS_RESUMEN, SNR_NIVELES, build_grid, cell_seed,
)


def test_todos_los_bloques_construyen():
    for bloque in BLOQUES:
        celdas = build_grid(bloque, reps=2)
        assert len(celdas) > 0, bloque


def test_recuentos_de_la_spec():
    esperado = {
        "principal":      135 * 2,
        "regularizacion": 135 * 2,
        "agregacion":     45 * 2 * len(METODOS_RESUMEN),
        "gcore":          45 * 2,
        "invariancia":    6 * 2,
        "computacional":  9 * 2,
    }
    assert len(build_grid("piloto", reps=3)) == len(SNR_NIVELES) * 3
    for bloque, n in esperado.items():
        assert len(build_grid(bloque, reps=2)) == n, bloque


def test_hay_cinco_metodos_de_resumen():
    etiquetas = [m.label for m in METODOS_RESUMEN]
    assert len(etiquetas) == 5
    assert len(set(etiquetas)) == 5
    assert "best_k:r2_median" in etiquetas
    assert "aggregate:tsp3:circmean" in etiquetas
    assert "aggregate:hods:circmean" in etiquetas


def test_principal_y_regularizacion_estan_emparejados():
    """
    Lo que hace pareada la comparacion: mismas matrices y mismas semillas,
    cambiando solo fmm_method.
    """
    p = build_grid("principal", reps=3)
    f = build_grid("regularizacion", reps=3)
    assert len(p) == len(f)
    for cp, cf in zip(p, f):
        assert cp.cfg == cf.cfg              # misma matriz exacta
        assert cp.cfg.seed == cf.cfg.seed
        assert cp.n_core_declared == cf.n_core_declared
        assert cp.fmm_method == "base"
        assert cf.fmm_method == "stabilized"


def test_las_semillas_son_deterministas_y_distintas_por_replica():
    a = cell_seed("principal", "x", 0)
    assert a == cell_seed("principal", "x", 0)
    assert a != cell_seed("principal", "x", 1)
    assert a != cell_seed("gcore", "x", 0)


def test_agregacion_barre_prop_rhythmic():
    celdas = build_grid("agregacion", reps=1)
    props = {c.cfg.prop_rhythmic for c in celdas}
    assert props == {0.01, 0.05, 0.20}


def test_invariancia_empareja_snr_con_escalas_distintas():
    celdas = build_grid("invariancia", reps=1)
    por_snr: dict[float, set[float]] = {}
    for c in celdas:
        por_snr.setdefault(round(c.cfg.snr, 6), set()).add(c.cfg.sigma)
    assert len(por_snr) == 3
    for sigmas in por_snr.values():
        assert len(sigmas) == 2      # dos escalas absolutas por cada SNR


def test_computacional_barre_el_tamano_y_usa_el_marco_estandar():
    celdas = build_grid("computacional", reps=1)
    assert {c.cfg.n_genes for c in celdas} == {1000, 5000, 20000}
    for c in celdas:
        assert c.cfg.sampling == "equiespaciado"
        assert c.cfg.prop_rhythmic == 0.05
        assert isinstance(c.cfg.omega, tuple)     # mezcla simetricas/asimetricas


def test_todas_las_celdas_son_generables():
    """
    Ninguna celda puede pedir mas genes core de los que caben en sus genes
    ritmicos. Con prop_rhythmic bajo y n_genes pequeno el generador falla, y
    ese fallo tiene que salir aqui y no a mitad de una ejecucion de horas.
    """
    for bloque in BLOQUES:
        for c in build_grid(bloque, reps=1):
            assert c.cfg.n_rhythmic >= c.cfg.n_core, (
                f"{bloque}: n_genes={c.cfg.n_genes} con "
                f"prop_rhythmic={c.cfg.prop_rhythmic} da "
                f"{c.cfg.n_rhythmic} ritmicos, menos que n_core={c.cfg.n_core}"
            )
            assert c.n_core_declared in c.cfg.core_sizes


def test_el_eje_de_senal_ruido_arranca_donde_hay_senal():
    """
    Medido en la Task 6: por debajo de A/sigma = 1 el orden no se reconstruye
    (oraculo 3,3 a 3,9 h) y gastar celdas ahi solo produce tasa de fallo. El
    1,0 se conserva como celda de frontera.
    """
    assert min(SNR_NIVELES) >= 1.0
    assert len(SNR_NIVELES) == 5
    # la transicion medida esta entre 1 y 2, y hay que muestrearla
    assert sum(1 for s in SNR_NIVELES if 1.0 <= s <= 2.0) >= 3


def test_bloque_sin_ruido_es_el_factorial_de_los_tres_factores():
    """
    Verificacion de que el pipeline funciona con datos perfectos. Es el
    factorial n x muestreo x omega, sin eje de senal-ruido porque aqui no
    aplica. Cualquier error residual no puede venir del ruido.
    """
    from simulation.grid import SIGMA_SIN_RUIDO
    celdas = build_grid("sin_ruido", reps=1)
    assert len(celdas) == 3 * 3 * 3
    assert {c.cfg.sigma for c in celdas} == {SIGMA_SIN_RUIDO}
    assert {c.cfg.n_samples for c in celdas} == {24, 48, 96}
    assert {c.cfg.sampling for c in celdas} == {"uniforme", "agrupado", "hueco"}
    assert {c.cfg.omega for c in celdas} == {0.1, 0.5, 0.9}
    # sigma no es 0 exacto: con 0 los genes de fondo quedan constantes y el
    # preprocesado los elimina, cambiando dos cosas a la vez
    assert SIGMA_SIN_RUIDO > 0.0
    assert all(c.cfg.snr > 1e5 for c in celdas)


def test_bloque_desconocido_falla():
    with pytest.raises(ValueError, match="Bloque desconocido"):
        build_grid("inventado", reps=1)
