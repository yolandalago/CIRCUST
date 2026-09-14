"""Tests de simulation/grid.py."""
import pytest

from simulation.grid import (
    BLOQUES, FMM_METODOS, GCORE_NIVELES, GENES_PROP, METODO_REFERENCIA,
    METODOS_RESUMEN, MUESTREOS, N_NIVELES, OMEGAS, SIGMA_SIN_RUIDO,
    SNR_NIVELES, build_grid, cell_seed, escenarios_rejilla,
)


def test_todos_los_bloques_construyen():
    for bloque in BLOQUES:
        celdas = build_grid(bloque, reps=1)
        assert len(celdas) > 0, bloque


def test_recuentos_de_la_rejilla():
    n_escenarios = (len(N_NIVELES) * len(MUESTREOS) * len(OMEGAS)
                    * len(SNR_NIVELES) * len(GENES_PROP))
    assert n_escenarios == 3 * 3 * 3 * 4 * 5 == 540
    assert len(escenarios_rejilla()) == n_escenarios
    por_escenario = len(GCORE_NIVELES) * len(FMM_METODOS) * len(METODOS_RESUMEN)
    assert por_escenario == 3 * 2 * 5
    assert len(build_grid("rejilla", reps=2)) == n_escenarios * por_escenario * 2


def test_recuentos_de_los_bloques_de_verificacion():
    esperado = {
        "sin_ruido":      27 * 2,
        "sin_ruido_stab": 27 * 2,
        "invariancia":    6 * 2,
        "computacional":  9 * 2,
    }
    for bloque, n in esperado.items():
        assert len(build_grid(bloque, reps=2)) == n, bloque


def test_hay_cinco_metodos_de_resumen_y_una_referencia():
    etiquetas = [m.label for m in METODOS_RESUMEN]
    assert len(etiquetas) == 5
    assert len(set(etiquetas)) == 5
    assert "best_k:r2_median" in etiquetas
    assert "aggregate:tsp3:circmean" in etiquetas
    assert "aggregate:hods:circmean" in etiquetas
    # la referencia es el comportamiento actual de CIRCUST antes de agregar
    assert METODO_REFERENCIA.label == "best_k:r2_median"


def test_la_rejilla_es_pareada_en_todas_las_variantes():
    """
    Lo que hace pareada la comparacion: para un escenario y una replica,
    las 30 variantes (3 core x 2 FMM x 5 metodos) llevan la MISMA semilla y
    por tanto corren sobre la misma matriz.
    """
    celdas = build_grid("rejilla", reps=2)
    grupos: dict[tuple, list] = {}
    for c in celdas:
        grupos.setdefault((c.cfg, c.rep), []).append(c)
    for (cfg, rep), variantes in grupos.items():
        assert len(variantes) == 30
        assert len({c.cfg.seed for c in variantes}) == 1
        assert {c.fmm_method for c in variantes} == set(FMM_METODOS)
        assert {c.n_core_declared for c in variantes} == set(GCORE_NIVELES)
        assert len({c.robust_cfg.label for c in variantes}) == 5


def test_replicas_distintas_llevan_semillas_distintas():
    celdas = build_grid("rejilla", reps=3)
    cfg0 = escenarios_rejilla()[0]
    semillas = {c.cfg.seed for c in celdas
                if c.cfg.n_samples == cfg0.n_samples
                and c.cfg.sampling == cfg0.sampling
                and c.cfg.omega == cfg0.omega
                and c.cfg.amplitude == cfg0.amplitude
                and c.cfg.n_genes == cfg0.n_genes
                and c.cfg.prop_rhythmic == cfg0.prop_rhythmic}
    assert len(semillas) == 3


def test_las_replicas_van_en_el_bucle_exterior():
    """
    Si la ejecucion se interrumpe a medias, lo que hay en disco tiene que ser
    la rejilla entera con menos replicas, no unas celdas completas y otras sin
    empezar. Eso exige que la replica sea el bucle mas externo.
    """
    celdas = build_grid("rejilla", reps=2)
    reps = [c.rep for c in celdas]
    mitad = len(reps) // 2
    assert set(reps[:mitad]) == {0}
    assert set(reps[mitad:]) == {1}


def test_las_semillas_son_deterministas_y_distintas_por_replica():
    a = cell_seed("rejilla", "x", 0)
    assert a == cell_seed("rejilla", "x", 0)
    assert a != cell_seed("rejilla", "x", 1)
    assert a != cell_seed("sin_ruido", "x", 0)


def test_no_hay_row_id_repetidos():
    celdas = build_grid("rejilla", reps=2)
    ids = [c.row_id for c in celdas]
    assert len(ids) == len(set(ids))


def test_genes_y_prop_cubren_los_regimenes_de_senal():
    """
    Cinco combinaciones que cubren desde un TOP que no llega a llenarse (50
    ritmicos) hasta uno saturado (1000), y dos pares con los mismos 50
    ritmicos y distinto fondo, que aislan el efecto del numero de genes.
    """
    ritmicos = sorted(int(round(g * p)) for g, p in GENES_PROP)
    assert ritmicos == [50, 50, 200, 250, 1000]
    assert (1000, 0.01) not in GENES_PROP     # 10 ritmicos no alojan 30 core


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
    Por debajo de A/sigma = 1 el orden no se reconstruye (oraculo 3,3 a 3,9 h)
    y gastar celdas ahi solo produce tasa de fallo. El 1,0 se conserva como
    celda de frontera. El 3 se quito tras medir que 3 y 4 son indistinguibles.
    """
    assert min(SNR_NIVELES) >= 1.0
    assert SNR_NIVELES == (1.0, 1.5, 2.0, 4.0)
    # la transicion medida esta entre 1 y 2, y hay que muestrearla
    assert sum(1 for s in SNR_NIVELES if 1.0 <= s <= 2.0) >= 3


def test_bloque_sin_ruido_es_el_factorial_de_los_tres_factores():
    """
    Verificacion de que el pipeline funciona con datos perfectos. Es el
    factorial n x muestreo x omega, sin eje de senal-ruido porque aqui no
    aplica. Cualquier error residual no puede venir del ruido.
    """
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
    with pytest.raises(ValueError, match="Bloque desconocido"):
        build_grid("principal", reps=1)      # bloque retirado


def test_sin_ruido_stab_esta_pareado_con_sin_ruido():
    """
    Los dos bloques generan las matrices IDENTICAS y solo cambia
    ``fmm_method``. Sin esto la comparacion seria entre celdas con ruido
    distinto y perderia toda la potencia del diseno pareado.
    """
    a = build_grid("sin_ruido", reps=3)
    b = build_grid("sin_ruido_stab", reps=3)
    assert len(a) == len(b)
    assert {c.fmm_method for c in a} == {"base"}
    assert {c.fmm_method for c in b} == {"stabilized"}

    def gen(celdas):
        return sorted((c.cfg.n_samples, c.cfg.sampling, c.cfg.omega,
                       c.cfg.sigma, c.cfg.seed, c.rep) for c in celdas)

    assert gen(a) == gen(b), (
        "las semillas de generacion no coinciden, la comparacion no seria pareada"
    )
