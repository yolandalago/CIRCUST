"""
Tests del criterio de orientacion del Paso 2.2 de CircularSynchronizer.

El criterio antiguo comparaba los picos con ``peak(g1) < peak(g2)``, que es
una comparacion lineal sobre datos circulares y por tanto depende de donde
este el corte. Con dos pares eso podia empatar y dejar pasar una inversion.
El criterio nuevo suma ``sin(peak(g2) - peak(g1))``, que es antisimetrico
bajo la reflexion y se autopondera por la separacion del par.
"""
import numpy as np
import pandas as pd
import pytest

from circust.synchronizer import CircularSynchronizer

TWO_PI = 2.0 * np.pi

ANCLA, OPUESTO, MEDIO = "A", "B", "M"
CORE = [ANCLA, OPUESTO, MEDIO]
PARES = [(OPUESTO, MEDIO), (MEDIO, ANCLA)]

# Picos verdaderos: B en 0, M en pi/3, A en pi. Es la disposicion de la
# simulacion, y los dos extremos caen en los puntos fijos de la reflexion.
PICOS_OK = np.array([np.pi, 0.0, np.pi / 3.0])
PICOS_REFLEJADOS = (TWO_PI - PICOS_OK) % TWO_PI


def _sincronizador(pares=None, r2_min=0.5):
    return CircularSynchronizer(
        anchor_gene       = ANCLA,
        direction_gene    = OPUESTO,
        consistency_pairs = PARES if pares is None else pares,
        consistency_r2_min= r2_min,
        verbose           = False,
    )


def _llamar(sync, picos, r2=None, n=8):
    """Invoca _basic_order con datos planos y devuelve si invirtio."""
    o   = np.arange(n)
    esc = np.linspace(0.0, TWO_PI, n, endpoint=False)
    mat = pd.DataFrame(
        np.random.default_rng(0).normal(size=(len(CORE), n)),
        index=CORE, columns=[f"S{i}" for i in range(n)],
    )
    par = np.zeros((len(CORE), 5))
    par[:, 2] = picos
    r2_fmm = np.full(len(CORE), 0.9) if r2 is None else np.asarray(r2, float)
    *_, flip = sync._basic_order(o, esc, mat, np.asarray(picos, float),
                                    par, r2_fmm, CORE)
    return flip


def test_no_invierte_cuando_la_orientacion_es_correcta():
    assert _llamar(_sincronizador(), PICOS_OK) is False


def test_invierte_cuando_la_orientacion_esta_reflejada():
    """
    La propiedad que fallaba. Con la regla lineal el voto salia 1 a 1 y el
    empate no invertia, asi que la reflexion sobrevivia.
    """
    assert _llamar(_sincronizador(), PICOS_REFLEJADOS) is True


def test_la_regla_lineal_habria_empatado_con_estos_mismos_picos():
    """Deja constancia de por que habia que cambiarla."""
    p = PICOS_REFLEJADOS
    idx = {g: i for i, g in enumerate(CORE)}
    a_favor = sum(p[idx[g1]] < p[idx[g2]] for g1, g2 in PARES)
    assert (a_favor, len(PARES) - a_favor) == (1, 1)


def test_un_par_antipodal_no_aporta_evidencia():
    """
    Dos genes separados 12 h no distinguen la direccion, porque en un ciclo
    de 24 h "12 horas despues" y "12 horas antes" son la misma posicion. Con
    la regla nueva aportan sin(pi) = 0 en vez de un voto de pleno derecho.
    """
    sync = _sincronizador(pares=[(OPUESTO, ANCLA)])   # B en 0, A en pi
    # ni con la orientacion correcta ni con la reflejada hay evidencia, asi
    # que se mantiene lo que decidio el Paso 2.1, que es no invertir
    assert _llamar(sync, PICOS_OK) is False
    assert _llamar(sync, PICOS_REFLEJADOS) is False


def test_los_pares_mal_ajustados_no_cuentan():
    """Un par con R2 por debajo del umbral se descarta, como antes."""
    r2_malo = [0.9, 0.9, 0.1]        # MEDIO mal ajustado, invalida los dos pares
    assert _llamar(_sincronizador(), PICOS_REFLEJADOS, r2=r2_malo) is False


def test_los_pares_con_genes_ausentes_no_cuentan():
    sync = _sincronizador(pares=[("NO_EXISTE", MEDIO), (MEDIO, "TAMPOCO")])
    assert _llamar(sync, PICOS_REFLEJADOS) is False


def test_una_cadena_larga_sigue_detectando_la_reflexion():
    """El criterio no depende de que la cadena tenga exactamente dos pares."""
    core  = ["A", "B", "M1", "M2", "M3"]
    picos = np.array([np.pi, 0.0, 0.4, 0.9, 1.8])
    pares = [("B", "M1"), ("M1", "M2"), ("M2", "M3"), ("M3", "A")]
    sync = CircularSynchronizer(anchor_gene="A", direction_gene="B",
                                consistency_pairs=pares, verbose=False)
    n = 8
    o, esc = np.arange(n), np.linspace(0.0, TWO_PI, n, endpoint=False)
    mat = pd.DataFrame(np.random.default_rng(0).normal(size=(len(core), n)),
                       index=core, columns=[f"S{i}" for i in range(n)])
    par = np.zeros((len(core), 5))
    r2  = np.full(len(core), 0.9)

    *_, flip_ok = sync._basic_order(o, esc, mat, picos, par, r2, core)
    reflejados = (TWO_PI - picos) % TWO_PI
    *_, flip_refl = sync._basic_order(o, esc, mat, reflejados, par, r2, core)
    assert flip_ok is False
    assert flip_refl is True


def test_el_signo_del_seno_es_antisimetrico_bajo_reflexion():
    """
    La garantia matematica del criterio, comprobada numericamente:
    sin((2pi - p2) - (2pi - p1)) = -sin(p2 - p1) para todo par.
    """
    rng = np.random.default_rng(3)
    p1, p2 = rng.uniform(0, TWO_PI, 500), rng.uniform(0, TWO_PI, 500)
    directo   = np.sin(p2 - p1)
    reflejado = np.sin(((TWO_PI - p2) % TWO_PI) - ((TWO_PI - p1) % TWO_PI))
    assert np.allclose(directo, -reflejado, atol=1e-12)


# ---------------------------------------------------------------------------
# El conocimiento a priori se consulta siempre, no solo si salta el aviso
# ---------------------------------------------------------------------------

# Gen de direccion lejos de 0 y de pi: 1-cos(0.75)=0.27 > 0.1, asi que la
# guarda ``aviso`` de R NO se dispara y con el codigo anterior el Paso 2.2 no
# llegaba a ejecutarse. Es la situacion medida en las celdas sinteticas que
# salian invertidas, donde dir_peak se estimaba entre 0,61 y 0,76 rad.
PICOS_SIN_AVISO      = np.array([np.pi, 0.75, 2.00])
PICOS_SIN_AVISO_REFL = (TWO_PI - PICOS_SIN_AVISO) % TWO_PI


def _aviso_saltaria(picos):
    """Reproduce la guarda de R para comprobar que efectivamente no salta."""
    dp = picos[CORE.index(OPUESTO)]
    p6am = 1.0 - np.cos(dp)
    p6pm = 1.0 - np.cos(dp - np.pi)
    peak_d = sum(1 for i, g in enumerate(CORE)
                 if g not in (ANCLA, OPUESTO) and 0 <= picos[i] < np.pi)
    mitad = int(np.floor((len(CORE) - 2) / 2))
    return (p6am <= 0.1) or (p6pm <= 0.1) or (peak_d < mitad)


def test_el_caso_de_prueba_efectivamente_no_dispara_el_aviso():
    """Sin esto el test siguiente no probaria lo que dice probar."""
    assert not _aviso_saltaria(PICOS_SIN_AVISO)
    assert not _aviso_saltaria(PICOS_SIN_AVISO_REFL)


def test_se_consultan_los_pares_aunque_no_salte_el_aviso():
    """
    La guarda de R solo consultaba el conocimiento a priori cuando el gen de
    direccion caia a menos de 0,45 rad de 0 o de pi. Fuera de esa ventana la
    decision del Paso 2.1 quedaba sin verificar, y esa decision se toma con
    un unico gen que, al estar en el antipoda del ancla, es invariante bajo
    reflexion y no aporta informacion.

    Teniendo conocimiento a priori disponible hay que usarlo siempre.
    """
    assert _llamar(_sincronizador(), PICOS_SIN_AVISO) is False
    assert _llamar(_sincronizador(), PICOS_SIN_AVISO_REFL) is True
