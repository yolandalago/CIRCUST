#!/usr/bin/env python3
"""
simulation/grid.py
==================
Definicion de la rejilla del estudio de simulacion: que celdas hay, con que
parametros, cuantas replicas y con que semillas.

Este modulo es solo datos y funciones puras. No genera nada ni evalua nada, de
modo que puede leerse como la especificacion del experimento y comprobarse sin
esperar horas.

Por que bloques y no factorial completo
---------------------------------------
Cruzar todos los factores daria miles de celdas y la mayoria no responderia a
ninguna pregunta. En su lugar hay un bloque por cada una de las tres mejoras
que este trabajo introduce, y cada bloque mueve solo los factores que afectan a
la suya, congelando el resto en la celda de referencia.

    P   eficiencia computacional     mueve el tamano del problema
    R   regularizacion del FMM       mueve muestreo y forma de onda
    A   orden representativo         mueve dificultad y cantidad de senal

Antes de los tres va una **fase 0 de calibracion**, que no es opcional. Una
celda con SNR 3 satura: todos los metodos dan rho 0,997 y ninguna aportacion
puede distinguirse, mientras que sobre babuino real el rho es 0,703. Sin
calibrar, el estudio saldria plano por construccion.

Calibracion y evaluacion van separadas
--------------------------------------
La fase 0 y el barrido de tau fijan parametros. Los bloques R y A los evaluan.
Usan bases de semillas distintas a proposito, de modo que ningun parametro se
elige y se juzga sobre los mismos datos.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from simulation.generate import SyntheticDataset, generate

# ===========================================================================
# Celda de referencia
# ===========================================================================
#
# Los tamanos estan anclados en los datos reales del propio repositorio:
# ``data/matrixIn.parquet``, el tejido de GTEx sobre el que corre CIRCUST, mide
# 56.200 genes por 480 muestras, y los tejidos de babuino de Mure et al. tienen
# 12 muestras.
#
# ``n_genes = 20.000`` es el orden de los genes efectivamente expresados. Bajar
# de ahi degenera las K repeticiones de la etapa robusta, que pasan a producir
# todas el mismo orden, y entonces el bloque A no mediria nada por falta de
# desacuerdo que resolver.
#
# El muestreo de referencia es ``aleatorio`` y no ``equiespaciado`` por dos
# razones. Es la condicion no equiespaciada del material suplementario de
# CIRCUST, de modo que la referencia coincide con su escenario, y no es
# degenerada, a diferencia del equiespaciado, que reparte las muestras con una
# regularidad que ningun estudio post mortem tiene.
#
# ``snr`` y ``n_samples`` son PROVISIONALES. Los fija la fase 0. Los valores de
# aqui son los del pilotaje del 14/08/2026, que sobre una celda de prueba dio
# rho entre 0,64 y 0,69 y separo a los seis metodos, o sea justo la ventana
# buscada.
REFERENCIA: dict = {
    "n_genes":         20_000,
    "n_samples":       48,
    "snr":             1.0,
    "prop_rhythmic":   0.05,
    "sampling":        "aleatorio",
    "waveform_mix":    "gtex",
    "beta_mode":       "aleatorio",
    "omega_jitter":    0.05,
    "core_phase_sd_h": 1.0,
}

# Rejilla de SNR de la fase 0. Recuerda que la desviacion tipica del ruido de
# un gen es su amplitud dividida por el SNR, luego valores mayores son
# problemas mas faciles. El suplemento de CIRCUST trabaja con el equivalente a
# SNR entre 4,5 y 6,3, todavia mas facil que el extremo derecho de esta
# rejilla.
GRID_SNR: tuple[float, ...] = (0.25, 0.5, 1.0, 2.0, 3.0)

# Ventana util de tau, medida sobre matrixIn_mini. Por debajo de 0,60 el
# criterio satura porque coincide con el umbral de ritmicidad del propio TOP, y
# por encima de 0,80 colapsa a cero genes.
GRID_TAU: tuple[float, ...] = (0.60, 0.65, 0.70, 0.75, 0.80)

TODOS_LOS_METODOS: tuple[str, ...] = (
    "median", "best_k:r2_median", "best_k:n_well_fitted", "best_k:medoid",
    "aggregate:tsp3", "aggregate:hods",
)

REPLICAS: dict[str, int] = {"fase0": 5, "tau": 5, "P": 3, "R": 20, "A": 20}

# Bases disjuntas por bloque. La separacion entre calibracion (fase0, tau) y
# evaluacion (R, A) no es cosmetica: garantiza que ningun parametro se elige y
# se juzga sobre los mismos datos.
SEED_BASE: dict[str, int] = {
    "fase0": 1_000_000,
    "tau":   1_500_000,
    "P":     2_000_000,
    "R":     3_000_000,
    "A":     4_000_000,
}

# Separacion entre celdas dentro de un bloque. Acota las replicas por celda,
# que es de sobra para las 20 previstas.
PASO_CELDA: int = 100


@dataclass(frozen=True)
class Celda:
    """Una condicion experimental de la rejilla.

    Attributes
    ----------
    bloque
        ``fase0``, ``tau``, ``P``, ``R`` o ``A``.
    nombre
        Identificador legible, valido como nombre de directorio.
    generate
        Argumentos de :func:`simulation.generate.generate`.
    fmm_methods
        Modelos de ajuste a evaluar sobre esta celda.
    metodos
        Metodos de resumen de las K repeticiones, por etiqueta o como pares
        ``(etiqueta, kwargs)``.
    robust_kwargs
        Ajustes extra del estimador robusto comunes a toda la celda.
    con_referencia
        Calcular ademas el error de picos contra la referencia estimada. En una
        celda sintetica los parametros verdaderos se conocen, de modo que solo
        interesa donde se quiera comparar de tu a tu con los datos reales.
    conservar_matriz, tambien_csv
        Solo el bloque P las necesita, porque el CIRCUST original en R tiene
        que leer exactamente los mismos datos.
    seed_base
        Semilla de la replica cero. La replica ``r`` usa ``seed_base + r``.
    """
    bloque:            str
    nombre:            str
    generate:          dict
    fmm_methods:       tuple = ("base",)
    metodos:           tuple = TODOS_LOS_METODOS
    robust_kwargs:     dict  = field(default_factory=dict)
    con_referencia:    bool  = False
    conservar_matriz:  bool  = False
    tambien_csv:       bool  = False
    seed_base:         int   = 0

    @property
    def replicas(self) -> int:
        return REPLICAS[self.bloque]


# ===========================================================================
# Utilidades
# ===========================================================================

def _num(x: float) -> str:
    """Formatea un numero para un nombre de directorio, sin punto flotante feo."""
    s = f"{x:g}"
    return s.replace("-", "m")


def _celdas(bloque: str, definiciones: list[tuple[str, dict]],
            **comunes) -> list[Celda]:
    """Construye las celdas de un bloque asignando semillas correlativas."""
    base = SEED_BASE[bloque]
    return [
        Celda(bloque=bloque, nombre=nombre, generate=params,
              seed_base=base + i * PASO_CELDA, **comunes)
        for i, (nombre, params) in enumerate(definiciones)
    ]


def vecinos_snr(snr: float) -> tuple[float, float, float]:
    """El SNR dado y sus dos vecinos en :data:`GRID_SNR`.

    En los extremos de la rejilla se desplaza hacia dentro en vez de devolver
    dos valores, de modo que el bloque A siempre tiene tres niveles.
    """
    if snr not in GRID_SNR:
        raise ValueError(
            f"snr {snr} no esta en la rejilla de calibracion {GRID_SNR}.")
    i = GRID_SNR.index(snr)
    i = min(max(i, 1), len(GRID_SNR) - 2)
    return GRID_SNR[i - 1:i + 2]


def ruta_de(celda: Celda, replica: int, raiz: Path | str = Path("output")
            ) -> Path:
    """Directorio de una replica concreta."""
    return Path(raiz) / celda.bloque / celda.nombre / f"rep{replica:02d}"


def materializar(celda: Celda, replica: int, **sobrescribir
                 ) -> SyntheticDataset:
    """Genera la replica pedida de una celda.

    ``sobrescribir`` existe para los tests, que necesitan matrices pequenas sin
    tener que duplicar la definicion de la rejilla.
    """
    params = {**celda.generate, **sobrescribir}
    nombre = f"{celda.bloque}/{celda.nombre}/rep{replica:02d}"
    return generate(nombre, seed=celda.seed_base + replica, **params)


# ===========================================================================
# Fase 0: calibracion de la celda de referencia
# ===========================================================================

def bloque_fase0() -> list[Celda]:
    """SNR por numero de muestras, para situar la celda de referencia.

    El eje de muestras recorre desde el limite inferior real, los 12 de un
    tejido de babuino, hasta un tejido mediano de GTEx.

    Se evaluan solo dos configuraciones, la que representa el resumen robusto
    habitual y la que representa la agregacion, porque lo que se busca aqui no
    es comparar metodos sino localizar el regimen. El criterio de eleccion esta
    fijado de antemano: rho en torno a 0,8, MedAE del orden entre una y dos
    horas, y sd de rho entre las K repeticiones por encima de 0,03.
    """
    definiciones = [
        (f"snr{_num(snr)}__n{n}", {**REFERENCIA, "snr": snr, "n_samples": n})
        for snr in GRID_SNR
        for n in (12, 48, 192)
    ]
    return _celdas("fase0", definiciones,
                   fmm_methods=("base",),
                   metodos=("best_k:r2_median", "aggregate:tsp3"))


def bloque_tau(referencia: dict | None = None) -> list[Celda]:
    """Barrido del umbral tau sobre la referencia y los extremos de SNR.

    ``tau`` es ``selection_r2_threshold``, y solo interviene en el criterio
    ``best_k:n_well_fitted``, una de las seis combinaciones. Los bloques P y R
    son insensibles a el.

    Se calibra aqui, con semillas de calibracion, y no dentro del bloque A,
    porque elegir el valor y reportarlo en la misma tabla que lo evalua seria
    sobreajuste.

    Aviso sobre el desenlace. El criterio desempata por mediana de R2 cuando
    varias repeticiones cuentan los mismos genes, de modo que con tau demasiado
    alto colapsa exactamente en ``best_k:r2_median``. Si resulta que solo se
    distingue de el en una ventana estrecha, la conclusion legitima no es un
    valor de tau sino que el criterio no aporta y procede retirarlo.
    """
    ref = referencia or REFERENCIA
    metodos = tuple(
        (f"best_k:n_well_fitted@{_num(tau)}",
         {"method": "best_k", "selection_criterion": "n_well_fitted",
          "selection_r2_threshold": tau})
        for tau in GRID_TAU
    )
    definiciones = [
        (f"snr{_num(snr)}", {**ref, "snr": snr})
        for snr in vecinos_snr(ref["snr"])
    ]
    return _celdas("tau", definiciones,
                   fmm_methods=("base",), metodos=metodos)


# ===========================================================================
# Bloque P: eficiencia computacional
# ===========================================================================

def bloque_p() -> list[Celda]:
    """Tamano del problema, que es de lo unico que depende el tiempo.

    No se cruzan factores de dificultad porque no afectan al coste. Los dos
    ejes recorren desde un tejido de babuino hasta la escala real de GTEx,
    56.200 genes por 480 muestras, que es la dimension exacta de
    ``data/matrixIn.parquet``.

    Se comparan la reimplementacion configurada como el original
    (``base`` mas ``median``) y la reimplementacion con las mejoras activas,
    sobre los mismos datos. El CIRCUST original en R corre aparte sobre las
    matrices que estas celdas dejan en disco, y de ahi que sean las unicas que
    las conservan y que escriban tambien CSV.

    Aviso de coste. La celda mayor son unos tres minutos por configuracion en
    Python, pero con el factor 22 ya medido para R se va a unas tres horas por
    replica.
    """
    definiciones = [
        (f"g{g}__n{n}", {**REFERENCIA, "n_genes": g, "n_samples": n})
        for g in (5_000, 20_000, 56_200)
        for n in (12, 48, 192, 480)
    ]
    return _celdas("P", definiciones,
                   fmm_methods=("base", "stabilized"),
                   metodos=("median", "aggregate:tsp3"),
                   conservar_matriz=True, tambien_csv=True)


# ===========================================================================
# Bloque R: regularizacion del FMM
# ===========================================================================

def bloque_r(referencia: dict | None = None) -> list[Celda]:
    """Muestreo por forma de onda, con el resumen robusto congelado.

    Lo unico que cambia entre las dos configuraciones de cada celda es el
    modelo de ajuste, ``base`` frente a ``stabilized``, de modo que cualquier
    diferencia es atribuible a la regularizacion y no al resumen de las K
    repeticiones.

    Hipotesis contrastable: la regularizacion mejora en los esquemas de huecos
    y concentrado y con formas apuntadas, y no empeora ni en el equiespaciado
    ni en la mezcla simetrica.

    Evidencia preliminar a favor, medida sobre los atlas ya calculados: en
    babuino, con doce muestras, el FMM base situa el 53 por ciento de los genes
    en el extremo apuntado, y el estabilizado baja esa fraccion al 8 por ciento
    acercandola a la forma de GTEx. Queda por ver, y para eso esta este bloque,
    si al conocer la verdad el estabilizado acierta o se pasa de suavizado.
    """
    ref = referencia or REFERENCIA
    definiciones = [
        (f"samp-{s}__mix-{m}", {**ref, "sampling": s, "waveform_mix": m})
        for s in ("equiespaciado", "aleatorio", "huecos", "concentrado")
        for m in ("simetrica", "gtex", "apuntada", "extrema")
    ]
    return _celdas("R", definiciones,
                   fmm_methods=("base", "stabilized"),
                   metodos=("best_k:r2_median",))


# ===========================================================================
# Bloque A: orden representativo
# ===========================================================================

def bloque_a(referencia: dict | None = None) -> list[Celda]:
    """Dificultad y cantidad de senal, con el modelo FMM congelado en base.

    Los tres factores estan elegidos porque son los que gobiernan cuanto
    discrepan entre si las K repeticiones, que es la variable que decide si
    agregar sirve de algo. Bajar el SNR, reducir el numero de muestras y
    reducir la proporcion de ritmicos aumentan ese desacuerdo.

    El experimento sobre babuino ya establecio que la agregacion no mejora
    cuando las repeticiones coinciden, con un mecanismo identificado: si no hay
    desacuerdo no hay nada que resolver. Por eso el analisis de este bloque
    estratifica por ``sd_rho_reps`` en lugar de promediar sobre todas las
    celdas, que es lo que enmascararia el efecto.
    """
    ref = referencia or REFERENCIA
    definiciones = [
        (f"snr{_num(snr)}__n{n}__prop{_num(p)}",
         {**ref, "snr": snr, "n_samples": n, "prop_rhythmic": p})
        for snr in vecinos_snr(ref["snr"])
        for n in (12, 48, 192)
        for p in (0.01, 0.05, 0.15)
    ]
    return _celdas("A", definiciones,
                   fmm_methods=("base",), metodos=TODOS_LOS_METODOS)


# ===========================================================================
# La rejilla completa
# ===========================================================================

BLOQUES = {
    "fase0": bloque_fase0,
    "tau":   bloque_tau,
    "P":     bloque_p,
    "R":     bloque_r,
    "A":     bloque_a,
}


def todas_las_celdas(referencia: dict | None = None) -> list[Celda]:
    """Todas las celdas de la rejilla, en orden de ejecucion."""
    ref = referencia or REFERENCIA
    celdas: list[Celda] = []
    for nombre, constructor in BLOQUES.items():
        celdas += (constructor() if nombre in ("fase0", "P")
                   else constructor(ref))
    return celdas
