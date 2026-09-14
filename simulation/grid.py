"""
Definicion de las celdas del estudio de simulacion.

Un bloque es una lista de celdas, y una celda es una configuracion de
generacion mas una configuracion de algoritmo mas un numero de replica. Las
semillas se derivan de forma determinista de la celda, asi que cualquier fila
del CSV se reproduce sola sin depender del orden de ejecucion.

La rejilla unificada
--------------------
El bloque ``rejilla`` es el que evalua las tres aportaciones del TFG. Es un
factorial de los factores de generacion, y sobre CADA celda y replica se
corren las dos variantes del FMM (base y estabilizado) y los cinco metodos de
resumen de las K repeticiones, todos sobre la misma matriz. Asi:

- el contraste base contra estabilizado existe en todas las celdas, y ademas
  por metodo de resumen,
- el contraste entre metodos de resumen existe en todas las celdas, y ademas
  por variante del FMM,
- el tiempo de cada configuracion sale gratis, por etapa.

Los antiguos bloques ``principal``, ``regularizacion``, ``agregacion`` y
``gcore`` evaluaban cada mejora en un subconjunto de ejes distinto, con lo que
no se podia responder, por ejemplo, si la agregacion ayuda con muestreo
agrupado o si el estabilizado ayuda con pocos genes ritmicos. Esta rejilla
los sustituye. ``sin_ruido`` e ``invariancia`` se mantienen como verificacion,
y ``computacional`` es el bloque de R, que va aparte porque R es unas 20
veces mas lento y necesita un proceso por celda.
"""
from __future__ import annotations

import hashlib
from dataclasses import dataclass, replace

from simulation.runner import RobustConfig
from simulation.truth import SimConfig

# La escala absoluta no la ve el pipeline (normalizacion min-max por gen), asi
# que el eje de senal-ruido se realiza fijando sigma y moviendo la amplitud.
SIGMA_BASE = 1.0

# Niveles de A/sigma. Por debajo de 1 el orden no se reconstruye (oraculo 3,3
# a 3,9 h, medido) y gastar replicas ahi solo produce tasa de fallo. La calidad
# del orden satura en 1,5 y la orientacion en 2, asi que la transicion esta
# entre 1 y 2 y se muestrea con tres niveles. El 4 es el regimen saturado. El
# 3 se quito tras medir en la rejilla anterior que 3 y 4 son indistinguibles
# en todas las medidas (medae 1,24 contra 1,17 h, orientacion 1,0 en ambos).
SNR_NIVELES = (1.0, 1.5, 2.0, 4.0)

N_NIVELES     = (24, 48, 96)
MUESTREOS     = ("uniforme", "agrupado", "hueco")
OMEGAS        = (0.1, 0.5, 0.9)
GCORE_NIVELES = (10, 20, 30)
FMM_METODOS   = ("base", "stabilized")

# Tamano de la matriz y proporcion de ritmicos, como eje conjunto. Lo que el
# pipeline ve de verdad es cuantos genes ritmicos hay (el TOP satura en 100
# genes) y cuanto fondo. El factorial completo de los dos ejes tendria celdas
# redundantes, asi que se toman cinco combinaciones que cubren ritmicos en
# {50, 200, 250, 1000} y, dentro, dos pares con los mismos 50 ritmicos y
# distinto fondo (1000 x 0,05 y 5000 x 0,01), que aislan el efecto del numero
# de genes por si solo. (1000, 0,01) es imposible: 10 ritmicos no alojan los
# 30 genes core.
GENES_PROP = (
    (1000, 0.05),
    (1000, 0.20),
    (5000, 0.01),
    (5000, 0.05),
    (5000, 0.20),
)

# Celda central, referencia de los bloques de verificacion.
N_CENTRAL        = 48
MUESTREO_CENTRAL = "uniforme"
OMEGA_CENTRAL    = 0.5
GCORE_CENTRAL    = 20

# Marco estandar del bloque computacional, segun el parrafo de apertura del
# capitulo: tiempos equiespaciados, ruido medio, ondas simetricas y
# asimetricas, 5 % de ritmicos.
OMEGA_MEZCLA        = (0.1, 0.5, 0.9)
GENES_COMPUTACIONAL = (1000, 5000, 20000)

METODOS_RESUMEN: list[RobustConfig] = [
    RobustConfig(robust_method="best_k",    selection_criterion="r2_median"),
    RobustConfig(robust_method="best_k",    selection_criterion="medoid"),
    RobustConfig(robust_method="aggregate", aggregation="tsp3", consensus_scale="circmean"),
    RobustConfig(robust_method="aggregate", aggregation="hods", consensus_scale="circmean"),
    RobustConfig(robust_method="aggregate", aggregation="tsp3", consensus_scale="cire"),
]

# El metodo de referencia: es lo que hace CIRCUST tal cual esta en el codigo
# antes de esta aportacion, y contra el se miden los demas.
METODO_REFERENCIA = METODOS_RESUMEN[0]

# Ruido del bloque de verificacion. No es 0 exacto a proposito: con sigma=0
# los genes de fondo quedan constantes y el preprocesado los elimina, con lo
# que se cambiarian dos cosas a la vez, el ruido y la composicion de la
# matriz. Con 1e-6 el ruido es despreciable y el fondo se conserva. Medido:
# las dos variantes dan los mismos resultados en las 27 celdas.
SIGMA_SIN_RUIDO = 1e-6

BLOQUES = (
    "rejilla", "sin_ruido", "sin_ruido_stab", "invariancia", "computacional",
)


@dataclass(frozen=True)
class Cell:
    """Una ejecucion del estudio: que generar, como correrlo y que replica."""

    bloque:          str
    cfg:             SimConfig
    n_core_declared: int
    fmm_method:      str
    robust_cfg:      RobustConfig
    rep:             int

    @property
    def key(self) -> str:
        """Identidad de la celda, sin la replica ni la semilla."""
        return (
            f"{self.bloque}|n{self.cfg.n_samples}|{self.cfg.sampling}"
            f"|w{self.cfg.omega}|A{self.cfg.amplitude}|s{self.cfg.sigma}"
            f"|g{self.cfg.n_genes}|p{self.cfg.prop_rhythmic}"
            f"|c{self.n_core_declared}|{self.fmm_method}"
            f"|{self.robust_cfg.label}"
        )

    @property
    def row_id(self) -> str:
        """Identidad de la fila del CSV, para poder reanudar."""
        return f"{self.key}|r{self.rep}"


def cell_seed(bloque: str, cfg_key: str, rep: int) -> int:
    """
    Semilla determinista de una celda.

    Se deriva de un hash estable del bloque, la identidad de la celda y la
    replica. No depende del orden de ejecucion, asi que se puede reanudar o
    reordenar la rejilla sin cambiar ningun resultado.
    """
    crudo = f"{bloque}|{cfg_key}|{rep}".encode("utf-8")
    return int.from_bytes(hashlib.sha256(crudo).digest()[:4], "big")


def _cfg(
    n_samples: int,
    sampling:  str,
    omega,
    snr:       float,
    *,
    sigma:         float = SIGMA_BASE,
    n_genes:       int = 1000,
    prop_rhythmic: float = 0.20,
) -> SimConfig:
    """Construye una ``SimConfig`` desde la razon senal-ruido."""
    return SimConfig(
        n_samples     = n_samples,
        sampling      = sampling,
        omega         = omega,
        amplitude     = snr * sigma,
        sigma         = sigma,
        n_genes       = n_genes,
        prop_rhythmic = prop_rhythmic,
    )


def _con_semilla(bloque: str, cfg: SimConfig, rep: int) -> SimConfig:
    """
    Fija la semilla de la celda a partir de su identidad de GENERACION.

    La clave no lleva ni el FMM, ni el metodo de resumen, ni los core
    declarados, porque ninguno de los tres cambia la matriz. Es lo que hace
    que todas las variantes de una celda corran sobre datos identicos y las
    comparaciones sean pareadas.
    """
    clave = (
        f"n{cfg.n_samples}|{cfg.sampling}|w{cfg.omega}|A{cfg.amplitude}"
        f"|s{cfg.sigma}|g{cfg.n_genes}|p{cfg.prop_rhythmic}"
    )
    return replace(cfg, seed=cell_seed(bloque, clave, rep))


def escenarios_rejilla() -> list[SimConfig]:
    """Las configuraciones de generacion de la rejilla unificada, sin semilla."""
    salida = []
    for n in N_NIVELES:
        for muestreo in MUESTREOS:
            for omega in OMEGAS:
                for snr in SNR_NIVELES:
                    for n_genes, prop in GENES_PROP:
                        salida.append(_cfg(n, muestreo, omega, snr,
                                           n_genes=n_genes, prop_rhythmic=prop))
    return salida


def build_grid(bloque: str, reps: int) -> list[Cell]:
    """Todas las celdas de un bloque, con sus replicas."""
    por_defecto = METODO_REFERENCIA
    celdas: list[Cell] = []

    if bloque == "rejilla":
        # Las replicas van en el bucle EXTERIOR a proposito: si la ejecucion
        # se interrumpe a medias, lo que hay en disco es la rejilla entera con
        # menos replicas, que sigue siendo un diseno equilibrado y analizable.
        # Con las replicas dentro quedarian unas celdas completas y otras sin
        # empezar.
        escenarios = escenarios_rejilla()
        for rep in range(reps):
            for cfg in escenarios:
                sembrada = _con_semilla(bloque, cfg, rep)
                for gcore in GCORE_NIVELES:
                    for fmm in FMM_METODOS:
                        for metodo in METODOS_RESUMEN:
                            celdas.append(Cell(bloque, sembrada, gcore, fmm,
                                               metodo, rep))
        return celdas

    if bloque in ("sin_ruido", "sin_ruido_stab"):
        # Verificacion de que el pipeline reconstruye el orden cuando los
        # datos son perfectos. Factorial de los tres factores de forma, sin
        # el eje de senal-ruido, que aqui no aplica. Cualquier error residual
        # no puede venir del ruido, asi que aisla los sesgos del algoritmo.
        #
        # La variante _stab corre las MISMAS matrices con el FMM estabilizado
        # y aisla el sesgo propio de la regularizacion.
        fmm = "base" if bloque == "sin_ruido" else "stabilized"
        for n in N_NIVELES:
            for muestreo in MUESTREOS:
                for omega in OMEGAS:
                    # amplitud 1 y sigma 1e-6, o sea A/sigma = 1e6. Ojo:
                    # _cfg calcula amplitude = snr * sigma, asi que el
                    # primer argumento es el cociente, no la amplitud.
                    cfg = _cfg(n, muestreo, omega, 1.0 / SIGMA_SIN_RUIDO,
                               sigma=SIGMA_SIN_RUIDO)
                    for rep in range(reps):
                        celdas.append(Cell(
                            bloque, _con_semilla("sin_ruido", cfg, rep),
                            GCORE_CENTRAL, fmm, por_defecto, rep))
        return celdas

    if bloque == "invariancia":
        # Tres pares con el mismo SNR y distinta escala absoluta. Si el
        # pipeline solo ve A/sigma, los dos miembros de cada par tienen que
        # dar metricas identicas.
        for snr in (1.0, 2.0, 4.0):
            for sigma in (SIGMA_BASE, SIGMA_BASE * 4.0):
                cfg = _cfg(N_CENTRAL, MUESTREO_CENTRAL, OMEGA_CENTRAL, snr,
                           sigma=sigma)
                for rep in range(reps):
                    celdas.append(Cell(bloque, _con_semilla(bloque, cfg, rep),
                                       GCORE_CENTRAL, "base", por_defecto, rep))
        return celdas

    if bloque == "computacional":
        for n_genes in GENES_COMPUTACIONAL:
            for n in N_NIVELES:
                cfg = _cfg(n, "equiespaciado", OMEGA_MEZCLA, 2.0,
                           n_genes=n_genes, prop_rhythmic=0.05)
                for rep in range(reps):
                    celdas.append(Cell(bloque, _con_semilla(bloque, cfg, rep),
                                       GCORE_CENTRAL, "base", por_defecto, rep))
        return celdas

    raise ValueError(
        f"Bloque desconocido: {bloque!r}. Validos: {', '.join(BLOQUES)}"
    )
