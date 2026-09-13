"""
Definicion de las celdas del estudio de simulacion.

Un bloque es una lista de celdas, y una celda es una configuracion de
generacion mas una configuracion de algoritmo mas un numero de replica. Las
semillas se derivan de forma determinista de la celda, asi que cualquier fila
del CSV se reproduce sola sin depender del orden de ejecucion.
"""
from __future__ import annotations

import hashlib
from dataclasses import dataclass, replace

from simulation.runner import RobustConfig
from simulation.truth import SimConfig

# La escala absoluta no la ve el pipeline (normalizacion min-max por gen), asi
# que el eje de senal-ruido se realiza fijando sigma y moviendo la amplitud.
SIGMA_BASE = 1.0

# Niveles medidos en la Task 6 sobre 27 celdas piloto. Por debajo de
# A/sigma = 1 el orden no se reconstruye (oraculo 3,3 a 3,9 h) y gastar
# replicas ahi solo produce tasa de fallo. La calidad del orden satura en 1,5
# (oraculo 0,44 h y ya no mejora) y la orientacion en 2, asi que la
# transicion interesante esta entre 1 y 2 y se muestrea con tres niveles.
# El 1,0 se conserva como celda de frontera de aplicabilidad.
SNR_NIVELES = (1.0, 1.5, 2.0, 3.0, 4.0)

N_NIVELES     = (24, 48, 96)
MUESTREOS     = ("uniforme", "agrupado", "hueco")
OMEGAS        = (0.1, 0.5, 0.9)
GCORE_NIVELES = (10, 20, 30)
PROP_NIVELES  = (0.01, 0.05, 0.20)

# Celda central, referencia del piloto y de los bloques con eje propio.
N_CENTRAL        = 48
MUESTREO_CENTRAL = "uniforme"
OMEGA_CENTRAL    = 0.5
GCORE_CENTRAL    = 20

# El bloque de agregacion usa 5000 genes y no 1000. Con 1000 y prop 0,01 solo
# habria 10 genes ritmicos, menos que los 30 core, y el generador no podria
# alojarlos. Con 5000 el nivel mas bajo da 50 ritmicos, que si caben.
GENES_AGREGACION = 5000

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

# Ruido del bloque de verificacion. No es 0 exacto a proposito: con sigma=0
# los genes de fondo quedan constantes y el preprocesado los elimina, con lo
# que se cambiarian dos cosas a la vez, el ruido y la composicion de la
# matriz. Con 1e-6 el ruido es despreciable y el fondo se conserva. Medido:
# las dos variantes dan los mismos resultados en las 27 celdas.
SIGMA_SIN_RUIDO = 1e-6

BLOQUES = (
    "sin_ruido", "piloto", "principal", "regularizacion", "agregacion",
    "gcore", "invariancia", "computacional",
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
    """Fija la semilla de la celda a partir de su identidad de generacion."""
    clave = (
        f"n{cfg.n_samples}|{cfg.sampling}|w{cfg.omega}|A{cfg.amplitude}"
        f"|s{cfg.sigma}|g{cfg.n_genes}|p{cfg.prop_rhythmic}"
    )
    return replace(cfg, seed=cell_seed(bloque, clave, rep))


# El bloque de regularizacion tiene que generar exactamente las mismas
# matrices que el principal, asi que ambos usan esta semilla de generacion.
SEMILLA_COMPARTIDA = "principal"


def _base_principal() -> list[SimConfig]:
    """Las 135 configuraciones de generacion de los bloques P y F."""
    salida = []
    for n in N_NIVELES:
        for muestreo in MUESTREOS:
            for omega in OMEGAS:
                for snr in SNR_NIVELES:
                    salida.append(_cfg(n, muestreo, omega, snr))
    return salida


def build_grid(bloque: str, reps: int) -> list[Cell]:
    """Todas las celdas de un bloque, con sus replicas."""
    por_defecto = RobustConfig()
    celdas: list[Cell] = []

    if bloque == "sin_ruido":
        # Verificacion de que el pipeline reconstruye el orden cuando los
        # datos son perfectos. Factorial completo de los tres factores del
        # bloque principal, sin el eje de senal-ruido, que aqui no aplica.
        # Cualquier error residual no puede venir del ruido, asi que aisla
        # los sesgos del propio algoritmo.
        for n in N_NIVELES:
            for muestreo in MUESTREOS:
                for omega in OMEGAS:
                    # amplitud 1 y sigma 1e-6, o sea A/sigma = 1e6. Ojo:
                    # _cfg calcula amplitude = snr * sigma, asi que el
                    # primer argumento es el cociente, no la amplitud.
                    cfg = _cfg(n, muestreo, omega, 1.0 / SIGMA_SIN_RUIDO,
                               sigma=SIGMA_SIN_RUIDO)
                    for rep in range(reps):
                        celdas.append(Cell(bloque, _con_semilla(bloque, cfg, rep),
                                           GCORE_CENTRAL, "base", por_defecto, rep))
        return celdas

    if bloque == "piloto":
        for snr in SNR_NIVELES:
            for rep in range(reps):
                cfg = _cfg(N_CENTRAL, MUESTREO_CENTRAL, OMEGA_CENTRAL, snr)
                celdas.append(Cell(bloque, _con_semilla(bloque, cfg, rep),
                                   GCORE_CENTRAL, "base", por_defecto, rep))
        return celdas

    if bloque in ("principal", "regularizacion"):
        fmm = "base" if bloque == "principal" else "stabilized"
        for cfg in _base_principal():
            for rep in range(reps):
                # semilla del bloque principal en los dos casos: la
                # comparacion es pareada sobre la matriz identica
                celdas.append(Cell(
                    bloque, _con_semilla(SEMILLA_COMPARTIDA, cfg, rep),
                    GCORE_CENTRAL, fmm, por_defecto, rep,
                ))
        return celdas

    if bloque == "agregacion":
        for prop in PROP_NIVELES:
            for snr in SNR_NIVELES:
                for n in N_NIVELES:
                    cfg = _cfg(n, MUESTREO_CENTRAL, OMEGA_CENTRAL, snr,
                               n_genes=GENES_AGREGACION, prop_rhythmic=prop)
                    for rep in range(reps):
                        sembrada = _con_semilla(bloque, cfg, rep)
                        for metodo in METODOS_RESUMEN:
                            celdas.append(Cell(bloque, sembrada, GCORE_CENTRAL,
                                               "base", metodo, rep))
        return celdas

    if bloque == "gcore":
        for gcore in GCORE_NIVELES:
            for snr in SNR_NIVELES:
                for n in N_NIVELES:
                    cfg = _cfg(n, MUESTREO_CENTRAL, OMEGA_CENTRAL, snr)
                    for rep in range(reps):
                        celdas.append(Cell(bloque, _con_semilla(bloque, cfg, rep),
                                           gcore, "base", por_defecto, rep))
        return celdas

    if bloque == "invariancia":
        # Tres pares con el mismo SNR y distinta escala absoluta. Si la
        # decision 4 de la spec es correcta, los dos miembros de cada par
        # tienen que dar metricas identicas.
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
