"""
Generador de conjuntos de expresion sinteticos con tiempo circadiano conocido.

Produce una matriz genes x muestras y la ``SimulationTruth`` asociada. Los
tiempos verdaderos se usan para construir la senal y **no** se le entregan al
algoritmo.
"""
from __future__ import annotations

import math

import numpy as np
import pandas as pd

from simulation.truth import (
    TWO_PI, SimConfig, SimulationTruth,
    consistency_chain, core_layout, nested_core_sets,
)


def sample_times(cfg: SimConfig, rng: np.random.Generator) -> np.ndarray:
    """
    Tiempos circadianos verdaderos en ``[0, 24)``.

    Escenarios
    ----------
    uniforme
        ``t ~ U(0, 24)``, las muestras cubren el ciclo de forma homogenea.
    agrupado
        Mezcla a partes iguales de von Mises centradas en ``cluster_centers``.
        Se usa von Mises y no una normal truncada porque el grupo de las 20 h
        tiene masa que cruza las 24, y truncar ahi mete un artefacto justo en
        el punto donde el algoritmo tiene que cerrar el circulo.
    hueco
        Uniforme sobre ``[0, 24)`` menos la ventana ``gap``.
    equiespaciado
        Rejilla regular sin aleatoriedad. No es un escenario del diseno
        factorial, existe porque el marco estandar del bloque computacional
        lo pide.
    """
    n = cfg.n_samples

    if cfg.sampling == "uniforme":
        return rng.uniform(0.0, 24.0, n)

    if cfg.sampling == "equiespaciado":
        return np.linspace(0.0, 24.0, n, endpoint=False)

    if cfg.sampling == "agrupado":
        centros = np.asarray(cfg.cluster_centers, dtype=float)
        # kappa equivalente a una desviacion de cluster_sd horas
        kappa = (24.0 / (TWO_PI * cfg.cluster_sd)) ** 2
        cual  = rng.integers(0, len(centros), n)
        mu    = centros[cual] * TWO_PI / 24.0
        ang   = rng.vonmises(mu, kappa)          # en (-pi, pi]
        return (ang % TWO_PI) * 24.0 / TWO_PI

    if cfg.sampling == "hueco":
        lo, hi = cfg.gap
        ancho  = 24.0 - (hi - lo)
        u      = rng.uniform(0.0, ancho, n)
        return np.where(u < lo, u, u + (hi - lo))

    raise ValueError(
        f"Muestreo desconocido: {cfg.sampling!r}. "
        f"Validos: uniforme, agrupado, hueco, equiespaciado"
    )
