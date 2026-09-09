"""
Configuracion de una celda de simulacion y verdad subyacente.

Este modulo no genera datos. Solo define que se va a generar: los parametros
de la celda, la posicion circadiana de los genes core, los conjuntos anidados
que se le declaran al pipeline y la cadena de conocimiento a priori que
sustituye a ``DBP < CRY1 < ARNTL`` en los datos sinteticos.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
import pandas as pd

TWO_PI = 2.0 * math.pi

# Gen de anclaje, su pico se fija en pi y define el origen del ciclo.
ANCHOR_NAME = "CORE_A"
# Gen contrapuesto, su pico se fija en 0. Sirve para comprobar la oposicion,
# no para decidir la orientacion: a exactamente pi del ancla es el punto de
# indiferencia y su semicirculo es indistinguible bajo reflexion.
OPPOSITE_NAME = "CORE_B"


@dataclass(frozen=True)
class SimConfig:
    """Parametros de una celda del estudio de simulacion."""

    n_samples:       int
    sampling:        str = "uniforme"    # uniforme|agrupado|hueco|equiespaciado
    omega:           float | tuple[float, ...] = 0.5
    amplitude:       float = 1.0
    sigma:           float = 1.0
    n_genes:         int = 1000
    prop_rhythmic:   float = 0.20
    n_core:          int = 30
    core_sizes:      tuple[int, ...] = (10, 20, 30)
    day_fraction:    float = 0.60
    n_chain:         int = 3
    cluster_centers: tuple[float, ...] = (4.0, 12.0, 20.0)
    cluster_sd:      float = 2.0
    gap:             tuple[float, float] = (10.0, 16.0)
    seed:            int = 0

    @property
    def snr(self) -> float:
        """Razon senal-ruido. Es lo unico que ve el pipeline, porque la
        normalizacion min-max por gen es invariante a la escala absoluta."""
        return float(self.amplitude) / float(self.sigma)

    @property
    def n_rhythmic(self) -> int:
        return int(round(self.n_genes * self.prop_rhythmic))


def _spread_order(total: int) -> list[int]:
    """
    Permutacion de ``range(total)`` tal que cualquier prefijo queda repartido
    por todo el rango, por bisecion recursiva de intervalos.

    Es lo que hace que los conjuntos core anidados esten anidados **y** cubran
    el semicirculo. Tomar simplemente los primeros k amontonaria G_10 en un
    arco.
    """
    order: list[int] = []
    segments = [(0, total - 1)]
    while segments:
        siguientes: list[tuple[int, int]] = []
        for lo, hi in segments:
            if lo > hi:
                continue
            mid = (lo + hi) // 2
            order.append(mid)
            siguientes.append((lo, mid - 1))
            siguientes.append((mid + 1, hi))
        segments = siguientes
    return order


def _spread_indices(total: int, k: int) -> list[int]:
    """``k`` indices de ``range(total)`` a intervalos regulares, distintos."""
    if k <= 0:
        return []
    if k >= total:
        return list(range(total))
    return [int((2 * i + 1) * total // (2 * k)) for i in range(k)]


def core_layout(
    n_core:       int,
    day_fraction: float,
) -> tuple[list[str], np.ndarray]:
    """
    Nombres y posicion del pico de los ``n_core`` genes core.

    ``CORE_A`` va en ``pi`` y cuenta en el lado noche, ``CORE_B`` va en ``0``
    y cuenta en el lado dia. Los demas se reparten segun ``day_fraction`` y se
    colocan equiespaciados dentro de su semicirculo abierto.
    """
    if n_core < 2:
        raise ValueError(f"n_core debe ser al menos 2, se recibio {n_core}")

    n_day   = int(round(day_fraction * n_core))
    n_night = n_core - n_day
    if n_day < 1 or n_night < 1:
        raise ValueError(
            f"day_fraction={day_fraction} con n_core={n_core} deja un "
            f"semicirculo vacio (dia={n_day}, noche={n_night})"
        )

    n_day_extra   = n_day - 1
    n_night_extra = n_night - 1

    day_peaks = math.pi * np.arange(1, n_day_extra + 1) / (n_day_extra + 1)
    night_peaks = math.pi + math.pi * np.arange(1, n_night_extra + 1) / (
        n_night_extra + 1
    )

    names = [ANCHOR_NAME, OPPOSITE_NAME]
    peaks = [math.pi, 0.0]
    names += [f"CORE_D{i + 1:02d}" for i in range(n_day_extra)]
    peaks += list(day_peaks)
    names += [f"CORE_N{i + 1:02d}" for i in range(n_night_extra)]
    peaks += list(night_peaks)

    return names, np.asarray(peaks, dtype=float)


def nested_core_sets(
    names:        list[str],
    core_sizes:   tuple[int, ...],
    day_fraction: float,
) -> dict[int, list[str]]:
    """
    Conjuntos core anidados ``G_k``, uno por cada tamano de ``core_sizes``.

    Cada nivel conserva las dos anclas y respeta ``day_fraction``. El
    anidamiento sale de tomar prefijos de la misma ``_spread_order``.
    """
    day_extra   = [n for n in names if n.startswith("CORE_D")]
    night_extra = [n for n in names if n.startswith("CORE_N")]
    orden_dia   = _spread_order(len(day_extra))
    orden_noche = _spread_order(len(night_extra))

    sets: dict[int, list[str]] = {}
    for k in core_sizes:
        if k > len(names):
            raise ValueError(f"core_size {k} supera n_core={len(names)}")
        n_day = int(round(day_fraction * k))
        n_night = k - n_day
        if n_day < 1 or n_night < 1:
            raise ValueError(f"core_size {k} deja un semicirculo vacio")
        elegidos = [ANCHOR_NAME, OPPOSITE_NAME]
        elegidos += [day_extra[i] for i in orden_dia[: n_day - 1]]
        elegidos += [night_extra[i] for i in orden_noche[: n_night - 1]]
        sets[k] = elegidos
    return sets


def consistency_chain(
    names:     list[str],
    peaks:     np.ndarray,
    core_sets: dict[int, list[str]],
    n_chain:   int,
) -> list[str]:
    """
    Cadena de conocimiento a priori ``CORE_B < ... < CORE_A``.

    Los genes intermedios salen del conjunto core mas pequeno, de modo que la
    informacion a priori este disponible en los tres niveles de ``G_core`` y
    no crezca con el. Con ``n_chain=3`` produce los dos pares que replican
    ``DBP < CRY1 < ARNTL``.
    """
    if n_chain < 2:
        raise ValueError(f"n_chain debe ser al menos 2, se recibio {n_chain}")

    peak_of  = dict(zip(names, np.asarray(peaks, dtype=float)))
    menor    = core_sets[min(core_sets)]
    dia      = sorted(
        [g for g in menor if 0.0 < peak_of[g] < math.pi],
        key=lambda g: peak_of[g],
    )
    n_medios = n_chain - 2
    if n_medios > len(dia):
        raise ValueError(
            f"n_chain={n_chain} necesita {n_medios} genes intermedios pero el "
            f"conjunto core mas pequeno solo tiene {len(dia)} en (0, pi)"
        )
    medios = [dia[i] for i in _spread_indices(len(dia), n_medios)]
    return [OPPOSITE_NAME] + medios + [ANCHOR_NAME]


def consistency_pairs_from_chain(chain: list[str]) -> list[tuple[str, str]]:
    """Pares consecutivos ``(g1, g2)`` con ``peak(g1) < peak(g2)`` esperado."""
    return [(chain[i], chain[i + 1]) for i in range(len(chain) - 1)]


@dataclass(frozen=True)
class SimulationTruth:
    """Todo lo que se conoce del conjunto generado y no se le da al algoritmo."""

    cfg:          SimConfig
    t_true:       np.ndarray          # horas en [0, 24)
    theta_true:   np.ndarray          # radianes en [0, 2*pi)
    sample_names: list[str]
    core_names:   list[str]
    core_peaks:   np.ndarray
    core_sets:    dict[int, list[str]]
    chain:        list[str]
    gene_params:  pd.DataFrame        # gene, kind, M, A, alpha, beta, omega

    @property
    def anchor(self) -> str:
        return ANCHOR_NAME

    @property
    def direction(self) -> str:
        return OPPOSITE_NAME

    def core_genes(self, size: int) -> list[str]:
        """Nombres del conjunto core de tamano ``size`` que se declara."""
        if size not in self.core_sets:
            raise KeyError(
                f"No hay conjunto core de tamano {size}. "
                f"Disponibles: {sorted(self.core_sets)}"
            )
        return list(self.core_sets[size])

    def consistency_pairs(self) -> list[tuple[str, str]]:
        return consistency_pairs_from_chain(self.chain)

    def theta_for(self, sample_names: list[str]) -> np.ndarray:
        """
        Angulos verdaderos de las muestras nombradas, en ese orden.

        El pipeline elimina outliers y reordena columnas, asi que la
        comparacion contra la verdad se hace por nombre y nunca por posicion.
        """
        pos = {nm: i for i, nm in enumerate(self.sample_names)}
        return np.asarray([self.theta_true[pos[nm]] for nm in sample_names])
