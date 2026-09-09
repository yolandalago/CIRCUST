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


def fmm_wave(
    theta: np.ndarray,
    M:     float,
    A:     float,
    alpha: float,
    beta:  float,
    omega: float,
) -> np.ndarray:
    """
    Onda FMM monocomponente evaluada en ``theta``.

        f(theta) = M + A*cos( beta + 2*arctan( omega*tan((theta - alpha)/2) ) )
    """
    mob = 2.0 * np.arctan(omega * np.tan((np.asarray(theta) - alpha) / 2.0))
    return M + A * np.cos(beta + mob)


def alpha_for_peak(theta_peak: float, beta: float, omega: float) -> float:
    """
    ``alpha`` tal que el maximo de la onda caiga en ``theta_peak``.

    ``FMMModel.peak_time`` es lineal en ``alpha``:

        t_U = alpha + 2*atan2( sin(-beta/2)/omega , cos(-beta/2) )   mod 2*pi

    asi que la inversion es exacta y no necesita optimizacion numerica.
    """
    desplazamiento = 2.0 * np.arctan2(
        np.sin(-beta / 2.0) / omega,
        np.cos(-beta / 2.0),
    )
    return float((theta_peak - desplazamiento) % TWO_PI)


def assign_omega(omega: float | tuple[float, ...], n: int) -> np.ndarray:
    """
    Vector de ``n`` valores de omega.

    Con un escalar todos los genes comparten omega, que es como funciona el
    diseno factorial. Con una lista se reparte a partes iguales por **resto
    mayor** y no por sorteo multinomial, para que la composicion de la mezcla
    sea exacta y reproducible.
    """
    if np.isscalar(omega):
        return np.full(n, float(omega))

    valores = np.asarray(omega, dtype=float)
    exacto  = np.full(len(valores), n / len(valores))
    cuenta  = np.floor(exacto).astype(int)
    resto   = n - int(cuenta.sum())
    orden   = np.argsort(-(exacto - cuenta), kind="stable")
    cuenta[orden[:resto]] += 1
    return np.repeat(valores, cuenta)


def generate_dataset(cfg: SimConfig) -> tuple[pd.DataFrame, SimulationTruth]:
    """
    Genera la matriz de expresion y la verdad de una celda.

    Todos los genes tienen ``M_g = 0`` y ruido ``N(0, sigma^2)``. Los
    ritmicos anaden una onda FMM. Los core se distinguen solo en que su
    ``alpha`` sale de invertir una posicion de pico conocida.
    """
    rng = np.random.default_rng(cfg.seed)

    t     = sample_times(cfg, rng)
    theta = t * TWO_PI / 24.0

    n_ritmicos = cfg.n_rhythmic
    if n_ritmicos < cfg.n_core:
        raise ValueError(
            f"Hacen falta al menos {cfg.n_core} genes ritmicos para alojar los "
            f"core, pero prop_rhythmic={cfg.prop_rhythmic} sobre "
            f"n_genes={cfg.n_genes} solo da {n_ritmicos}"
        )
    n_no_core = n_ritmicos - cfg.n_core
    n_fondo   = cfg.n_genes - n_ritmicos

    core_names, core_peaks = core_layout(cfg.n_core, cfg.day_fraction)
    core_sets = nested_core_sets(core_names, cfg.core_sizes, cfg.day_fraction)
    chain     = consistency_chain(core_names, core_peaks, core_sets,
                                  cfg.n_chain)

    # Los omegas se permutan antes de partir en core y no core para que, si
    # omega es una mezcla, ambos grupos reciban toda la variedad.
    omegas = rng.permutation(assign_omega(cfg.omega, n_ritmicos))
    betas  = rng.uniform(0.0, TWO_PI, n_ritmicos)
    om_core, om_libre = omegas[: cfg.n_core], omegas[cfg.n_core :]
    be_core, be_libre = betas[: cfg.n_core], betas[cfg.n_core :]

    al_core = np.asarray([
        alpha_for_peak(p, b, w)
        for p, b, w in zip(core_peaks, be_core, om_core)
    ])
    al_libre = rng.uniform(0.0, TWO_PI, n_no_core)

    filas:   list[np.ndarray] = []
    nombres: list[str] = []
    params:  list[dict] = []

    for nm, a, b, w in zip(core_names, al_core, be_core, om_core):
        filas.append(fmm_wave(theta, 0.0, cfg.amplitude, a, b, w))
        nombres.append(nm)
        params.append(dict(gene=nm, kind="core", M=0.0, A=cfg.amplitude,
                           alpha=a, beta=b, omega=w))

    for i, (a, b, w) in enumerate(zip(al_libre, be_libre, om_libre)):
        nm = f"RIT{i + 1:04d}"
        filas.append(fmm_wave(theta, 0.0, cfg.amplitude, a, b, w))
        nombres.append(nm)
        params.append(dict(gene=nm, kind="ritmico", M=0.0, A=cfg.amplitude,
                           alpha=a, beta=b, omega=w))

    for i in range(n_fondo):
        nm = f"GEN{i + 1:04d}"
        filas.append(np.zeros(cfg.n_samples))
        nombres.append(nm)
        params.append(dict(gene=nm, kind="fondo", M=0.0, A=0.0,
                           alpha=np.nan, beta=np.nan, omega=np.nan))

    senal = np.vstack(filas)
    ruido = rng.normal(0.0, cfg.sigma, senal.shape) if cfg.sigma > 0 else 0.0

    sample_names = [f"S{i + 1:04d}" for i in range(cfg.n_samples)]
    matrix = pd.DataFrame(senal + ruido, index=nombres, columns=sample_names)

    # Barajar las filas: ningun indice debe filtrar el tipo de gen.
    matrix = matrix.iloc[rng.permutation(cfg.n_genes)]

    truth = SimulationTruth(
        cfg          = cfg,
        t_true       = t,
        theta_true   = theta,
        sample_names = sample_names,
        core_names   = core_names,
        core_peaks   = core_peaks,
        core_sets    = core_sets,
        chain        = chain,
        gene_params  = pd.DataFrame(params).set_index("gene"),
    )
    return matrix, truth
