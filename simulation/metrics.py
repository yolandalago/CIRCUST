"""
Medidas de validacion de la reconstruccion temporal.

Todas trabajan con angulos en ``[0, 2*pi)`` y devuelven los errores en horas,
que es lo interpretable. La comparacion no necesita alineamiento previo: el
gen ancla tiene su pico en ``pi`` por construccion y el sincronizador coloca
el pico estimado en ``pi``, asi que las dos escalas viven en el mismo marco.
Las dos medidas ``oracle`` y ``best_cut`` existen solo como diagnostico.
"""
from __future__ import annotations

import numpy as np
from scipy.stats import kendalltau

from simulation.truth import TWO_PI

# Tolerancias en horas de la proporcion de muestras bien localizadas.
TOLERANCIAS_H = (1.0, 2.0, 4.0)


def circular_distance(a, b) -> np.ndarray:
    """Distancia circular en radianes, en ``[0, pi]``."""
    d = np.abs(np.asarray(a, dtype=float) - np.asarray(b, dtype=float)) % TWO_PI
    return np.minimum(d, TWO_PI - d)


def to_hours(rad) -> np.ndarray:
    """Convierte radianes a horas suponiendo un periodo de 24 h."""
    return np.asarray(rad, dtype=float) * 24.0 / TWO_PI


# Resultante minima para que la media circular tenga sentido. Solo salta con
# muestras exactamente equiespaciadas, donde R es del orden de 1e-17. Una
# muestra uniforme aleatoria de n=24 ya da R ~ 0,18 y pasa de sobra.
_R_MINIMA = 1e-8


def _mean_direction(x: np.ndarray) -> float:
    """
    Direccion media circular, o ``nan`` si la resultante es despreciable.

    Con muestras exactamente equiespaciadas las sumas de senos y cosenos se
    cancelan a cero y ``arctan2`` devuelve un angulo arbitrario que ademas no
    respeta la equivarianza bajo rotacion. Eso hacia que ``rho`` de 0,26 en
    vez de 1 al comparar una escala consigo misma rotada. Es preferible decir
    que no esta definida a colar un numero sin sentido en una tabla.
    """
    c, s = np.cos(x).mean(), np.sin(x).mean()
    if np.hypot(c, s) < _R_MINIMA:
        return float("nan")
    return float(np.arctan2(s, c))


def circular_correlation(a, b) -> float:
    """
    Correlacion circular-circular de Jammalamadaka.

        r = sum sin(a - a_) sin(b - b_) / sqrt( sum sin^2(a - a_) sum sin^2(b - b_) )

    con ``a_`` y ``b_`` las medias circulares. Es invariante a rotar
    cualquiera de las dos escalas, porque la media circular es equivariante.

    Devuelve ``nan`` si alguna de las dos escalas esta tan repartida por el
    circulo que su direccion media no esta definida. Eso ocurre con muestreo
    exactamente equiespaciado, que solo usa el bloque computacional, donde la
    medida de interes es el tiempo y no la calidad.
    """
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    ma = _mean_direction(a)
    mb = _mean_direction(b)
    if np.isnan(ma) or np.isnan(mb):
        return float("nan")
    sa = np.sin(a - ma)
    sb = np.sin(b - mb)
    denom = np.sqrt((sa ** 2).sum() * (sb ** 2).sum())
    if denom == 0.0:
        return float("nan")
    return float((sa * sb).sum() / denom)


def kendall_tau_fixed_cut(a, b) -> float:
    """
    Tau de Kendall sobre los rangos, con el corte comun en 0.

    Es la medida que pide la spec. Castiga mucho un desfase pequeno, porque
    las muestras que cruzan el corte saltan al extremo opuesto del ranking.
    Por eso existe ``kendall_tau_best_cut`` al lado.
    """
    a = np.asarray(a, dtype=float) % TWO_PI
    b = np.asarray(b, dtype=float) % TWO_PI
    resultado = kendalltau(a, b)
    return float(resultado.statistic)


def kendall_tau_best_cut(theta_true, theta_hat) -> float:
    """
    Maximo de ``kendall_tau_fixed_cut`` sobre los cortes de la escala estimada.

    El corte se mueve **solo en la estimada**, no en las dos a la vez. La
    escala verdadera tiene origen propio, el tiempo circadiano real, mientras
    que la estimada lo tiene arbitrario, y es ese el que se deja libre. Rotar
    las dos juntas no puede deshacer un desfase entre ellas: con la estimada
    adelantada 0,4 rad sobre 48 muestras, barrer el corte en ambas deja tau en
    0,761 mientras que barrerlo solo en la estimada lo lleva a 1,000.

    Con ``n <= 96`` el coste es despreciable.
    """
    a = np.asarray(theta_true, dtype=float) % TWO_PI
    b = np.asarray(theta_hat,  dtype=float) % TWO_PI
    mejor = -1.0
    for corte in np.unique(b):
        t = kendall_tau_fixed_cut(a, (b - corte) % TWO_PI)
        if t > mejor:
            mejor = t
    return float(mejor)


def oracle_mae_h(theta_true, theta_hat, n_grid: int = 720) -> float:
    """
    MAE en horas bajo la mejor rotacion y orientacion posibles.

    No es una medida de calidad del metodo, es un diagnostico: la diferencia
    con ``mae_h`` es el error atribuible a la sincronizacion y no a la
    reconstruccion del orden.
    """
    tt = np.asarray(theta_true, dtype=float)
    th = np.asarray(theta_hat, dtype=float)
    mejor = np.inf
    for s in (1.0, -1.0):
        base = (s * th) % TWO_PI
        for delta in np.linspace(0.0, TWO_PI, n_grid, endpoint=False):
            m = to_hours(circular_distance(tt, (base + delta) % TWO_PI)).mean()
            if m < mejor:
                mejor = m
    return float(mejor)


def evaluate(theta_true, theta_hat) -> dict[str, float]:
    """Las siete medidas de la spec sobre un par de escalas circulares."""
    tt = np.asarray(theta_true, dtype=float)
    th = np.asarray(theta_hat, dtype=float)
    if tt.shape != th.shape:
        raise ValueError(
            f"Las dos escalas deben tener la misma longitud, se recibio "
            f"{tt.shape} y {th.shape}"
        )

    dh = to_hours(circular_distance(tt, th))
    salida = {
        "mae_h":        float(dh.mean()),
        "medae_h":      float(np.median(dh)),
        "rho_circ":     circular_correlation(tt, th),
        "tau":          kendall_tau_fixed_cut(tt, th),
        "tau_best_cut": kendall_tau_best_cut(tt, th),
    }
    for tol in TOLERANCIAS_H:
        salida[f"p{int(tol)}"] = float((dh <= tol).mean())
    salida["mae_h_oracle"] = oracle_mae_h(tt, th)
    return salida
