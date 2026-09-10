"""
Medidas de validacion de la reconstruccion temporal.

Todas trabajan con angulos en ``[0, 2*pi)`` y devuelven los errores en horas,
que es lo interpretable. La comparacion no necesita alineamiento previo: el
gen ancla tiene su pico en ``pi`` por construccion y el sincronizador coloca
el pico estimado en ``pi``, asi que las dos escalas viven en el mismo marco.
Las dos medidas ``mae_h_oracle`` y ``tau_circ`` existen solo como
diagnostico: dicen si el fallo fue del orden o solo del anclaje.
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
    Por eso existe ``tau_circular`` al lado, que no necesita corte.
    """
    a = np.asarray(a, dtype=float) % TWO_PI
    b = np.asarray(b, dtype=float) % TWO_PI
    resultado = kendalltau(a, b)
    return float(resultado.statistic)


def _orientacion(x: np.ndarray, i: int, J: np.ndarray, K: np.ndarray) -> np.ndarray:
    """Sentido de giro del triplete ``(i, J, K)`` sobre el circulo, en +-1."""
    return np.sign(
        np.sin(x[J] - x[i]) + np.sin(x[K] - x[J]) + np.sin(x[i] - x[K])
    )


def tau_circular(theta_true, theta_hat) -> float:
    """
    Concordancia circular por tripletes, sin corte y sin barrido.

    Dos puntos de un circulo no tienen orden mientras no se elija un corte, y
    de ahi salia todo el problema de ``tau``: la muestra que se estima pasado
    2*pi vuelve al principio y arrastra el ranking entero. Tres puntos si
    tienen orden: recorren el circulo en un sentido o en el otro, y eso no
    depende de ningun corte. Se cuenta en cuantos tripletes coinciden las dos
    escalas.

    Propiedades, todas medidas sobre 40 muestras:

    - invariante a rotar cualquiera de las dos escalas. Un desfase de 3 rad
      deja ``tau`` con corte fijo en 0,015, que se lee como si no hubiera
      relacion ninguna, mientras que esta da 1,000, que es la verdad.
    - el signo separa las dos patologias. ``+1`` con rotacion significa que el
      orden esta bien y fallo el origen, ``-1`` que fallo la orientacion.
      Reflejado y rotado a la vez deja ``tau`` en -0,015, otra vez como si no
      hubiera relacion, y esta lo llama -1,000.
    - vale ``0`` sin relacion y degrada de forma suave con el ruido.

    Coste ``O(n^3)`` en tiempo pero ``O(n^2)`` en memoria, y con ``n <= 96``
    son unos milisegundos.
    """
    a = np.asarray(theta_true, dtype=float) % TWO_PI
    b = np.asarray(theta_hat,  dtype=float) % TWO_PI
    n = a.size
    if n < 3:
        return float("nan")

    suma  = 0.0
    total = 0
    for i in range(n - 2):
        resto = np.arange(i + 1, n)
        jj, kk = np.triu_indices(resto.size, 1)
        J, K = resto[jj], resto[kk]
        suma  += float(np.sum(_orientacion(a, i, J, K)
                              * _orientacion(b, i, J, K)))
        total += J.size
    return suma / total


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
        # en el marco tal cual, es lo que se reporta
        "mae_h":    float(dh.mean()),
        "medae_h":  float(np.median(dh)),
        "rho_circ": circular_correlation(tt, th),
        "tau":      kendall_tau_fixed_cut(tt, th),
        # invariante al marco, es el diagnostico
        "tau_circ": tau_circular(tt, th),
    }
    # El sentido de giro, como columna propia. Va aparte del error porque son
    # cosas distintas: el orden puede estar bien recuperado y el ciclo
    # reportarse al reves, y entonces mae_h mide la suerte de la orientacion y
    # no la calidad de la reconstruccion. Promediada sobre replicas da la tasa
    # de acierto de orientacion.
    salida["orientacion_ok"] = float(salida["tau_circ"] > 0.0)
    for tol in TOLERANCIAS_H:
        salida[f"p{int(tol)}"] = float((dh <= tol).mean())
    salida["mae_h_oracle"] = oracle_mae_h(tt, th)
    return salida
