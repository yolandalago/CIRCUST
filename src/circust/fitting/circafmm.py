"""
circust/fitting/circafmm.py
==================================
CircaFMM (Abston et al., 2026): version mejorada del ajuste FMM.

A diferencia del FMM base (``fmm.py``), este enfoque:

1. Usa una formulacion **numericamente estable** de la fase Mobius,
   sin singularidades en ``t = alpha +/- pi`` donde se tiende a inf.
2. Aplica **pesos por densidad temporal** (kernels KDE von Mises) que compensan
   distribuciones no uniformes de las muestras en el circulo.
3. Reemplaza la minimizacion del RSS por un **criterio de potencia
   penalizado**: J(alpha, omega) = P(alpha, omega) - lambda * g(omega) * R(alpha, omega).
4. Opcionalmente refina (alpha, omega) con busqueda de patron 2D (Algoritmo
   S2) + WLS cerrado para los parametros lineales, en vez de Nelder-Mead 5D.
   Desactivado por defecto, como en el pipeline publicado de CircaFMM
   (``runCircaFMM.R``, ``refine = FALSE``) y como en la calibracion de lambda,
   que nunca refina.
5. Calibra ``lambda`` automaticamente bajo H0 mediante Monte Carlo + codo.

El resultado es un ajuste mas robusto al ruido y mas eficiente para
datasets grandes (la rejilla y el cache se precomputan una sola vez y
se reutilizan para todos los genes).

Valores por defecto
-------------------
Donde el articulo (Abston, Larriba y Fernandez, 2026) no fija un valor se
toman los del script de analisis publicado por los autores,
``runCircaFMM.R`` del repositorio ``FMMGroupVa/CircaFMM``: 48 valores de
alpha, 50 valores de omega equiespaciados en ``[0.02, 1]``, rugosidad
``R_dense`` con 512 puntos, ``omega0 = 0.1``, ``p = 3``, sin refinamiento,
y calibracion con B = 1000, R = 5, q = 0.975 y mediana movil de ventana 3.
"""
from __future__ import annotations

import math
import threading
from dataclasses import dataclass, field
from typing import Optional

import numpy as np
from scipy.special import i0e as _bessel_i0e

from circust.fitting.rhythm_model import FitResult, RhythmModel

# ── Aceleracion opcional con Numba (mismo patron que fmm.py) ───────────────
try:
    from numba import njit
    _HAS_NUMBA = True
except ImportError:
    _HAS_NUMBA = False
    def njit(*args, **kwargs):                                  # type: ignore
        def deco(fn): return fn
        if len(args) == 1 and callable(args[0]):
            return args[0]
        return deco


# ═══════════════════════════════════════════════════════════════════════════
# Seccion 1 — Helpers numericos (privados)
# ═══════════════════════════════════════════════════════════════════════════

def _normalize_weights_to_n(w: np.ndarray) -> np.ndarray:
    """Reescala pesos para que ``sum(w) == len(w)``.

    Usado en multiples sitios (preparacion de pesos de densidad, calibracion
    bajo H0). Mantener como funcion separada porque se llama desde varios
    contextos con misma logica de validacion.
    """
    w = np.asarray(w, dtype=float)
    if not np.all(np.isfinite(w)):
        raise ValueError("Pesos no finitos.")
    if np.any(w < 0):
        raise ValueError("Los pesos deben ser >= 0.")
    sw = float(w.sum())
    if sw <= 0.0:
        raise ValueError("Suma de pesos debe ser > 0.")
    return w.shape[0] * w / sw


def _time_density_weights(
    times:      np.ndarray,
    kappa:      Optional[float]      = None,
    kappa_grid: Optional[np.ndarray] = None,
    eps:        float                = 1e-12,
) -> np.ndarray:
    """Pesos de densidad temporal von Mises:  ``w_i ~ 1 / f_hat(t_i)``,  ``sum(w) = n``.

    Construye una KDE von Mises sobre los ``times`` y devuelve pesos
    inversamente proporcionales a la densidad estimada las muestras en
    zonas temporales dispersas reciben mas peso.

    Si ``kappa`` no se especifica, se selecciona automaticamente por
    Leave-One-Out cross-validation sobre ``kappa_grid``.
    """
    times = np.asarray(times, dtype=float) % (2.0 * np.pi)
    n = times.shape[0]
    if n < 2:
        raise ValueError("Se necesitan al menos 2 muestras.")

    # Matriz n x n de cos(delta) — compartida por todos los candidatos kappa
    cos_delta = np.cos(times[:, None] - times[None, :])

    # Seleccion automatica de kappa por LOO-CV
    if kappa is None:
        if kappa_grid is None:
            kappa_grid = np.exp(np.linspace(np.log(1.0), np.log(64.0), 10))
        kappa_grid = np.asarray(kappa_grid, dtype=float)
        if np.any(~np.isfinite(kappa_grid)) or np.any(kappa_grid <= 0):
            raise ValueError("kappa_grid debe ser finito y > 0.")

        lcv = np.empty(kappa_grid.shape[0])
        for j, k in enumerate(kappa_grid):
            # Nucleo escalado: exp(k cos d)/I0(k) = exp(k (cos d - 1))/i0e(k).
            # Es la misma cantidad, pero ni el numerador ni I0 desbordan con
            # k grande (exp e i0 dan inf por encima de k ~ 709).
            # K[i,i] = exp(0)/(2π·i0e) es constante; restamos diagonal del row sum
            row_sums = np.exp(k * (cos_delta - 1.0)).sum(axis=1)
            f_loo = (row_sums - 1.0) / ((n - 1) * 2.0 * math.pi * float(_bessel_i0e(k)))
            lcv[j] = np.log(np.maximum(f_loo, eps)).sum()

        kappa = float(kappa_grid[int(np.argmax(lcv))])

    if not (np.isfinite(kappa) and kappa > 0):
        raise ValueError("kappa debe ser > 0.")

    # Densidad final (con todos los puntos, sin LOO)
    fhat = (np.exp(kappa * (cos_delta - 1.0)).sum(axis=1)
            / (n * 2.0 * math.pi * float(_bessel_i0e(kappa))))
    return _normalize_weights_to_n(1.0 / np.maximum(fhat, eps))


# ═══════════════════════════════════════════════════════════════════════════
# Seccion 2 — Fase Mobius estable
# ═══════════════════════════════════════════════════════════════════════════

@njit(cache=True, fastmath=True)
def _mobius_cos_sin_jit(
    t: np.ndarray,
    alpha: float,
    omega: float,
    den_eps: float,
) -> tuple:
    """
    Calcula ``cos(phi)`` y ``sin(phi)`` para la fase Mobius sin singularidades.

    phi(t; alpha, omega) = 2 * atan(omega * tan((t - alpha) / 2))

    Pero en vez de calcular ``phi`` explicitamente (que diverge cuando
    ``t - alpha = pi``), se calculan ``cos(phi)`` y ``sin(phi)`` directamente
    via identidades del angulo doble:

        u = (t - alpha) / 2
        den = cos^2(u) + omega^2 * sin^2(u)
        C = (cos^2(u) - omega^2 * sin^2(u)) / den
        S = (2 * omega * sin(u) * cos(u)) / den

    El denominador es siempre positivo (suma de cuadrados), nunca diverge.

    Renormaliza ``C^2 + S^2 = 1`` por seguridad numerica.
    """
    n = t.shape[0]
    C = np.empty(n)
    S = np.empty(n)

    if not np.isfinite(omega) or omega <= 0.0:
        # Limite exacto: phi = 0  =>  cos = 1, sin = 0
        for i in range(n):
            C[i] = 1.0
            S[i] = 0.0
        return C, S

    # Acotar omega a [0, 1]
    if omega > 1.0:
        omega = 1.0

    omega2 = omega * omega
    two_pi = 2.0 * math.pi

    for i in range(n):
        # Envolver t y alpha
        ti = t[i] % two_pi
        ai = alpha % two_pi
        u = 0.5 * (ti - ai)
        s = math.sin(u)
        c = math.cos(u)

        den = c * c + omega2 * s * s + den_eps
        C_i = (c * c - omega2 * s * s) / den
        S_i = (2.0 * omega * s * c) / den

        # Renormalizar a modulo 1
        r = math.sqrt(C_i * C_i + S_i * S_i) + den_eps
        C[i] = C_i / r
        S[i] = S_i / r

    return C, S


def _mobius_cos_sin(
    t: np.ndarray,
    alpha: float,
    omega: float,
    den_eps: float = np.finfo(float).eps,
) -> tuple[np.ndarray, np.ndarray]:
    """Wrapper publico de ``_mobius_cos_sin_jit``.

    Devuelve ``(C, S)`` con ``C[i] = cos(phi(t[i]; alpha, omega))`` y analogamente ``S``.
    """
    t_arr = np.ascontiguousarray(t, dtype=np.float64)
    return _mobius_cos_sin_jit(t_arr, float(alpha), float(omega), float(den_eps))


# ═══════════════════════════════════════════════════════════════════════════
# Seccion 3 — Estructuras de datos publicas
# ═══════════════════════════════════════════════════════════════════════════

@dataclass
class FMMGrid:
    """
    Rejilla precomputada de cantidades trigonometricas y de rugosidad.

    Se construye una sola vez por dataset (via ``precompute_grid``) y se
    reutiliza para el ajuste de todos los genes y para la calibracion de
    lambda.

    Attributes
    ----------
    times : np.ndarray, shape (n,)
        Puntos temporales originales (en radianes, [0, 2pi)).

    alpha_grid : np.ndarray, shape (nA,)
    omega_grid : np.ndarray, shape (nW,)
        Rejillas de valores para los parametros no lineales.

    alpha_vec : np.ndarray, shape (nG,)
    omega_vec : np.ndarray, shape (nG,)
        Aplanados con ``nG = nA * nW``: valor (alpha, omega) en cada punto.

    C, S : np.ndarray, shape (nG, n)
        Cantidades ``cos(phi)`` y ``sin(phi)`` para cada punto de rejilla y
        cada muestra. Calculadas via ``_mobius_cos_sin``.

    C2, S2, CS : np.ndarray, shape (nG, n)
        Productos elementales ``C*C``, ``S*S``, ``C*S``. Precomputados para
        evitar recalcularlos en cada ajuste.

    roughness : np.ndarray, shape (nG,)
        Rugosidad de curvatura ``R_dense(alpha, omega)`` por punto de rejilla,
        evaluada sobre la rejilla temporal densa y uniforme de ``n_dense``
        puntos que define la Seccion S3.4.2 de CircaFMM.

    n_dense : int
        Tamano de la rejilla temporal densa. Por defecto 512.

    order : np.ndarray, shape (n,)
        Permutacion que ordena los tiempos circularmente.

    rough_om_log, rough_val_log : np.ndarray
        Tabla 1-D log-log de ``R_dense`` frente a omega, para el refinamiento
        local fuera de rejilla. Ver la nota sobre independencia de alpha en
        :func:`_refine_jit`.
    """
    times:           np.ndarray
    alpha_grid:      np.ndarray
    omega_grid:      np.ndarray
    alpha_vec:       np.ndarray
    omega_vec:       np.ndarray
    C:               np.ndarray
    S:               np.ndarray
    C2:              np.ndarray
    S2:              np.ndarray
    CS:              np.ndarray
    roughness:       np.ndarray
    order:           np.ndarray
    n_dense:         int
    rough_om_log:    np.ndarray
    rough_val_log:   np.ndarray


@dataclass
class FMMCache:
    """
    Cache WLS dependiente de los pesos.

    Para cada punto de la rejilla, almacena la inversa ``(X'WX)^{-1}`` con
    ``X = [1, C, S]`` y los pesos fijos. Esto permite resolver la WLS exacta
    en cualquier punto de la rejilla sin recalcular Cholesky cada vez.

    Attributes
    ----------
    weights : np.ndarray, shape (n,)
        Pesos normalizados ``sum(w) = n``.

    inv_XWX : np.ndarray, shape (nG, 3, 3)
        Inversa de la matriz normal ponderada para cada punto de rejilla.

    den : np.ndarray, shape (nG,)
        Denominador del criterio de potencia,
        ``den[g] = sum_i w_i (C[g,i]^2 + S[g,i]^2)``. No depende del gen, de
        modo que se calcula una sola vez aqui en lugar de en cada ajuste.
        Nota: como la base ``(C, S)`` se renormaliza a modulo 1, este vector es
        constante e igual a ``sum(w) = n`` salvo error de redondeo. Se almacena
        calculado y no como escalar para no depender de esa propiedad.
    """
    weights: np.ndarray
    inv_XWX: np.ndarray
    den:     np.ndarray


@dataclass
class FMMPreparedDataset:
    """
    Bundle completo para ajuste con CircaFMM.

    Se crea **una sola vez** por dataset (via ``prepare_dataset``) y se
    reutiliza para todos los ajustes posteriores de genes.

    Attributes
    ----------
    times : np.ndarray, shape (n,)
        Puntos temporales (radianes en [0, 2pi)).

    weights : np.ndarray, shape (n,)
        Pesos normalizados a ``sum(w) = n``. Pueden venir del usuario o
        haberse calculado por densidad von Mises.

    grid : FMMGrid
        Rejilla precomputada.

    cache : FMMCache
        Cache WLS para los pesos.

    lambda_star : float
        Fuerza de la penalizacion de rugosidad. Calibrada por Monte Carlo
        bajo H0, o proporcionada por el usuario.

    omega0 : float
        Parametro del factor de penalizacion ``g(omega) = (omega0/(omega+omega0))^p``.

    penalty_power : float
        Exponente ``p`` del factor de penalizacion.
    """
    times:         np.ndarray
    weights:       np.ndarray
    grid:          FMMGrid
    cache:         FMMCache
    lambda_star:   float
    omega0:        float = 0.10
    penalty_power: float = 3.0
    # Vector ``lambda* g(omega) R(alpha, omega)``, que no depende del gen.
    penalty:       Optional[np.ndarray] = None


# ═══════════════════════════════════════════════════════════════════════════
# Seccion 4 — Precomputos (rejilla + cache)
# ═══════════════════════════════════════════════════════════════════════════

@njit(cache=True, fastmath=True)
def _fill_basis_jit(
    times:     np.ndarray,
    alpha_vec: np.ndarray,
    omega_vec: np.ndarray,
    den_eps:   float,
):
    """Rellena la base ``C, S`` de forma ``(nGenes, nMuestras)`` para todos los puntos de la rejilla.

    Bucle plano JIT sobre ``(gen, i)``. Equivale a llamar ``_mobius_cos_sin_jit``
    nGenes veces, pero sin overhead de funcion y permitiendo SIMD/loop fusion.
    """
    nG = alpha_vec.shape[0]
    n  = times.shape[0]
    C = np.empty((nG, n))
    S = np.empty((nG, n))
    two_pi = 2.0 * math.pi

    for g in range(nG):
        alpha = alpha_vec[g]
        omega = omega_vec[g]

        if not np.isfinite(omega) or omega <= 0.0:
            for i in range(n):
                C[g, i] = 1.0
                S[g, i] = 0.0
            continue
        if omega > 1.0:
            omega = 1.0

        omega2 = omega * omega
        ai = alpha % two_pi

        for i in range(n):
            ti = times[i] % two_pi
            u = 0.5 * (ti - ai)
            s = math.sin(u)
            c = math.cos(u)

            den = c * c + omega2 * s * s + den_eps
            C_i = (c * c - omega2 * s * s) / den
            S_i = (2.0 * omega * s * c) / den

            r = math.sqrt(C_i * C_i + S_i * S_i) + den_eps
            C[g, i] = C_i / r
            S[g, i] = S_i / r

    return C, S


def _curvature_roughness(
    C:     np.ndarray,
    S:     np.ndarray,
    dt:    np.ndarray,
    order: np.ndarray,
) -> np.ndarray:
    """Rugosidad de curvatura por punto de rejilla:  ``Sum(ddC^2) + Sum(ddS^2)``.

    Discretiza ``int |C''(t)|^2 dt + int |S''(t)|^2 dt`` mediante diferencias
    finitas circulares: pendiente local entre vecinos consecutivos, y luego
    diferencia de pendientes (segunda derivada discreta).

    Vectorizado con ``np.roll`` para los desplazamientos circulares.
    """
    C_ord = C[:, order]
    S_ord = S[:, order]
    dt_b  = dt[np.newaxis, :]

    slope_C = (np.roll(C_ord, -1, axis=1) - C_ord) / dt_b
    slope_S = (np.roll(S_ord, -1, axis=1) - S_ord) / dt_b

    ddC = slope_C - np.roll(slope_C, 1, axis=1)
    ddS = slope_S - np.roll(slope_S, 1, axis=1)

    return (ddC * ddC).sum(axis=1) + (ddS * ddS).sum(axis=1)


def precompute_grid(
    times:      np.ndarray,
    alpha_grid: Optional[np.ndarray] = None,
    omega_grid: Optional[np.ndarray] = None,
    n_dense:    int = 512,
) -> FMMGrid:
    """Precomputa la base trigonometrica y la rugosidad para la rejilla.

    Para cada par ``(alpha_g, omega_g)`` de la rejilla:

    - ``C[g, i], S[g, i]``: coseno y seno de la fase Mobius en cada muestra.
    - ``C2, S2, CS``: productos elementales (precomputados para WLS rapida).
    - ``roughness[g]``: ``R_dense(alpha, omega)`` de la Seccion S3.4.2.

    Sobre el modo de rugosidad
    --------------------------
    CircaFMM define la rugosidad de dos maneras. ``R_obs`` la evalua sobre los
    propios instantes de observacion y sus huecos circulares, de modo que
    depende del diseno de muestreo. ``R_dense`` la evalua sobre una rejilla
    uniforme de ``n_dense`` puntos con espaciado constante, de modo que depende
    unicamente de ``(alpha, omega)``.

    Aqui se implementa solo ``R_dense``. El motivo es el que da el propio
    articulo. Con huecos amplios en el muestreo la transicion de la fase en
    torno a ``t = alpha + pi`` puede caer integramente dentro de un intervalo
    no observado, con lo que ``R_obs`` no la ve y la penalizacion no se
    dispara. A la inversa, instantes casi duplicados dan incrementos
    ``Delta t`` minusculos que inflan las pendientes de forma erratica. Ambas
    situaciones son las de CIRCUST, donde los instantes son la escala circular
    inducida por el CPCA y por tanto irregulares por construccion. Es tambien
    la eleccion que hacen los autores para GTEx-13T, cuyos tiempos proceden
    igualmente de un orden estimado y no de un diseno regular.

    Parameters
    ----------
    times : np.ndarray, shape (n,)
        Puntos temporales en radianes (se envuelve a ``[0, 2*pi)``).
    alpha_grid : np.ndarray, opcional
        Rejilla de alpha. Por defecto: 48 valores en ``[0, 2*pi)``.
    omega_grid : np.ndarray, opcional
        Rejilla de omega. Por defecto: 50 valores equiespaciados en
        ``[0.02, 1]``, que es la rejilla del pipeline publicado de CircaFMM
        (``runCircaFMM.R``). La funcion R ``fmm_grid_precompute`` trae como
        valor por defecto 30 valores log-espaciados en ``[1e-4, 1]``; el
        articulo no fija ninguna de las dos.
    n_dense : int
        Tamano de la rejilla temporal densa. Por defecto 512.
    """
    times = np.ascontiguousarray(times, dtype=np.float64) % (2.0 * np.pi)
    n = times.shape[0]
    if n < 4:
        raise ValueError("Se necesitan al menos 4 muestras.")
    if n_dense < 16:
        raise ValueError("n_dense debe ser >= 16.")

    if alpha_grid is None:
        alpha_grid = np.linspace(0.0, 2.0 * np.pi, 48, endpoint=False)
    if omega_grid is None:
        omega_grid = np.linspace(0.02, 1.0, 50)

    alpha_grid = np.asarray(alpha_grid, dtype=np.float64) % (2.0 * np.pi)
    omega_grid = np.clip(np.asarray(omega_grid, dtype=np.float64), 0.0, 1.0)
    if np.any(omega_grid <= 0):
        raise ValueError("omega_grid debe contener valores en (0, 1].")

    # Aplanado de la rejilla a vectores de longitud nG = nA * nW
    A_mesh, W_mesh = np.meshgrid(alpha_grid, omega_grid, indexing="ij")
    alpha_vec = A_mesh.ravel()
    omega_vec = W_mesh.ravel()

    eps = float(np.finfo(float).eps)

    # Base trigonometrica sobre los tiempos observados (kernel JIT)
    C, S = _fill_basis_jit(times, alpha_vec, omega_vec, eps)

    # Productos elementales (se reusan para la WLS y para la potencia)
    C2 = C * C
    S2 = S * S
    CS = C * S

    # Orden circular de las muestras (lo usa el diagnostico de salto maximo)
    order = np.argsort(times)

    # ── R_dense sobre la rejilla temporal uniforme ───────────────────────
    t_dense  = np.linspace(0.0, 2.0 * np.pi, n_dense, endpoint=False)
    dt_dense = np.full(n_dense, 2.0 * np.pi / n_dense)
    ord_dense = np.arange(n_dense)

    C_d, S_d = _fill_basis_jit(t_dense, alpha_vec, omega_vec, eps)
    roughness = _curvature_roughness(C_d, S_d, dt_dense, ord_dense)

    # Tabla 1-D R_dense(omega) para el refinamiento local fuera de rejilla.
    # Basta con un valor de alpha porque R_dense no depende de el (rejilla
    # temporal uniforme sobre toda la circunferencia).
    n_table = 512
    om_tab  = np.exp(np.linspace(np.log(1e-4), np.log(1.0), n_table))
    C_t, S_t = _fill_basis_jit(t_dense, np.zeros(n_table), om_tab, eps)
    R_tab = _curvature_roughness(C_t, S_t, dt_dense, ord_dense)

    return FMMGrid(
        times         = times,
        alpha_grid    = alpha_grid,
        omega_grid    = omega_grid,
        alpha_vec     = alpha_vec,
        omega_vec     = omega_vec,
        C             = C,
        S             = S,
        C2            = C2,
        S2            = S2,
        CS            = CS,
        roughness     = roughness,
        order         = order,
        n_dense       = int(n_dense),
        rough_om_log  = np.log(om_tab),
        rough_val_log = np.log(np.maximum(R_tab, 1e-300)),
    )


def prepare_cache(
    weights: np.ndarray,
    grid:    FMMGrid,
    ridge:   float = 1e-10,
) -> FMMCache:
    """Calcula la inversa ``(X'WX)^{-1}`` para cada punto de rejilla.

    Construye, para cada ``g``, la matriz normal ponderada::

        X = [1 | C[g] | S[g]]   (n x 3)
        X'WX = [ sum(w)        sum(w*C)      sum(w*S)     ]
               [ sum(w*C)      sum(w*C^2)    sum(w*C*S)   ]
               [ sum(w*S)      sum(w*C*S)    sum(w*S^2)   ]

    e invierte con ``np.linalg.inv`` el tensor ``(nG, 3, 3)`` completo de
    una sola pasada. Una pequena regularizacion ``ridge`` se anade a la
    diagonal para estabilidad numerica.

    Parameters
    ----------
    weights : np.ndarray, shape (n,)
        Pesos (se normalizan internamente a ``sum(w) = n``).
    grid : FMMGrid
        Rejilla precomputada.
    ridge : float
        Regularizacion diagonal para estabilidad. Por defecto 1e-10, el mismo
        valor que ``fmm_cache_prepare`` en el codigo R de CircaFMM.
    """
    w = _normalize_weights_to_n(weights)
    if w.shape[0] != grid.times.shape[0]:
        raise ValueError("weights debe tener la misma longitud que grid.times.")

    nG = grid.C.shape[0]

    # Productos elementales con los pesos -> sumas escalares y por punto de rejilla
    a11 = float(w.sum())
    a12 = grid.C  @ w     # (nG,)
    a13 = grid.S  @ w     # (nG,)
    a22 = grid.C2 @ w     # (nG,)
    a23 = grid.CS @ w     # (nG,)
    a33 = grid.S2 @ w     # (nG,)

    # Tensor (nG, 3, 3) de matrices normales ponderadas + ridge
    A = np.empty((nG, 3, 3), dtype=np.float64)
    A[:, 0, 0] = a11 + ridge
    A[:, 0, 1] = a12; A[:, 1, 0] = a12
    A[:, 0, 2] = a13; A[:, 2, 0] = a13
    A[:, 1, 1] = a22 + ridge
    A[:, 1, 2] = a23; A[:, 2, 1] = a23
    A[:, 2, 2] = a33 + ridge

    # Invertir todo el tensor de una vez (NumPy/LAPACK vectorizado)
    inv_XWX = np.linalg.inv(A)

    # Denominador del criterio de potencia: a22 + a33 = sum_i w_i (C^2 + S^2)
    den = a22 + a33

    return FMMCache(weights=w, inv_XWX=inv_XWX, den=den)


# ═══════════════════════════════════════════════════════════════════════════
# Seccion 5 — Ajuste por gen (criterio de potencia + penalizacion)
# ═══════════════════════════════════════════════════════════════════════════

def _penalty_vector(
    grid:          FMMGrid,
    lambda_:       float,
    omega0:        float,
    penalty_power: float,
) -> np.ndarray:
    """Termino de penalizacion ``lambda * g(omega) * R(alpha, omega)`` por punto de rejilla.

    No depende de los datos del gen, solo de la rejilla y de los parametros de
    regularizacion, de modo que se calcula una sola vez por conjunto de datos.
    """
    if lambda_ <= 0:
        return np.zeros(grid.alpha_vec.shape[0], dtype=np.float64)
    g_omega = (omega0 / (grid.omega_vec + omega0)) ** penalty_power
    return lambda_ * g_omega * grid.roughness


def _compute_power_scores(
    data:          np.ndarray,
    grid:          FMMGrid,
    weights:       np.ndarray,
    lambda_:       float,
    omega0:        float,
    penalty_power: float,
    den:           Optional[np.ndarray] = None,
    penalty:       Optional[np.ndarray] = None,
) -> np.ndarray:
    """Calcula el criterio  ``J(alpha, omega) = P - lambda*g(omega)*R``  para toda la rejilla.

    - ``P(alpha, omega) = (u^2 + v^2) / den``   (potencia normalizada Lomb-Scargle)
      con  ``u = sum(w*x_centered*C)``,  ``v = sum(w*x_centered*S)``,
      ``den = sum(w*C^2) + sum(w*S^2)``.
    - ``g(omega) = (omega0 / (omega + omega0))^p``
    - ``R(alpha, omega) = R_dense(alpha, omega)``, en ``grid.roughness``.

    Operaciones todas vectorizadas matriz-vector — un solo gen cuesta O(nG*n).

    ``den`` y ``penalty`` son cantidades que no dependen del gen. Si se pasan
    precomputadas se ahorran, por cada gen, los dos productos matriz-vector que
    calculan el denominador y el vector de penalizacion completo.
    """
    # Centrar los datos por la media ponderada (asi el intercepto desaparece de la ecuacion)
    mu_w = float((weights * data).sum() / weights.sum())
    wx0 = weights * (data - mu_w)

    # Potencia: u, v para todos los puntos a la vez
    u = grid.C @ wx0          # (nG,)
    v = grid.S @ wx0          # (nG,)
    pow_ = u * u + v * v

    # Normalizacion Lomb-Scargle: dividir por la norma de la base ponderada
    if den is None:
        den = grid.C2 @ weights + grid.S2 @ weights
    den = np.where(den > 0, den, np.nan)
    pow_ = pow_ / den

    if lambda_ <= 0:
        return pow_

    if penalty is None:
        penalty = _penalty_vector(grid, lambda_, omega0, penalty_power)
    return pow_ - penalty


def _wls_fit_at(
    data:    np.ndarray,
    alpha:   float,
    omega:   float,
    weights: np.ndarray,
    times:   np.ndarray,
    ridge:   float = 1e-10,
) -> dict:
    """Resuelve la WLS exacta 3x3 en un punto ``(alpha, omega)`` arbitrario.

    Construye la base ``[1, C, S]``, resuelve  ``(X'WX) coef = X'Wx``, y devuelve
    los parametros del FMM con  ``A = sqrt(b^2 + g^2)``,  ``beta = atan2(-g, b)``.

    Devuelve dict con ``M, A, alpha, beta, omega, fitted, rss``.
    """
    C, S = _mobius_cos_sin(times, alpha, omega)

    # Sistema normal ponderado 3x3
    a11 = weights.sum() + ridge
    a12 = (weights * C).sum()
    a13 = (weights * S).sum()
    a22 = (weights * C * C).sum() + ridge
    a23 = (weights * C * S).sum()
    a33 = (weights * S * S).sum() + ridge

    A = np.array([[a11, a12, a13],
                  [a12, a22, a23],
                  [a13, a23, a33]])
    b = np.array([(weights * data).sum(),
                  (weights * data * C).sum(),
                  (weights * data * S).sum()])

    coef = np.linalg.solve(A, b)
    M, b_coef, g_coef = coef

    fitted = M + b_coef * C + g_coef * S
    rss = float((weights * (data - fitted) ** 2).sum())

    A_amp = float(np.sqrt(b_coef ** 2 + g_coef ** 2))
    beta  = float(np.arctan2(-g_coef, b_coef) % (2.0 * np.pi))

    return dict(
        M = float(M), A = A_amp, alpha = float(alpha) % (2.0 * np.pi),
        beta = beta, omega = float(omega),
        fitted = fitted, rss = rss,
    )


@njit(cache=True, fastmath=True)
def _interp_loglog(x: float, xp_log: np.ndarray, fp_log: np.ndarray) -> float:
    """Interpolacion lineal en escala log-log, con extrapolacion constante.

    ``xp_log`` debe estar ordenado de forma creciente. Se usa para leer
    ``R_dense(omega)`` de una tabla precomputada.
    """
    lx = math.log(x)
    m  = xp_log.shape[0]
    if lx <= xp_log[0]:
        return math.exp(fp_log[0])
    if lx >= xp_log[m - 1]:
        return math.exp(fp_log[m - 1])
    lo, hi = 0, m - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if xp_log[mid] <= lx:
            lo = mid
        else:
            hi = mid
    w = (lx - xp_log[lo]) / (xp_log[hi] - xp_log[lo])
    return math.exp(fp_log[lo] + w * (fp_log[hi] - fp_log[lo]))


@njit(cache=True, fastmath=True)
def _refine_jit(
    times:      np.ndarray,
    data:       np.ndarray,
    weights:    np.ndarray,
    wx0:        np.ndarray,
    om_log:     np.ndarray,
    rough_log:  np.ndarray,
    alpha0:     float,
    omega_ini:  float,
    score_init: float,
    lambda_:    float,
    omega0:     float,
    pen_power:  float,
    n_steps:    int,
    alpha_step: float,
    omega_mult: float,
    omega_min:  float,
    omega_max:  float,
    tau:        float,
    ridge:      float,
    den_eps:    float,
):
    """Refinamiento local de ``(alpha, omega)``, Algoritmo S2 de CircaFMM.

    Tras el maximizador grueso sobre la rejilla se aplica una busqueda de
    patron sobre un entorno discreto de ocho vecinos, dos movimientos en
    ``alpha``, dos en ``omega`` y cuatro combinados, con ``alpha`` envuelto
    modulo ``2*pi`` y ``omega`` recortado al intervalo admisible.

    La eleccion entre los ocho candidatos no es el simple maximo de
    ``J(alpha, omega)``. El articulo define tres etapas.

    Primero se retienen los candidatos cuyo criterio queda a menos de una
    tolerancia relativa del mejor,

        J(alpha, omega) >= Jmax - tau * max(1, |Jmax|),   tau = 1e-4.

    Despues, entre los retenidos, se elige el que minimiza el objetivo
    secundario ``WRSS(alpha, omega) + lambda g(omega) R(alpha, omega)``, donde
    ``WRSS`` procede del ajuste por minimos cuadrados ponderados condicionado
    al par candidato. Por ultimo, si ese objetivo empata, se prefiere el
    ``omega`` mayor, con lo que se evita el desplazamiento innecesario hacia la
    frontera de ``omega`` pequenos, que es justamente la region que toda la
    regularizacion trata de evitar.

    El candidato asi elegido se acepta solo si mejora estrictamente el
    criterio actual. En caso contrario se contrae el entorno, ``alpha_step``
    a la mitad y ``omega_mult`` a su raiz cuadrada.

    La rugosidad se lee por interpolacion log-log de una tabla en ``omega``.
    Es licito porque ``R_dense`` no depende de ``alpha``: la rejilla temporal
    densa es uniforme y cubre toda la circunferencia, de modo que trasladar
    ``alpha`` solo rota que punto de la rejilla evalua cada valor de la base y
    la suma de curvaturas discretas queda invariante. Se comprueba
    numericamente que la variacion con ``alpha`` es del orden de 1e-16 en
    terminos relativos, y el error de la interpolacion con 512 nodos es como
    mucho de 3e-4 relativo.

    Devuelve ``(alpha, omega)``.
    """
    n      = times.shape[0]
    two_pi = 2.0 * math.pi

    da      = alpha_step
    om_step = omega_mult
    a_cur   = alpha0
    o_cur   = omega_ini
    sc_cur  = score_init

    cand_a   = np.empty(8)
    cand_o   = np.empty(8)
    cand_J   = np.empty(8)
    cand_sec = np.empty(8)
    Cb       = np.empty(n)
    Sb       = np.empty(n)

    for _ in range(n_steps):

        for c in range(8):
            # 2 movimientos en alpha, 2 en omega y 4 combinados
            if c == 0:
                a = a_cur + da; o = o_cur
            elif c == 1:
                a = a_cur - da; o = o_cur
            elif c == 2:
                a = a_cur;      o = o_cur * om_step
            elif c == 3:
                a = a_cur;      o = o_cur / om_step
            elif c == 4:
                a = a_cur + da; o = o_cur * om_step
            elif c == 5:
                a = a_cur - da; o = o_cur * om_step
            elif c == 6:
                a = a_cur + da; o = o_cur / om_step
            else:
                a = a_cur - da; o = o_cur / om_step

            a = a % two_pi
            if o < omega_min:
                o = omega_min
            elif o > omega_max:
                o = omega_max

            cand_a[c] = a
            cand_o[c] = o

            # ── base de la fase y acumuladores ────────────────────────────
            o2  = o * o
            u   = 0.0
            v   = 0.0
            den = 0.0
            a11 = 0.0; a12 = 0.0; a13 = 0.0
            a22 = 0.0; a23 = 0.0; a33 = 0.0
            b1  = 0.0; b2  = 0.0; b3  = 0.0
            for i in range(n):
                uu = 0.5 * ((times[i] % two_pi) - a)
                s  = math.sin(uu)
                cc = math.cos(uu)
                d  = cc * cc + o2 * s * s + den_eps
                C  = (cc * cc - o2 * s * s) / d
                S  = (2.0 * o * s * cc) / d
                r  = math.sqrt(C * C + S * S) + den_eps
                C /= r
                S /= r
                Cb[i] = C
                Sb[i] = S

                wi = weights[i]
                xi = data[i]
                u   += wx0[i] * C
                v   += wx0[i] * S
                den += wi * (C * C + S * S)

                a11 += wi
                a12 += wi * C
                a13 += wi * S
                a22 += wi * C * C
                a23 += wi * C * S
                a33 += wi * S * S
                b1  += wi * xi
                b2  += wi * xi * C
                b3  += wi * xi * S

            if den <= 0.0:
                cand_J[c]   = -1.0e308
                cand_sec[c] =  1.0e308
                continue

            # ── criterio J = P - lambda g(omega) R ────────────────────────
            P = (u * u + v * v) / den
            pen = 0.0
            if lambda_ > 0.0:
                R = _interp_loglog(o, om_log, rough_log)
                g = (omega0 / (o + omega0)) ** pen_power
                pen = lambda_ * g * R
            cand_J[c] = P - pen

            # ── objetivo secundario: WRSS + lambda g(omega) R ─────────────
            a11 += ridge; a22 += ridge; a33 += ridge
            m11 = a22 * a33 - a23 * a23
            m12 = a12 * a33 - a23 * a13
            m13 = a12 * a23 - a22 * a13
            det = a11 * m11 - a12 * m12 + a13 * m13
            if abs(det) < 1e-300:
                cand_sec[c] = 1.0e308
                continue
            M   = (b1 * m11 - a12 * (b2 * a33 - a23 * b3)
                   + a13 * (b2 * a23 - a22 * b3)) / det
            dlt = (a11 * (b2 * a33 - a23 * b3) - b1 * m12
                   + a13 * (a12 * b3 - b2 * a13)) / det
            gam = (a11 * (a22 * b3 - b2 * a23) - a12 * (a12 * b3 - b2 * a13)
                   + b1 * m13) / det
            wrss = 0.0
            for i in range(n):
                e = data[i] - (M + dlt * Cb[i] + gam * Sb[i])
                wrss += weights[i] * e * e
            cand_sec[c] = wrss + pen

        # ── etapa 1: candidatos dentro de la tolerancia relativa ──────────
        Jmax = -1.0e308
        for c in range(8):
            if cand_J[c] > Jmax:
                Jmax = cand_J[c]
        if Jmax <= -1.0e307:
            da      *= 0.5
            om_step  = math.sqrt(om_step)
            continue
        thr = Jmax - tau * max(1.0, abs(Jmax))

        # ── etapas 2 y 3: minimo del objetivo secundario, empate a omega mayor
        best = -1
        for c in range(8):
            if cand_J[c] < thr:
                continue
            if best < 0:
                best = c
            elif cand_sec[c] < cand_sec[best] - 1e-12:
                best = c
            elif abs(cand_sec[c] - cand_sec[best]) <= 1e-12 and cand_o[c] > cand_o[best]:
                best = c

        # ── aceptacion solo si mejora estrictamente ───────────────────────
        if best >= 0 and cand_J[best] > sc_cur:
            a_cur  = cand_a[best]
            o_cur  = cand_o[best]
            sc_cur = cand_J[best]
        else:
            da      *= 0.5
            om_step  = math.sqrt(om_step)

    return a_cur, o_cur


def fit_gene(
    data:         np.ndarray,
    prepared:     FMMPreparedDataset,
    refine:       bool = False,
    refine_steps: int  = 3,
) -> FitResult:
    """Ajuste CircaFMM para un gen.

    Algoritmo:
      1. Calcula ``J(alpha, omega) = P - lambda*g(omega)*R`` para toda la rejilla.
      2. Selecciona el ``(alpha*, omega*)`` que maximiza ``J``.
      3. (Opcional) Refinamiento local, Algoritmo S2. Ver :func:`_refine_jit`.
      4. WLS exacta 3x3 en ``(alpha*, omega*)`` para obtener ``M, A, beta``.
      5. Construye ``FitResult`` compatible con ``FMMModel``.

    En el paso 2 se toma el maximo sin tolerancias ni desempates, que es como
    lo define la Seccion S3.4.3 del articulo. La regla de tolerancia relativa y
    el criterio secundario aparecen unicamente en el refinamiento del paso 3.

    ``refine`` esta desactivado por defecto: el articulo lo describe como
    opcional, el pipeline publicado (``runCircaFMM.R``) lo ejecuta con
    ``refine = FALSE`` y la calibracion de ``lambda`` (Algoritmo S3) nunca
    refina, de modo que con ``refine=True`` el ajuste por gen no seria el
    mismo procedimiento con el que se calibro ``lambda*``. En ``omega = 1`` el
    Algoritmo S2 es ademas degenerado (la direccion ``alpha`` es plana porque
    ``beta`` la absorbe) y su resultado depende del redondeo.

    ``refine_steps`` es el numero ``S`` de pasos de la busqueda de patron. El
    valor por defecto, 3, es el del articulo.

    El ``model_name`` del resultado es ``"circafmm"`` (vs ``"fmm"`` en el base),
    por si se quiere distinguir mas adelante. Las claves de ``params`` son las
    mismas (M, A, alpha, beta, omega) y la formula de ``peak_time`` es identica.

    El ``r2`` es el ponderado con los pesos de densidad temporal,
    ``1 - sum(w e^2) / sum(w (x - mean_w(x))^2)``, coherente con la WLS del
    paso 4. ``sse`` sigue siendo la suma de cuadrados sin pesos.
    """
    data = np.asarray(data, dtype=float)
    if data.shape[0] != prepared.times.shape[0]:
        raise ValueError("data debe tener la misma longitud que prepared.times.")
    if not np.all(np.isfinite(data)):
        raise ValueError("data contiene valores no finitos.")

    grid    = prepared.grid
    weights = prepared.cache.weights
    lambda_ = float(prepared.lambda_star)
    omega0  = float(prepared.omega0)
    p       = float(prepared.penalty_power)

    # 1. Score sobre toda la rejilla. ``den`` y ``penalty`` van precomputados
    # porque no dependen del gen.
    scores = _compute_power_scores(
        data, grid, weights, lambda_, omega0, p,
        den=prepared.cache.den, penalty=prepared.penalty,
    )
    scores_safe = np.where(np.isfinite(scores), scores, -np.inf)
    if not np.any(np.isfinite(scores_safe)):
        raise RuntimeError("Ningun score finito en la rejilla.")

    # 2. Seleccion sobre la rejilla: maximizador de J(alpha, omega), sin
    # tolerancias ni desempates (Seccion S3.4.3).
    idx = int(np.argmax(scores_safe))
    alpha = float(grid.alpha_vec[idx])
    omega = float(grid.omega_vec[idx])

    # 3. Refinamiento local opcional (Algoritmo S2)
    if refine:
        mu_w = float((weights * data).sum() / weights.sum())
        wx0  = np.ascontiguousarray(weights * (data - mu_w))
        alpha, omega = _refine_jit(
            grid.times, np.ascontiguousarray(data), weights, wx0,
            grid.rough_om_log, grid.rough_val_log,
            float(alpha), float(omega), float(scores_safe[idx]),
            lambda_, omega0, p,
            int(refine_steps),
            (2.0 * np.pi) / 48.0,                            # Delta_alpha (S3.4.4)
            1.25,                                            # eta
            1e-4, 1.0,                                       # omega_min, omega_max
            1e-4,                                            # tau
            1e-10,                                           # ridge de la WLS
            float(np.finfo(float).eps),
        )
        alpha, omega = float(alpha), float(omega)

    # 4. WLS final exacta en el punto seleccionado
    fit = _wls_fit_at(data, alpha, omega, weights, prepared.times)

    # 5. Construir FitResult (mismo formato que FMMModel)
    fitted = fit["fitted"]
    residuals = data - fitted
    # ddof=5: parametros libres del modelo (M, A, alpha, beta, omega), igual
    # que en FMMModel. El anexo "Deteccion de muestras atipicas" usa ddof=1,
    # pero eso solo es insesgado si el unico parametro estimado es la media;
    # aqui se estimaron 5, asi que el denominador correcto es n-5.
    res_std = residuals.std(ddof=5)
    residuals_std = (residuals - residuals.mean()) / res_std if res_std > 0 \
                    else np.zeros_like(residuals)
    sse = float(np.sum(residuals ** 2))
    mu_w = float((weights * data).sum() / weights.sum())
    ss_tot_w = float((weights * (data - mu_w) ** 2).sum())
    r2 = float(1.0 - (weights * residuals ** 2).sum() / ss_tot_w) if ss_tot_w > 0 else 0.0

    peak_t = float(
        (fit["alpha"] + 2.0 * np.arctan2(
            (1.0 / fit["omega"]) * np.sin(-fit["beta"] / 2.0),
            np.cos(-fit["beta"] / 2.0),
        )) % (2.0 * np.pi)
    )

    return FitResult(
        fitted        = fitted,
        params        = {
            "M":     fit["M"],
            "A":     fit["A"],
            "alpha": fit["alpha"],
            "beta":  fit["beta"],
            "omega": fit["omega"],
        },
        peak_time     = peak_t,
        r2            = r2,
        residuals     = residuals,
        residuals_std = residuals_std,
        sse           = sse,
        model_name    = "circafmm",
    )


# ═══════════════════════════════════════════════════════════════════════════
# Seccion 6 — Calibracion de lambda
# ═══════════════════════════════════════════════════════════════════════════

def _spike_max_jump(y: np.ndarray, order: np.ndarray) -> float:
    """Salto maximo entre valores consecutivos en orden circular.

    Diagnostico de sobreajuste para calibrar lambda. Bajo H0 (datos planos),
    si el modelo ajusta picos artificiales, los valores ``fitted`` tendran
    saltos grandes entre vecinos consecutivos.
    """
    yo = y[order]
    diffs = np.abs(np.concatenate([np.diff(yo), np.array([yo[0] - yo[-1]])]))
    return float(diffs.max())


def _runmed_r(y: np.ndarray, k: int) -> np.ndarray:
    """``stats::runmed(y, k, endrule = "median")`` de R.

    Mediana movil de ventana impar ``k`` en el interior y regla de Tukey en
    los extremos (``stats::smoothEnds``). Para ``k = 3`` la regla es

        y[0]  = med3(y[0],  s[1],  3 s[1]  - 2 s[2])
        y[-1] = med3(y[-1], s[-2], 3 s[-2] - 2 s[-3])

    con ``s`` la serie suavizada. Para ``k > 3`` R suaviza primero los
    ``(k-1)/2`` puntos de cada extremo con ventanas crecientes ``3, 5, ...``
    aplicadas sobre la serie que sale de la mediana movil (interior suavizado,
    extremos sin tocar), y despues aplica la regla anterior al primer y
    ultimo punto. Se replica esa misma secuencia.
    """
    y = np.asarray(y, dtype=float)
    n = y.shape[0]
    if k <= 1 or n < 3:
        return y.copy()
    k = min(k, n if n % 2 == 1 else n - 1)
    h = k // 2

    # runmed con iend = 0 ("keep"): interior suavizado, h puntos de cada
    # extremo sin tocar
    res = y.copy()
    for i in range(h, n - h):
        res[i] = np.median(y[i - h:i + h + 1])

    # smoothEnds(res, k)
    def med3(a, b, c):
        return float(np.median([a, b, c]))

    sm = res.copy()
    if h >= 2:
        sm[1]     = med3(res[0], res[1], res[2])
        sm[n - 2] = med3(res[n - 1], res[n - 2], res[n - 3])
        for i in range(3, h + 1):           # i = 3..h, ventana j = 2i-1
            j = 2 * i - 1
            sm[i - 1] = np.median(res[0:j])
            sm[n - i] = np.median(res[n - j:n])
    sm[0]     = med3(res[0],     sm[1],     3.0 * sm[1]     - 2.0 * sm[2])
    sm[n - 1] = med3(res[n - 1], sm[n - 2], 3.0 * sm[n - 2] - 2.0 * sm[n - 3])
    return sm


def _pick_lambda_elbow(
    lambda_grid: np.ndarray,
    q_spike:     np.ndarray,
    smooth_k:    int  = 3,
    use_log10:   bool = True,
) -> float:
    """Detecta el codo de ``q_spike(lambda)`` por distancia a la cuerda.

    En la curva ``(log10(lambda), q_spike)``:
      1. Suaviza ``q_spike`` con mediana movil (ventana ``smooth_k``).
      2. Traza la cuerda del primer al ultimo punto.
      3. El codo es el punto interno mas alejado de la cuerda.

    La mediana movil replica ``stats::runmed(y, k, endrule = "median")`` de
    R, que es lo que usa ``pick_lambda_elbow`` en el codigo de CircaFMM. Los
    puntos interiores se sustituyen por la mediana de su ventana y los dos
    extremos por la regla de Tukey (``smoothEnds``): para ``k = 3``

        y[0]  = mediana(y[0],  s[1],  3 s[1]  - 2 s[2]),
        y[-1] = mediana(y[-1], s[-2], 3 s[-2] - 2 s[-3]),

    con ``s`` la serie ya suavizada en el interior. Importa porque los dos
    extremos definen la cuerda del criterio de distancia. Para ventanas
    mayores que 3 se aplica la recursion de ``smoothEnds`` de R.
    """
    lam = np.asarray(lambda_grid, dtype=float)
    y_raw = np.asarray(q_spike,   dtype=float)
    if lam.shape != y_raw.shape:
        raise ValueError("lambda_grid y q_spike deben tener la misma forma.")

    # Ordenar por lambda crecientes (defensa)
    o = np.argsort(lam)
    lam = lam[o]; y_raw = y_raw[o]

    if smooth_k and smooth_k > 1:
        k = smooth_k if smooth_k % 2 == 1 else smooth_k + 1
        y = _runmed_r(y_raw, k)
    else:
        y = y_raw

    x = np.log10(lam) if use_log10 else lam

    # Distancia perpendicular del punto (x_i, y_i) a la recta P1-PN
    x1, y1 = x[0], y[0]
    xN, yN = x[-1], y[-1]
    a = yN - y1
    b = -(xN - x1)
    c = (xN - x1) * y1 - (yN - y1) * x1
    denom = math.sqrt(a * a + b * b)
    if not (denom > 0):
        raise ValueError("Cuerda degenerada en deteccion de codo.")
    d = np.abs(a * x + b * y + c) / denom
    d[0] = -np.inf
    d[-1] = -np.inf
    return float(lam[int(np.argmax(d))])


def _simulate_spikes_vectorized(
    Xnull:      np.ndarray,
    weights:    np.ndarray,
    grid:       FMMGrid,
    cache:      FMMCache,
    lambda_:    float,
    omega0:     float,
    penalty_power: float,
) -> np.ndarray:
    """Simula B fits + spike_max_jump de una sola pasada, totalmente vectorizado.

    Replica B replicas nulas en bucle Python tarda ~4 ms cada una; esta version
    hace todas las B simultaneamente usando productos matriciales y argmax
    sobre tensores ``(B, nG)``, eliminando completamente el bucle.

    Shape conventions:
        Xnull       (B, n)
        weights     (n,)
        grid.C/S    (nG, n)
        cache.inv_XWX (nG, 3, 3)

    Returns
    -------
    spikes : np.ndarray, shape (B,)
        Diagnostico spike_max_jump para cada simulacion.
    """
    B = Xnull.shape[0]
    sw = float(weights.sum())

    # 1. Centrar las B replicas (sustraer media ponderada)
    mu_b = (Xnull @ weights) / sw                    # (B,)
    wx0  = weights[None, :] * (Xnull - mu_b[:, None])  # (B, n)

    # 2. Potencia: U, V para todas las simulaciones y todos los puntos de rejilla
    U = wx0 @ grid.C.T                               # (B, nG)
    V = wx0 @ grid.S.T                               # (B, nG)

    # 3. Denominador y rugosidad — solo dependen de la rejilla, son constantes B
    den = grid.C2 @ weights + grid.S2 @ weights      # (nG,)
    den_safe = np.where(den > 0, den, np.nan)
    pow_ = (U * U + V * V) / den_safe[None, :]       # (B, nG)

    if lambda_ > 0:
        g_omega = (omega0 / (grid.omega_vec + omega0)) ** penalty_power
        score = pow_ - lambda_ * g_omega[None, :] * grid.roughness[None, :]
    else:
        score = pow_

    # 4. Argmax de J por replica. Misma regla de seleccion que ``fit_gene``:
    # la calibracion tiene que reproducir exactamente el ajuste que se usara
    # despues, o el lambda* calibrado no corresponderia al procedimiento real.
    score_safe = np.where(np.isfinite(score), score, -np.inf)
    idx_b = np.argmax(score_safe, axis=1)            # (B,)

    # 5. WLS final usando cache: coef_b = inv_XWX[idx_b] @ b_vec_b
    C_at = grid.C[idx_b]                             # (B, n)
    S_at = grid.S[idx_b]                             # (B, n)
    b_vec = np.stack([
        (Xnull * weights).sum(axis=1),
        (Xnull * weights * C_at).sum(axis=1),
        (Xnull * weights * S_at).sum(axis=1),
    ], axis=1)                                        # (B, 3)
    inv_at = cache.inv_XWX[idx_b]                    # (B, 3, 3)
    coef   = np.einsum("bij,bj->bi", inv_at, b_vec)  # (B, 3)

    fitted = (
        coef[:, 0:1]
        + coef[:, 1:2] * C_at
        + coef[:, 2:3] * S_at
    )                                                 # (B, n)

    # 6. spike_max_jump vectorizado en orden circular
    fitted_ord = fitted[:, grid.order]
    diffs = np.abs(np.diff(fitted_ord, axis=1, append=fitted_ord[:, :1]))
    return diffs.max(axis=1)


def calibrate_lambda(
    times:         np.ndarray,
    weights:       np.ndarray,
    grid:          FMMGrid,
    cache:         FMMCache,
    lambda_grid:   Optional[np.ndarray] = None,
    n_simulations: int    = 1000,
    quantile:      float  = 0.975,
    sigma:         float  = 1.0,
    omega0:        float  = 0.10,
    penalty_power: float  = 3.0,
    seed:          int    = 1,
    smooth_k:      int    = 3,
    center_null:   bool   = True,
) -> dict:
    """Calibracion de lambda por Monte Carlo bajo H0 + codo.

    Para cada lambda candidato simula B replicas bajo H0 (Common Random
    Numbers) y mide el cuantil ``q_spike`` del diagnostico de salto maximo.
    El ``lambda_star`` es el codo de la curva ``q_spike(lambda)``.

    Implementacion vectorizada: las B simulaciones para un lambda dado se
    procesan en NumPy puro mediante un solo conjunto de productos matriciales
    sobre tensores ``(B, nG)``, sin bucle Python. Esto es ~30x mas rapido que
    paralelizar el bucle con threading porque elimina completamente el
    overhead per-tarea de joblib (cada simulacion individual cuesta ~1 ms).

    Devuelve dict con ``lambda_star``, ``lambda_grid``, ``q_spike``.
    """
    if lambda_grid is None:
        # Rejilla logaritmica equivalente a la del codigo R (24 valores en [2^-1, 2^4])
        lambda_grid = 2.0 ** np.linspace(-1, 4, 24)
    lambda_grid = np.asarray(lambda_grid, dtype=float)

    rng = np.random.default_rng(seed)
    n = times.shape[0]

    # CRN: las B replicas se generan UNA vez y se reusan en todos los lambdas
    Xnull = sigma * rng.standard_normal((n_simulations, n))
    if center_null:
        sw = float(weights.sum())
        mu = (Xnull * weights).sum(axis=1) / sw
        Xnull = Xnull - mu[:, np.newaxis]

    q_spike = np.empty(lambda_grid.shape[0])
    for i, lam in enumerate(lambda_grid):
        spikes = _simulate_spikes_vectorized(
            Xnull, weights, grid, cache,
            float(lam), omega0, penalty_power,
        )
        q_spike[i] = float(np.quantile(spikes, quantile))

    lambda_star = _pick_lambda_elbow(lambda_grid, q_spike, smooth_k=smooth_k)

    return dict(
        lambda_star = lambda_star,
        lambda_grid = lambda_grid,
        q_spike     = q_spike,
    )


def calibrate_lambda_repeat(
    n_repeats: int   = 5,
    seed0:     int   = 1,
    agg:       str   = "median",
    smooth_k:  int   = 3,
    **kwargs,
) -> float:
    """Repite calibracion ``n_repeats`` veces y agrega q_spike via mediana/media.

    Reduce la varianza Monte Carlo del Q-cuantil. Cada repeticion usa una
    semilla distinta (CRN dentro de cada repeticion). La curva final
    ``q_spike_agg(lambda)`` es la mediana (default) o media entre repeticiones.
    """
    if agg not in ("median", "mean"):
        raise ValueError("agg debe ser 'median' o 'mean'.")

    Q = []
    lambda_grid = None
    for r in range(n_repeats):
        out = calibrate_lambda(seed=seed0 + r, smooth_k=smooth_k, **kwargs)
        Q.append(out['q_spike'])
        lambda_grid = out['lambda_grid']

    Q = np.array(Q)
    Q_agg = Q.mean(axis=0) if agg == "mean" else np.median(Q, axis=0)
    return _pick_lambda_elbow(lambda_grid, Q_agg, smooth_k=smooth_k)


# ═══════════════════════════════════════════════════════════════════════════
# Seccion 7 — Funcion publica de preparacion
# ═══════════════════════════════════════════════════════════════════════════

def prepare_dataset(
    times:             np.ndarray,
    weights:           Optional[np.ndarray] = None,
    lambda_star:       Optional[float]      = None,
    length_alpha_grid: int   = 48,
    length_omega_grid: int   = 50,
    omega_min:         float = 0.02,
    omega_max:         float = 1.0,
    omega0:            float = 0.10,
    penalty_power:     float = 3.0,
    n_dense:           int   = 512,
    n_simulations:     int   = 1000,
    n_repeats:         int   = 5,
    quantile:          float = 0.975,
    seed:              int   = 1,
) -> FMMPreparedDataset:
    """Preparacion completa del dataset para ajuste con CircaFMM.

    Llamar UNA SOLA VEZ por dataset. El ``FMMPreparedDataset`` resultante
    se reutiliza para todos los genes posteriores y para cualquier test
    estadistico (LRT, etc.).

    Pasos:

    1. Si ``weights`` es None: pesos de densidad temporal von Mises
       (kappa seleccionado por LOO-CV).
    2. Construir rejillas alpha/omega (defaults del pipeline publicado
       ``runCircaFMM.R``: 48 alphas y 50 omegas equiespaciados en [0.02, 1]).
    3. ``precompute_grid`` (incluye la rugosidad densa).
    4. ``prepare_cache`` con los pesos.
    5. Si ``lambda_star`` es None: ``calibrate_lambda_repeat`` con
       ``n_repeats`` repeticiones agregadas por mediana.

    Parameters
    ----------
    times : np.ndarray, shape (n,)
        Tiempos circulares en radianes (se envolvera a [0, 2*pi)).
    weights : np.ndarray opcional
        Si se proporciona, se usa directamente (normalizado a sum=n).
        Si es None, se calcula por densidad temporal von Mises.
    lambda_star : float opcional
        Si se proporciona, salta la calibracion Monte Carlo.
    length_alpha_grid, length_omega_grid : int
        Tamano de las rejillas alpha y omega.
    omega_min, omega_max : float
        Extremos de la rejilla de omega, equiespaciada (defaults 0.02 y 1.0).
    omega0, penalty_power : float
        Parametros de la funcion de penalizacion ``g(omega)``.
    n_dense : int
        Tamano de la rejilla densa para rugosidad estable.
    n_simulations, n_repeats, quantile, seed : float/int
        Parametros del Monte Carlo de calibracion (solo si lambda_star=None).
        Los valores por defecto son los del articulo: B = 1000 replicas nulas,
        R = 5 repeticiones y cuantil q = 0,975.
    """
    times = np.ascontiguousarray(times, dtype=np.float64) % (2.0 * np.pi)

    # 1. Pesos
    if weights is None:
        weights = _time_density_weights(times)
    else:
        weights = _normalize_weights_to_n(weights)

    # 2-3. Rejillas + grid + rugosidad densa
    alpha_grid = np.linspace(0.0, 2.0 * np.pi, length_alpha_grid, endpoint=False)
    omega_grid = np.linspace(omega_min, omega_max, length_omega_grid)

    grid = precompute_grid(times, alpha_grid=alpha_grid, omega_grid=omega_grid,
                           n_dense=n_dense)

    # 4. Cache WLS
    cache = prepare_cache(weights, grid)

    # 5. Calibracion de lambda (si no se proporciona)
    if lambda_star is None:
        lambda_star = calibrate_lambda_repeat(
            n_repeats     = n_repeats,
            seed0         = seed,
            agg           = "median",
            times         = times,
            weights       = weights,
            grid          = grid,
            cache         = cache,
            n_simulations = n_simulations,
            quantile      = quantile,
            sigma         = 1.0,
            omega0        = omega0,
            penalty_power = penalty_power,
        )

    penalty = _penalty_vector(
        grid, float(lambda_star), float(omega0), float(penalty_power),
    )

    return FMMPreparedDataset(
        times         = times,
        weights       = weights,
        grid          = grid,
        cache         = cache,
        lambda_star   = float(lambda_star),
        omega0        = float(omega0),
        penalty_power = float(penalty_power),
        penalty       = penalty,
    )


# ═══════════════════════════════════════════════════════════════════════════
# Seccion 8 — Clase modelo (interfaz RhythmModel)
# ═══════════════════════════════════════════════════════════════════════════

class CircaFMMModel(RhythmModel):
    """
    CircaFMM con criterio de potencia + penalizacion de rugosidad.

    A diferencia del FMM base (``FMMModel``), este modelo:

    - Usa pesos por densidad temporal von Mises
    - Selecciona (alpha, omega) por criterio de potencia penalizado
    - Calibra lambda automaticamente bajo H0
    - Es sustancialmente mas rapido en batch (preparacion compartida)

    El ``FitResult`` resultante es **compatible** con el de ``FMMModel``:
    misma estructura de ``params`` (M, A, alpha, beta, omega) y misma formula
    de ``peak_time``. Es drop-in en el resto del pipeline.

    Modos de uso:

    1. **Eficiente (recomendado)**: el usuario crea el ``FMMPreparedDataset``
       una vez y lo pasa al constructor::

           prepared = prepare_dataset(times)
           model = CircaFMMModel(prepared)
           for gene in genes:
               result = model.fit(data[gene], times)

    2. **Lazy (drop-in)**: la clase prepara automaticamente al primer fit
       y cachea el resultado mientras los tiempos no cambien::

           model = CircaFMMModel()
           for gene in genes:
               result = model.fit(data[gene], times)   # primera prep + cache

    Parameters
    ----------
    prepared : FMMPreparedDataset, opcional
        Si se proporciona, se usa directamente. Si es None, se prepara
        de forma lazy en el primer ``fit()``.
    refine : bool
        Si aplicar el refinamiento local del Algoritmo S2 tras la seleccion
        sobre la rejilla. Por defecto ``False`` (ver :func:`fit_gene`).
    refine_steps : int
        Numero ``S`` de pasos del refinamiento. Por defecto 3, el del articulo.
    **prep_kwargs
        Parametros extra para ``prepare_dataset`` (solo en modo lazy).
    """

    def __init__(
        self,
        prepared:     Optional[FMMPreparedDataset] = None,
        refine:       bool = False,
        refine_steps: int  = 3,
        **prep_kwargs,
    ) -> None:
        self.prepared     = prepared
        self.refine       = refine
        self.refine_steps = refine_steps
        self._prep_kwargs = prep_kwargs
        # La preparacion lazy se protege con un lock: en los caminos con hilos
        # (TopGeneSelector y RobustOrderEstimator comparten una sola instancia
        # entre workers) sin el, cada hilo dispararia su propia calibracion de
        # lambda antes de que el primero terminase de cachear el resultado.
        self._prep_lock = threading.Lock()

    def _needs_prepare(self, time_points: np.ndarray) -> bool:
        """¿Hay que (re)preparar el dataset para estos ``time_points``?

        Solo si no hay ninguno preparado o si la escala circular ha cambiado.
        La comparacion se hace contra ``prepared.times``, que ya esta envuelto
        a ``[0, 2*pi)``, de modo que un ``prepared`` recibido en el constructor
        se reutiliza tal cual en lugar de recalibrarse.
        """
        if self.prepared is None:
            return True
        t = np.asarray(time_points, dtype=np.float64) % (2.0 * np.pi)
        ref = self.prepared.times
        return t.shape != ref.shape or not np.allclose(t, ref, rtol=0.0, atol=1e-12)

    @staticmethod
    def peak_time(alpha: float, beta: float, omega: float) -> float:
        """
        Tiempo de pico — formula identica a ``FMMModel.peak_time``.

        t_U = alpha + 2*atan2((1/omega)*sin(-beta/2), cos(-beta/2))  mod 2*pi
        """
        return float(
            (alpha + 2.0 * np.arctan2(
                (1.0 / omega) * np.sin(-beta / 2.0),
                np.cos(-beta / 2.0),
            )) % (2.0 * np.pi)
        )

    def fit(
        self,
        data:        np.ndarray,
        time_points: np.ndarray,
    ) -> FitResult:
        """
        Ajusta CircaFMM a ``data``.

        Si no se ha proporcionado ``prepared`` en el constructor, se prepara
        en la primera llamada y se cachea. Si los ``time_points`` cambian
        respecto a la preparacion previa, se vuelve a preparar.
        """
        # Preparar (lazy) solo si no hay ``prepared`` o cambian los tiempos
        if self._needs_prepare(time_points):
            with self._prep_lock:
                if self._needs_prepare(time_points):      # doble comprobacion
                    self.prepared = prepare_dataset(time_points, **self._prep_kwargs)

        return fit_gene(
            data         = data,
            prepared     = self.prepared,
            refine       = self.refine,
            refine_steps = self.refine_steps,
        )


# ═══════════════════════════════════════════════════════════════════════════
# Factory simple para selector base/CircaFMM en el pipeline
# ═══════════════════════════════════════════════════════════════════════════

def make_fmm_model(method: str = "base", **kwargs) -> RhythmModel:
    """
    Factory para seleccionar entre FMM base y CircaFMM.

    Parameters
    ----------
    method : str
        - ``"base"``       : retorna ``FMMModel`` (fmm.py, RSS + Nelder-Mead)
        - ``"circafmm"`` : retorna ``CircaFMMModel`` (circafmm.py)
    **kwargs
        Argumentos del constructor del modelo elegido.
    """
    if method == "base":
        from circust.fitting.fmm import FMMModel
        return FMMModel(**kwargs)
    elif method == "circafmm":
        return CircaFMMModel(**kwargs)
    else:
        raise ValueError(
            f"Unknown fmm_method='{method}'. Use 'base' or 'circafmm'."
        )
