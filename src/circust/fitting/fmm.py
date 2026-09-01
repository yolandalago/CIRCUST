"""
circust/fitting/fmm.py
=======================
Frecuency Modulated Mobius (FMM) — ajuste ritmico de un solo componente.

Modelo matematico
-----------------
    y(t) = M + A * cos(beta + 2*atan(omega*tan((t - alpha)/2)))

El termino interno  Phi(t) = beta + 2*atan(omega*tan((t - alpha)/2))  es una
transformacion de Mobius del circulo.  Cuando omega = 1 se reduce a un
simple desfase. Cuando omega -> 1 la forma de onda se aproxima a un coseno
(Cosinor). La ventaja clave sobre Cosinor es que omega controla la *asimetria
de la forma de onda*: el modelo puede representar picos asimetricos presentes
en datos circadianos reales.

Parametros
----------
M     mesor       — nivel medio
A     amplitud    — siempre > 0
alpha             — localizacion de la region del pico en [0, 2*pi)
beta              — parametro de forma en [0, 2*pi)
omega             — asimetria en (0, omegaMax]

Algoritmo de estimacion
-----------------------

Paso 1 — Busqueda en rejilla
    Para cada combinacion (alpha, omega) en una rejilla, se fijan esos dos
    parametros y se estiman M, A, beta por OLS sobre el eje temporal
    transformado por Mobius. Se selecciona el candidato con menor RSS que
    tambien satisfaga las restricciones de estabilidad de amplitud.

Paso 2 — Refinamiento Nelder-Mead
    Se pule la estimacion del Paso 1 minimizando el RSS sobre los cinco
    parametros simultaneamente, sujeto a las mismas restricciones de
    estabilidad (puntos infactibles devuelven infinito).

Los Pasos 1+2 se repiten ``num_reps`` veces, cada vez estrechando la
rejilla alrededor de la mejor estimacion anterior.

Devuelve
--------
FitResult
    Claves de ``params``: ``M``, ``A``, ``alpha``, ``beta``, ``omega``
"""
import numpy as np
from scipy.optimize import minimize

from circust.fitting.rhythm_model import FitResult, RhythmModel

# ── Aceleracion opcional con Numba (CPU JIT) ───────────────────────────────
# Si Numba esta disponible, compilamos la evaluacion del modelo y la funcion
# objetivo que Nelder-Mead llama miles de veces por gen.
try:
    from numba import njit
    _HAS_NUMBA = True
except ImportError:
    _HAS_NUMBA = False
    def njit(*args, **kwargs):                              # type: ignore
        def deco(fn): return fn
        if len(args) == 1 and callable(args[0]):
            return args[0]
        return deco

# ── Aceleracion opcional con Numba CUDA (GPU) ──────────────────────────────
# Opcion A: el Paso 1 (busqueda en rejilla) se ejecuta en GPU si esta
# disponible.  Para datasets grandes (muchas muestras o rejillas densas) el
# grid search es el cuello de botella principal — se presta bien a GPU porque
# cada punto (alpha, omega) es independiente.
#
# Opcion B (trabajo futuro): paralelizar el ajuste de N genes en GPU
# requeriria reimplementar Nelder-Mead en CUDA, lo que esta fuera del alcance
# actual del TFG.  Ver comentario en ``_step1_grid_cuda`` para detalles.
try:
    from numba import cuda as _numba_cuda
    _HAS_CUDA = _numba_cuda.is_available()
except (ImportError, Exception):
    _HAS_CUDA = False


@njit(cache=True, fastmath=True)
def _fmm_predict_jit(
    t: np.ndarray,
    M: float, A: float, alpha: float, beta: float, omega: float,
    out: np.ndarray,
) -> None:
    """Escribe la forma de onda FMM en *out* (sin asignacion)."""
    n = t.shape[0]
    for i in range(n):
        mob  = 2.0 * np.arctan(omega * np.tan((t[i] - alpha) / 2.0))
        out[i] = M + A * np.cos(beta + mob)


@njit(cache=True, fastmath=True)
def _fmm_objective_jit(
    params: np.ndarray,
    data:   np.ndarray,
    t:      np.ndarray,
    upper:  float,
    lower:  float,
    omega_max: float,
    buf:    np.ndarray,
) -> float:
    """
    Objetivo Nelder-Mead compilado. Devuelve +inf cuando los parametros
    violan las restricciones de estabilidad de amplitud o rangos validos.
    """
    M     = params[0]
    A     = params[1]
    alpha = params[2]
    beta  = params[3]
    omega = params[4]

    if (M + A) > upper: return np.inf
    if (M - A) < lower: return np.inf
    if A <= 0.0:        return np.inf
    if omega <= 0.0:    return np.inf
    if omega > omega_max: return np.inf

    n = t.shape[0]
    rss = 0.0
    for i in range(n):
        mob  = 2.0 * np.arctan(omega * np.tan((t[i] - alpha) / 2.0))
        pred = M + A * np.cos(beta + mob)
        d    = pred - data[i]
        rss += d * d
    return rss / n


# ---------------------------------------------------------------------------
# Funciones internas auxiliares (replican funciones ocultas de R)
# ---------------------------------------------------------------------------



def _mobius(t: np.ndarray, alpha: float, omega: float) -> np.ndarray:
    """Transformacion de Mobius del eje temporal.

    Phi(t) = 2*atan(omega*tan((t - alpha)/2))
    """
    return 2.0 * np.arctan(omega * np.tan((t - alpha) / 2.0))


def _circular_abs_diff(a: float, b: float) -> float:
    """Distancia angular minima entre dos angulos en radianes (resultado en [0, pi]).

    Necesario para early-stopping del bucle ``num_reps``: si entre iteraciones
    alpha pasa de 6.27 a 0.01 (cruzando 2*pi), la diferencia "lineal" es 6.26
    pero la distancia circular real es solo 0.02. Sin esto el early stopping
    fallaria justo cuando deberia activarse.
    """
    d = abs(a - b) % (2.0 * np.pi)
    return float(min(d, 2.0 * np.pi - d))


def _fmm_predict(
    t: np.ndarray,
    M: float, A: float, alpha: float, beta: float, omega: float,
) -> np.ndarray:
    """Evalua la forma de onda FMM en los puntos temporales ``t``."""
    phi = beta + _mobius(t, alpha, omega)
    return M + A * np.cos(phi)


def _step1(
    alpha: float,
    omega: float,
    data: np.ndarray,
    t: np.ndarray,
) -> np.ndarray:
    """
    Paso 1 (escalar): fija alpha y omega, estima M, A, beta por OLS Cosinor.

    Devuelve [M, A, alpha, beta, omega, RSS] — mismo layout que step1FMM de R.
    Se usa solo como fallback; la version vectorizada es preferible.
    """
    mob      = _mobius(t, alpha, omega)
    t_star   = alpha + mob
    xx       = np.cos(t_star)
    zz       = np.sin(t_star)
    X        = np.column_stack([np.ones(len(data)), xx, zz])
    coeffs, *_ = np.linalg.lstsq(X, data, rcond=None)
    M_est    = coeffs[0]
    b        = coeffs[1]
    g        = coeffs[2]
    phi_est  = np.arctan2(-g, b)
    A_est    = np.sqrt(b**2 + g**2)
    beta_est = (phi_est + alpha) % (2 * np.pi)

    fitted   = M_est + A_est * np.cos(beta_est + mob)
    rss      = np.sum((data - fitted)**2) / len(t)
    return np.array([M_est, A_est, alpha, beta_est, omega, rss])


@njit(cache=True, fastmath=True)
def _step1_grid_jit(
    alpha_grid: np.ndarray,
    omega_grid: np.ndarray,
    data:       np.ndarray,
    t:          np.ndarray,
    out:        np.ndarray,
) -> None:
    """Paso 1 vectorizado, compilado y sin tensores intermedios.

    Escribe en ``out`` una fila ``[M, A, alpha, beta, omega, RSS]`` por punto
    de la rejilla, con el mismo layout que ``_step1`` escalar.

    Dos ideas evitan casi todo el trabajo de la version con NumPy puro.

    La primera es que ``tan((t - alpha)/2)`` depende solo de ``(alpha, t)`` y no
    de ``omega``, de modo que basta con calcular ``sin`` y ``cos`` de
    ``u = (t - alpha)/2`` una vez por cada ``alpha``, en lugar de una vez por
    cada par ``(alpha, omega)``. La segunda es que la fase de Mobius no hace
    falta evaluarla: por la identidad del angulo doble, con
    ``phi = 2 arctan(omega tan u)`` se tiene

        cos(phi) = (cos^2 u - omega^2 sin^2 u) / (cos^2 u + omega^2 sin^2 u),
        sin(phi) = 2 omega sin u cos u       / (cos^2 u + omega^2 sin^2 u),

    de manera que los regresores ``cos(alpha + phi)`` y ``sin(alpha + phi)`` se
    obtienen con las formulas de adicion y sin una sola llamada a ``tan``,
    ``arctan``, ``cos`` o ``sin`` dentro del bucle sobre ``omega``. El
    denominador es una suma de cuadrados, siempre positivo, asi que ademas
    desaparece la singularidad en ``t = alpha + pi``.

    Los valores ajustados se obtienen como ``M + b*xx + g*zz``, que es
    algebraicamente identico a ``M + A cos(beta + phi)`` por construccion de la
    reparametrizacion lineal, y ahorra otra pasada trigonometrica.
    """
    An = alpha_grid.shape[0]
    Wn = omega_grid.shape[0]
    N  = t.shape[0]
    two_pi = 2.0 * np.pi

    ybar = 0.0
    for i in range(N):
        ybar += data[i]
    ybar /= N

    s_u  = np.empty(N)
    c_u  = np.empty(N)
    xbuf = np.empty(N)
    zbuf = np.empty(N)

    for ia in range(An):
        alpha = alpha_grid[ia]
        ca = np.cos(alpha)
        sa = np.sin(alpha)
        # Un unico seno y coseno por (alpha, muestra), reutilizados en todos los omega
        for i in range(N):
            u = 0.5 * (t[i] - alpha)
            s_u[i] = np.sin(u)
            c_u[i] = np.cos(u)

        for iw in range(Wn):
            omega = omega_grid[iw]
            o2 = omega * omega

            sx = 0.0
            sz = 0.0
            for i in range(N):
                s = s_u[i]
                c = c_u[i]
                cc = c * c
                ss = s * s
                den = cc + o2 * ss
                Cm  = (cc - o2 * ss) / den
                Sm  = (2.0 * omega * s * c) / den
                x = ca * Cm - sa * Sm          # cos(alpha + phi)
                z = sa * Cm + ca * Sm          # sin(alpha + phi)
                xbuf[i] = x
                zbuf[i] = z
                sx += x
                sz += z
            xbar = sx / N
            zbar = sz / N

            sxx = 0.0; sxz = 0.0; szz = 0.0; sxy = 0.0; szy = 0.0
            for i in range(N):
                xc = xbuf[i] - xbar
                zc = zbuf[i] - zbar
                yc = data[i] - ybar
                sxx += xc * xc
                sxz += xc * zc
                szz += zc * zc
                sxy += xc * yc
                szy += zc * yc

            det = sxx * szz - sxz * sxz
            if abs(det) < 1e-12:
                det = 1e-12
            b = (szz * sxy - sxz * szy) / det
            g = (sxx * szy - sxz * sxy) / det
            M = ybar - b * xbar - g * zbar

            rss = 0.0
            for i in range(N):
                d = M + b * xbuf[i] + g * zbuf[i] - data[i]
                rss += d * d
            rss /= N

            A_est = np.sqrt(b * b + g * g)
            beta  = (np.arctan2(-g, b) + alpha) % two_pi

            row = ia * Wn + iw
            out[row, 0] = M
            out[row, 1] = A_est
            out[row, 2] = alpha
            out[row, 3] = beta
            out[row, 4] = omega
            out[row, 5] = rss


def _step1_grid(
    alpha_grid: np.ndarray,
    omega_grid: np.ndarray,
    data: np.ndarray,
    t: np.ndarray,
) -> np.ndarray:
    """
    Paso 1 vectorizado: evalua TODOS los puntos (alpha, omega) simultaneamente.

    Si Numba esta disponible se delega en ``_step1_grid_jit``, que evita tanto
    los tensores intermedios de forma ``(A, W, N)`` como las llamadas
    trigonometricas repetidas para cada omega. Si no, se usa la version con
    NumPy puro que hay a continuacion, de semantica identica.

    Estrategia
    ----------
    Para cada par (alpha, omega) el problema OLS se reduce a un sistema de
    ecuaciones normales 2x2 sobre los regresores centrados [cos(t*), sin(t*)].
    El centrado hace que el sistema este bien condicionado y evita las matrices
    3x3 casi singulares que aparecen cuando omega esta cerca de cero. El
    intercepto M se recupera analiticamente de la ecuacion de medias.

    Parameters
    ----------
    alpha_grid : array (A,) de valores alpha
    omega_grid : array (W,) de valores omega
    data       : array (N,) senal observada
    t          : array (N,) puntos temporales en [0, 2*pi)

    Returns
    -------
    array (A*W, 6) — columnas: [M, A, alpha, beta, omega, RSS]
    Mismo layout de columnas que ``_step1`` escalar.
    """
    An, Wn, N = len(alpha_grid), len(omega_grid), len(t)
    AW = An * Wn

    if _HAS_NUMBA:
        out = np.empty((AW, 6), dtype=np.float64)
        _step1_grid_jit(
            np.ascontiguousarray(alpha_grid, dtype=np.float64),
            np.ascontiguousarray(omega_grid, dtype=np.float64),
            np.ascontiguousarray(data,       dtype=np.float64),
            np.ascontiguousarray(t,          dtype=np.float64),
            out,
        )
        return out

    # ── Broadcasting de toda la trigonometria: formas (A, W, N) ─────────
    alpha_b = alpha_grid[:, None, None]   # (A, 1, 1)
    omega_b = omega_grid[None, :, None]   # (1, W, 1)
    t_b     = t[None, None, :]            # (1, 1, N)

    mob    = 2.0 * np.arctan(omega_b * np.tan((t_b - alpha_b) / 2.0))  # (A,W,N)
    t_star = alpha_b + mob                                               # (A,W,N)
    xx     = np.cos(t_star).reshape(AW, N)   # (AW, N)
    zz     = np.sin(t_star).reshape(AW, N)   # (AW, N)
    mob_2d = mob.reshape(AW, N)              # (AW, N)

    # ── OLS 2x2 centrado — evita mal condicionamiento cerca de omega=0 ──
    xx_c = xx - xx.mean(axis=1, keepdims=True)   # (AW, N)
    zz_c = zz - zz.mean(axis=1, keepdims=True)   # (AW, N)
    y_c  = data - data.mean()                     # (N,) broadcast escalar

    sxx = (xx_c ** 2).sum(1)          # (AW,)
    sxz = (xx_c * zz_c).sum(1)       # (AW,)
    szz = (zz_c ** 2).sum(1)         # (AW,)
    sxy = (xx_c * y_c).sum(1)        # (AW,)
    szy = (zz_c * y_c).sum(1)        # (AW,)

    det = sxx * szz - sxz ** 2
    det = np.where(np.abs(det) < 1e-12, 1e-12, det)   # proteger casos singulares

    b     = (szz * sxy - sxz * szy) / det    # coeficiente coseno
    g     = (sxx * szy - sxz * sxy) / det    # coeficiente seno
    M_est = data.mean() - b * xx.mean(1) - g * zz.mean(1)

    # ── Recuperacion de parametros ───────────────────────────────────────
    A_est   = np.sqrt(b ** 2 + g ** 2)
    phi_est = np.arctan2(-g, b)

    alpha_flat = np.repeat(alpha_grid, Wn)     # (AW,)
    omega_flat = np.tile(omega_grid, An)       # (AW,)
    beta_est   = (phi_est + alpha_flat) % (2 * np.pi)

    # ── RSS ──────────────────────────────────────────────────────────────
    fitted = M_est[:, None] + A_est[:, None] * np.cos(beta_est[:, None] + mob_2d)
    rss    = np.sum((data[None, :] - fitted) ** 2, axis=1) / N

    return np.column_stack([M_est, A_est, alpha_flat, beta_est, omega_flat, rss])


def _step1_grid_cuda(
    alpha_grid: np.ndarray,
    omega_grid: np.ndarray,
    data: np.ndarray,
    t: np.ndarray,
) -> np.ndarray:
    """
    Paso 1 vectorizado con Numba CUDA — Opcion A de aceleracion GPU.

    Cada hilo GPU procesa un punto (alpha, omega) de la rejilla de forma
    independiente, calculando los coeficientes OLS y el RSS correspondiente.

    La GPU es especialmente util cuando ``len(alpha_grid) * len(omega_grid)``
    es grande (>= varios miles de combinaciones) o cuando ``len(t)`` es
    grande (>= cientos de muestras).

    Requisito: ``numba.cuda.is_available() == True``.

    Opcion B — trabajo futuro
    -------------------------
    Paralelizar el ajuste de N genes simultaneamente en GPU requeriria
    reimplementar el optimizador Nelder-Mead como kernel CUDA, ya que
    ``scipy.optimize.minimize`` no puede ejecutarse en GPU. Las referencias
    de implementacion son:
      - Luersen & Le Riche (2004), «Globalized Nelder-Mead method for
        engineering optimization», Computers & Structures 82(23-26), 2253-2260.
      - Implementaciones existentes de CUDA NM: Vehtari et al. / PyMC3 GPU.
    Esta opcion queda documentada como trabajo futuro una vez el pipeline
    principal este validado.

    Returns
    -------
    array (A*W, 6) — mismo layout que ``_step1_grid`` (CPU).
    """
    from numba import cuda as _cuda
    import math

    An, Wn, N = len(alpha_grid), len(omega_grid), len(t)
    AW = An * Wn

    alpha_flat = np.repeat(alpha_grid, Wn).astype(np.float64)
    omega_flat = np.tile(omega_grid, An).astype(np.float64)
    data_d     = _cuda.to_device(data.astype(np.float64))
    t_d        = _cuda.to_device(t.astype(np.float64))
    alpha_d    = _cuda.to_device(alpha_flat)
    omega_d    = _cuda.to_device(omega_flat)
    out_d      = _cuda.device_array((AW, 6), dtype=np.float64)

    @_cuda.jit
    def _grid_kernel(alpha_arr, omega_arr, data_arr, t_arr, out, n_samples):
        idx = _cuda.grid(1)
        if idx >= alpha_arr.shape[0]:
            return
        alpha = alpha_arr[idx]
        omega = omega_arr[idx]
        N     = n_samples

        # Acumula sumas para OLS centrado
        sum_xx = 0.0; sum_zz = 0.0; sum_xz = 0.0
        sum_xy = 0.0; sum_zy = 0.0
        sum_x  = 0.0; sum_z  = 0.0
        data_mean = 0.0
        for i in range(N):
            data_mean += data_arr[i]
        data_mean /= N

        for i in range(N):
            mob    = 2.0 * math.atan(omega * math.tan((t_arr[i] - alpha) / 2.0))
            t_star = alpha + mob
            x_i    = math.cos(t_star)
            z_i    = math.sin(t_star)
            sum_x  += x_i
            sum_z  += z_i

        x_mean = sum_x / N
        z_mean = sum_z / N

        for i in range(N):
            mob    = 2.0 * math.atan(omega * math.tan((t_arr[i] - alpha) / 2.0))
            t_star = alpha + mob
            x_c    = math.cos(t_star) - x_mean
            z_c    = math.sin(t_star) - z_mean
            y_c    = data_arr[i] - data_mean
            sum_xx += x_c * x_c
            sum_zz += z_c * z_c
            sum_xz += x_c * z_c
            sum_xy += x_c * y_c
            sum_zy += z_c * y_c

        det = sum_xx * sum_zz - sum_xz * sum_xz
        if abs(det) < 1e-12:
            det = 1e-12

        b     = (sum_zz * sum_xy - sum_xz * sum_zy) / det
        g     = (sum_xx * sum_zy - sum_xz * sum_xy) / det

        # Calcular medias de x, z para recuperar M
        sum_x2 = 0.0; sum_z2 = 0.0
        for i in range(N):
            mob    = 2.0 * math.atan(omega * math.tan((t_arr[i] - alpha) / 2.0))
            t_star = alpha + mob
            sum_x2 += math.cos(t_star)
            sum_z2 += math.sin(t_star)
        M_est = data_mean - b * (sum_x2 / N) - g * (sum_z2 / N)
        A_est = math.sqrt(b * b + g * g)
        phi_est = math.atan2(-g, b)
        beta_est = math.fmod(phi_est + alpha, 2.0 * math.pi)
        if beta_est < 0.0:
            beta_est += 2.0 * math.pi

        # RSS
        rss = 0.0
        for i in range(N):
            mob    = 2.0 * math.atan(omega * math.tan((t_arr[i] - alpha) / 2.0))
            fitted = M_est + A_est * math.cos(beta_est + mob)
            d      = fitted - data_arr[i]
            rss   += d * d
        rss /= N

        out[idx, 0] = M_est
        out[idx, 1] = A_est
        out[idx, 2] = alpha
        out[idx, 3] = beta_est
        out[idx, 4] = omega
        out[idx, 5] = rss

    threads_per_block = 128
    blocks = (AW + threads_per_block - 1) // threads_per_block
    _grid_kernel[blocks, threads_per_block](
        alpha_d, omega_d, data_d, t_d, out_d, N
    )
    return out_d.copy_to_host()


def _best_step1(
    data: np.ndarray,
    grid_results: np.ndarray,
) -> np.ndarray | None:
    """
    Selecciona el candidato de la rejilla con menor RSS que satisfaga
    las restricciones de estabilidad de amplitud.

    Estabilidad: M+A <= max(datos) + 10% rango  Y  M-A >= min(datos) - 10% rango.
    """
    order    = np.argsort(grid_results[:, 5])   # ordenar por RSS
    data_min = data.min()
    data_max = data.max()
    slack    = 0.1 * (data_max - data_min)

    for idx in order:
        row   = grid_results[idx]
        M, A  = row[0], row[1]
        if np.isnan(M) or np.isnan(A):
            continue
        if (M + A) <= (data_max + slack) and (M - A) >= (data_min - slack):
            return row
    return None


def _make_step2_objective(
    data: np.ndarray,
    t: np.ndarray,
    omega_max: float,
) -> callable:
    """
    Construye la funcion objetivo para Nelder-Mead con limites precalculados.

    Evita recalcular min/max/margen en cada evaluacion (se llama miles de
    veces por ajuste).

    Si Numba esta disponible, la funcion objetivo es un wrapper delgado
    sobre un nucleo JIT-compilado. Si no, cae en la implementacion pura
    NumPy con semantica identica.
    """
    data_min = float(data.min())
    data_max = float(data.max())
    slack    = 0.1 * (data_max - data_min)
    upper    = data_max + slack
    lower    = data_min - slack
    n        = len(t)

    if _HAS_NUMBA:
        buf = np.empty(n, dtype=np.float64)
        def objective(params: np.ndarray) -> float:
            return float(_fmm_objective_jit(
                np.ascontiguousarray(params, dtype=np.float64),
                data, t, upper, lower, omega_max, buf,
            ))
        return objective

    def objective(params: np.ndarray) -> float:
        M, A, alpha, beta, omega = params
        if (M + A) > upper:    return np.inf
        if (M - A) < lower:    return np.inf
        if A <= 0:             return np.inf
        if omega <= 0:         return np.inf
        if omega > omega_max:  return np.inf
        fitted = _fmm_predict(t, M, A, alpha, beta, omega)
        return float(np.sum((fitted - data)**2) / n)

    return objective


# ---------------------------------------------------------------------------
# Clase FMMModel
# ---------------------------------------------------------------------------

class FMMModel(RhythmModel):
    """
    Ajustador ritmico FMM (Frecuency Modulated Mobius).

    Ajusta  y(t) = M + A*cos(beta + 2*atan(omega*tan((t - alpha)/2)))
    usando un algoritmo de dos pasos: busqueda en rejilla + Nelder-Mead.

    Parameters
    ----------
    length_alpha_grid : int
        Numero de valores de alpha en la rejilla. Valor R por defecto: 48.
    length_omega_grid : int
        Numero de valores de omega en la rejilla. Valor R por defecto: 24.
    omega_max : float
        Cota superior de omega. Valor R por defecto: 1.0.
    num_reps : int
        Numero de iteraciones de refinamiento de rejilla. Valor R por defecto: 3.

    Example
    -------
    >>> import numpy as np
    >>> from circust.fitting.fmm import FMMModel
    >>>
    >>> t    = np.linspace(0, 2*np.pi, 50, endpoint=False)
    >>> data = 0.3 + 0.6 * np.cos(t - 1.2) + np.random.normal(0, 0.05, 50)
    >>> result = FMMModel().fit(data, t)
    >>> print(result.summary())
    """

    def __init__(
        self,
        length_alpha_grid: int   = 48,
        length_omega_grid: int   = 24,
        omega_max:         float = 1.0,
        num_reps:          int   = 3,
        skip_nm_r2:        float = 0.95,
        early_stop_tol:    float = 1e-3,
    ) -> None:
        """
        Parameters
        ----------
        skip_nm_r2 : float
            Si tras el grid search del paso 1 el R² ya supera este umbral, se
            omite el refinamiento Nelder-Mead (rara vez aporta mas que decimales
            cuando el grid ya da buen ajuste). Por defecto 0.95. Ponerlo a >1.0
            (ej. 1.1) desactiva esta optimizacion (NM siempre se ejecuta).

        early_stop_tol : float
            Tolerancia para parar el bucle ``num_reps`` cuando los parametros
            convergen entre iteraciones consecutivas. Si todos los cambios
            absolutos en (alpha, beta, omega) son menores que esta tolerancia,
            se omite la siguiente iteracion. Por defecto 1e-3. Ponerlo a 0
            desactiva esta optimizacion.
        """
        self.length_alpha_grid = length_alpha_grid
        self.length_omega_grid = length_omega_grid
        self.omega_max         = omega_max
        self.num_reps          = num_reps
        self.skip_nm_r2        = skip_nm_r2
        self.early_stop_tol    = early_stop_tol

    # ------------------------------------------------------------------
    # API publica
    # ------------------------------------------------------------------

    @staticmethod
    def peak_time(alpha: float, beta: float, omega: float) -> float:
        """
        Tiempo de pico del modelo FMM — formula compUU de R.

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
        Ajusta el modelo FMM a ``data``.

        Parameters
        ----------
        data : np.ndarray, forma (n_muestras,)
            Valores de expresion observados (normalizados a [-1, 1]).
        time_points : np.ndarray, forma (n_muestras,)
            Eje temporal circular en [0, 2*pi), tipicamente ``CPCAResult.circular_scale``.

        Returns
        -------
        FitResult
            Claves de ``params``: ``M``, ``A``, ``alpha``, ``beta``, ``omega``
        """
        data        = np.asarray(data,        dtype=float)
        time_points = np.asarray(time_points, dtype=float)

        # Rejillas iniciales (mismos valores que R por defecto)
        alpha_grid = np.linspace(0, 2 * np.pi, self.length_alpha_grid, endpoint=False)
        omega_grid = np.exp(
            np.linspace(
                np.log(1e-4),
                np.log(self.omega_max),
                self.length_omega_grid,
            )
        )
        # Pasos de las rejillas iniciales. Se propagan de una iteracion a la
        # siguiente para poder estrechar las rejillas sin tener que deducir el
        # paso de una rejilla ya envuelta a [0, 2*pi) (ver _refine_alpha_grid).
        alpha_step = (2 * np.pi) / self.length_alpha_grid
        omega_step = (float(np.mean(np.diff(omega_grid)))
                      if self.length_omega_grid > 1 else self.omega_max)

        best_par  = None
        par_final = None
        prev_par  = None         # parametros de la iteracion anterior (early stopping)
        objective = _make_step2_objective(data, time_points, self.omega_max)
        ss_tot    = float(np.sum((data - data.mean()) ** 2))   # para R² del grid

        for rep in range(self.num_reps):

            # ── Paso 1: busqueda en rejilla ─────────────────────────────
            # Si hay GPU disponible (Opcion A), la rejilla se evalua en CUDA
            # con un hilo por combinacion (alpha, omega).  Si no, se usa la
            # version CPU vectorizada con NumPy (~3x vs bucle for puro).
            if _HAS_CUDA:
                grid_results = _step1_grid_cuda(
                    alpha_grid, omega_grid, data, time_points
                )
            else:
                grid_results = _step1_grid(alpha_grid, omega_grid, data, time_points)

            prev_best = best_par
            best_par  = _best_step1(data, grid_results)

            if best_par is None:
                # Sin candidato estable — revertir al anterior si existe
                best_par = prev_best
                break

            # ── Optimizacion 1: skip Nelder-Mead si el grid ya da buen R²
            # Calculamos el R² del mejor punto del grid; si supera el umbral,
            # NM aporta poco mas que ruido y nos lo saltamos.
            M_g, A_g, alpha_g, beta_g, omega_g = best_par[:5]
            fitted_g = _fmm_predict(time_points, M_g, A_g, alpha_g, beta_g, omega_g)
            r2_grid  = 1.0 - float(np.sum((data - fitted_g) ** 2)) / ss_tot if ss_tot > 0 else 0.0

            if r2_grid >= self.skip_nm_r2:
                par_final = best_par[:5].copy()
            else:
                # ── Paso 2: refinamiento Nelder-Mead ────────────────────
                result = minimize(
                    objective,
                    x0      = best_par[:5],
                    method  = "Nelder-Mead",
                    options = {"xatol": 1e-6, "fatol": 1e-6, "maxiter": 5000},
                )
                par_final = result.x.copy()

            # Envolver alpha y beta en [0, 2*pi)
            par_final[2] = par_final[2] % (2 * np.pi)
            par_final[3] = par_final[3] % (2 * np.pi)

            # ── Optimizacion 2: early stopping si los parametros convergen
            # Si entre la iteracion anterior y la actual los cambios en
            # (alpha, beta, omega) son menores que la tolerancia, parar.
            # alpha y beta son angulares, asi que comparamos su distancia
            # circular minima.
            if self.early_stop_tol > 0 and prev_par is not None:
                d_alpha = _circular_abs_diff(par_final[2], prev_par[2])
                d_beta  = _circular_abs_diff(par_final[3], prev_par[3])
                d_omega = abs(par_final[4] - prev_par[4])
                if d_alpha < self.early_stop_tol and d_beta < self.early_stop_tol \
                        and d_omega < self.early_stop_tol:
                    break
            prev_par = par_final.copy()

            # Estrechar rejillas alrededor del mejor actual para siguiente iteracion
            if rep < self.num_reps - 1:
                alpha_grid, alpha_step = self._refine_alpha_grid(
                    par_final[2], self.length_alpha_grid, alpha_step,
                )
                omega_grid, omega_step = self._refine_omega_grid(
                    par_final[4], self.length_omega_grid, omega_step,
                )

        if par_final is None:
            # Fallback: modelo plano
            par_final = np.array([data.mean(), 0.0, 0.0, 0.0, 0.5])

        M, A, alpha, beta, omega = par_final

        # ── Valores ajustados finales y estadisticos ───────────────────
        fitted        = _fmm_predict(time_points, M, A, alpha, beta, omega)
        residuals     = data - fitted
        # ddof=5: parametros libres del modelo (M, A, alpha, beta, omega).
        # El anexo "Deteccion de muestras atipicas" usa ddof=1 (via sd() de R),
        # pero esa formula solo es insesgada cuando el unico parametro estimado
        # es la media; aqui se estimaron 5, asi que el denominador correcto de
        # la varianza residual es n-5, no n-1.
        residuals_std = self._standardise_residuals(residuals, ddof=5)
        r2            = self._r2(data, fitted)
        sse           = float(np.sum(residuals**2))

        return FitResult(
            fitted        = fitted,
            params        = {
                "M":     float(M),
                "A":     float(A),
                "alpha": float(alpha),
                "beta":  float(beta),
                "omega": float(omega),
            },
            peak_time     = FMMModel.peak_time(float(alpha), float(beta), float(omega)),
            r2            = r2,
            residuals     = residuals,
            residuals_std = residuals_std,
            sse           = sse,
            model_name    = "fmm",
        )

    # ------------------------------------------------------------------
    # Auxiliares para refinamiento de rejilla
    # ------------------------------------------------------------------

    def _refine_alpha_grid(
        self,
        centre:   float,
        n_points: int,
        step:     float,
    ) -> tuple[np.ndarray, float]:
        """
        Estrecha la rejilla de alpha alrededor de ``centre``.

        El paso de la rejilla anterior se recibe explicitamente en ``step`` en
        vez de reconstruirse con ``mean(diff(prev_grid))``. La reconstruccion
        falla porque la rejilla ya viene envuelta a ``[0, 2*pi)``: si ``centre``
        cae cerca del origen, la rejilla contiene valores proximos a 0 y a
        2*pi a la vez, de modo que ``mean(diff(sort(grid)))`` devuelve
        aproximadamente ``2*pi/(n-1)`` y la rejilla no se estrecha en ninguna
        de las iteraciones siguientes.

        Returns
        -------
        (nueva_rejilla, nuevo_paso)
        """
        amplitude = 1.5 * step
        new_grid  = np.linspace(centre - amplitude, centre + amplitude, n_points)
        new_step  = (2.0 * amplitude / (n_points - 1)) if n_points > 1 else step
        return new_grid % (2 * np.pi), new_step

    def _refine_omega_grid(
        self,
        centre:   float,
        n_points: int,
        step:     float,
    ) -> tuple[np.ndarray, float]:
        """
        Estrecha la rejilla de omega alrededor de ``centre``.

        Igual que ``_refine_alpha_grid`` pero sin envoltura circular. Los
        extremos se recortan a ``(0, omega_max]``, asi que el paso efectivo se
        mide sobre la rejilla generada.

        Returns
        -------
        (nueva_rejilla, nuevo_paso)
        """
        amplitude = 1.5 * step
        lo = max(centre - amplitude, 1e-6)
        hi = min(centre + amplitude, self.omega_max)
        if hi <= lo:
            lo, hi = max(self.omega_max - 1e-6, 1e-6), self.omega_max
        new_grid = np.linspace(lo, hi, n_points)
        new_step = float(np.mean(np.diff(new_grid))) if n_points > 1 else step
        return new_grid, new_step
