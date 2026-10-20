"""
circust/fitting/ori.py
=======================
Modelo ritmico ORI — regresion isotonica circular unimodal (no parametrica).

Referencia
----------
Larriba, Rueda & Fernandez (2016). «Identification of rhythmic genes
in complex biological datasets and application to the study of human
circadian system gene expression.» *Statistical Applications in
Genetics and Molecular Biology*, 15(3), 205-223.

Modelo matematico
-----------------
Dado un perfil de expresion genica ordenado por tiempo circular, el
modelo ORI busca el mejor ajuste de forma subida-bajada (un pico, un
valle) mediante el Algoritmo Pool-Adjacent-Violators (PAVA).

Para cada par candidato (L=valle, U=pico) de extremos locales se ajusta:
  - Un PAVA creciente de L a U (circularmente).
  - Un PAVA decreciente de U a L (circularmente).

Se selecciona el par (L_opt, U_opt) que minimiza el MSE del ajuste
unimodal resultante. La puntuacion de ritmicidad es:

    R2_ORI = 1 - MSE_unimodal / MSE_plano

donde MSE_plano = Var(datos) es el error del modelo nulo (media constante).

Algoritmo PAVA (Pool-Adjacent-Violators)
-----------------------------------------
Complejidad: O(n) tiempo, O(n) espacio.

Dado un vector de observaciones *y* con pesos *w*, produce el vector
y_hat que minimiza  sum w_i (y_i - y_hat_i)^2  sujeto a restricciones
de monotonia (creciente o decreciente).

Aceleracion
-----------
Cuando Numba esta disponible, las funciones criticas (PAVA, busqueda
de pares L-U) se compilan JIT a codigo maquina, obteniendo ~20-50x de
aceleracion dado que el bucle (L, U) x PAVA es O(n^3) por gen.
"""
from __future__ import annotations

import numpy as np

from circust.fitting.rhythm_model import FitResult, RhythmModel

# ── Aceleracion opcional con Numba ────────────────────────────────────────
try:
    from numba import njit
    _HAS_NUMBA = True
except ImportError:
    _HAS_NUMBA = False
    def njit(*args, **kwargs):                          # type: ignore
        def deco(fn): return fn
        if len(args) == 1 and callable(args[0]):
            return args[0]
        return deco


# ═══════════════════════════════════════════════════════════════════════════
# PAVA — Pool-Adjacent-Violators Algorithm
# ═══════════════════════════════════════════════════════════════════════════

@njit(cache=True, fastmath=True)
def _pava_inc_numba(y: np.ndarray, w: np.ndarray) -> np.ndarray:
    """
    PAVA creciente compilado JIT.

    Minimiza  sum w_i (y_i - y_hat_i)^2  s.a.  y_hat_1 <= ... <= y_hat_n.
    """
    n = y.shape[0]
    if n <= 1:
        out = np.empty(n, dtype=np.float64)
        for i in range(n):
            out[i] = y[i]
        return out

    block_sum = np.empty(n, dtype=np.float64)
    block_w   = np.empty(n, dtype=np.float64)
    block_end = np.empty(n, dtype=np.int64)
    nb = 0

    for i in range(n):
        block_sum[nb] = y[i] * w[i]
        block_w[nb]   = w[i]
        block_end[nb] = i
        nb += 1
        while nb >= 2:
            avg_prev = block_sum[nb - 2] / block_w[nb - 2]
            avg_curr = block_sum[nb - 1] / block_w[nb - 1]
            if avg_prev <= avg_curr:
                break
            block_sum[nb - 2] += block_sum[nb - 1]
            block_w[nb - 2]   += block_w[nb - 1]
            block_end[nb - 2]  = block_end[nb - 1]
            nb -= 1

    out = np.empty(n, dtype=np.float64)
    start = 0
    for b in range(nb):
        end = block_end[b] + 1
        avg = block_sum[b] / block_w[b]
        for k in range(start, end):
            out[k] = avg
        start = end
    return out


@njit(cache=True, fastmath=True)
def _pava_inc_range(buf: np.ndarray, out: np.ndarray, m: int) -> None:
    """PAVA creciente sobre buf[:m], escribe resultado en out[:m]. Sin alloc."""
    if m <= 0:
        return
    block_sum = np.empty(m, dtype=np.float64)
    block_w   = np.empty(m, dtype=np.float64)
    block_end = np.empty(m, dtype=np.int64)
    nb = 0
    for i in range(m):
        block_sum[nb] = buf[i]
        block_w[nb]   = 1.0
        block_end[nb] = i
        nb += 1
        while nb >= 2:
            ap = block_sum[nb - 2] / block_w[nb - 2]
            ac = block_sum[nb - 1] / block_w[nb - 1]
            if ap <= ac:
                break
            block_sum[nb - 2] += block_sum[nb - 1]
            block_w[nb - 2]   += block_w[nb - 1]
            block_end[nb - 2]  = block_end[nb - 1]
            nb -= 1
    start = 0
    for b in range(nb):
        end = block_end[b] + 1
        avg = block_sum[b] / block_w[b]
        for k in range(start, end):
            out[k] = avg
        start = end


def pava_increasing(y: np.ndarray, w: np.ndarray | None = None) -> np.ndarray:
    """
    Regresion isotonica creciente mediante PAVA.

    Complejidad: O(n) tiempo, O(n) espacio.

    Parameters
    ----------
    y : array (n,)
        Valores observados.
    w : array (n,) o None
        Pesos positivos. Si es None, pesos uniformes.

    Returns
    -------
    array (n,) — valores ajustados monotonicamente no decrecientes.
    """
    y = np.ascontiguousarray(np.asarray(y, dtype=np.float64))
    n = len(y)
    if n <= 1:
        return y.copy()

    w_arr = (np.ones(n, dtype=np.float64) if w is None
             else np.ascontiguousarray(np.asarray(w, dtype=np.float64)))

    if _HAS_NUMBA:
        return _pava_inc_numba(y, w_arr)

    # Fallback Python puro
    block_sum = np.empty(n, dtype=np.float64)
    block_w   = np.empty(n, dtype=np.float64)
    block_end = np.empty(n, dtype=np.intp)
    nb = 0

    for i in range(n):
        block_sum[nb] = y[i] * w_arr[i]
        block_w[nb]   = w_arr[i]
        block_end[nb] = i
        nb += 1

        while nb >= 2:
            avg_prev = block_sum[nb - 2] / block_w[nb - 2]
            avg_curr = block_sum[nb - 1] / block_w[nb - 1]
            if avg_prev <= avg_curr:
                break
            block_sum[nb - 2] += block_sum[nb - 1]
            block_w[nb - 2]   += block_w[nb - 1]
            block_end[nb - 2]  = block_end[nb - 1]
            nb -= 1

    result = np.empty(n, dtype=np.float64)
    start = 0
    for b in range(nb):
        end = block_end[b] + 1
        result[start:end] = block_sum[b] / block_w[b]
        start = end

    return result


def pava_decreasing(y: np.ndarray, w: np.ndarray | None = None) -> np.ndarray:
    """
    Regresion isotonica decreciente.

    Equivalente a negar *y*, ejecutar ``pava_increasing``, y negar de vuelta.
    """
    return -pava_increasing(-np.asarray(y, dtype=np.float64), w)


# ═══════════════════════════════════════════════════════════════════════════
# Deteccion de extremos locales circulares
# ═══════════════════════════════════════════════════════════════════════════

def find_local_extrema(v: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    Encuentra minimos locales (valles) y maximos locales (picos) usando
    vecinos circulares.

    Un punto i es un maximo local si  v[i-1] <= v[i]  Y  v[i] >= v[i+1],
    con indices envolventes circulares. Minimos definidos simetricamente.

    Parameters
    ----------
    v : array (n,)

    Returns
    -------
    candL : array int — indices de minimos locales (valles).
    candU : array int — indices de maximos locales (picos).
    """
    v = np.asarray(v, dtype=np.float64)
    n = len(v)
    if n < 3:
        return np.array([], dtype=int), np.array([], dtype=int)

    prev  = np.roll(v, 1)
    next_ = np.roll(v, -1)

    candU = np.where((prev <= v) & (v >= next_))[0]    # local maxima
    candL = np.where((prev >= v) & (v <= next_))[0]    # local minima

    return candL, candU


# ═══════════════════════════════════════════════════════════════════════════
# Ajuste isotonico circular unimodal — nucleo Numba
# ═══════════════════════════════════════════════════════════════════════════

@njit(cache=True, fastmath=True)
def _circular_unimodal_fit_numba(
    v: np.ndarray, candL: np.ndarray, candU: np.ndarray,
):
    """
    Version JIT del bucle (L, U) x PAVA. Devuelve:
        (fitted, mse, L_opt, U_opt, ok)
    ``ok`` es True si se encontro un ajuste valido.
    """
    n = v.shape[0]
    nL = candL.shape[0]
    nU = candU.shape[0]

    best_mse  = 1e300
    best_L    = 0
    best_U    = 0
    found     = False
    best_fit  = np.empty(n, dtype=np.float64)

    # Buffers reutilizables
    ordered_U = np.empty(nU, dtype=np.int64)
    inc_buf   = np.empty(n + 2, dtype=np.float64)
    inc_out   = np.empty(n + 2, dtype=np.float64)
    full_inc  = np.empty(n + 2, dtype=np.float64)
    dec_buf   = np.empty(n + 2, dtype=np.float64)
    dec_out   = np.empty(n + 2, dtype=np.float64)
    valid_Us      = np.empty(nU, dtype=np.int64)
    valid_inc_len = np.empty(nU, dtype=np.int64)
    valid_inc_mat = np.empty((nU, n + 2), dtype=np.float64)
    fitted        = np.empty(n, dtype=np.float64)

    for li in range(nL):
        indL = candL[li]

        # Ordenar candU: primero >= indL, luego < indL
        k = 0
        for j in range(nU):
            if candU[j] >= indL:
                ordered_U[k] = candU[j]; k += 1
        for j in range(nU):
            if candU[j] < indL:
                ordered_U[k] = candU[j]; k += 1

        n_valid = 0

        for uj in range(nU):
            indU = ordered_U[uj]

            # ── _increasing_segment ──
            if indL == indU:
                mean_v = 0.0
                for t in range(n):
                    mean_v += v[t]
                mean_v /= n
                mse = 0.0
                for t in range(n):
                    d = v[t] - mean_v
                    mse += d * d
                mse /= n
                if mse < best_mse:
                    best_mse = mse
                    best_L = indL
                    best_U = indU
                    for t in range(n):
                        best_fit[t] = mean_v
                    found = True
                continue

            if indL < indU:
                k_interior = indU - indL - 1
                if k_interior > 0:
                    for t in range(k_interior):
                        inc_buf[t] = v[indL + 1 + t]
            else:
                k_interior = (n - indL - 1) + indU
                if k_interior > 0:
                    p = 0
                    for t in range(indL + 1, n):
                        inc_buf[p] = v[t]; p += 1
                    for t in range(0, indU):
                        inc_buf[p] = v[t]; p += 1

            inc_len = k_interior + 2
            if k_interior > 0:
                _pava_inc_range(inc_buf, inc_out, k_interior)
                full_inc[0] = v[indL]
                for t in range(k_interior):
                    full_inc[1 + t] = inc_out[t]
                full_inc[k_interior + 1] = v[indU]
            else:
                full_inc[0] = v[indL]
                full_inc[1] = v[indU]

            # Validez
            if inc_len > 2:
                if full_inc[1] <= v[indL]:
                    break
                is_valid_U = full_inc[inc_len - 2] < v[indU]
            else:
                is_valid_U = True

            if is_valid_U:
                valid_Us[n_valid] = indU
                valid_inc_len[n_valid] = inc_len
                for t in range(inc_len):
                    valid_inc_mat[n_valid, t] = full_inc[t]
                n_valid += 1

        # Pasada hacia atras
        for rev in range(n_valid - 1, -1, -1):
            indU    = valid_Us[rev]
            inc_len = valid_inc_len[rev]

            if indL == indU:
                continue

            # ── _decreasing_segment ──
            if indU < indL:
                k_int = indL - indU - 1
                if k_int > 0:
                    for t in range(k_int):
                        dec_buf[t] = v[indU + 1 + t]
            else:
                k_int = (n - indU - 1) + indL
                if k_int > 0:
                    p = 0
                    for t in range(indU + 1, n):
                        dec_buf[p] = v[t]; p += 1
                    for t in range(0, indL):
                        dec_buf[p] = v[t]; p += 1

            if k_int == 0:
                is_valid_dec = True
            else:
                for t in range(k_int):
                    dec_buf[t] = -dec_buf[t]
                _pava_inc_range(dec_buf, dec_out, k_int)
                for t in range(k_int):
                    dec_out[t] = -dec_out[t]

                if dec_out[k_int - 1] < v[indL]:
                    break
                is_valid_dec = dec_out[0] <= v[indU]

            if not is_valid_dec:
                continue

            # ── _assemble ──
            if indL <= indU:
                for t in range(inc_len):
                    fitted[indL + t] = valid_inc_mat[rev, t]
            else:
                p = 0
                for t in range(indL, n):
                    fitted[t] = valid_inc_mat[rev, p]; p += 1
                for t in range(0, indU + 1):
                    fitted[t] = valid_inc_mat[rev, p]; p += 1

            if k_int > 0:
                if indU < indL:
                    for t in range(k_int):
                        fitted[indU + 1 + t] = dec_out[t]
                else:
                    p = 0
                    for t in range(indU + 1, n):
                        fitted[t] = dec_out[p]; p += 1
                    for t in range(0, indL):
                        fitted[t] = dec_out[p]; p += 1

            # MSE
            mse = 0.0
            for t in range(n):
                d = v[t] - fitted[t]
                mse += d * d
            mse /= n

            if mse < best_mse:
                best_mse = mse
                best_L = indL
                best_U = indU
                for t in range(n):
                    best_fit[t] = fitted[t]
                found = True

        if best_mse == 0.0:
            break

    return best_fit, best_mse, best_L, best_U, found


# ═══════════════════════════════════════════════════════════════════════════
# Ajuste isotonico circular unimodal — API Python
# ═══════════════════════════════════════════════════════════════════════════

def _increasing_segment(
    v: np.ndarray, v2: np.ndarray,
    indL: int, indU: int, n: int,
) -> tuple[np.ndarray, bool] | None:
    """
    PAVA creciente de indL a indU (circular).

    Devuelve (full_inc, is_valid_U) o None si no hay mas U's validos.
    """
    if indL == indU:
        return np.full(n, v.mean()), True

    if indL < indU:
        k_interior = indU - indL - 1
        interior_vals = v[indL + 1 : indU] if k_interior > 0 else np.array([], dtype=np.float64)
    else:
        k_interior = (n - indL - 1) + indU
        interior_vals = v2[indL + 1 : indL + 1 + k_interior] if k_interior > 0 else np.array([], dtype=np.float64)

    if len(interior_vals) > 0:
        inc_fit  = pava_increasing(interior_vals)
        full_inc = np.empty(len(interior_vals) + 2, dtype=np.float64)
        full_inc[0]    = v[indL]
        full_inc[1:-1] = inc_fit
        full_inc[-1]   = v[indU]
    else:
        full_inc = np.array([v[indL], v[indU]], dtype=np.float64)

    if len(full_inc) > 2:
        if full_inc[1] <= v[indL]:
            return None
        is_valid = full_inc[-2] < v[indU]
    else:
        is_valid = True

    return full_inc, is_valid


def _decreasing_segment(
    v: np.ndarray, v2: np.ndarray,
    indL: int, indU: int, n: int,
) -> tuple[np.ndarray, bool] | None:
    """
    PAVA decreciente de indU a indL (circular, interior).

    Devuelve (dec_fit, is_valid) o None si no hay mas U's validos.
    """
    if indU < indL:
        k_interior = indL - indU - 1
        interior_vals = v[indU + 1 : indL] if k_interior > 0 else np.array([], dtype=np.float64)
    else:
        k_interior = (n - indU - 1) + indL
        interior_vals = v2[indU + 1 : indU + 1 + k_interior] if k_interior > 0 else np.array([], dtype=np.float64)

    if len(interior_vals) == 0:
        return np.array([], dtype=np.float64), True

    dec_fit = pava_decreasing(interior_vals)

    if dec_fit[-1] < v[indL]:
        return None

    is_valid = dec_fit[0] <= v[indU]
    return dec_fit, is_valid


def _assemble(
    v: np.ndarray,
    indL: int, indU: int,
    full_inc: np.ndarray,
    dec_fit: np.ndarray,
    n: int,
) -> np.ndarray:
    """Coloca arcos creciente y decreciente en un array de longitud n."""
    fitted = np.empty(n, dtype=np.float64)

    if indL <= indU:
        inc_pos = np.arange(indL, indU + 1)
    else:
        inc_pos = np.concatenate([np.arange(indL, n), np.arange(0, indU + 1)])
    fitted[inc_pos] = full_inc

    if len(dec_fit) > 0:
        if indU < indL:
            dec_pos = np.arange(indU + 1, indL)
        else:
            dec_pos = np.concatenate([np.arange(indU + 1, n), np.arange(0, indL)])
        fitted[dec_pos] = dec_fit

    return fitted


def circular_unimodal_fit(
    v: np.ndarray,
) -> tuple[np.ndarray, float, int, int] | None:
    """
    Ajusta una regresion isotonica circular unimodal (subida-bajada).

    Busca sobre todos los pares candidatos (valle, pico) — minimos y
    maximos locales de *v* — para encontrar el par (L, U) que produce
    el menor MSE al ajustar un PAVA creciente de L->U y un PAVA
    decreciente de U->L (circularmente).

    Parameters
    ----------
    v : array (n,)
        Valores de expresion genica ordenados por tiempo circular.

    Returns
    -------
    fitted : array (n,) — valores ajustados.
    mse    : float      — error cuadratico medio.
    L_opt  : int        — indice del valle optimo (base 0).
    U_opt  : int        — indice del pico optimo (base 0).
    None si no se encuentra un ajuste valido.
    """
    v = np.asarray(v, dtype=np.float64)
    n = len(v)

    if n == 0:
        return None

    candL, candU = find_local_extrema(v)

    if len(candL) == 0 or len(candU) == 0:
        fitted = np.full(n, v.mean())
        mse = float(np.mean((v - fitted) ** 2))
        return fitted, mse, 0, 0

    # ── Ruta rapida Numba ────────────────────────────────────────────────
    if _HAS_NUMBA:
        fitted, mse, L_opt, U_opt, ok = _circular_unimodal_fit_numba(
            v,
            np.ascontiguousarray(candL.astype(np.int64)),
            np.ascontiguousarray(candU.astype(np.int64)),
        )
        if ok:
            return fitted, float(mse), int(L_opt), int(U_opt)
        return None

    # ── Fallback Python puro ─────────────────────────────────────────────
    v2 = np.concatenate([v, v])

    best_mse    = np.inf
    best_fitted = None
    best_L      = 0
    best_U      = 0

    for indL in candL:
        mask      = candU >= indL
        ordered_U = np.concatenate([candU[mask], candU[~mask]])

        valid_Us:   list[int]        = []
        valid_incs: list[np.ndarray] = []

        for indU in ordered_U:
            inc_result = _increasing_segment(v, v2, indL, indU, n)

            if inc_result is None:
                break

            full_inc, is_valid_U = inc_result
            if is_valid_U:
                valid_Us.append(indU)
                valid_incs.append(full_inc)

        for idx_rev, indU in enumerate(reversed(valid_Us)):
            j = len(valid_Us) - 1 - idx_rev
            full_inc = valid_incs[j]

            if indL == indU:
                fitted = full_inc
                mse = float(np.mean((v - fitted) ** 2))
                if mse < best_mse:
                    best_mse, best_fitted = mse, fitted
                    best_L, best_U = indL, indU
                continue

            dec_result = _decreasing_segment(v, v2, indL, indU, n)

            if dec_result is None:
                break

            dec_fit, is_valid = dec_result
            if not is_valid:
                continue

            fitted = _assemble(v, indL, indU, full_inc, dec_fit, n)
            mse = float(np.mean((v - fitted) ** 2))

            if mse < best_mse:
                best_mse, best_fitted = mse, fitted
                best_L, best_U = indL, indU

        if best_mse == 0.0:
            break

    if best_fitted is None:
        return None

    return best_fitted, best_mse, best_L, best_U


# ═══════════════════════════════════════════════════════════════════════════
# ORIModel — implementacion de RhythmModel
# ═══════════════════════════════════════════════════════════════════════════

class ORIModel(RhythmModel):
    """
    Modelo de ritmicidad no parametrica basado en regresion isotonica
    circular unimodal (ORI / ORIOS).

    A diferencia de Cosinor y FMM, este modelo no asume una forma
    de onda parametrica. Busca la mejor forma subida-bajada (unimodal)
    compatible con restricciones de monotonia impuestas por PAVA.

    Los "parametros" del modelo son los indices del valle (L_opt)
    y pico (U_opt) optimos; el R2 compara contra el modelo plano.

    Example
    -------
    >>> import numpy as np
    >>> from circust.fitting.ori import ORIModel
    >>>
    >>> t    = np.linspace(0, 2*np.pi, 50, endpoint=False)
    >>> data = 0.3 + 0.6 * np.cos(t - 1.2) + np.random.normal(0, 0.05, 50)
    >>> result = ORIModel().fit(data, t)
    >>> print(result.r2)
    >>> print(result.summary())
    """

    @staticmethod
    def peak_time(U_opt: int, time_points: np.ndarray) -> float:
        """
        Tiempo de pico del modelo ORI.

        El pico se localiza en la posicion U_opt del vector ordenado
        circularmente. El tiempo de pico es ``time_points[U_opt]``.
        """
        if len(time_points) == 0:
            return 0.0
        return float(time_points[U_opt % len(time_points)])

    def fit(
        self,
        data:        np.ndarray,
        time_points: np.ndarray,
    ) -> FitResult:
        """
        Ajusta el modelo ORI a ``data``.

        Parameters
        ----------
        data : np.ndarray, forma (n_muestras,)
            Valores de expresion ordenados por el tiempo circular
            preliminar.
        time_points : np.ndarray, forma (n_muestras,)
            Eje temporal circular en [0, 2*pi), producido por CPCA.

        Returns
        -------
        FitResult
            Claves de ``params``: ``L_opt``, ``U_opt``.

            - L_opt : indice del valle optimo (base 0)
            - U_opt : indice del pico optimo (base 0)

            ``peak_time`` = ``time_points[U_opt]``.
            ``r2`` = R2_ORI = 1 - MSE_unimodal / MSE_plano.
        """
        data        = np.asarray(data,        dtype=np.float64)
        time_points = np.asarray(time_points, dtype=np.float64)
        n           = len(data)

        # ── Ajuste unimodal circular (PAVA) ───────────────────────────
        result = circular_unimodal_fit(data)

        if result is not None:
            fitted, _, L_opt, U_opt = result
        else:
            fitted = np.full(n, data.mean())
            L_opt  = 0
            U_opt  = 0

        # ── Residuos y estadisticos (metodos heredados de RhythmModel) ─
        residuals     = data - fitted
        residuals_std = self._standardise_residuals(residuals)
        r2            = self._r2(data, fitted)
        sse           = float(np.sum(residuals ** 2))
        pk            = ORIModel.peak_time(U_opt, time_points)

        return FitResult(
            fitted        = fitted,
            params        = {"L_opt": int(L_opt), "U_opt": int(U_opt)},
            peak_time     = pk,
            r2            = r2,
            residuals     = residuals,
            residuals_std = residuals_std,
            sse           = sse,
            model_name    = "ori",
        )
