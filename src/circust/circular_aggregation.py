"""
circust/circular_aggregation.py
================================
Agregacion de ordenes circulares — paper Barragan, Rueda y Fernandez (2015).

Toma K ordenes circulares observados sobre las mismas n items (cada orden es
un vector de angulos en [0, 2*pi)) y produce un orden consenso. Util para
combinar resultados de distintas repeticiones del pipeline RobustOrder.

Alcance respecto al paper
-------------------------
El paper organiza sus metodos en dos pasos, "a first step [where] an initial
solution denoted as O_0 is provided, which is refined in step 2".

Aqui se implementa **solo el paso 1**. Para el TSP el paper no fija una
heuristica concreta y prescribe en su lugar la regla de decision: "a better
approximation to the optimum is obtained by repeatedly running different
heuristics and selecting the best solution **in terms of the objective
function**". Es lo que hace :func:`_solve_tsp_multi_start`.

El paso 2, la busqueda local CLMA sobre ternas de elementos consecutivos, se
ha omitido deliberadamente. El motivo es que CIRCUST opera fuera del regimen
en que compensa. Los autores ya lo acotan, "this second step does not improve
too much the results for the methods whose step 1 already yielded good results
(such as TSP3 or HODs under MSCE) [...] CLMA is suggested just for those
situations where it is convenient to refine the initial solution", y la medida
sobre datos reales lo confirma. Las K repeticiones de CIRCUST son subconjuntos
del mismo conjunto TOP sobre las mismas muestras, de modo que su discrepancia
angular es muy pequena, con una desviacion tipica circular mediana de 0,044
radianes sobre matrixIn_mini. En ese regimen CLMA acepta del orden de una
permutacion y mejora el orden un 8 % en desplazamiento medio de rango, sobre
un error que ya es de 0,17 posiciones, mientras que su coste es O(K n^4) por
pasada porque evalua el MSCE 5n veces. Medido sobre el propio pipeline con
n=114 y K=5, el MSCE no cambio en cuatro decimales y el paso costo 9 min 44 s.

Metodos implementados
---------------------
- **TSP3**: TSP con distancia direccional asimetrica con penalizacion alpha=3
  para movimientos en la direccion equivocada. Recomendado por el paper para
  el criterio MSCE (Mean Sum of Circular Errors). Los costes de arista siguen
  la Tabla 1 y la ecuacion (7) del paper.

- **HODs**: Hodge score-based triplewise (Tabla 2 y Teorema 4.1). El flujo
  solo es no nulo en los tripletes que contienen el ancla ``l``, donde vale la
  diferencia de angulos rotados a ``theta_l = 0``. Los angulos se agregan entre
  experimentos por media circular y el ancla se elige por ``gamma_l``. Es la
  variante ``nat`` de ``isocir::ACO(method="CH")`` (que fija el ancla en el
  primer item), verificada contra ella en ``tests/test_circular_aggregation.py``.
  Una version anterior de este modulo implementaba por error la variante
  ``ave`` de isocir, que no aparece en el paper y da un MSCE mucho peor.

Funcion objetivo (MSCE)
-----------------------
El MSCE se calcula con el **CIRE** (Circular Isotonic Regression Estimator) de
Rueda, Fernandez y Peddada (2009), igual que el paquete R ``isocir``. Para cada
repeticion se proyecta su vector angular sobre el conjunto de vectores que
verifican el orden candidato dando una sola vuelta al circulo (regresion
isotonica circular), y el MSCE de esa repeticion es ``(1/n) * SCE`` con
``SCE = sum(1 - cos(theta_i - theta_tilde_i))``. El MSCE global es la media
(ponderada) de los MSCE por repeticion.

La implementacion del CIRE aqui es una version numba casi-exacta validada
contra ``isocir`` (coincide salvo desviaciones <= 0.024 en el SCE en casos con
bloques muy anchos, irrelevantes para ordenar metodos). Se usa como criterio de
seleccion entre TSP3 y HODs, no como cifra a reportar.

API publica
-----------
- ``aggregate_circular_orders(Theta, weights, method="tsp3")``: devuelve un
  orden consenso (permutacion de items) y su MSCE.
- ``circular_isotonic_regression(theta, order)``: CIRE de un vector angular
  bajo un orden dado (diagnostico / tests).

Equivalente en R: paquete ``isocir`` (funciones ``CIRE`` y ``ACO``).
"""
from __future__ import annotations

from typing import Optional
import numpy as np

# ── Aceleracion opcional con Numba (mismo patron que circust.fitting.ori) ──
try:
    from numba import njit
    _HAS_NUMBA = True
except ImportError:                                   # pragma: no cover
    _HAS_NUMBA = False
    def njit(*args, **kwargs):                         # type: ignore
        def deco(fn): return fn
        if len(args) == 1 and callable(args[0]):
            return args[0]
        return deco

TWO_PI = 2.0 * np.pi


# ═══════════════════════════════════════════════════════════════════════════
# Seccion 1 — Helpers numericos
# ═══════════════════════════════════════════════════════════════════════════

def _circular_dist(a: float, b: float) -> float:
    """Distancia angular minima entre dos angulos (resultado en ``[0, pi]``).

    Equivalente a la metrica `min(d, 2π−d)` con d=|a−b| mod 2π.
    """
    d = abs(a - b) % (2.0 * np.pi)
    return float(min(d, 2.0 * np.pi - d))


def _circular_dist_vec(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Version vectorizada de ``_circular_dist`` (broadcasting)."""
    d = np.abs(a - b) % (2.0 * np.pi)
    return np.minimum(d, 2.0 * np.pi - d)


# ═══════════════════════════════════════════════════════════════════════════
# Seccion 2 — CIRE y MSCE (funcion objetivo cardinal del paper)
# ═══════════════════════════════════════════════════════════════════════════

@njit(cache=True, fastmath=True)
def _cire_core(d: np.ndarray, S: np.ndarray, C: np.ndarray):
    """Nucleo del CIRE circular sobre datos YA ordenados ``d`` (con S=sin, C=cos).

    Enumera el arco que forma el primer bloque (corte ``c`` x longitud ``m``),
    fija su media circular ``mu0`` como referencia de una sola vuelta, resuelve
    el resto con una PAVA en la ventana ``[mu0, mu0+2pi)`` (pooling por media
    circular, violacion cuando la media-ventana decrece) y se queda con el
    minimo SCE. Devuelve ``(fitted_posiciones, SCE)``.

    Es la version eficiente (O(n^3)) validada contra ``isocir``.
    """
    q = d.shape[0]
    best_sce = 1e18
    best_fit = np.empty(q)
    bS = np.empty(q); bC = np.empty(q)
    bN = np.empty(q, dtype=np.int64); bW = np.empty(q)
    cur_fit = np.empty(q)
    for c in range(q):
        for m in range(1, q + 1):
            # primer bloque = posiciones (c..c+m-1) mod q
            S0 = 0.0; C0 = 0.0
            for t in range(m):
                idx = (c + t) % q
                S0 += S[idx]; C0 += C[idx]
            mu0 = np.arctan2(S0, C0) % TWO_PI
            cost = 0.0
            for t in range(m):
                idx = (c + t) % q
                cur_fit[idx] = mu0
                cost += 1.0 - np.cos(d[idx] - mu0)
            # PAVA del resto en ventana [mu0, mu0+2pi)
            nb = 0
            for t in range(m, q):
                idx = (c + t) % q
                cs = S[idx]; cc = C[idx]; cn = 1
                win = (np.arctan2(cs, cc) - mu0) % TWO_PI
                while nb > 0 and bW[nb - 1] > win + 1e-12:
                    nb -= 1
                    cs += bS[nb]; cc += bC[nb]; cn += bN[nb]
                    win = (np.arctan2(cs, cc) - mu0) % TWO_PI
                bS[nb] = cs; bC[nb] = cc; bN[nb] = cn; bW[nb] = win; nb += 1
            # expandir bloques del resto -> cur_fit y acumular coste
            pos = m
            for bi in range(nb):
                mu = np.arctan2(bS[bi], bC[bi]) % TWO_PI
                for _ in range(bN[bi]):
                    idx = (c + pos) % q
                    cur_fit[idx] = mu
                    cost += 1.0 - np.cos(d[idx] - mu)
                    pos += 1
            if cost < best_sce - 1e-12:
                best_sce = cost
                for i in range(q):
                    best_fit[i] = cur_fit[i]
    return best_fit, best_sce


def circular_isotonic_regression(
    theta: np.ndarray, order: np.ndarray,
) -> tuple[np.ndarray, float]:
    """CIRE de ``theta`` bajo el orden circular dado por ``order``.

    Parameters
    ----------
    theta : np.ndarray (n,)
        Angulos observados en ``[0, 2*pi)`` (indice por item).
    order : np.ndarray (n,)
        Permutacion de items: ``order[pos]`` es el item en la posicion ``pos``
        del orden circular candidato.

    Returns
    -------
    (fitted, sce) : np.ndarray (n,), float
        ``fitted`` es el CIRE alineado al indice de ``theta``; ``sce`` la suma
        de errores circulares ``sum(1 - cos(theta_i - fitted_i))``.
    """
    theta = np.asarray(theta, dtype=float) % TWO_PI
    order = np.asarray(order, dtype=np.int64)
    n = theta.shape[0]
    d = np.ascontiguousarray(theta[order])
    fit_pos, sce = _cire_core(d, np.sin(d), np.cos(d))
    fitted = np.empty(n)
    fitted[order] = fit_pos
    return fitted, float(sce)


def _msce_cire(
    Theta: np.ndarray,
    target_order: np.ndarray,
    weights: Optional[np.ndarray] = None,
) -> float:
    """MSCE del orden ``target_order`` frente a las K repeticiones de ``Theta``.

    Para cada repeticion k:  ``msce_k = SCE(CIRE(Theta[k] | target_order)) / n``.
    El MSCE global es la media (ponderada por ``weights``) de los ``msce_k``.
    Es el criterio cardinal del paper (eqs. 2-4) usando el CIRE de Rueda et al.

    Parameters
    ----------
    Theta : np.ndarray (K, n)
        ``Theta[k, i]`` angulo asignado al item ``i`` en la repeticion ``k``.
    target_order : np.ndarray (n,)
        Permutacion de items que define el orden consenso a evaluar.
    weights : np.ndarray opcional (K,)
        Pesos por repeticion. Default uniforme.
    """
    Theta = np.asarray(Theta, dtype=float) % TWO_PI
    K, n = Theta.shape
    order = np.asarray(target_order, dtype=np.int64)
    msces = np.empty(K)
    for k in range(K):
        d = np.ascontiguousarray(Theta[k][order])
        _, sce = _cire_core(d, np.sin(d), np.cos(d))
        msces[k] = sce / n
    if weights is None:
        return float(msces.mean())
    w = np.asarray(weights, dtype=float)
    if w.shape[0] != K:
        raise ValueError(f"weights debe tener longitud {K}, got {w.shape[0]}")
    return float(np.average(msces, weights=w))


# ═══════════════════════════════════════════════════════════════════════════
# Seccion 3 — TSP3 (alpha=3, recomendado para MSCE)
# ═══════════════════════════════════════════════════════════════════════════

def _distance_matrix_tsp3_single(theta: np.ndarray) -> np.ndarray:
    """Matriz de distancias TSP3 ``(n, n)`` para una sola repeticion.

    Para cada par (h, k):

      d_R(h, k) = 1 - cos(theta_k - theta_h)             si  diff <= pi
                  3 - cos(theta_k - theta_h - pi)        si  diff > pi
      d_C(h, k) = d_R(k, h)
      E[h, k]   = min(d_R(h, k), 3 * d_C(h, k))

    donde ``diff = (theta_k - theta_h) mod 2*pi``. La matriz ``E`` es
    asimetrica: capturar la direccion circular es la clave del enfoque TSP3.
    """
    theta = np.asarray(theta, dtype=float)
    diff  = (theta[None, :] - theta[:, None]) % (2.0 * np.pi)   # (n, n)
    # d_R: distancia en direccion "directa" (corto arco) o penalizada si hay que dar la vuelta
    d_R = np.where(
        diff <= np.pi,
        1.0 - np.cos(diff),
        3.0 - np.cos(diff - np.pi),
    )
    d_C = d_R.T                                                  # d_C(h,k) = d_R(k,h)
    E = np.minimum(d_R, 3.0 * d_C)
    np.fill_diagonal(E, 0.0)
    return E


def _build_distance_matrix_tsp3(
    Theta:   np.ndarray,
    weights: Optional[np.ndarray] = None,
) -> np.ndarray:
    """Matriz de distancias TSP3 agregada sobre las K repeticiones.

    ``E[h, k] = sum_j w_j * E^j[h, k]`` donde ``E^j`` es la matriz TSP3 de la
    repeticion j (formula del paper, eq. tras Sec. 3.1).
    """
    Theta = np.asarray(Theta, dtype=float)
    K, n = Theta.shape

    if weights is None:
        weights = np.ones(K) / K
    else:
        weights = np.asarray(weights, dtype=float)
        if weights.shape[0] != K:
            raise ValueError(f"weights debe tener longitud {K}, got {weights.shape[0]}")
        s = weights.sum()
        if s <= 0:
            raise ValueError("Suma de weights debe ser > 0.")
        weights = weights / s

    E = np.zeros((n, n), dtype=float)
    for j in range(K):
        E += weights[j] * _distance_matrix_tsp3_single(Theta[j])
    return E


def _solve_tsp_nearest_neighbor(D: np.ndarray, start: int) -> tuple[np.ndarray, float]:
    """TSP nearest-neighbor desde un nodo de inicio.

    Heuristica simple pero robusta: empezar en ``start`` y, en cada paso,
    moverse al nodo no visitado mas cercano. Devuelve ``(tour, length)``.
    """
    n = D.shape[0]
    visited = np.zeros(n, dtype=bool)
    tour    = np.empty(n, dtype=int)

    tour[0] = start
    visited[start] = True

    for i in range(1, n):
        last  = tour[i - 1]
        dists = D[last].copy()
        dists[visited] = np.inf
        next_node = int(np.argmin(dists))
        tour[i] = next_node
        visited[next_node] = True

    # Longitud incluyendo el cierre del tour (vuelta al inicio)
    length = float(sum(D[tour[i], tour[(i + 1) % n]] for i in range(n)))
    return tour, length


def _tour_length(D: np.ndarray, tour: np.ndarray) -> float:
    """Longitud del ciclo hamiltoniano ``tour`` sobre la matriz de costes ``D``."""
    nxt = np.roll(tour, -1)
    return float(D[tour, nxt].sum())


def _or_opt(
    D:          np.ndarray,
    tour:       np.ndarray,
    max_seg:    int = 3,
    max_passes: int = 20,
) -> np.ndarray:
    """Mejora local Or-opt sobre un tour, valida para TSP **asimetrico**.

    Un movimiento Or-opt extrae un segmento de 1 a ``max_seg`` nodos
    consecutivos y lo reinserta en otra posicion **sin invertirlo**. Es la
    familia de movimientos apropiada aqui: el 2-opt clasico invierte un tramo
    del recorrido, lo que en una matriz asimetrica como la de TSP3 cambia el
    coste de todas las aristas interiores y deja de poder evaluarse en O(1).

    Para un segmento ``S = (s_0, ..., s_{L-1})`` extraido entre ``p`` y ``q``, y
    reinsertado entre ``a`` y ``b``, la variacion de coste es

        delta = [D(a,s_0) + D(s_{L-1},b) + D(p,q)]
              - [D(p,s_0) + D(s_{L-1},q) + D(a,b)],

    de modo que el interior del segmento no interviene. El bucle sobre las
    posiciones de reinsercion esta vectorizado, asi que cada pasada cuesta
    O(max_seg * n^2) operaciones NumPy.

    Se aplica el mejor movimiento de cada barrido (steepest descent) hasta que
    ninguno mejora o se agotan las pasadas.
    """
    t = np.asarray(tour, dtype=np.int64).copy()
    n = t.shape[0]
    if n < 5:
        return t

    for _ in range(max_passes):
        best_delta = -1e-12
        best_move  = None

        for L in range(1, min(max_seg, n - 3) + 1):
            for i in range(n):
                idx_seg = (np.arange(i, i + L)) % n
                seg     = t[idx_seg]
                p       = t[(i - 1) % n]           # nodo anterior al segmento
                q       = t[(i + L) % n]           # nodo posterior al segmento
                s0, sl  = seg[0], seg[-1]

                # Resto del ciclo tras extraer el segmento, empezando en q.
                rest = t[(np.arange(i + L, i + n)) % n]      # (n - L,)
                a    = rest
                b    = np.roll(rest, -1)

                delta = (D[a, s0] + D[sl, b] + D[p, q]) - (
                         D[p, s0] + D[sl, q] + D[a, b])
                m = int(np.argmin(delta))
                if delta[m] < best_delta:
                    best_delta = float(delta[m])
                    best_move  = (seg.copy(), rest.copy(), m)

        if best_move is None:
            break

        seg, rest, m = best_move
        t = np.concatenate([rest[: m + 1], seg, rest[m + 1:]])

    return t


def _solve_tsp_multi_start(
    D:           np.ndarray,
    n_starts:    Optional[int] = None,
    or_opt:      bool = True,
    warm_starts: Optional[list] = None,
    score_fn:    Optional[callable] = None,
) -> tuple[np.ndarray, float]:
    """Paso 1 de Barragan et al. (2015): solucion inicial ``O_0`` del TSP.

    El paper no fija una heuristica concreta para resolver (6). Lo que
    prescribe, al final de su Seccion 3, es la regla de decision:

        "there is not an absolute winner and that a better approximation to
         the optimum is obtained by repeatedly running different heuristics
         and selecting the best solution in terms of the objective function."

    De ahi las dos decisiones de esta funcion.

    La primera es generar una coleccion amplia de candidatos en lugar de uno.
    Se usan dos familias de recorridos, los ``n_starts`` del vecino mas cercano
    tomando sucesivamente como origen distintas muestras, y los de
    ``warm_starts``, tipicamente el orden que induce la media circular de las K
    repeticiones. Cada uno de ellos se pule despues con ``_or_opt``. La version
    anterior se quedaba solo con el recorrido de vecino mas cercano de menor
    longitud y pulia unicamente ese, de modo que descartaba 29 de los 30
    candidatos con un criterio, la longitud del tour, que no es la funcion
    objetivo. Medido sobre la familia del vecino mas cercano, entre 14 y 27 de
    cada 30 arranques pulen a un MSCE mejor que el que la longitud selecciona.

    La segunda es que la eleccion entre los candidatos pulidos la hace
    ``score_fn``, que ``aggregate_circular_orders`` fija en el MSCE, es decir en
    la funcion objetivo ``d*(Theta, O)`` de las ecuaciones (3) y (4) del paper.
    Este es el punto exacto que la cita anterior exige y la razon de que no
    valga elegir por longitud del tour: con ruido por replica grande frente a
    la separacion angular entre muestras el objetivo del TSP se sobreajusta al
    ruido, y existen ciclos de coste menor que el del orden verdadero
    construidos encadenando aristas cuyo coste bajo es fruto del azar.

    Devuelve ``(tour, longitud)``.
    """
    n = D.shape[0]
    if n_starts is None:
        n_starts = min(n, 30)               # 30 starts es suficiente en practica
    n_starts = max(1, min(n_starts, n))

    # Familia 1: vecino mas cercano desde distintos origenes. Se conservan
    # TODOS, no solo el de menor longitud.
    seeds: list[np.ndarray] = [
        _solve_tsp_nearest_neighbor(D, start=int(s))[0]
        for s in np.linspace(0, n - 1, n_starts, dtype=int)
    ]
    # Familia 2: arranques en caliente proporcionados por el llamante.
    seeds += [np.asarray(w, dtype=np.int64) for w in (warm_starts or [])]

    score = score_fn if score_fn is not None else (lambda t: _tour_length(D, t))

    best_tour, best_score = None, np.inf
    for seed in seeds:
        tour = _or_opt(D, seed) if or_opt else seed
        s    = float(score(tour))
        if s < best_score:
            best_tour, best_score = tour, s

    return best_tour, _tour_length(D, best_tour)


def _circular_mean_order(
    Theta: np.ndarray, weights: Optional[np.ndarray] = None,
) -> np.ndarray:
    """Orden que induce la media circular, muestra a muestra, de las K fases.

    Se usa como arranque en caliente del TSP. Las K repeticiones ya estan
    sincronizadas (mismo origen y sentido), de modo que promediar sus fases
    muestra a muestra es licito.
    """
    Theta = np.asarray(Theta, dtype=float)
    K = Theta.shape[0]
    w = (np.ones(K) if weights is None else np.asarray(weights, dtype=float))
    w = w / w.sum()
    mu = np.arctan2(
        (w[:, None] * np.sin(Theta)).sum(axis=0),
        (w[:, None] * np.cos(Theta)).sum(axis=0),
    ) % TWO_PI
    return np.argsort(mu).astype(np.int64)


# ═══════════════════════════════════════════════════════════════════════════
# Seccion 3b — HODs (Hodge score-based triplewise; Barragan et al. 2015)
# ═══════════════════════════════════════════════════════════════════════════

def _hodge_anchored_means(
    Theta: np.ndarray, w: np.ndarray,
) -> np.ndarray:
    """Medias circulares ponderadas de los angulos rotados a cada ancla.

    ``A[l, i]`` es la media circular, ponderada por ``w`` sobre los K
    experimentos, de ``(theta_i^j - theta_l^j) mod 2pi``. Es el paso con el que
    el paper pasa de angulos a scores en ``V^(l)``: eliminar ``l`` y tomarlo como
    origen del circulo, con lo que los angulos rotados son ya valores lineales
    en ``[0, 2pi)`` y ``psi_lhk = a_k - a_h``. La agregacion entre experimentos
    por media circular es la que el paper senala para ``Psi^j`` y la que usa
    ``isocir`` (variante ``nat``). ``A[l, l] = 0`` por construccion.
    """
    rot = (Theta[:, None, :] - Theta[:, :, None]) % TWO_PI      # (K, n, n)
    S = np.einsum("k,kli->li", w, np.sin(rot))
    C = np.einsum("k,kli->li", w, np.cos(rot))
    return np.arctan2(S, C) % TWO_PI


def _hodge_order(
    Theta:   np.ndarray,
    weights: Optional[np.ndarray] = None,
    anchor:  Optional[int] = None,
) -> np.ndarray:
    """HODs: orden consenso por Hodge con flujo score-based (Barragan et al.).

    Es la variante ``HODs`` de la Tabla 2 del paper. Para el ancla ``l`` el
    flujo triplewise vale ``psi_lhk = sign(l,h,k) * |theta_h - theta_k|`` sobre
    los angulos rotados a ``theta_l = 0``, que se reduce a ``a_k - a_h``, y es
    cero en los tripletes que no contienen a ``l``. Con esa hipermatriz el
    Teorema 4.1 da ``Y_ih = (a_h - a_i)/n`` para ``i, h != l`` y scores
    ``s_i = -(1/(n-1)) sum_{h != l} Y_ih = (a_i - media(a))/n``, de modo que el
    orden circular ``l <= i_1 <= ... <= i_{n-1} <= l`` es el de ``a`` creciente.
    Esto es lo que se calcula aqui de forma directa, en O(K n^2), sin construir
    la hipermatriz. El test ``test_hods_coincide_con_teorema_4_1_literal``
    verifica esta reduccion contra la aplicacion literal del teorema.

    El ancla se elige como propone el paper para HODp y HODs,
    ``l0 = argmax_l gamma_l`` con ``gamma_l = sum_{h,k in V^(l)} psi_lhk^2``.
    ``anchor`` permite fijarla, util para reproducir ``isocir``, cuya variante
    ``nat`` usa siempre el primer item.

    En el regimen de CIRCUST, con las K repeticiones ya sincronizadas, este
    metodo coincide en la practica con ordenar por la media circular de las
    fases de cada muestra, lo que explica que TSP3 y HODs empaten en MSCE.

    Returns
    -------
    np.ndarray (n,) : permutacion de items con el ancla en la primera posicion.
    """
    Theta = np.asarray(Theta, dtype=float) % TWO_PI
    K, n = Theta.shape
    if weights is None:
        w = np.ones(K) / K
    else:
        w = np.asarray(weights, dtype=float)
        w = w / w.sum()
    A = _hodge_anchored_means(Theta, w)                         # (n, n)

    if anchor is None:
        # gamma_l = sum_{h,k != l} (a_k - a_h)^2 = 2 (n-1) sum a^2 - 2 (sum a)^2
        # (las sumas excluyen l, cuyo a_l es 0).
        gamma = 2.0 * (n - 1) * (A ** 2).sum(axis=1) - 2.0 * A.sum(axis=1) ** 2
        l0 = int(np.argmax(gamma))
    else:
        l0 = int(anchor)
        if not 0 <= l0 < n:
            raise ValueError(f"anchor debe estar en [0, {n}), got {l0}")

    a = A[l0]
    keep = np.array([i for i in range(n) if i != l0], dtype=np.int64)
    order_keep = keep[np.argsort(a[keep], kind="stable")]
    return np.concatenate([[l0], order_keep]).astype(np.int64)


# ═══════════════════════════════════════════════════════════════════════════
# Seccion 4 — API publica
# ═══════════════════════════════════════════════════════════════════════════

def aggregate_circular_orders(
    Theta:      np.ndarray,
    weights:    Optional[np.ndarray] = None,
    method:     str = "auto",
    n_starts:   Optional[int] = None,
) -> dict:
    """Agrega K ordenes circulares en un orden consenso.

    Implementa el paso 1 de Barragan et al. (2015) con la tecnica elegida,
    TSP3 o HODs, y arbitra entre ellas por el MSCE. El paso 2 del paper no se
    implementa, por los motivos documentados en la cabecera del modulo.

    Parameters
    ----------
    Theta : np.ndarray, shape (K, n)
        ``Theta[k, i]`` es el angulo en ``[0, 2*pi)`` asignado al item ``i``
        en la repeticion ``k``.
    weights : np.ndarray opcional, shape (K,)
        Pesos por repeticion. Por defecto uniforme. El paper los define como
        ``w_j = kappa_j / sum(kappa)`` para recoger diferencias de variabilidad
        entre experimentos. En CIRCUST las K repeticiones proceden del mismo
        procedimiento sobre distintos subconjuntos de genes TOP, de modo que
        son intercambiables y el peso uniforme es el caso natural.
    method : str
        - ``"auto"`` (default): corre TSP3 y HODs y devuelve el de menor MSCE.
        - ``"tsp3"``: TSP con alpha=3.
        - ``"hods"``: Hodge score-based triplewise (Teorema 4.1).
    n_starts : int opcional
        Numero de puntos de inicio para multi-start NN (solo TSP3).

    Returns
    -------
    dict con:
        - ``order``: np.ndarray (n,), permutacion de items en orden consenso
        - ``method``: str, metodo solicitado ("auto"/"tsp3"/"hods")
        - ``selected_method``: str, metodo ganador (== method si no es "auto")
        - ``tour_length``: float o None, longitud del tour (solo si gana TSP3)
        - ``msce``: float, MSCE (CIRE) del orden ganador
        - ``msce_by_method``: dict, MSCE de cada metodo evaluado
    """
    Theta = np.asarray(Theta, dtype=float) % (2.0 * np.pi)
    K, n = Theta.shape
    if K < 1:
        raise ValueError("Se requiere al menos 1 repeticion en Theta.")
    if n < 3:
        raise ValueError("Se requieren al menos 3 items para agregacion circular.")

    method = method.lower()
    if method == "auto":
        candidates = ("tsp3", "hods")
    elif method in ("tsp3", "hods"):
        candidates = (method,)
    else:
        raise ValueError(f"method desconocido: {method!r}")

    # Resolver cada tecnica candidata y evaluar su MSCE (CIRE)
    solved = {}
    for m in candidates:
        if m == "tsp3":
            D = _build_distance_matrix_tsp3(Theta, weights)
            order_m, tour_m = _solve_tsp_multi_start(
                D, n_starts=n_starts,
                warm_starts=[_circular_mean_order(Theta, weights)],
                score_fn=lambda t: _msce_cire(Theta, t, weights=weights),
            )
        else:  # "hods"
            order_m, tour_m = _hodge_order(Theta, weights), None

        solved[m] = {
            "order":       order_m,
            "tour_length": tour_m,
            "msce":        _msce_cire(Theta, order_m, weights=weights),
        }

    # Auto-seleccion: el metodo con menor MSCE (criterio cardinal del paper)
    selected = min(solved, key=lambda m: solved[m]["msce"])
    best = solved[selected]

    return {
        "order":           best["order"],
        "method":          method,
        "selected_method": selected,
        "tour_length":     (float(best["tour_length"])
                            if best["tour_length"] is not None else None),
        "msce":            best["msce"],
        "msce_by_method":  {m: solved[m]["msce"] for m in solved},
    }
