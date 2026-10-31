"""
circust/circular_aggregation.py
================================
Consenso de repeticiones circulares y agregacion de ordenes de Barragan,
Rueda y Fernandez (2015).

Toma K repeticiones circulares sobre los mismos n items (cada una es un vector
de angulos en [0, 2*pi)) y produce un orden de consenso y una fase por item.
Combina las repeticiones de la Etapa 4 de CIRCUST. El metodo por defecto es la
media circular ponderada; TSP y Hodge son la agregacion de ordenes del paper.

Alcance respecto al paper
-------------------------
El paper organiza sus metodos en dos pasos, "a first step [where] an initial
solution denoted as O_0 is provided, which is refined in step 2".

Aqui se implementa **solo el paso 1**. Para el TSP el paper no fija una
heuristica concreta, "a better approximation to the optimum is obtained by
repeatedly running different heuristics and selecting the best solution in
terms of the objective function". :func:`_solve_tsp_multi_start` sigue el
esquema de ``isocir``, los recorridos mas cortos y desempate por MSCE, sin
la media circular entre los candidatos porque es un metodo aparte (ver su
docstring).

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
- **Media circular** (``method="circmean"``, defecto): media circular de las K
  fases de cada item ponderada por kappa. Ver ``aggregate_circular_orders``.

- **TSP** (``method="tsp"``): TSP con distancia direccional asimetrica y
  penalizacion ``tsp_alpha`` para movimientos en la direccion equivocada. Por
  defecto alpha=3 (TSP3), recomendado por el paper para el criterio MSCE
  (Mean Sum of Circular Errors). Los costes de arista siguen la Tabla 1 y la
  ecuacion (7) del paper.

- **Hodge** (``method="hod"``, variante ``hodge_variant``, por defecto
  HODs): Hodge score-based triplewise (Tabla 2 y Teorema 4.1). El flujo
  solo es no nulo en los tripletes que contienen el ancla ``l``, donde vale la
  diferencia de angulos rotados a ``theta_l = 0``. Los angulos se agregan entre
  experimentos por media circular y el ancla se elige por ``gamma_l``. Es la
  variante ``nat`` de ``isocir::ACO(method="CH")`` (que fija el ancla en el
  primer item). Las cuatro variantes coinciden con las de ``isocir`` y con el
  Teorema 4.1 aplicado literalmente (``tests/test_hods.py``).
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
bloques muy anchos, irrelevantes para ordenar metodos). Se usa para elegir el
recorrido del TSP, para reconstruir las fases bajo el orden de TSP o Hodge y
como diagnostico, no como cifra a reportar.

API publica
-----------
- ``aggregate_circular_orders(Theta, weights, method="circmean")``: orden de
  consenso y fase de cada item. Metodos ``"circmean"`` (defecto, media
  circular ponderada), ``"tsp"`` (con ``tsp_alpha``) y ``"hod"`` (con
  ``hodge_variant``).
- ``estimate_kappa_pairwise(Theta)`` y ``kappa_weights(Theta)``: concentracion
  de cada repeticion y pesos ``w_j = kappa_j / sum(kappa)`` del paper.
- ``circular_isotonic_regression(theta, order, weights)``: CIRE (ponderado) de
  un vector angular bajo un orden dado.

Equivalente en R: paquete ``isocir`` (funciones ``CIRE`` y ``ACO``).
"""
from __future__ import annotations

import itertools
from typing import Optional
import numpy as np
from scipy.optimize import brentq
from scipy.special import i0e, i1e

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
def _cire_core(d: np.ndarray, S: np.ndarray, C: np.ndarray, r: np.ndarray):
    """Nucleo del CIRE circular sobre datos YA ordenados ``d``.

    Enumera el arco que forma el primer bloque (corte ``c`` x longitud ``m``),
    fija su media circular ``mu0`` como referencia de una sola vuelta, resuelve
    el resto con una PAVA en la ventana ``[mu0, mu0+2pi)`` (pooling por media
    circular, violacion cuando la media-ventana decrece) y se queda con el
    minimo SCE. Devuelve ``(fitted_posiciones, SCE)``.

    ``r`` son pesos por item: el coste es ``sum r_i (1 - cos(d_i - mu))`` y
    los bloques se agrupan por la direccion de ``sum S``, ``sum C``, que deben
    venir ya ponderados (``S_i = r_i sin d_i``). Con ``r = 1`` es el CIRE de
    Rueda et al. (2009).

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
                cost += r[idx] * (1.0 - np.cos(d[idx] - mu0))
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
                    cost += r[idx] * (1.0 - np.cos(d[idx] - mu))
                    pos += 1
            if cost < best_sce - 1e-12:
                best_sce = cost
                for i in range(q):
                    best_fit[i] = cur_fit[i]
    return best_fit, best_sce


def circular_isotonic_regression(
    theta: np.ndarray, order: np.ndarray, weights: Optional[np.ndarray] = None,
) -> tuple[np.ndarray, float]:
    """CIRE de ``theta`` bajo el orden circular dado por ``order``.

    Parameters
    ----------
    theta : np.ndarray (n,)
        Angulos observados en ``[0, 2*pi)`` (indice por item).
    order : np.ndarray (n,)
        Permutacion de items: ``order[pos]`` es el item en la posicion ``pos``
        del orden circular candidato.
    weights : np.ndarray (n,) opcional
        Peso de cada item en el coste, ``sum w_i (1 - cos(theta_i - fit_i))``.
        Por defecto 1, el CIRE de Rueda et al. (2009).

    Returns
    -------
    (fitted, sce) : np.ndarray (n,), float
        ``fitted`` es el CIRE alineado al indice de ``theta``; ``sce`` la suma
        (ponderada) de errores circulares.
    """
    theta = np.asarray(theta, dtype=float) % TWO_PI
    order = np.asarray(order, dtype=np.int64)
    n = theta.shape[0]
    d = np.ascontiguousarray(theta[order])
    r = (np.ones(n) if weights is None
         else np.ascontiguousarray(np.asarray(weights, dtype=float)[order]))
    fit_pos, sce = _cire_core(d, r * np.sin(d), r * np.cos(d), r)
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
        _, sce = _cire_core(d, np.sin(d), np.cos(d), np.ones(n))
        msces[k] = sce / n
    if weights is None:
        return float(msces.mean())
    w = np.asarray(weights, dtype=float)
    if w.shape[0] != K:
        raise ValueError(f"weights debe tener longitud {K}, got {w.shape[0]}")
    return float(np.average(msces, weights=w))


# ═══════════════════════════════════════════════════════════════════════════
# Seccion 3 — TSPalpha (por defecto alpha=3, recomendado para MSCE)
# ═══════════════════════════════════════════════════════════════════════════

def _check_tsp_alpha(alpha: float) -> float:
    """El paper define la penalizacion en ``alpha >= 1`` (``inf`` incluido)."""
    alpha = float(alpha)
    if np.isnan(alpha) or alpha < 1.0:
        raise ValueError(f"tsp_alpha debe ser >= 1 (puede ser inf), recibido {alpha!r}")
    return alpha


def _distance_matrix_tsp_single(theta: np.ndarray, alpha: float = 3.0) -> np.ndarray:
    """Matriz de distancias TSPalpha ``(n, n)`` para una sola repeticion.

    Para cada par (h, k), ecuacion (7) y Tabla 1 del paper:

      d_R(h, k) = 1 - cos(theta_k - theta_h)             si  diff <= pi
                  3 - cos(theta_k - theta_h - pi)        si  diff > pi
      d_C(h, k) = d_R(k, h)
      E[h, k]   = min(d_R(h, k), alpha * d_C(h, k))

    donde ``diff = (theta_k - theta_h) mod 2*pi``. Con ``alpha = inf`` queda
    ``E = d_R``. Por defecto ``alpha = 3`` (TSP3), la variante que recomienda
    el paper. En simulacion TSP2, TSP3 y TSP4 empatan, ``alpha = 1`` (matriz
    simetrica) pierde el sentido de giro y ``alpha = inf`` falla en cuanto el
    orden tiene un error, porque corregirlo obliga a dar casi una vuelta.
    """
    alpha = _check_tsp_alpha(alpha)
    theta = np.asarray(theta, dtype=float)
    diff  = (theta[None, :] - theta[:, None]) % (2.0 * np.pi)   # (n, n)
    # d_R: distancia en direccion "directa" (corto arco) o penalizada si hay que dar la vuelta
    d_R = np.where(
        diff <= np.pi,
        1.0 - np.cos(diff),
        3.0 - np.cos(diff - np.pi),
    )
    d_C = d_R.T                                                  # d_C(h,k) = d_R(k,h)
    E = d_R.copy() if np.isinf(alpha) else np.minimum(d_R, alpha * d_C)
    np.fill_diagonal(E, 0.0)
    return E


def _build_distance_matrix_tsp(
    Theta:   np.ndarray,
    weights: Optional[np.ndarray] = None,
    alpha:   float = 3.0,
) -> np.ndarray:
    """Matriz de distancias TSPalpha agregada sobre las K repeticiones.

    ``E[h, k] = sum_j w_j * E^j[h, k]`` donde ``E^j`` es la matriz de la
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
        E += weights[j] * _distance_matrix_tsp_single(Theta[j], alpha)
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


def _canonical_cycle(tour: np.ndarray) -> tuple:
    """El mismo ciclo dirigido empezado por el item 0, para detectar repetidos."""
    i = int(np.flatnonzero(tour == 0)[0])
    return tuple(np.roll(tour, -i).tolist())


def _solve_tsp_multi_start(
    D:           np.ndarray,
    n_starts:    Optional[int] = None,
    or_opt:      bool = True,
    score_fn:    Optional[callable] = None,
    n_shortest:  int = 5,
) -> tuple[np.ndarray, float]:
    """Paso 1 de Barragan et al. (2015): solucion inicial ``O_0`` del TSP.

    En el paper ``O_0`` es la solucion de (6), el recorrido mas corto sobre la
    matriz agregada, y para aproximarla recomiendan ejecutar varias heuristicas
    y quedarse con la mejor. ``isocir`` (``TSPmethod`` en ``ACO.R``) lo hace en
    dos fases: se queda con los ``n`` recorridos de longitud distinta mas
    cortos y entre ellos elige el de menor MSCE. Aqui se sigue ese esquema.

    Los candidatos son los ``n_starts`` recorridos del vecino mas cercano desde
    distintos origenes, pulidos con ``_or_opt``. De ellos pasan a la lista los
    ``n_shortest`` de menor longitud, contando una sola vez los que son el
    mismo ciclo, y entre los de la lista gana el de menor ``score_fn``, que
    ``aggregate_circular_orders`` fija en el MSCE.

    El orden de la media circular de las K repeticiones no entra entre los
    candidatos. Hasta octubre de 2026 entraba siempre, porque con ruido alto
    la longitud del recorrido se sobreajusta (el optimo de (6) encadena
    aristas cuyo coste bajo es fruto del azar) y ese orden, mas largo, era
    mucho mejor. Pero la media circular es ahora un metodo de consenso propio
    (``method="circmean"``), y dejarla dentro del TSP mezclaba los dos metodos
    y hacia imposible compararlos. Lo que cuesta quitarla, medido con K=5 y
    12 replicas (tau circular frente a la verdad, media circular contra TSP):

        n=48,  sigma=0,05   0,97 contra 0,97
        n=48,  sigma=0,15   0,93 contra 0,84
        n=100, sigma=0,15   0,93 contra 0,45
        n=100, sigma=0,30   0,86 contra 0,29

    Por que ``n_shortest=5`` y no los ``n`` de ``isocir``: en esas mismas
    simulaciones 10 candidatos dan lo mismo que 5 (diferencias de 0,01 de
    tau como mucho) y 30 dan bastante menos (0,14 a 0,61), porque el MSCE
    tambien acaba eligiendo recorridos sobreajustados. Ademas cada candidato
    cuesta K evaluaciones del CIRE, de coste cubico en el numero de muestras.

    Sin ``score_fn`` devuelve simplemente el recorrido mas corto, que es la
    solucion de (6). Devuelve ``(tour, longitud)``.
    """
    n = D.shape[0]
    if n_starts is None:
        n_starts = min(n, 30)               # 30 starts es suficiente en practica
    n_starts = max(1, min(n_starts, n))
    pulir = (lambda t: _or_opt(D, t)) if or_opt else (lambda t: t)

    # Vecino mas cercano desde distintos origenes, sin ciclos repetidos.
    nn: dict[tuple, np.ndarray] = {}
    for s in np.linspace(0, n - 1, n_starts, dtype=int):
        tour = pulir(_solve_tsp_nearest_neighbor(D, start=int(s))[0])
        nn.setdefault(_canonical_cycle(tour), tour)
    nn_tours = sorted(nn.values(), key=lambda t: _tour_length(D, t))

    if score_fn is None:
        best = nn_tours[0]
        return best, _tour_length(D, best)

    lista = nn_tours[:max(1, n_shortest)]
    best = min(lista, key=lambda t: float(score_fn(t)))
    return best, _tour_length(D, best)


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


HODGE_VARIANTS = ("s", "p", "b", "r")


def _position_angles(Theta: np.ndarray) -> np.ndarray:
    """Angulo de la posicion de cada item en cada repeticion, ``2 pi tau_ij / n``.

    Es la entrada de HODp, que solo usa el orden de cada repeticion: Tabla 2,
    ``lambda = (2 pi / n) |tau_h - tau_k|``, y la variante ``natp`` de ``isocir``.
    """
    K, n = Theta.shape
    tau = np.argsort(np.argsort(Theta, axis=1, kind="stable"), axis=1)
    return TWO_PI * tau / n


def _hodge_Y_all_triplets(Theta: np.ndarray, w: np.ndarray, variant: str) -> np.ndarray:
    """``Ybar = (1/n) delta_1*(Psibar)`` para HODb y HODr.

    En estas variantes ``psi_ihk = sign(i,h,k) * lambda_ihk`` esta definido en
    todos los tripletes (ecuacion 8 y Tabla 2), con ``lambda = 1`` en HODb y
    ``lambda = 1 - Rbar_ihk``, la longitud resultante media de los tres
    angulos, en HODr. ``sign(i,h,k) = sign(t_h - t_i) + sign(t_k - t_h) +
    sign(t_i - t_k)`` vale +1 o -1 segun el sentido de giro del triplete y no
    depende del origen del circulo. Las repeticiones se agregan con la media
    aritmetica ponderada de las ``Psi^j``, como en ``hodgefusion`` de
    ``isocir``. Se recorre fila a fila para no guardar la hipermatriz ``n^3``.
    """
    K, n = Theta.shape
    Y = np.zeros((n, n))
    for j in range(K):
        th = Theta[j]
        sg = np.sign(th[None, :] - th[:, None])              # sg[a, b] = sign(t_b - t_a)
        e = np.exp(1j * th)
        for i in range(n):
            # psi[h, k] = sign(t_h - t_i) + sign(t_k - t_h) + sign(t_i - t_k)
            psi = sg[i][:, None] + sg + sg[:, i][None, :]
            if variant == "r":
                psi = psi * (1.0 - np.abs(e[i] + e[:, None] + e[None, :]) / 3.0)
            Y[i] += w[j] * psi.sum(axis=1)
    return Y / n


def _hodge_order(
    Theta:   np.ndarray,
    weights: Optional[np.ndarray] = None,
    anchor:  Optional[int] = None,
    variant: str = "s",
) -> np.ndarray:
    """Orden consenso por Hodge (Barragan et al., Seccion 4 y Tabla 2).

    ``variant`` elige el flujo triplewise ``lambda``:

    - ``"s"`` (HODs, defecto): diferencias de fase con el ancla. Ver abajo.
    - ``"p"`` (HODp): igual que HODs sobre los angulos de posicion
      ``2 pi tau / n``, de modo que solo usa el orden de cada repeticion.
    - ``"b"`` (HODb): ``lambda = 1``, solo el sentido de giro de cada triplete.
    - ``"r"`` (HODr): ``lambda = 1 - Rbar``, mas peso a los tripletes cuyos
      tres angulos estan mas separados.

    La eleccion del ancla depende de la familia. En HODb y HODr la hipermatriz
    no depende de ningun ancla y ``l0`` sale del Teorema 4.1,
    ``l0 = argmax_l sum_h Ybar_lh^2``. En HODs y HODp ``lambda`` solo es
    distinta de cero en los tripletes que contienen al ancla, de modo que hay
    que fijarla antes, y el paper propone ``l0 = argmax_l gamma_l``.

    Lo que sigue describe HODs, la variante por defecto.

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
    variant = variant.lower()
    if variant not in HODGE_VARIANTS:
        raise ValueError(f"hodge_variant debe ser uno de {HODGE_VARIANTS}, recibido {variant!r}")
    Theta = np.asarray(Theta, dtype=float) % TWO_PI
    K, n = Theta.shape
    if weights is None:
        w = np.ones(K) / K
    else:
        w = np.asarray(weights, dtype=float)
        w = w / w.sum()

    if variant in ("b", "r"):
        Y = _hodge_Y_all_triplets(Theta, w, variant)
        if anchor is None:
            l0 = int(np.argmax((Y ** 2).sum(axis=1)))           # Teorema 4.1
        else:
            l0 = int(anchor)
            if not 0 <= l0 < n:
                raise ValueError(f"anchor debe estar en [0, {n}), got {l0}")
        keep = np.array([i for i in range(n) if i != l0], dtype=np.int64)
        s = -Y[np.ix_(keep, keep)].sum(axis=1) / (n - 1)        # s_i del Teorema 4.1
        return np.concatenate([[l0], keep[np.argsort(s, kind="stable")]]).astype(np.int64)

    if variant == "p":
        Theta = _position_angles(Theta)
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
# Seccion 3b — Pesos por repeticion: kappa por pares
# ═══════════════════════════════════════════════════════════════════════════

# Tope de kappa. A(kappa) se acerca a 1 tan deprisa que por encima de este
# valor la inversa deja de estar bien condicionada en doble precision.
KAPPA_MAX = 1e5


def _a1(k: float) -> float:
    """``A(kappa) = I1(kappa) / I0(kappa)``.

    Se evalua con las funciones de Bessel escaladas ``i1e`` e ``i0e``: el
    factor ``exp(-kappa)`` se cancela en el cociente y ninguna desborda, al
    contrario que ``i0`` e ``i1`` por encima de kappa ~ 700.
    """
    return float(i1e(k) / i0e(k)) if k > 0 else 0.0


def _a1_inv(r: float) -> float:
    """Inversa de ``A`` en ``[0, KAPPA_MAX]``. Devuelve 0 si ``r <= 0``."""
    if r <= 0.0:
        return 0.0
    if r >= _a1(KAPPA_MAX):
        return KAPPA_MAX
    return float(brentq(lambda k: _a1(k) - r, 1e-12, KAPPA_MAX))


def estimate_kappa_pairwise(Theta: np.ndarray) -> np.ndarray:
    """Concentracion ``kappa_j`` de cada repeticion, por comparacion por pares.

    Modelo de Barragan et al. (2015), con un giro libre por repeticion:
    ``theta_ij = phi_i + beta_j + e_ij`` con ``e_ij ~ vM(0, kappa_j)``
    independientes. Al restar dos repeticiones sobre el mismo item la fase
    ``phi_i`` se cancela, y como el primer momento trigonometrico de una suma
    de variables independientes es el producto de los de cada una,

        |E exp(i (theta_ij - theta_il))| = A(kappa_j) A(kappa_l).

    El modulo elimina ademas el giro ``beta_j - beta_l``. Con logaritmos cada
    par da una ecuacion lineal, ``log rho_jl = a_j + a_l`` con
    ``a_j = log A(kappa_j)``, y las ``K (K-1) / 2`` ecuaciones se resuelven por
    minimos cuadrados. Hace falta ``K >= 3``. Con ``K = 3`` la solucion es
    cerrada: ``A(kappa_1) = sqrt(rho_12 rho_13 / rho_23)``.

    Por que asi y no como Fernandez et al. (2012, Supp. 1.3), de donde
    Barragan toma los kappa: esa estimacion maximiza la verosimilitud del
    modelo con un ``alpha_i`` libre por item, y esa verosimilitud no esta
    acotada (si ``alpha_i`` copia una repeticion, su residuo es cero y su kappa
    crece sin limite). Iterada hasta converger, el kappa de la repeticion mas
    concentrada diverge. Aqui ``phi_i`` no se estima nunca, de modo que el
    problema no aparece y el estimador es consistente al crecer ``n``.

    Parameters
    ----------
    Theta : np.ndarray, shape (K, n)
        Angulo del item ``i`` en la repeticion ``j``.

    Returns
    -------
    np.ndarray (K,) con los kappa estimados, en ``[0, KAPPA_MAX]``.
    """
    Theta = np.asarray(Theta, dtype=float)
    K = Theta.shape[0]
    if K < 3:
        raise ValueError(
            f"La estimacion por pares necesita al menos 3 repeticiones, hay {K}."
        )
    pares = list(itertools.combinations(range(K), 2))
    X = np.zeros((len(pares), K))
    y = np.empty(len(pares))
    for p, (j, l) in enumerate(pares):
        d = Theta[j] - Theta[l]
        rho = np.hypot(np.cos(d).mean(), np.sin(d).mean())
        X[p, j] = X[p, l] = 1.0
        y[p] = np.log(max(rho, 1e-12))
    a, *_ = np.linalg.lstsq(X, y, rcond=None)
    # a_j = log A(kappa_j) <= 0. El ruido de muestreo puede darlo positivo en
    # repeticiones casi perfectas, y entonces kappa sale al tope.
    return np.array([_a1_inv(float(np.exp(min(aj, 0.0)))) for aj in a])


def kappa_weights(Theta: np.ndarray) -> np.ndarray:
    """Pesos del paper, ``w_j = kappa_j / sum(kappa)``, con kappa por pares.

    Si ningun kappa es positivo (repeticiones sin relacion entre si) no hay
    informacion para ponderar y se devuelven pesos iguales.
    """
    kappa = estimate_kappa_pairwise(Theta)
    total = kappa.sum()
    if total <= 0.0:
        return np.full(kappa.shape[0], 1.0 / kappa.shape[0])
    return kappa / total


def _resolve_weights(Theta: np.ndarray, weights: Optional[np.ndarray]) -> np.ndarray:
    """Pesos normalizados a suma 1: los dados, o los de kappa por pares si None.

    Con menos de 3 repeticiones kappa no es estimable por pares y, sin otra
    informacion, se usan pesos iguales.
    """
    K = Theta.shape[0]
    if weights is None:
        return kappa_weights(Theta) if K >= 3 else np.full(K, 1.0 / K)
    w = np.asarray(weights, dtype=float)
    if w.shape != (K,):
        raise ValueError(f"weights debe tener longitud {K}, got {w.shape}")
    if np.any(w < 0) or w.sum() <= 0:
        raise ValueError("weights debe ser no negativo y de suma > 0.")
    return w / w.sum()


# ═══════════════════════════════════════════════════════════════════════════
# Seccion 4 — API publica
# ═══════════════════════════════════════════════════════════════════════════

AGGREGATION_METHODS = ("circmean", "tsp", "hod")
_ALIAS_METODO = {"tsp3": "tsp", "hods": "hod"}


def _weighted_circular_mean(Theta: np.ndarray, w: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Media circular ponderada de cada item y su longitud resultante ``R_i``."""
    S = (w[:, None] * np.sin(Theta)).sum(axis=0)
    C = (w[:, None] * np.cos(Theta)).sum(axis=0)
    return np.arctan2(S, C) % TWO_PI, np.hypot(S, C)


def aggregate_circular_orders(
    Theta:      np.ndarray,
    weights:    Optional[np.ndarray] = None,
    method:     str = "circmean",
    n_starts:   Optional[int] = None,
    tsp_alpha:  float = 3.0,
    hodge_variant: str = "s",
) -> dict:
    """Consenso de K repeticiones circulares: orden de los items y fase de cada uno.

    Tres metodos:

    - ``"circmean"`` (defecto): la fase de cada item es la media circular de
      sus K angulos ponderada por ``weights``, y el orden es el que inducen esas
      fases. Si ``theta_ij ~ vM(phi_i, kappa_j)`` con una fase comun por item,
      es el estimador de maxima verosimilitud de ``phi_i``. Requiere que las
      repeticiones compartan origen y sentido, cosa que en CIRCUST garantiza
      la sincronizacion.
    - ``"tsp"``: agregacion de ordenes de Barragan et al. (2015) por TSP, con
      penalizacion ``tsp_alpha``. ``"tsp3"`` se acepta como sinonimo.
    - ``"hod"``: agregacion de ordenes de Barragan et al. por Hodge, con la
      variante ``hodge_variant``. ``"hods"`` se acepta como sinonimo.

    TSP y Hodge solo producen un orden, porque el paper no reconstruye fases.
    Aqui la fase de cada item es la de maxima verosimilitud bajo ese orden: el
    CIRE de la media circular ponderada, con cada item pesado por su longitud
    resultante ``R_i``. Desarrollando ``sum_j w_j sum_i (1 - cos(theta_ij -
    phi_i))`` queda ``sum_i (1 - R_i cos(theta_bar_i - phi_i))``, de modo que
    los items en que las repeticiones coinciden apenas se mueven.

    Por que ``circmean`` es el defecto: en simulacion, con repeticiones
    sincronizadas, ninguna forma de imponer el orden de TSP o Hodge a las fases
    mejoro el error de fase de la media circular, con ningun patron de
    muestreo (equiespaciado, uniforme, concentrado, bimodal, con hueco o con
    tiempos repetidos), y el orden que induce la media circular era mas exacto
    que el de consenso. La agregacion de ordenes esta pensada para
    experimentos cuyas fases no son comparables entre si.

    El paso 2 del paper (CLMA) no se implementa, por los motivos documentados
    en la cabecera del modulo.

    Parameters
    ----------
    Theta : np.ndarray, shape (K, n)
        ``Theta[k, i]`` es el angulo en ``[0, 2*pi)`` asignado al item ``i``
        en la repeticion ``k``.
    weights : np.ndarray opcional, shape (K,)
        Pesos por repeticion, que se normalizan a suma 1. Si no se dan, se usan
        los del paper, ``w_j = kappa_j / sum(kappa)``, con kappa estimado por
        pares (:func:`estimate_kappa_pairwise`).
    method : str
        ``"circmean"``, ``"tsp"`` o ``"hod"``.
    n_starts : int opcional
        Numero de puntos de inicio para multi-start NN (solo TSP).
    tsp_alpha : float
        Penalizacion ``alpha >= 1`` de la ecuacion (7) del paper, que fija la
        variante TSPalpha (``inf`` admitido). Por defecto 3, TSP3.
    hodge_variant : str
        Flujo triplewise de la Tabla 2 del paper: ``"s"`` (HODs, defecto),
        ``"p"``, ``"b"`` o ``"r"``. Ver :func:`_hodge_order`.

    Returns
    -------
    dict con:
        - ``order``: np.ndarray (n,), permutacion de items en orden consenso
        - ``phases``: np.ndarray (n,), fase de consenso de cada item (indice
          por item), que verifica ``order`` salvo en ``"circmean"``, donde el
          orden se deriva de ella
        - ``method``: str, metodo usado
        - ``tour_length``: float o None, longitud del recorrido (solo TSP)
        - ``msce``: float, MSCE (CIRE) del orden, como diagnostico
        - ``weights``: np.ndarray (K,), pesos usados, normalizados a suma 1
    """
    Theta = np.asarray(Theta, dtype=float) % (2.0 * np.pi)
    K, n = Theta.shape
    if K < 1:
        raise ValueError("Se requiere al menos 1 repeticion en Theta.")
    if n < 3:
        raise ValueError("Se requieren al menos 3 items para agregacion circular.")
    method = _ALIAS_METODO.get(method.lower(), method.lower())
    if method not in AGGREGATION_METHODS:
        raise ValueError(f"method debe ser uno de {AGGREGATION_METHODS}, recibido {method!r}")
    weights = _resolve_weights(Theta, weights)
    tsp_alpha = _check_tsp_alpha(tsp_alpha)
    hodge_variant = hodge_variant.lower()
    if hodge_variant not in HODGE_VARIANTS:
        raise ValueError(f"hodge_variant debe ser uno de {HODGE_VARIANTS}, recibido {hodge_variant!r}")

    media, R = _weighted_circular_mean(Theta, weights)
    tour_length = None
    if method == "circmean":
        order = np.argsort(media, kind="stable").astype(np.int64)
        phases = media
    else:
        if method == "tsp":
            D = _build_distance_matrix_tsp(Theta, weights, tsp_alpha)
            order, tour_length = _solve_tsp_multi_start(
                D, n_starts=n_starts,
                score_fn=lambda t: _msce_cire(Theta, t, weights=weights),
            )
            tour_length = float(tour_length)
        else:
            order = _hodge_order(Theta, weights, variant=hodge_variant)
        phases, _ = circular_isotonic_regression(media, order, weights=R)

    return {
        "order":       order,
        "phases":      phases,
        "method":      method,
        "tour_length": tour_length,
        "msce":        _msce_cire(Theta, order, weights=weights),
        "weights":     weights,
    }
