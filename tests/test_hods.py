"""
Tests de HODs, la variante score-based del metodo de Hodge de Barragan et al.

``_hodge_order`` no construye la hipermatriz n x n x n: usa que, para HODs,
el Teorema 4.1 se reduce a ordenar las medias circulares de los angulos
rotados al ancla. Aqui se aplica el teorema tal cual y se comprueba que el
orden coincide.
"""
import numpy as np

from circust.circular_aggregation import _hodge_anchored_means, _hodge_order

TWO_PI = 2.0 * np.pi


def _hods_literal(Theta, w):
    """Seccion 4 del paper: hipermatriz (ec. 8 y Tabla 2), delta_1* y Teorema 4.1."""
    n = Theta.shape[1]
    A = _hodge_anchored_means(Theta, w)

    def psi(l):
        a = A[l]
        P = np.zeros((n, n, n))
        for h in range(n):
            for k in range(n):
                if len({l, h, k}) < 3:
                    continue
                signo = np.sign(a[h] - a[l]) + np.sign(a[k] - a[h]) + np.sign(a[l] - a[k])
                v = signo * abs(a[h] - a[k])
                for (x, y, z), s in (((l, h, k), 1), ((h, k, l), 1), ((k, l, h), 1),
                                     ((h, l, k), -1), ((k, h, l), -1), ((l, k, h), -1)):
                    P[x, y, z] = s * v
        return P

    gamma = [(psi(l)[l] ** 2).sum() for l in range(n)]
    l0 = int(np.argmax(gamma))
    Y = psi(l0).sum(axis=2) / n
    resto = [i for i in range(n) if i != l0]
    s = {i: -sum(Y[i, h] for h in resto) / (n - 1) for i in resto}
    return np.array([l0] + sorted(resto, key=lambda i: s[i]))


def test_hods_coincide_con_teorema_4_1_literal():
    rng = np.random.default_rng(1)
    for _ in range(25):
        n = int(rng.integers(4, 10))
        K = int(rng.integers(2, 6))
        phi = rng.uniform(0, TWO_PI, n)
        Theta = (phi + rng.normal(0, rng.uniform(0, 1.5), (K, n))) % TWO_PI
        w = rng.uniform(0.2, 1, K)
        w /= w.sum()
        assert np.array_equal(_hods_literal(Theta, w), _hodge_order(Theta, w))


def test_sin_ruido_recupera_el_orden_verdadero():
    rng = np.random.default_rng(2)
    phi = rng.uniform(0, TWO_PI, 30)
    orden = _hodge_order(np.tile(phi, (4, 1)))
    verdad = np.argsort(phi)
    i = int(np.flatnonzero(verdad == orden[0])[0])
    assert np.array_equal(np.roll(verdad, -i), orden)


# ── Variantes de la Tabla 2: HODb, HODr y HODp ─────────────────────────────

def _hod_todos_los_tripletes_literal(Theta, w, variante):
    """HODb y HODr: hipermatriz completa (ec. 8), media aritmetica y Teorema 4.1."""
    K, n = Theta.shape
    Psi = np.zeros((n, n, n))
    for j in range(K):
        t = Theta[j]
        for i in range(n):
            for h in range(n):
                for k in range(n):
                    signo = np.sign(t[h] - t[i]) + np.sign(t[k] - t[h]) + np.sign(t[i] - t[k])
                    lam = 1.0 if variante == "b" else 1.0 - abs(np.exp(1j * t[[i, h, k]]).sum()) / 3
                    Psi[i, h, k] += w[j] * signo * lam
    Y = Psi.sum(axis=2) / n
    l0 = int(np.argmax((Y ** 2).sum(axis=1)))
    resto = [i for i in range(n) if i != l0]
    s = {i: -sum(Y[i, h] for h in resto) / (n - 1) for i in resto}
    return np.array([l0] + sorted(resto, key=lambda i: s[i]))


def test_hodb_y_hodr_coinciden_con_el_teorema_literal():
    """Se compara el ciclo y no la permutacion: en HODb las puntuaciones son
    sumas de signos y dos anclas pueden empatar exactamente en gamma. Cada
    implementacion rompe el empate segun su redondeo, pero el orden circular
    es el mismo, solo cambia la muestra por la que se empieza a leer."""
    from circust.circular_aggregation import _canonical_cycle
    rng = np.random.default_rng(3)
    for variante in ("b", "r"):
        for _ in range(10):
            n = int(rng.integers(4, 9))
            K = int(rng.integers(2, 5))
            Theta = (rng.uniform(0, TWO_PI, n) + rng.normal(0, 0.8, (K, n))) % TWO_PI
            w = rng.uniform(0.2, 1, K)
            w /= w.sum()
            assert _canonical_cycle(_hod_todos_los_tripletes_literal(Theta, w, variante)) == \
                   _canonical_cycle(_hodge_order(Theta, w, variant=variante))


def test_hodp_es_hods_sobre_las_posiciones():
    from circust.circular_aggregation import _position_angles
    rng = np.random.default_rng(4)
    Theta = (rng.uniform(0, TWO_PI, 15) + rng.normal(0, 0.5, (4, 15))) % TWO_PI
    assert np.array_equal(_hodge_order(Theta, variant="p"),
                          _hodge_order(_position_angles(Theta), variant="s"))


def test_alias_y_variante_invalida():
    import pytest
    from circust.circular_aggregation import aggregate_circular_orders
    rng = np.random.default_rng(6)
    Theta = (rng.uniform(0, TWO_PI, 12) + rng.normal(0, 0.3, (4, 12))) % TWO_PI
    a = aggregate_circular_orders(Theta, method="hods")
    b = aggregate_circular_orders(Theta, method="hod")
    assert a["method"] == b["method"] == "hod"
    assert np.array_equal(a["order"], b["order"])
    with pytest.raises(ValueError):
        aggregate_circular_orders(Theta, method="hod", hodge_variant="x")
