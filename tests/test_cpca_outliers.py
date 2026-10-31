"""
Prioridad de las muestras atipicas cuando hay mas candidatas que el tope.

``max_outlier_fraction`` limita cuantas muestras se eliminan. Las candidatas
del spline se recorrian en el orden circular, de modo que, si superaban el
tope, se eliminaban las mas proximas al origen del circulo y no las mas
atipicas. Ahora se ordenan por la fraccion de genes que las marcan y, a
igualdad, por el residuo medio absoluto. Las radiales, por norma creciente.
"""
import numpy as np
import pandas as pd

from circust.cpca import CPCA


def _residuos(n_genes, n_muestras, marcas):
    """Residuos estandarizados nulos salvo en ``marcas`` = {posicion: valores por gen}."""
    res = np.zeros((n_genes, n_muestras))
    for pos, valores in marcas.items():
        res[:, pos] = valores
    genes = [f"G{g}" for g in range(n_genes)]
    return pd.DataFrame(res, index=genes), genes


def test_con_tope_se_eliminan_las_mas_atipicas_y_no_las_primeras_del_circulo():
    # Tres candidatas y tope de 2. La de la posicion 1 es la primera del
    # circulo pero la menos atipica: la marcan 2 de 4 genes con residuo 4.
    # La de la posicion 4 tiene la misma fraccion pero mas residuo, y la de
    # la 7 la marcan los 4 genes.
    df, genes = _residuos(4, 10, {
        1: [4.0, 4.0, 0.0, 0.0],
        4: [9.0, 9.0, 0.0, 0.0],
        7: [4.0, 4.0, 4.0, 4.0],
    })
    orden = np.arange(10)[::-1]          # posicion -> indice original de muestra
    cpca = CPCA(max_outlier_fraction=0.2, spline_outlier_fraction=0.5, verbose=False)
    _, spline_outs, eliminadas = cpca._detect_outliers(
        np.array([], dtype=int), orden, df, genes,
    )
    assert eliminadas == [int(orden[7]), int(orden[4])]
    assert spline_outs == eliminadas


def test_sin_tope_se_eliminan_todas_las_candidatas():
    df, genes = _residuos(4, 10, {1: [4.0, 4.0, 0, 0], 7: [4.0] * 4})
    cpca = CPCA(max_outlier_fraction=1.0, spline_outlier_fraction=0.5, verbose=False)
    _, _, eliminadas = cpca._detect_outliers(np.array([], dtype=int), np.arange(10), df, genes)
    assert sorted(eliminadas) == [1, 7]


def test_las_candidatas_radiales_salen_de_menor_a_mayor_norma():
    rng = np.random.default_rng(0)
    t = np.sort(rng.uniform(0, 2 * np.pi, 40))
    fases = rng.uniform(0, 2 * np.pi, 6)
    core = pd.DataFrame(
        [np.cos(t - f) + rng.normal(0, 0.3, t.size) for f in fases],
        index=[f"G{g}" for g in range(6)],
    )
    # Con radios enormes todas las muestras son candidatas radiales.
    cpca = CPCA(tight_radius=10.0, verbose=False)
    *_, pc1, pc2, _, _, candidatas = cpca._cpca_step(core)
    norma = np.sqrt(pc1 ** 2 + pc2 ** 2)
    assert len(candidatas) == t.size
    assert np.all(np.diff(norma[candidatas]) >= 0)
