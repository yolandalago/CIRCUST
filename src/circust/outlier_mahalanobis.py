r"""
circust/outlier_mahalanobis.py
================================
Criterio de deteccion de outliers via distancia de Mahalanobis robusta
(MinCovDet) sobre los residuos estandarizados de los genes core.

NO esta conectado al pipeline de ``CPCA`` (que usa el criterio de
spline periodico + porcentaje de genes, ver ``fitting/cyclic_spline.py``
y ``CPCA._detect_outliers``). Se conserva aqui como alternativa evaluada
y descartada como metodo principal del TFG, por si se quiere retomar la
comparacion en las simulaciones.

Motivo de no ser el metodo principal (memoria de proyecto):
- Necesita n>=5p (idealmente 10p) para estimar la covarianza robusta con
  garantias; con p=12 genes core eso son 60-120 muestras, inviable en
  los casos reales mas importantes del TFG (babuino, n=11-12).
- El umbral chi2 teorico esta mal calibrado aqui (los residuos FMM no
  son realmente normales multivariantes, colas mucho mas pesadas de lo
  que predice la teoria); el umbral empirico alternativo evita eso, pero
  siempre recorta un porcentaje fijo de muestras aunque el lote sea
  limpio, cosa que el criterio final (spline + porcentaje de genes) no
  hace.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import chi2
from sklearn.covariance import MinCovDet


def mahalanobis_empirical_flags(
    std_res_df: pd.DataFrame,
    genes: list[str],
    alpha: float = 0.01,
    seed: int | None = None,
) -> np.ndarray:
    """MinCovDet + distancia de Mahalanobis sobre los residuos
    estandarizados de los genes core (sin reducir dimension). Umbral:
    percentil empirico (1-alpha) de las distancias observadas en la
    propia muestra (no un cuantil chi2 teorico, mal calibrado aqui).
    """
    core_res = std_res_df.loc[genes].values.T  # (n_samples, n_genes)
    mcd = MinCovDet(random_state=seed).fit(core_res)
    dist2 = mcd.mahalanobis(core_res)
    threshold = np.percentile(dist2, 100 * (1 - alpha))
    return dist2 > threshold


def mahalanobis_chi2_flags(
    std_res_df: pd.DataFrame,
    genes: list[str],
    alpha: float = 0.01,
    seed: int | None = None,
) -> np.ndarray:
    """Igual que ``mahalanobis_empirical_flags`` pero con el umbral
    chi2 teorico (df=p). Se conserva solo como referencia: en la
    validacion disparaba muchos falsos positivos porque los residuos no
    son realmente normales multivariantes aqui.
    """
    core_res = std_res_df.loc[genes].values.T
    mcd = MinCovDet(random_state=seed).fit(core_res)
    dist2 = mcd.mahalanobis(core_res)
    threshold = chi2.ppf(1 - alpha, df=core_res.shape[1])
    return dist2 > threshold
