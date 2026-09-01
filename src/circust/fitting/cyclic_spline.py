r"""
circust/fitting/cyclic_spline.py
==================================
Spline cubico periodico con nodos uniformes, usado para la deteccion de
outliers en CPCA (independiente del modelo FMM que se usa despues en el
pipeline, por peticion explicita de diseno: no conviene detectar
outliers con el mismo modelo que luego se ajustara a los datos).

Construccion: ``sklearn.preprocessing.SplineTransformer`` genera la
matriz de base del spline con nodos explicitos uniformes y
``extrapolation="periodic"`` (la propia documentacion de sklearn
recomienda fijar los nodos a mano en el caso periodico, en vez de
dejar que se calculen del rango muestral). El ajuste a los datos es una
regresion lineal ordinaria sobre esa base. Con pocos nodos (5 por
defecto) el spline tiene muchos menos grados de libertad que muestras,
asi que no puede pasar por cada punto (a diferencia de un spline
interpolante), y una desviacion real deja un residuo visible.

Validado en ``validation/spline_outlier_prototype.py`` con contaminacion
sintetica controlada: mantiene correctamente muestras con una desviacion
aislada de un solo gen y detecta desviaciones coordinadas entre varios
genes, a partir de n>=24 aprox. (por debajo de eso, ver limitacion
conocida en la memoria de proyecto "fmm-sobreajuste-oculta-outliers-n-bajo",
que aplica igual de mal al spline que a FMM en n tan bajo).
"""
from __future__ import annotations

import numpy as np
from sklearn.preprocessing import SplineTransformer

PERIOD = 2 * np.pi


def cyclic_spline_basis(x: np.ndarray, n_knots: int = 5) -> np.ndarray:
    """Matriz de diseño (n_muestras, n_knots-1) del spline cubico
    periodico con nodos explicitos, uniformemente repartidos en el
    ciclo completo [0, 2*pi] (no en el rango muestral de ``x``).
    """
    explicit_knots = np.linspace(0.0, PERIOD, n_knots).reshape(-1, 1)
    st = SplineTransformer(degree=3, knots=explicit_knots, extrapolation="periodic")
    x_mod = np.mod(x, PERIOD).reshape(-1, 1)
    st.fit(x_mod)
    return st.transform(x_mod)


def fit_cyclic_spline_residuals(
    y: np.ndarray, circular_scale: np.ndarray, n_knots: int = 5,
) -> np.ndarray:
    """Ajusta el spline por minimos cuadrados y devuelve los residuos
    estandarizados (residuo / desviacion tipica de los residuos).
    """
    basis = cyclic_spline_basis(circular_scale, n_knots)
    coef, *_ = np.linalg.lstsq(basis, y, rcond=None)
    fitted = basis @ coef
    residuals = y - fitted
    sd = residuals.std(ddof=1) if len(residuals) > 1 else 1.0
    return residuals / sd if sd > 0 else residuals
