"""
Interfaz compartida para todos los modelos de ajuste ritmico del pipeline CIRCUST.

Nota de diseno
--------------
``FitResult`` es un contenedor de datos puro: almacena todo lo que ``fit()``
produce, incluido el tiempo de pico precalculado por el modelo. Asi el
resultado es autosuficiente — nadie necesita llamar al modelo despues de
obtenerlo.

La logica especifica de cada modelo vive exclusivamente en ``RhythmModel`` y
sus subclases:

    ``fit()`` — abstracto, produce un FitResult con peak_time ya calculado.

Para el calculo del pico con floats sueltos (reparametrizacion manual) cada
modelo expone un metodo estatico auxiliar, p.ej. ``FMMModel.peak_time()``.
"""
from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass

import numpy as np


# ---------------------------------------------------------------------------
# Contenedor de resultados — datos puros, sin logica de modelo
# ---------------------------------------------------------------------------

@dataclass
class FitResult:
    """
    Resultado de cualquier modelo de ajuste ritmico.

    Attributes
    ----------
    fitted : np.ndarray, dim (n_muestras,)
        Valores predichos por el modelo en cada punto temporal.

    params : dict
        Parametros estimados del modelo.

        Claves Cosinor : ``M``, ``A``, ``phi``
            - M   : mesor (nivel medio)
            - A   : amplitud
            - phi : acrofase en [0, 2*pi)

        Claves FMM : ``M``, ``A``, ``alpha``, ``beta``, ``omega``
            - M     : mesor
            - A     : amplitud
            - alpha : parametro de localizacion en [0, 2*pi)
            - beta  : parametro de forma en [0, 2*pi)
            - omega : parametro de asimetria en (0, 1]

    peak_time : float
        Tiempo de pico en [0, 2*pi), precalculado por el modelo al hacer fit().
        Cosinor: phi mod 2*pi  (phi es directamente el tiempo de pico).
        FMM: formula compUU — alpha + 2*atan2(...) mod 2*pi.

    r2 : float
        Coeficiente de determinacion  R-cuadrado = 1 - SS_res / SS_tot.

    residuals : np.ndarray, forma (n_muestras,)
        Residuos brutos: datos - ajustado.

    residuals_std : np.ndarray, forma (n_muestras,)
        Residuos estandarizados: (residuos - media) / desv_std.

    sse : float
        Suma de errores al cuadrado.

    model_name : str
        ``"cosinor"`` o ``"fmm"``.
    """

    fitted:        np.ndarray
    params:        dict
    peak_time:     float
    r2:            float
    residuals:     np.ndarray
    residuals_std: np.ndarray
    sse:           float
    model_name:    str

    def summary(self) -> str:
        lines = [
            f"=== Ajuste {self.model_name.upper()} ===",
            f"  R-cuadrado : {self.r2:.4f}",
            f"  peak_time  : {self.peak_time:.4f} rad",
            f"  SSE        : {self.sse:.6f}",
            f"  parametros : { {k: round(float(v), 5) for k, v in self.params.items()} }",
        ]
        return "\n".join(lines)


# ---------------------------------------------------------------------------
# Modelo base abstracto
# ---------------------------------------------------------------------------

class RhythmModel(ABC):
    """
    Interfaz comun para todos los modelos de ajuste de ritmos circadianos.

    Las subclases deben implementar :meth:`fit`.

    Los parametros pasados al constructor son hiperparametros que permanecen
    constantes entre llamadas a ``fit()`` (p.ej. tamanos de rejilla para FMM).
    """

    @abstractmethod
    def fit(
        self,
        data:        np.ndarray,
        time_points: np.ndarray,
    ) -> FitResult:
        """
        Ajusta el modelo a ``data`` observados en ``time_points``.

        El ``FitResult`` devuelto incluye ``peak_time`` ya calculado.

        Parameters
        ----------
        data : np.ndarray, forma (n_muestras,)
            Valores de expresion normalizados (tipicamente en [-1, 1]).
        time_points : np.ndarray, forma (n_muestras,)
            Eje temporal circular en [0, 2*pi), producido por CPCA.

        Returns
        -------
        FitResult
        """

    # ------------------------------------------------------------------
    # Utilidades compartidas por todos los modelos
    # ------------------------------------------------------------------

    @staticmethod
    def _r2(data: np.ndarray, fitted: np.ndarray) -> float:
        """R-cuadrado = 1 - SS_res / SS_tot."""
        ss_res = np.sum((data - fitted) ** 2)
        ss_tot = np.sum((data - data.mean()) ** 2)
        if ss_tot == 0:
            return 0.0
        return float(1.0 - ss_res / ss_tot)

    @staticmethod
    def _standardise_residuals(residuals: np.ndarray, ddof: int = 1) -> np.ndarray:
        """
        (res - media) / desv_std  con ``ddof`` grados de libertad.

        ``ddof`` debe ser el numero de parametros libres del modelo que
        genero ``residuals`` (p.ej. 3 para Cosinor, 5 para FMM), no un valor
        fijo de 1: SSE/(n-1) solo es insesgado para la varianza residual
        cuando el unico parametro estimado es la media; si el modelo estimo
        p parametros, el estimador insesgado es SSE/(n-p).
        """
        std = residuals.std(ddof=ddof)
        if std == 0:
            return np.zeros_like(residuals)
        return (residuals - residuals.mean()) / std
