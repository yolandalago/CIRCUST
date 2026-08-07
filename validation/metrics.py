"""
validation/metrics.py
=====================
Metricas de validacion para reconstruccion de orden temporal circular.

Cada estadistico es una funcion independiente con firma uniforme::

    fn(predicted, true, period=24.0) -> float

asi puede usarse tanto con datos sinteticos como con datos reales
(p.ej. tiempos GTEx, baboons con ZT conocido, etc.).

Contenido
---------
Helpers:
    - circular_abs_diff(a, b, period)        distancia circular elemento a elemento
    - mean_resultant_length(theta, period)   R_bar — concentracion (diagnostico)
    - circular_mean(theta, period)           MLE de mu bajo von Mises

Alineacion:
    - align_phases_by_mad(pred, true, period, try_reflection, grid_step)

Estadisticos individuales (sobre vectores de fases):
    - median_absolute_error(pred, true, period)        MedAE
    - std_absolute_error(pred, true, period)           SDAE
    - pct_within(pred, true, threshold_h, period)      %< Xh
    - auc_error_cdf(pred, true, period)                AUC del CDF
    - circular_correlation(pred, true, period)         CCC Jammalamadaka-Sarma
    - circular_correlation_trimmed(pred, true, period) rTrim (Mahmood 2022)

Estadisticos individuales (sobre matrices ajustadas FMM):
    - rre(data, fit_predicted)                         RRE Eq.(8) Larriba 2019
    - cre(fit_real, fit_predicted)                     CRE Eq.(9) Larriba 2019

Atajo:
    - evaluate_all(pred, true, period, ...)            todas las metricas en dict
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


PERIOD_HOURS: float = 24.0


# ===========================================================================
# Helpers numericos
# ===========================================================================

def circular_abs_diff(
    a:      np.ndarray,
    b:      np.ndarray,
    period: float = PERIOD_HOURS,
) -> np.ndarray:
    """
    Distancia circular absoluta elemento a elemento, en ``[0, period/2]``.

    Para puntos en un circulo, la distancia entre p.ej. 23 h y 1 h es 2 h
    (no 22 h). Esta funcion implementa ese atajo:

        d_i = min( |a_i - b_i| mod P,   P - |a_i - b_i| mod P )
    """
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    d = np.abs(a - b) % period
    return np.minimum(d, period - d)


def mean_resultant_length(
    theta:  np.ndarray,
    period: float = PERIOD_HOURS,
) -> float:
    """
    Longitud del vector resultante :math:`\\bar R` de un conjunto de fases.

    Definicion::

        R_bar = sqrt( mean(cos)^2 + mean(sin)^2 )    en [0, 1]

    Interpretacion
    --------------
    - ``R_bar -> 1`` : distribucion **unimodal y concentrada** en torno a
      una direccion preferida.
    - ``R_bar -> 0`` : distribucion **uniforme** en el circulo (sin direccion
      preferida).

    Por que importa
    ---------------
    El CCC de Jammalamadaka depende de las medias circulares ``mu_1`` y
    ``mu_2``, que son inestables cuando ``R_bar -> 0``. Usar esta funcion
    como **diagnostico** para saber si el CCC del dataset es fiable:

    - ``R_bar > 0.3``  el CCC suele ser fiable
    - ``R_bar < 0.1``  el CCC se vuelve inestable y puede dar valores
                       negativos aunque la prediccion sea casi perfecta
    """
    theta = np.asarray(theta, dtype=np.float64)
    rad   = 2.0 * np.pi * theta / period
    c     = float(np.mean(np.cos(rad)))
    s     = float(np.mean(np.sin(rad)))
    return float(np.sqrt(c * c + s * s))


def circular_mean(
    theta:  np.ndarray,
    period: float = PERIOD_HOURS,
) -> float:
    """
    Direccion media circular: MLE de :math:`\\mu` bajo von Mises.

    Formula::

        mu_hat = arctan2( sum sin(theta),  sum cos(theta) )

    Devuelve un escalar en ``[0, period)``. Es numericamente inestable
    cuando :func:`mean_resultant_length` esta cerca de 0.
    """
    theta  = np.asarray(theta, dtype=np.float64)
    rad    = 2.0 * np.pi * theta / period
    s      = float(np.sin(rad).sum())
    c      = float(np.cos(rad).sum())
    mu_rad = np.arctan2(s, c) % (2.0 * np.pi)
    return float(mu_rad * period / (2.0 * np.pi))


# ===========================================================================
# Alineacion: rotacion + reflexion optima por MAD-minimo
# ===========================================================================

@dataclass
class AlignmentResult:
    """
    Resultado de :func:`align_phases_by_mad`.

    Attributes
    ----------
    delta : float
        Desplazamiento optimo en horas (sumado a la prediccion).
    sign : int
        +1 si no se aplico reflexion, -1 si se invirtio el recorrido.
    aligned_pred : np.ndarray
        Prediccion ya alineada, en ``[0, period)``.
    mad : float
        MAD del error tras la alineacion (criterio de optimizacion).
    """
    delta:        float
    sign:         int
    aligned_pred: np.ndarray
    mad:          float


def align_phases_by_mad(
    predicted:      np.ndarray,
    true:           np.ndarray,
    period:         float = PERIOD_HOURS,
    try_reflection: bool  = True,
    grid_step:      float = 0.05,
) -> AlignmentResult:
    """
    Encuentra la rotacion (y opcionalmente reflexion) que minimiza la
    **mediana del error absoluto circular** entre prediccion y verdad.

    Necesario porque algoritmos *unsupervised* (CIRCUST, CYCLOPS, CHIRAL,
    DCPR) producen fases definidas salvo rotacion y/o reflexion del circulo.
    Sin esta alineacion, una prediccion perfecta podria reportar un error
    arbitrariamente grande solo por una eleccion distinta del "origen".

    Parameters
    ----------
    predicted, true
        Fases en horas, asumidas en ``[0, period)``.
    period
        Periodo del ritmo, en horas. Por defecto 24.
    try_reflection
        Si ``True``, prueba tambien ``-predicted`` (recorrido en sentido
        opuesto). Util porque CIRCUST puede invertir la direccion segun
        la asuncion biologica de ARNTL/DBP.
    grid_step
        Paso del grid de busqueda en horas. ``0.05`` da resolucion fina
        sin penalizar el tiempo de calculo.

    Returns
    -------
    AlignmentResult
    """
    predicted = np.asarray(predicted, dtype=np.float64) % period
    true      = np.asarray(true,      dtype=np.float64) % period

    grid  = np.arange(0.0, period, grid_step)
    signs = [+1, -1] if try_reflection else [+1]

    best_mad   = np.inf
    best_delta = 0.0
    best_sign  = +1

    for sign in signs:
        candidate = (sign * predicted) % period
        for delta in grid:
            shifted = (candidate + delta) % period
            mad     = float(np.median(circular_abs_diff(shifted, true, period)))
            if mad < best_mad:
                best_mad   = mad
                best_delta = float(delta)
                best_sign  = sign

    aligned = (best_sign * predicted + best_delta) % period
    return AlignmentResult(
        delta        = best_delta,
        sign         = best_sign,
        aligned_pred = aligned,
        mad          = best_mad,
    )


# ===========================================================================
# Estadisticos individuales (una funcion por metrica)
# ===========================================================================

def median_absolute_error(
    predicted: np.ndarray,
    true:      np.ndarray,
    period:    float = PERIOD_HOURS,
) -> float:
    """
    Mediana del error absoluto circular (MedAE), en horas.
    """
    return float(np.median(circular_abs_diff(predicted, true, period)))


def std_absolute_error(
    predicted: np.ndarray,
    true:      np.ndarray,
    period:    float = PERIOD_HOURS,
) -> float:
    """Desviacion estandar del error absoluto circular (SDAE), en horas."""
    return float(np.std(circular_abs_diff(predicted, true, period), ddof=0))


def pct_within(
    predicted:   np.ndarray,
    true:        np.ndarray,
    threshold_h: float,
    period:      float = PERIOD_HOURS,
) -> float:
    """
    Porcentaje de muestras con error absoluto **estrictamente menor** que
    ``threshold_h`` (en horas). Devuelve un valor en ``[0, 100]``.

    Tipicos usados en cronobiologia: ``threshold_h = 1`` y ``threshold_h = 2``.
    """
    errors = circular_abs_diff(predicted, true, period)
    if errors.size == 0:
        return float("nan")
    return 100.0 * float(np.mean(errors < threshold_h))


def auc_error_cdf(
    predicted: np.ndarray,
    true:      np.ndarray,
    period:    float = PERIOD_HOURS,
) -> float:
    """
    Area bajo el CDF empirico del error absoluto, normalizada a ``[0, 1]``.

    Interpretacion:
        - ``AUC = 1.0``: predictor perfecto (todos los errores = 0).
        - ``AUC ~ 0.5``: predictor aleatorio (errores distribuidos
          uniformemente en ``[0, period/2]``).
    """
    errors = circular_abs_diff(predicted, true, period)
    n      = errors.size
    if n == 0:
        return float("nan")

    max_err  = period / 2.0
    sorted_e = np.sort(errors)
    cdf_x    = np.concatenate([[0.0], sorted_e, [max_err]])
    cdf_y    = np.concatenate([[0.0], np.arange(1, n + 1) / n, [1.0]])
    return float(np.trapezoid(cdf_y, cdf_x) / max_err)


def circular_correlation(
    predicted: np.ndarray,
    true:      np.ndarray,
    period:    float = PERIOD_HOURS,
) -> float:
    """
    Coeficiente de correlacion circular **rho** de Jammalamadaka & SenGupta
    (referencia [44] del paper CIRCUST de Larriba et al. 2023):

    .. math::

        \\rho =
        \\frac{\\sum_i \\sin(\\theta_{1i} - \\mu_1)\\, \\sin(\\theta_{2i} - \\mu_2)}
             {\\sqrt{\\sum_i \\sin^2(\\theta_{1i} - \\mu_1)\\,
                     \\sum_i \\sin^2(\\theta_{2i} - \\mu_2)}}

    con :math:`\\mu_1, \\mu_2` estimados por **maxima verosimilitud bajo von Mises**,
    que coincide con la media circular muestral
    :math:`\\hat\\mu = \\arctan2(\\sum \\sin\\theta, \\sum \\cos\\theta)`.

    Interpretacion:
        - ``rho = +1``: prediccion perfecta tras un shift constante.
        - ``rho =  0``: no hay correlacion lineal-circular.
        - ``rho = -1``: prediccion anticorrelacionada (sentido invertido).

    Limitacion (importante)
    -----------------------
    Esta metrica depende de las medias circulares. Cuando los datos cubren
    uniformemente el circulo (``mean_resultant_length -> 0``), las medias
    saltan a posiciones arbitrarias y rho se vuelve numericamente
    inestable, dando valores negativos incluso para predicciones casi
    perfectas. Para diagnosticar la fiabilidad:

        R_bar = mean_resultant_length(true, period)

    - ``R_bar > 0.3``: rho generalmente fiable.
    - ``R_bar < 0.1``: rho inestable; prefiere MedAE / AUC / %within Xh.
    """
    predicted = np.asarray(predicted, dtype=np.float64)
    true      = np.asarray(true,      dtype=np.float64)

    p_rad = 2.0 * np.pi * predicted / period
    t_rad = 2.0 * np.pi * true      / period

    mu_p = np.arctan2(np.sin(p_rad).sum(), np.cos(p_rad).sum())
    mu_t = np.arctan2(np.sin(t_rad).sum(), np.cos(t_rad).sum())

    sp = np.sin(p_rad - mu_p)
    st = np.sin(t_rad - mu_t)

    num = float((sp * st).sum())
    den = float(np.sqrt((sp ** 2).sum() * (st ** 2).sum()))
    return num / den if den > 0.0 else 0.0


def _circular_trimmed_mean_with_mask(
    theta_rad:     np.ndarray,
    trim_fraction: float,
) -> tuple[float, np.ndarray]:
    """
    Media circular *trimmed*: la media circular estandar tras eliminar la
    proporcion ``trim_fraction`` de observaciones cuya **distancia angular**
    a la media circular inicial sea mayor (el equivalente circular del
    trimmed mean lineal).

    En datos circulares no existen "valores maximos" o "minimos" porque
    ``0 = 2pi``, asi que Mahmood (2022) propone podar las observaciones mas
    alejadas del grueso de los datos, medido por distancia angular a la
    media circular inicial.

    Devuelve la media trimmed (en radianes, en ``(-pi, pi]``) y la mascara
    booleana de observaciones conservadas, util para sincronizar la poda
    entre dos variables emparejadas.
    """
    n = theta_rad.size
    p = int(np.floor(trim_fraction * n))

    mu_init = float(
        np.arctan2(np.sin(theta_rad).sum(), np.cos(theta_rad).sum())
    )

    if p <= 0:
        return mu_init, np.ones(n, dtype=bool)

    diff = theta_rad - mu_init
    diff = np.arctan2(np.sin(diff), np.cos(diff))   # wrap a (-pi, pi]
    dist = np.abs(diff)

    keep_idx        = np.argsort(dist)[: n - p]
    kept            = np.zeros(n, dtype=bool)
    kept[keep_idx]  = True

    theta_kept = theta_rad[kept]
    mu_trim    = float(
        np.arctan2(np.sin(theta_kept).sum(), np.cos(theta_kept).sum())
    )
    return mu_trim, kept


def circular_correlation_trimmed(
    predicted:     np.ndarray,
    true:          np.ndarray,
    period:        float = PERIOD_HOURS,
    trim_fraction: float = 0.10,
) -> float:
    """
    Coeficiente de correlacion circular-circular **robusto** basado en la
    **media circular trimmed** (rTrim, Mahmood 2022).

    Adaptacion robusta del CCC de Jammalamadaka & SenGupta donde las medias
    circulares ``mu_1, mu_2`` se sustituyen por **medias trimmed**: en cada
    variable se descartan las observaciones cuya distancia angular a la
    media inicial sea mayor, y la suma se calcula sobre los pares en los
    que **ambas** componentes sobreviven a la poda.

    .. math::

        r_{Trim} =
        \\frac{\\sum_{i \\in K} \\sin(\\theta_i - \\bar\\theta_{Trim})
                              \\, \\sin(\\varphi_i - \\bar\\varphi_{Trim})}
             {\\sqrt{\\sum_{i \\in K} \\sin^2(\\theta_i - \\bar\\theta_{Trim})
                     \\, \\sum_{i \\in K} \\sin^2(\\varphi_i - \\bar\\varphi_{Trim})}}

    donde ``K`` es la interseccion de las mascaras de inliers de cada variable
    (un par se conserva solo si ambas componentes lo son).

    Parameters
    ----------
    predicted, true
        Fases en horas, en ``[0, period)``.
    period
        Periodo del ritmo, en horas (24 por defecto).
    trim_fraction
        Proporcion :math:`\\delta \\in [0, 0.5]` de observaciones a podar en
        cada variable. ``0.10`` es un valor tipico (Mahmood 2022 reporta
        resultados estables hasta 0.20-0.25).

    Notas
    -----
    - El paper original (Eq. 12) tiene un aparente erratum: el limite de la
      segunda suma del denominador aparece como ``n`` pero por consistencia
      dimensional debe ser ``n - p``. Aqui se usa el mismo conjunto ``K``
      en numerador y denominador (lo coherente con la idea de "trimmed").
    - El criterio de poda es por distancia angular a la media circular
      inicial de cada variable (Mahmood 2022, sec. 4.2). Se descarta el par
      cuando alguna de las dos componentes es outlier.
    - Para ``trim_fraction = 0`` se recupera el CCC clasico de
      :func:`circular_correlation`.
    """
    if not 0.0 <= trim_fraction < 0.5:
        raise ValueError("trim_fraction debe estar en [0, 0.5)")

    predicted = np.asarray(predicted, dtype=np.float64)
    true      = np.asarray(true,      dtype=np.float64)

    p_rad = 2.0 * np.pi * predicted / period
    t_rad = 2.0 * np.pi * true      / period

    mu_p, mask_p = _circular_trimmed_mean_with_mask(p_rad, trim_fraction)
    mu_t, mask_t = _circular_trimmed_mean_with_mask(t_rad, trim_fraction)

    kept = mask_p & mask_t
    if not np.any(kept):
        return 0.0

    sp = np.sin(p_rad[kept] - mu_p)
    st = np.sin(t_rad[kept] - mu_t)

    num = float((sp * st).sum())
    den = float(np.sqrt((sp ** 2).sum() * (st ** 2).sum()))
    return num / den if den > 0.0 else 0.0


def circular_regression_rho(
    predicted:      np.ndarray,
    true:           np.ndarray,
    period:         float = PERIOD_HOURS,
    try_reflection: bool  = True,
) -> float:
    """
    Bondad de ajuste **rho** de una **regresion circular-circular**, la
    metrica que usa el paper CIRCUST (Larriba et al. 2023) para validar el
    babuino (main, p.12; Supp Tabla A y Fig C). Alli rho se define como "un
    analogo de la suma de cuadrados residual de una regresion lineal", calculado
    con el paquete R ``circular`` (``lm.circular(type="c-c")``, Fisher 1993 [43],
    Adzhar 2011 [44]). Esta es una reimplementacion **en Python sin dependencia
    de R**.

    Modelo
    ------
    Se regresa la fase estimada ``y`` (respuesta) sobre el tiempo real ``x``
    (covariable), ambas angulares, con el modelo de **rotacion + reflexion**
    que es el unico grado de libertad que deja CPCA (el orden circular queda
    definido salvo origen y sentido):

        y_hat_i = mu_0 + s * x_i   (mod 2*pi),   s in {+1, -1}

    Para cada ``s``, el intercepto optimo ``mu_0`` es la **media circular** de
    los residuos ``y_i - s*x_i`` (minimiza la desviacion circular residual). Se
    elige el ``s`` con menor desviacion residual.

    Definicion de rho (coeficiente de determinacion circular)
    --------------------------------------------------------
    Con la desviacion circular ``D(a, b) = sum_i (1 - cos(a_i - b_i))`` como
    analogo circular de la suma de cuadrados:

        rho = 1 - D_res / D_tot
            = 1 - sum_i (1 - cos(y_i - y_hat_i)) / sum_i (1 - cos(y_i - ybar))

    donde ``ybar`` es la media circular de ``y``. Equivalentemente, con
    ``Rbar`` la longitud resultante media (:func:`mean_resultant_length`):

        rho = 1 - (1 - Rbar_residuos) / (1 - Rbar_y)

    Interpretacion y robustez
    -------------------------
    - ``rho -> 1``: los puntos caen sobre la diagonal de la regresion (orden
      recuperado). ``rho ~ 0``: no mejor que predecir la fase media constante.
    - **Es robusta al muestreo equiespaciado** (a diferencia del CCC de
      Jammalamadaka, :func:`circular_correlation`): ``D_res`` solo depende de
      los **residuos** del modelo ajustado, no de la concentracion marginal de
      ``x`` ni de ``y``. Por eso el paper la usa justamente para el babuino,
      cuyo ZT es uniforme. En ese regimen (``y`` tambien ~uniforme,
      ``Rbar_y ~ 0``) se reduce a ``rho ~ Rbar_residuos``, que es lo que
      ``lm.circular`` reporta como bondad de ajuste, asi que los valores son
      comparables a su Tabla A.

    Parameters
    ----------
    predicted, true
        Fases en horas. ``predicted`` = fase estimada (respuesta ``y``),
        ``true`` = tiempo real / ZT (covariable ``x``).
    period
        Periodo del ritmo, en horas (24 por defecto).
    try_reflection
        Si ``True`` (defecto), prueba ambos sentidos de giro y se queda con el
        mejor ajuste. CIRCUST puede invertir la direccion segun ARNTL/DBP.

    Returns
    -------
    float
        ``rho`` en ``(-inf, 1]`` (tipicamente ``[0, 1]``; valores negativos
        indican un orden peor que la fase media constante).
    """
    predicted = np.asarray(predicted, dtype=np.float64)
    true      = np.asarray(true,      dtype=np.float64)

    two_pi = 2.0 * np.pi
    y = (two_pi * predicted / period) % two_pi   # respuesta
    x = (two_pi * true      / period) % two_pi   # covariable

    # Desviacion total: respuesta alrededor de su media circular
    ybar  = np.arctan2(np.sin(y).sum(), np.cos(y).sum())
    d_tot = float((1.0 - np.cos(y - ybar)).sum())
    if d_tot <= 0.0:
        return float("nan")

    signs = (1.0, -1.0) if try_reflection else (1.0,)
    d_res = np.inf
    for s in signs:
        resid = y - s * x
        # intercepto optimo = media circular de los residuos
        mu0   = np.arctan2(np.sin(resid).sum(), np.cos(resid).sum())
        d_s   = float((1.0 - np.cos(resid - mu0)).sum())
        if d_s < d_res:
            d_res = d_s

    return 1.0 - d_res / d_tot


# ===========================================================================
# RRE / CRE — metricas a nivel de matrices ajustadas FMM (Larriba 2019)
# ===========================================================================
#
# Estas dos metricas NO comparan vectores de fases sino **matrices completas
# de expresion ajustada**. Validan el orden de muestras de forma indirecta:
# si el orden predicho es bueno, los ajustes FMM de los genes top que se
# obtienen con ese orden reconstruyen bien (a) los datos crudos (RRE) o
# (b) los ajustes FMM bajo el orden real (CRE).
#
# Por eso la firma es distinta: requieren matrices ``(G, n)`` de genes top,
# no los vectores ``(n,)`` de fases por sample. El re-fit del FMM bajo cada
# orden lo hace el orquestador antes de llamar a estas funciones.

def _sre_relative(
    reference: np.ndarray,
    fitted:    np.ndarray,
    eps:       float = 1e-12,
) -> np.ndarray:
    """
    Suma de errores cuadraticos *relativos*, por gen.

    Para cada gen ``k``::

        SRE^k = sum_i ((reference[k,i] - fitted[k,i]) / reference[k,i])^2

    Devuelve un array ``(G,)`` con un valor por gen.

    Notas
    -----
    La definicion del paper (Larriba 2019, Eqs. 8-9) usa errores normalizados
    por el valor de referencia, lo cual asume ``reference[k,i] != 0``. En
    expresion genica log-transformada (lo habitual) esto se cumple, pero
    se anade un ``eps`` minimo por seguridad numerica.
    """
    denom = np.where(np.abs(reference) < eps, eps, reference)
    return np.sum(((reference - fitted) / denom) ** 2, axis=1)


def rre(
    data:          np.ndarray,
    fit_predicted: np.ndarray,
) -> float:
    """
    **Relative Reconstruction Error** (Larriba 2019, Eq. 8).

    .. math::

        RRE_T(o) = \\frac{\\sum_k \\sum_i \\left(\\frac{Y_{k,i} - \\hat Y^k_{o,i}}{Y_{k,i}}\\right)^2}
                        {\\sum_k \\sum_i \\left(\\frac{Y_{k,i} - \\overline{Y_k}}{Y_{k,i}}\\right)^2}

    Mide cuanto se aleja el ajuste FMM **bajo el orden predicho** de los
    datos crudos, normalizado por la dispersion total. Es esencialmente
    1 - R^2 agregado sobre todos los genes top, en escala relativa.

    Interpretacion
    --------------
    - ``RRE -> 0``: el FMM bajo orden predicho reconstruye casi perfecto los
      datos crudos. Senal: el orden es consistente con la estructura ritmica
      de los genes top.
    - ``RRE -> 1``: el ajuste no es mejor que una constante (el orden
      predicho destruye la estructura ritmica).
    - ``RRE > 1``: el ajuste es **peor** que una constante (orden patologico).

    Util incluso **sin ground truth**: solo requiere los datos crudos y el
    fit FMM bajo el orden predicho. Aplicable a datos reales (GTEx, baboons,
    skeletal muscle).

    Parameters
    ----------
    data
        Matriz de expresion observada de los genes top, forma ``(G, n)``.
    fit_predicted
        Matriz de valores ajustados FMM bajo el orden predicho, ``(G, n)``,
        evaluada en las mismas posiciones-muestra que ``data``.

    Returns
    -------
    float
        RRE agregado sobre los ``G`` genes.

    Referencias
    -----------
    Larriba Y, Rueda C, Fernandez MA, Peddada SD (2019). Order Restricted
    Inference in Chronobiology. *Statistics in Medicine*, Eq. (8). La
    sustitucion IR -> FMM es indicacion de la tutora.
    """
    data          = np.asarray(data,          dtype=np.float64)
    fit_predicted = np.asarray(fit_predicted, dtype=np.float64)
    if data.shape != fit_predicted.shape:
        raise ValueError(
            f"data{data.shape} y fit_predicted{fit_predicted.shape} "
            f"deben tener la misma forma (G, n)."
        )

    means = data.mean(axis=1, keepdims=True)
    sre_residual = _sre_relative(data, fit_predicted).sum()
    sre_total    = _sre_relative(data, np.broadcast_to(means, data.shape)).sum()
    return float(sre_residual / sre_total) if sre_total > 0.0 else float("nan")


def cre(
    fit_real:      np.ndarray,
    fit_predicted: np.ndarray,
) -> float:
    """
    **Concordance Relative Error** (Larriba 2019, Eq. 9).

    .. math::

        CRE_T(o, REAL) =
        \\frac{\\sum_k \\sum_i \\left(\\frac{\\hat Y^k_{REAL,i} - \\hat Y^k_{o,i}}{\\hat Y^k_{REAL,i}}\\right)^2}
             {\\sum_k \\sum_i \\left(\\frac{\\hat Y^k_{REAL,i} - \\overline{\\hat Y^k_{REAL}}}{\\hat Y^k_{REAL,i}}\\right)^2}

    Mide la concordancia entre **dos reconstrucciones FMM** del mismo gen:
    una usando el orden predicho de samples y otra usando el orden real.
    Solo aplicable cuando se conoce el ground truth (datos sinteticos,
    experimentos con tiempo de muerte conocido).

    Por que valida orden de muestras
    --------------------------------
    Las fases predichas no se evaluan en abstracto: si son correctas, los
    pares ``(phi_pred_i, y_{k,i})`` definen la misma curva FMM que los
    pares ``(phi_real_i, y_{k,i})`` -- solo cambia la etiqueta del eje t.
    Si las fases predichas estan barajadas, los dos FMM ajustan **a curvas
    distintas** sample-a-sample. CRE cuantifica esa discrepancia.

    Ventajas frente al CCC entre vectores de fases
    ----------------------------------------------
    - Invariante a rotacion global del orden (el FMM se ajusta igual).
    - Pondera los errores de fase por su **impacto en la reconstruccion**:
      errores en zonas donde los genes top tienen mucha pendiente cuentan
      mas que en zonas planas.
    - Agrega informacion de todos los genes top en un solo escalar.

    Interpretacion
    --------------
    - ``CRE -> 0``: orden predicho equivalente al real (predicccion perfecta).
    - ``CRE = 1``: no es mejor que una constante.
    - ``CRE > 1``: el orden predicho destruye mas senal de la que recupera.

    Parameters
    ----------
    fit_real
        Matriz ``(G, n)`` de FMM ajustado con las fases reales.
    fit_predicted
        Matriz ``(G, n)`` de FMM ajustado con las fases predichas,
        evaluada en las mismas posiciones-muestra.

    Returns
    -------
    float

    Referencias
    -----------
    Larriba Y, Rueda C, Fernandez MA, Peddada SD (2019). Order Restricted
    Inference in Chronobiology. *Statistics in Medicine*, Eq. (9). La
    sustitucion IR -> FMM es indicacion de la tutora.
    """
    fit_real      = np.asarray(fit_real,      dtype=np.float64)
    fit_predicted = np.asarray(fit_predicted, dtype=np.float64)
    if fit_real.shape != fit_predicted.shape:
        raise ValueError(
            f"fit_real{fit_real.shape} y fit_predicted{fit_predicted.shape} "
            f"deben tener la misma forma (G, n)."
        )

    means = fit_real.mean(axis=1, keepdims=True)
    sre_residual = _sre_relative(fit_real, fit_predicted).sum()
    sre_total    = _sre_relative(fit_real, np.broadcast_to(means, fit_real.shape)).sum()
    return float(sre_residual / sre_total) if sre_total > 0.0 else float("nan")


# ===========================================================================
# API de conveniencia
# ===========================================================================

def evaluate_all(
    predicted:    np.ndarray,
    true:         np.ndarray,
    period:       float          = PERIOD_HOURS,
    thresholds_h: tuple          = (1.0, 2.0),
    align:        bool           = True,
) -> dict:
    """
    Calcula todas las metricas en un solo dict. Util en scripts de validacion.

    Parameters
    ----------
    predicted, true
        Fases en horas.
    period
        Periodo del ritmo, en horas (24 por defecto).
    thresholds_h
        Umbrales (en horas) para ``pct_within``. Por defecto ``(1, 2)``.
    align
        Si ``True`` (recomendado para datos unsupervised), aplica
        :func:`align_phases_by_mad` antes de calcular metricas y anade
        ``delta``, ``sign`` y ``mad`` al dict.

    Returns
    -------
    dict con claves:
        - ``MedAE``, ``SDAE``, ``AUC``        (error, en horas)
        - ``rho``                             (concordancia, recomendada)
        - ``CCC``, ``rTrim``                  (concordancia, ver aviso)
        - ``pct_within_Xh`` para cada umbral
        - ``R_bar_true``, ``R_bar_pred``   (diagnostico de estabilidad)
        - ``n``                              (numero de muestras evaluadas)
        - ``delta``, ``sign``, ``mad``       (solo si ``align=True``)

    Aviso sobre la medida de concordancia
    -------------------------------------
    Usar ``rho`` (:func:`circular_regression_rho`) y no ``CCC``. El CCC de
    Jammalamadaka centra cada marginal en su direccion media, que queda
    **indeterminada** cuando las fases se reparten de forma uniforme sobre
    el ciclo. Con ``R_bar_true`` proxima a 0 su valor es arbitrario: se han
    observado valores en torno a ``-0,9`` sobre reconstrucciones con un
    error mediano de 30 minutos. ``rho`` solo depende de los residuos del
    ajuste de rotacion y reflexion, de modo que se mantiene valida en ese
    regimen, que es el de todo diseno de muestreo equiespaciado.
    """
    if align:
        alig = align_phases_by_mad(predicted, true, period=period)
        pred = alig.aligned_pred
    else:
        pred = np.asarray(predicted, dtype=np.float64) % period

    true_a = np.asarray(true, dtype=np.float64) % period

    result: dict = {
        "MedAE":      median_absolute_error(pred, true_a, period),
        "SDAE":       std_absolute_error(pred, true_a, period),
        "AUC":        auc_error_cdf(pred, true_a, period),
        # Medida de concordancia recomendada. A diferencia del CCC, rho no
        # se apoya en la direccion media de las marginales, de modo que
        # sigue siendo valida cuando las fases verdaderas son uniformes
        # sobre el ciclo (Rbar ~ 0), que es justo el caso de los disenos de
        # muestreo equiespaciado y del ZT del babuino. Ver la nota sobre
        # degeneracion del CCC mas abajo.
        "rho":        circular_regression_rho(pred, true_a, period),
        # CCC y rTrim se conservan por continuidad con analisis previos,
        # pero DEGENERAN cuando ``R_bar_true`` es proxima a 0: al anularse
        # la resultante media, la direccion media que centra el coeficiente
        # queda indeterminada y el valor resultante es arbitrario. En ese
        # regimen se han observado valores fuertemente negativos con
        # reconstrucciones casi perfectas. No usarlos como criterio sin
        # comprobar antes ``R_bar_true``.
        "CCC":        circular_correlation(pred, true_a, period),
        "rTrim":      circular_correlation_trimmed(pred, true_a, period),
        "R_bar_true": mean_resultant_length(true_a, period),
        "R_bar_pred": mean_resultant_length(pred,   period),
        "n":          int(true_a.size),
    }
    for th in thresholds_h:
        result[f"pct_within_{th:g}h"] = pct_within(pred, true_a, th, period)

    if align:
        result["delta"] = alig.delta
        result["sign"]  = alig.sign
        result["mad"]   = alig.mad

    return result
