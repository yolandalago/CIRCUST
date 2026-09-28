"""
circust/robust_order.py
========================
Etapa 4: Estimacion robusta del orden circular.

Pipeline en dos fases:

  **Fase A — Etapa comun** (no cambia entre metodos):
      1. Genera K subconjuntos aleatorios de genes TOP (tamano configurable,
         por defecto 2/3 del total), muestreados sobre el conjunto TOP
         completo: los genes TOP ya superaron los filtros de calidad de la
         etapa 3, forzados o no, asi que ninguno se excluye del muestreo.
      2. Para cada subconjunto k=1..K ejecuta:
         a) CPCA sobre la submatriz de genes seleccionada.
         b) Sincronizacion biologica (rotacion + orientacion).
         c) Re-ajuste FMM a cada gen del TOP para obtener
            la tabla de estadisticos de 12 columnas.
      3. Calcula la mediana R² (FMM) por repeticion sobre genes de evaluacion.

  **Fase B — Tarea 1: eleccion del orden robusto** (segun ``method``):
      - ``"best_k"``: elige UNA de las K repeticiones segun el criterio de
        seleccion indicado en ``selection_criterion``:
          * ``"r2_median"``   (default): maximiza la mediana de R² (FMM) sobre
            los genes de evaluacion.
          * ``"n_well_fitted"``: maximiza el numero de genes de evaluacion con
            R² >= ``selection_r2_threshold``.
          * ``"medoid"``: minimiza la distancia coseno media al resto de las K
            ordenaciones (orden mas representativo).
        Sencillo y rapido. En este modo se calcula ademas, como salida robusta
        extra del ORDEN, un consenso por media circular de las K ordenaciones
        (opcion A) y se contrasta frente al orden preliminar CPCA mediante
        correlacion circular.
      - ``"aggregate"``: agregacion circular de los K ordenes (Barragan et al.,
        2015) usando el modulo ``circular_aggregation``. Por defecto TSP3.

Posicion en el pipeline
-----------------------
    TopGeneSelector  ->  **RobustOrderEstimator**  ->  (Etapa 5: visualizacion)
"""
from __future__ import annotations

import os
import numpy as np
import pandas as pd
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass, field
from math import pi

from circust.circular_aggregation import (
    aggregate_circular_orders, circular_isotonic_regression,
)
from circust.cpca                 import CPCA
from circust.fitting.fmm          import FMMModel
from circust.fitting.fmm_stabilized import StabilizedFMMModel, prepare_dataset
from circust.synchronizer        import CircularSynchronizer, DEFAULT_CONSISTENCY_PAIRS
# top_genes solo depende de circust.fitting, asi que no hay ciclo.
from circust.top_genes           import circular_median


# ===========================================================================
# Helpers a nivel de modulo (serializables por pickle, reutilizables)
# ===========================================================================

# Tope de seguridad de remuestreos por repeticion (Fase A.2) cuando ninguna
# submuestra de genes produce un orden circular legitimo. Provisional, sin
# calibrar en el estudio de simulacion.
_MAX_LEGITIMACY_ATTEMPTS = 50

def _fit_gene(args):
    """Ajusta FMM para un unico gen en una repeticion.

    Funcion a nivel de modulo para serializacion por ``ProcessPoolExecutor``.
    El segundo elemento de ``args`` puede ser:
      - dict de kwargs para ``FMMModel`` (camino "base")
      - un ``StabilizedFMMModel`` ya preparado (camino "stabilized");
        en ese caso se reusa directamente sin construir nada.
    """
    (i, gene, vvv, esc_k, fmm_obj, arntl_peak_raw,
     orientation_changed, order_rot, n_samp) = args

    if isinstance(fmm_obj, StabilizedFMMModel):
        fmm_loc = fmm_obj                         # reusa modelo compartido
    else:
        fmm_loc = FMMModel(**fmm_obj)             # base, construye por gen

    two_pi = 2.0 * pi

    # ---- FMM ----
    fr = fmm_loc.fit(vvv, esc_k)
    M, A = fr.params["M"], fr.params["A"]
    alpha = fr.params["alpha"]
    beta  = fr.params["beta"]
    omega = fr.params["omega"]
    al = (alpha - arntl_peak_raw + pi) % two_pi
    if not orientation_changed:
        pars_fmm = (M, A, al, beta, omega)
    else:
        pars_fmm = (
            M, A,
            (two_pi - al) % two_pi,
            (two_pi - beta) % two_pi,
            omega,
        )
    pkU = FMMModel.peak_time(pars_fmm[2], pars_fmm[3], pars_fmm[4])
    # No es pkU + pi: eso supone que el valle cae a media vuelta del pico, lo
    # cual solo vale para una sinusoide. Con omega = 0,2 esa regla erraba el
    # minimo real en 7,6 h de mediana.
    pkL = FMMModel.trough_time(pars_fmm[2], pars_fmm[3], pars_fmm[4])
    peaks_fmm = (pkU, pkL, pkU / two_pi * 100, pkL / two_pi * 100)
    resid = vvv - fr.fitted
    sFMM   = float(np.sum(resid**2) / max(n_samp - 5, 1))
    mseFMM = float(np.sum(resid**2) / n_samp)
    r2FMM  = float(fr.r2)
    if not orientation_changed:
        fitted_fmm_reord = fr.fitted[order_rot]
    else:
        fitted_fmm_reord = fr.fitted[order_rot][::-1]
    stat_fmm = list(pars_fmm) + list(peaks_fmm) + [sFMM, mseFMM, r2FMM]

    return (i, gene, stat_fmm, fitted_fmm_reord)


def _build_theta(
    sample_orders: np.ndarray, circ_scales: np.ndarray, n_samp: int,
) -> np.ndarray:
    """Matriz ``Theta`` (K, n_samp) con la fase asignada a cada muestra por rep.

    ``Theta[k, i]`` es la fase de la muestra ``i`` (indice original, el mismo
    de las columnas en el orden preliminar) en la repeticion ``k``.
    ``sample_orders[k]`` es la permutacion de indices y ``circ_scales[k]`` las
    fases en ese orden, de modo que ``Theta[k, sample_orders[k]] = circ_scales[k]``.
    """
    K = sample_orders.shape[0]
    Theta = np.zeros((K, n_samp), dtype=np.float64)
    for k in range(K):
        Theta[k, sample_orders[k]] = circ_scales[k]
    return Theta


def _circular_correlation(a: np.ndarray, b: np.ndarray) -> float:
    """Correlacion circular-circular de Jammalamadaka-Sarma entre dos vectores
    de angulos (radianes) emparejados muestra a muestra. Resultado en
    ``[-1, 1]``: 1 = misma estructura angular, -1 = invertida, ~0 = sin relacion.

    Es invariante a rotaciones globales (la eleccion de origen del circulo no
    afecta), justo lo deseable para comparar dos ordenaciones circulares. Es el
    analogo a la ``circular correlation`` del paper (Rao & Sengupta).
    """
    two_pi = 2.0 * np.pi
    a = np.asarray(a, dtype=np.float64) % two_pi
    b = np.asarray(b, dtype=np.float64) % two_pi
    a_bar = np.arctan2(np.sin(a).mean(), np.cos(a).mean())
    b_bar = np.arctan2(np.sin(b).mean(), np.cos(b).mean())
    sa = np.sin(a - a_bar)
    sb = np.sin(b - b_bar)
    den = np.sqrt(np.sum(sa ** 2) * np.sum(sb ** 2))
    if den == 0.0:
        return 0.0
    return float(np.sum(sa * sb) / den)


def _fisher_lee(a: np.ndarray, b: np.ndarray) -> float:
    """Correlacion circular de Fisher y Lee (1983) entre dos vectores de angulos.

        rho = sum_{i<j} sin(a_i-a_j) sin(b_i-b_j)
              / sqrt[ sum_{i<j} sin^2(a_i-a_j) * sum_{i<j} sin^2(b_i-b_j) ]

    Devuelve ``+1`` si las dos escalas recorren el circulo en el mismo
    sentido, ``-1`` si en sentidos opuestos y ``~0`` si no guardan relacion.

    Sustituye a ``_circular_correlation`` para **decidir orientacion**. La
    diferencia esta en de donde saca cada una su referencia: aquella resta a
    cada serie su propia direccion media, que con muestras repartidas por
    todo el ciclo queda mal determinada, mientras que esta trabaja con
    diferencias entre pares y no calcula ninguna media.

    Medido sobre 720 escenarios de inversion conocida, cubriendo n de 24 a 96
    y ruido de 0,05 a 0,8 rad: Jammalamadaka acierta el signo el 95,3 % de las
    veces y esta el 100 %. Sobre repeticiones genuinamente invertidas del
    pipeline, aquella llegaba a devolver valores POSITIVOS de +0,11 a +0,36,
    de modo que la condicion ``corr < 0`` no se cumplia y no se corregia
    ninguna repeticion.

    Es invariante a rotar cualquiera de las dos escalas, asi que el origen de
    cada repeticion no influye, y antisimetrica bajo reflexion, que es la
    propiedad que hace falta aqui.

    Se evalua por la forma computacional de Fisher (1993), que descompone la
    doble suma en sumas simples:

        A = sum cos(a) cos(b)    B = sum sin(a) sin(b)
        C = sum cos(a) sin(b)    D = sum sin(a) cos(b)
        E = sum cos(2a)          F = sum sin(2a)
        G = sum cos(2b)          H = sum sin(2b)

        rho = 4(AB - CD) / sqrt[ (n^2 - E^2 - F^2)(n^2 - G^2 - H^2) ]

    Es algebraicamente identica: el numerador sale de expandir el producto de
    senos y reagrupar la doble suma en productos de sumas, y el denominador
    de ``sin^2 x = (1 - cos 2x)/2``, que convierte
    ``sum_{i<j} sin^2(a_i - a_j)`` en ``(n^2 - E^2 - F^2)/4``. Pero cuesta
    **O(n) en vez de O(n^2)**: con n=480, el tamano del tejido real de GTEx,
    son 108 us contra 11,6 ms. Verificado contra la definicion literal sobre
    2000 pares al azar, con diferencia maxima 2e-16.
    """
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    n = a.size
    if n < 2 or b.size != n:
        return float("nan")
    ca, sa = np.cos(a), np.sin(a)
    cb, sb = np.cos(b), np.sin(b)
    A = float((ca * cb).sum())
    B = float((sa * sb).sum())
    C = float((ca * sb).sum())
    D = float((sa * cb).sum())
    E = float(np.cos(2.0 * a).sum())
    F = float(np.sin(2.0 * a).sum())
    G = float(np.cos(2.0 * b).sum())
    H = float(np.sin(2.0 * b).sum())
    den = np.sqrt((n * n - E * E - F * F) * (n * n - G * G - H * H))
    if den <= 0.0:
        return float("nan")
    return float(4.0 * (A * B - C * D) / den)


def _max_circular_gap(theta: np.ndarray) -> float:
    """Mayor hueco angular entre valores consecutivos de ``theta`` (radianes).

    Ordena los angulos y calcula las diferencias consecutivas, incluyendo el
    hueco de cierre del circulo (del ultimo valor de vuelta al primero, mas
    ``2*pi``). Es la base de las dos condiciones de legitimidad de un orden
    circular candidato que exige la Fase A.2 (ver ``_run_repetitions``):
    cobertura de mas de media circunferencia (hueco maximo < pi) y hueco
    maximo no mayor que el del orden preliminar de referencia.
    """
    two_pi = 2.0 * np.pi
    theta_sorted = np.sort(np.asarray(theta, dtype=np.float64) % two_pi)
    if theta_sorted.size < 2:
        return 0.0
    gaps = np.diff(theta_sorted)
    wrap = theta_sorted[0] + two_pi - theta_sorted[-1]
    return float(max(gaps.max(), wrap))


def _circular_mean_consensus(
    Theta: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Consenso por **media circular** de las K fases por muestra.

    1. Para cada muestra, media circular de las K fases:
       ``theta_bar_i = atan2(sum sin, sum cos)``.
    2. Longitud resultante media ``R_i in [0, 1]`` por muestra: concentracion de
       las K estimaciones (1 = la muestra cae siempre en la misma fase, ~0 =
       posicion inestable). Sirve como confianza por muestra.

    Promediar directamente es licito porque las K repeticiones ya han pasado por
    la etapa de sincronizacion, que les fija un origen y un sentido de recorrido
    comunes, y porque las que hubieran quedado invertidas ya se han corregido en
    :meth:`RobustOrderEstimator._fix_orientations`.

    Devuelve ``(theta_bar, resultant)``, ambos indexados por muestra.
    """
    two_pi = 2.0 * np.pi
    K = Theta.shape[0]
    C = np.cos(Theta).sum(axis=0)
    S = np.sin(Theta).sum(axis=0)
    theta_bar = np.arctan2(S, C) % two_pi
    resultant = np.sqrt(C ** 2 + S ** 2) / K
    return theta_bar, resultant


# ===========================================================================
# Resultado
# ===========================================================================

@dataclass
class RobustOrderResult:
    """Salida de :class:`RobustOrderEstimator`.

    Atributos -- orden final
    ------------------------
    sample_order : np.ndarray int, shape (n_muestras,)
        Orden de muestras final (best_k o agregacion). **Vacio con
        ``method="median"``**: CIRCUST original no define ningun orden de
        muestras consenso, calcula K y los descarta tras ajustar. Comprueba
        ``sample_order.size`` antes de consumirlo.
    circular_scale : np.ndarray, shape (n_muestras,)
        Escala circular del orden final. Vacio con ``method="median"``.
    method_used : str
        ``"median"``, ``"best_k"`` o ``"aggregate"``.
    aggregation_method : str | None
        Metodo de agregacion ganador (``"tsp3"`` / ``"hods"``) si
        method=aggregate, None si best_k.
    selected_rep : int | None
        Indice (0-based) de la repeticion seleccionada (solo en best_k).
    selection_criterion_used : str | None
        Criterio con el que se eligio la repeticion (solo en best_k).
    selection_score : float | None
        Valor del criterio en la rep seleccionada (solo en best_k): mediana R²,
        numero de genes bien ajustados o distancia coseno media, segun el caso.
    aggregation_msce : float | None
        MSCE (CIRE) del consenso (solo en aggregate).

    Atributos -- por repeticion (diagnostico)
    -----------------------------------------
    sample_orders : np.ndarray int, shape (K, n_muestras)
    circular_scales : np.ndarray, shape (K, n_muestras)
    direction_flipped : np.ndarray bool, shape (K,)
    stats_table : np.ndarray, shape (K * n_top, 12)
    per_rep_fmm_fit : list de pd.DataFrame
    median_r2_per_rep : np.ndarray, shape (K,)
    n_well_fitted_per_rep : np.ndarray, shape (K,)

    Atributos -- seleccion aleatoria
    --------------------------------
    selection_names : np.ndarray object, shape (K, tam)
    top_gene_names : list[str]
    """

    # ── Orden final ──────────────────────────────────────────────────────
    sample_order:        np.ndarray
    circular_scale:      np.ndarray
    method_used:         str
    aggregation_method:  str | None        = None
    selected_rep:        int | None        = None
    selection_criterion_used: str | None   = None
    selection_score:     float | None      = None
    aggregation_msce:    float | None      = None

    # ── Por repeticion ────────────────────────────────────────────────────
    sample_orders:     np.ndarray = field(default_factory=lambda: np.array([]))
    circular_scales:   np.ndarray = field(default_factory=lambda: np.array([]))
    direction_flipped: np.ndarray = field(default_factory=lambda: np.array([]))
    stats_table:       np.ndarray = field(default_factory=lambda: np.zeros((0, 12)))
    per_rep_fmm_fit:   list       = field(default_factory=list)
    median_r2_per_rep: np.ndarray = field(default_factory=lambda: np.array([]))
    n_well_fitted_per_rep: np.ndarray = field(default_factory=lambda: np.array([]))

    # ── Seleccion aleatoria ──────────────────────────────────────────────
    selection_names: np.ndarray = field(default_factory=lambda: np.empty((0, 0), dtype=object))
    top_gene_names:  list       = field(default_factory=list)
    # Atlas de best_k (Larriba 2023, Step 4): parametros FMM de la REPETICION
    # SELECCIONADA (sel_rep), consistentes con el orden que se devuelve.
    # Columnas: gene_name, M, amplitude, alpha, beta, omega, peak_phase_rad,
    # trough_phase_rad, r2_fmm. Solo se rellena con method="best_k" (en aggregate
    # el atlas se re-ajusta sobre el orden consenso, ver _build_atlas).
    selected_fmm:    "pd.DataFrame" = field(default_factory=pd.DataFrame)

    # ── Opcion A: consenso del ORDEN por media circular (solo best_k) ──────
    # Agregacion de las K ordenaciones por media circular de las fases por
    # muestra, alineadas (Procrustes circular) al orden preliminar. Es la
    # estimacion robusta del ORDEN de muestras por agregacion de las K
    # ordenaciones. Solo se rellena con method="best_k".
    #   circmean_order     : permutacion consenso (indices originales)
    #   circmean_scale     : fases del consenso, en el orden de circmean_order
    #   circmean_resultant : R_i in [0,1] por muestra (confianza), mismo orden
    circmean_order:     np.ndarray | None = None
    circmean_scale:     np.ndarray | None = None
    circmean_resultant: np.ndarray | None = None
    # Contraste frente al orden preliminar CPCA: correlacion circular de
    # Jammalamadaka. Se incluye tambien best_k vs preliminar para comparar si
    # la agregacion (opcion A) se acerca mas o menos al preliminar que la
    # repeticion unica elegida.
    circmean_vs_prelim_corr:  float | None = None
    bestk_vs_prelim_corr:     float | None = None

    def summary(self) -> str:
        K = self.sample_orders.shape[0] if self.sample_orders.ndim == 2 else 0
        n_top = len(self.top_gene_names)
        flipped = int(self.direction_flipped.sum()) if len(self.direction_flipped) > 0 else 0
        lines = [
            "=== Resumen de Orden Robusto ===",
            f"  Repeticiones (K)         : {K}",
            f"  Genes en TOP             : {n_top}",
            f"  Metodo de seleccion      : {self.method_used}",
        ]
        if self.method_used == "aggregate":
            lines.append(f"  Metodo de agregacion     : {self.aggregation_method}")
            if self.aggregation_msce is not None:
                lines.append(f"  MSCE consenso            : {self.aggregation_msce:.4f}")
        elif self.method_used == "best_k" and self.selected_rep is not None:
            lines.append(f"  Repeticion seleccionada  : {self.selected_rep + 1}")
            lines.append(f"  Criterio de seleccion    : {self.selection_criterion_used}")
            lines.append(f"  Puntuacion del criterio  : {self.selection_score:.4f}")
            if self.circmean_vs_prelim_corr is not None:
                lines.append(
                    f"  Consenso A vs prelim     : corr={self.circmean_vs_prelim_corr:.4f}"
                )
                lines.append(
                    f"  best_k   vs prelim       : corr={self.bestk_vs_prelim_corr:.4f}"
                )
        lines.append(f"  Repeticiones con flip    : {flipped}/{K}")
        if len(self.median_r2_per_rep) > 0:
            lines.append(
                f"  Mediana R² por rep       : "
                f"{np.round(self.median_r2_per_rep, 4).tolist()}"
            )
        return "\n".join(lines)


# ===========================================================================
# RobustOrderEstimator
# ===========================================================================

class RobustOrderEstimator:
    """Etapa 4 de CIRCUST: estimacion robusta del orden circular.

    Parameters
    ----------
    n_reps : int
        Numero K de repeticiones aleatorias. Por defecto 5.
    sample_size_fraction : float
        Fraccion de genes TOP por extraccion. Por defecto 2/3.
    method : str
        Tarea 1 — eleccion del orden final:
        - ``"best_k"``: seleccionar una de las K repeticiones (ver
          ``selection_criterion``).
        - ``"aggregate"``: agregar los K ordenes via ``circular_aggregation``.
    selection_criterion : str
        Solo si method=best_k. Criterio de calidad con el que se elige la
        repeticion:
        - ``"r2_median"`` (default): mediana de R² (FMM) sobre los genes de
          evaluacion, maxima.
        - ``"n_well_fitted"``: numero de genes de evaluacion con
          R² >= ``selection_r2_threshold``, maximo. Ante empate se desempata
          por mediana de R².
        - ``"medoid"``: distancia coseno media al resto de ordenaciones,
          minima (orden mas representativo del conjunto).
    selection_r2_threshold : float
        Umbral tau del criterio ``"n_well_fitted"``. Valor PROVISIONAL 0.65,
        pendiente de calibrar en el estudio de simulacion. No usar 0.5: es el
        mismo umbral con el que se filtra la ritmicidad al construir el
        conjunto TOP, de modo que sobre ese conjunto el criterio satura por
        construccion (97-100 % de genes lo superan en todas las reps).

        TODO(tau): cuando las simulaciones fijen el valor definitivo hay que
        actualizarlo aqui, en el default de ``--selection-r2-threshold`` de
        scripts/run_pipeline.py y en docs/main.tex (subseccion "Criterios de
        seleccion", donde hay un bloque de comentario PENDIENTE con el
        detalle de lo que falta escribir).
    aggregation_method : str
        Solo si method=aggregate. ``"auto"`` (default; corre TSP3 y HODs y se
        queda con el de menor MSCE), ``"tsp3"`` o ``"hods"``.
    consensus_scale_mode : str
        Solo si method=aggregate. Como se construye la escala circular asociada
        al orden consenso.

        - ``"circmean"`` (default): media circular por muestra de las K fases.
          Es el mejor estimador de la fase, pero no garantiza que
          ``argsort(circular_scale)`` reproduzca ``sample_order``.
        - ``"cire"``: proyeccion isotonica circular de esa media sobre el orden
          consenso. Verifica ``phi_{o_1} <= ... <= phi_{o_n}`` por construccion,
          como pide el anexo, pero agrupa muestras en bloques y en simulacion
          multiplica por 2 a 5 el error mediano de fase.
    anchor_gene, direction_gene : str
        Genes ancla para sincronizacion.
    consistency_pairs : list[tuple[str, str]] | None
        Pares de orden a priori para la verificacion de consistencia del
        Paso 2.2 de :class:`CircularSynchronizer`. Ver su docstring.
    fmm_length_alpha_grid, fmm_length_omega_grid, fmm_num_reps : int
        Hiperparametros FMM base.
    fmm_method : str
        ``"base"`` (FMMModel) o ``"stabilized"`` (StabilizedFMMModel).
    fmm_stab_kwargs : dict | None
        Kwargs adicionales para ``prepare_dataset`` cuando stab.
    lambda_star : float | None
        Lambda heredado de CPCA (modo stab) — evita recalibrar.
    fmm_refine : bool
        Refinamiento local opcional de ``(alpha, omega)`` del FMM estabilizado
        (Seccion S3.4.4 del suplemento de CircaFMM) en los ajustes de cada
        repeticion. Por defecto False. Se ignora con ``fmm_method="base"``.
    recalibrate_lambda_per_rep : bool
        Si True, recalibra lambda del FMM estabilizado en cada una de las K
        repeticiones. Si False (default), reusa el ``lambda_star`` heredado
        del CPCA.

        El default es False porque recalibrar cuesta mucho y no aporta nada,
        las dos cosas medidas:

        **No aporta.** La calibracion por Monte Carlo no determina lambda de
        forma estable. Sobre la MISMA escala, cambiando solo la semilla,
        devuelve 1,23, 1,23, 4,12, 4,12 y 8,76: un factor 7. Esa dispersion
        (193 %) es MAYOR que la que hay entre repeticiones distintas (150 %),
        de modo que lo que varia entre repeticiones es ruido del estimador y
        no informacion sobre los datos. Y el valor tampoco cambia el
        resultado donde el metodo funciona: con omega=0,5 y A/sigma=2, pasar
        de lambda 1,23 a 8,76 da el mismo orden bit a bit. Solo mueve algo
        con omega=0,1, donde el conjunto TOP ya ha colapsado y el resultado
        es incorrecto con los tres valores.

        **Cuesta.** Cada calibracion son 27 s y se hacian seis por celda, una
        en el CPCA y una por repeticion: el 99 % del sobrecoste del FMM
        estabilizado. Perfilado de la Etapa 4: 132,7 s de 135,7 s en
        prepare_dataset, con 600 llamadas a _simulate_spikes_vectorized.
    n_jobs : int
        Workers para ajuste paralelo de genes (-1 = todos).
    seed : int | None
        Semilla para reproducibilidad.
    verbose : bool
        Mensajes de progreso.
    """

    def __init__(
        self,
        n_reps:                int   = 5,
        sample_size_fraction:  float = 2.0 / 3.0,
        method:                str   = "aggregate",
        selection_criterion:   str   = "r2_median",
        selection_r2_threshold: float = 0.65,
        aggregation_method:    str   = "auto",
        consensus_scale_mode:  str   = "circmean",
        anchor_gene:           str   = "ARNTL",
        direction_gene:        str   = "DBP",
        consistency_pairs:     list[tuple[str, str]] | None = None,
        fmm_length_alpha_grid: int   = 48,
        fmm_length_omega_grid: int   = 24,
        fmm_num_reps:          int   = 3,
        fmm_method:            str   = "base",
        fmm_stab_kwargs:       dict|None  = None,
        lambda_star:           float|None = None,
        recalibrate_lambda_per_rep: bool   = False,
        fmm_refine:            bool  = False,
        n_jobs:                int   = 1,
        seed:                  int|None   = None,
        verbose:               bool  = True,
    ) -> None:
        if method not in ("median", "best_k", "aggregate"):
            raise ValueError(
                "method debe ser 'median', 'best_k' o 'aggregate', "
                f"recibido: {method!r}"
            )
        if selection_criterion not in ("r2_median", "n_well_fitted", "medoid"):
            raise ValueError(
                "selection_criterion debe ser 'r2_median', 'n_well_fitted' o "
                f"'medoid', recibido: {selection_criterion!r}"
            )
        if aggregation_method not in ("auto", "tsp3", "hods"):
            raise ValueError(f"aggregation_method debe ser 'auto', 'tsp3' o 'hods', recibido: {aggregation_method!r}")
        if consensus_scale_mode not in ("circmean", "cire"):
            raise ValueError(
                "consensus_scale_mode debe ser 'circmean' o 'cire', "
                f"recibido: {consensus_scale_mode!r}"
            )
        if fmm_method not in ("base", "stabilized"):
            raise ValueError(f"fmm_method debe ser 'base' o 'stabilized', recibido: {fmm_method!r}")

        self.n_reps               = n_reps
        self.sample_size_fraction = sample_size_fraction
        self.method               = method
        self.selection_criterion  = selection_criterion
        self.selection_r2_threshold = selection_r2_threshold
        self.aggregation_method   = aggregation_method
        self.consensus_scale_mode = consensus_scale_mode
        self.anchor_gene          = anchor_gene
        self.direction_gene       = direction_gene
        self.consistency_pairs    = consistency_pairs
        self._fmm_kwargs = dict(
            length_alpha_grid=fmm_length_alpha_grid,
            length_omega_grid=fmm_length_omega_grid,
            num_reps=fmm_num_reps,
        )
        self.fmm_method      = fmm_method
        self.fmm_stab_kwargs = dict(fmm_stab_kwargs or {})
        self.lambda_star_inherited      = lambda_star
        self.recalibrate_lambda_per_rep = recalibrate_lambda_per_rep
        self.fmm_refine = fmm_refine
        self.n_jobs  = n_jobs
        self.seed    = seed
        self.verbose = verbose

    # ------------------------------------------------------------------
    # API publica
    # ------------------------------------------------------------------

    def run(
        self,
        top_matrix:     pd.DataFrame,
        theta_pre:      np.ndarray,
        expr_full_norm: pd.DataFrame,
        core_genes:     list[str],
        eval_genes:     list[str] | None = None,
    ) -> RobustOrderResult:
        """Ejecuta la estimacion robusta completa.

        Parameters
        ----------
        top_matrix : pd.DataFrame (n_top, n_muestras)
            Matriz de los genes TOP, en el orden circular preliminar. Sus
            indices de fila son el conjunto TOP sobre el que se muestrea.
        theta_pre : np.ndarray (n_muestras,)
            Escala circular preliminar. Sirve de referencia en tres puntos:
            la condicion de legitimidad de A.2, la correccion de
            orientaciones de A.2b y el contraste del consenso de la opcion A.
        expr_full_norm : pd.DataFrame (n_genes, n_muestras)
            Matriz normalizada completa, en el orden circular preliminar. De
            aqui se extraen las filas de cada submuestra para su CPCA.
        core_genes : list[str]
            Genes reloj centrales, usados para sincronizar cada repeticion.
        eval_genes : list[str] | None
            Genes sobre los que se evaluan los criterios de seleccion.
            Default: el conjunto TOP completo.
        """
        self._log("=== Etapa 4: Estimacion de Orden Robusto ===")
        self._log(f"  Metodo: {self.method}  |  K={self.n_reps}  |  "
                  f"fraccion={self.sample_size_fraction:.2f}")
        if self.method == "aggregate":
            self._log(f"  Agregacion: {self.aggregation_method}")

        names_top  = list(top_matrix.index)
        n_top      = len(names_top)
        n_samp     = top_matrix.shape[1]

        # Conjunto de evaluacion de los criterios de seleccion. Por defecto el
        # conjunto TOP completo y no los genes core: con ~12 genes core la
        # mediana de R² es muy inestable y el conteo de genes bien ajustados
        # satura (12/12 en casi todas las reps), mientras que sobre el TOP
        # ambos criterios discriminan. Ninguno de los dos conjuntos es
        # independiente del procedimiento (el TOP se selecciona con el orden
        # preliminar, que a su vez sale de los core), de modo que la eleccion
        # se hace por potencia estadistica.
        if eval_genes is None:
            eval_genes = list(names_top)
        self._log(f"  Genes de evaluacion: {len(eval_genes)}")

        # ── Fase A.1: Seleccion aleatoria ───────────────────────────────
        rng, tam = self._random_selection(names_top)
        K = self.n_reps
        theta_pre = np.asarray(theta_pre, dtype=np.float64)

        # ── Fase A.2: CPCA + Sincronizacion + ajustes por rep ─────────
        (sample_orders, circ_scales, flipped_arr,
         stats_table, per_rep_fmm, sel_names) = \
            self._run_repetitions(
                rng, tam, top_matrix, expr_full_norm,
                core_genes, names_top, n_top, n_samp, K, theta_pre,
            )

        # ── Fase A.2b: correccion de orientaciones ────────────────────
        sample_orders, circ_scales, stats_table, per_rep_fmm = \
            self._fix_orientations(
                sample_orders, circ_scales, stats_table, per_rep_fmm,
                theta_pre = theta_pre,
                n_samp    = n_samp,
            )

        # ── Fase A.3: Mediana R² por rep ──────────────────────────────
        median_r2 = self._compute_median_r2(
            stats_table, names_top, eval_genes, K, n_top,
        )
        self._log(f"  Mediana R² por rep: {np.round(median_r2, 4).tolist()}")

        n_well_fitted = self._compute_n_well_fitted(
            stats_table, names_top, eval_genes, K, n_top,
        )
        self._log(
            f"  Genes con R² >= {self.selection_r2_threshold:g} por rep: "
            f"{n_well_fitted.astype(int).tolist()}"
        )

        # ── Fase B: Tarea 1 — eleccion del orden ──────────────────────
        if self.method == "median":
            # CIRCUST original (Larriba et al. 2023, paso 4). Agrega los
            # parametros FMM por gen y NO define ningun orden de muestras
            # consenso: los K ordenes se calculan, se usan para ajustar y se
            # descartan. Por eso sample_order y circular_scale salen vacios,
            # y es la razon de que este metodo no pueda entrar en las
            # metricas de orden (rho, MAD), solo en las de atlas.
            self._log("  Metodo 'median': mediana de los parametros FMM "
                      "sobre las K repeticiones; sin orden de muestras.")
            return RobustOrderResult(
                sample_order       = np.array([], dtype=int),
                circular_scale     = np.array([], dtype=np.float64),
                method_used        = "median",
                aggregation_method = None,
                selected_rep       = None,
                selection_criterion_used = None,
                selection_score    = None,
                sample_orders      = sample_orders,
                circular_scales    = circ_scales,
                direction_flipped  = flipped_arr,
                stats_table        = stats_table,
                per_rep_fmm_fit    = per_rep_fmm,
                median_r2_per_rep  = median_r2,
                n_well_fitted_per_rep = n_well_fitted,
                selection_names    = sel_names,
                top_gene_names     = names_top,
                selected_fmm       = self._median_fmm_atlas(
                    stats_table, names_top, K, n_top,
                ),
            )

        if self.method == "best_k":
            final_order, final_scale, sel_rep, sel_score = \
                self._select_best_k(
                    sample_orders, circ_scales, median_r2, n_well_fitted, n_samp,
                )
            # Atlas de best_k: parametros FMM de la rep elegida, consistentes con
            # su orden. No se usa la mediana sobre las K (promediaria ordenes ya
            # descartados por peores segun el R² mediano).
            selected_fmm = self._rep_fmm_atlas(
                stats_table, names_top, sel_rep, n_top,
            )
            # Opcion A (extra robusto del orden): consenso por media circular de
            # las K ordenaciones + contraste vs preliminar.
            circ_extra = self._circmean_consensus_and_contrast(
                sample_orders, circ_scales,
                best_order = final_order,
                best_scale = final_scale,
                theta_pre  = theta_pre,
                n_samp     = n_samp,
            )
            return RobustOrderResult(
                sample_order       = final_order,
                circular_scale     = final_scale,
                method_used        = "best_k",
                aggregation_method = None,
                selected_rep       = sel_rep,
                selection_criterion_used = self.selection_criterion,
                selection_score    = sel_score,
                sample_orders      = sample_orders,
                circular_scales    = circ_scales,
                direction_flipped  = flipped_arr,
                stats_table        = stats_table,
                per_rep_fmm_fit    = per_rep_fmm,
                median_r2_per_rep  = median_r2,
                n_well_fitted_per_rep = n_well_fitted,
                selection_names    = sel_names,
                top_gene_names     = names_top,
                selected_fmm       = selected_fmm,
                **circ_extra,
            )
        else:
            final_order, final_scale, agg_msce, agg_selected = \
                self._aggregate_orders_via_consensus(
                    sample_orders, circ_scales, n_samp,
                )
            return RobustOrderResult(
                sample_order       = final_order,
                circular_scale     = final_scale,
                method_used        = "aggregate",
                aggregation_method = agg_selected,
                selected_rep       = None,
                selection_criterion_used = None,
                selection_score    = None,
                aggregation_msce   = agg_msce,
                sample_orders      = sample_orders,
                circular_scales    = circ_scales,
                direction_flipped  = flipped_arr,
                stats_table        = stats_table,
                per_rep_fmm_fit    = per_rep_fmm,
                median_r2_per_rep  = median_r2,
                n_well_fitted_per_rep = n_well_fitted,
                selection_names    = sel_names,
                top_gene_names     = names_top,
            )

    # ==================================================================
    # Fase A.1 — Seleccion aleatoria controlada
    # ==================================================================

    def _random_selection(
        self,
        names_top: list[str],
    ) -> tuple[np.random.Generator, int]:
        """Prepara el generador aleatorio y el tamano de cada submuestra.

        Se muestrea sobre el conjunto TOP completo, sin excluir ningun gen:
        los genes TOP (forzados como semilla o no) ya superaron los filtros de
        calidad de la etapa 3, asi que se consideran igualmente validos para
        entrar en la submatriz de cada repeticion.

        El sorteo en si de cada subconjunto de tamano ``tam`` ocurre dentro de
        :meth:`_run_repetitions`, porque una submuestra puede tener que
        redibujarse si el orden circular que produce no resulta legitimo (ver
        el docstring de ese metodo).

        Returns
        -------
        (rng, tam)
        """
        self._log("  --- Seleccion aleatoria ---")

        rng = np.random.default_rng(self.seed)

        n_top = len(names_top)
        if n_top < 2:
            raise ValueError(
                f"El conjunto TOP tiene {n_top} genes; se necesitan al menos 2."
            )

        tam = max(1, min(n_top, int(np.ceil(self.sample_size_fraction * n_top))))

        self._log(f"    Genes TOP: {n_top}  |  tam={tam}  |  K={self.n_reps}")

        return rng, tam

    # ==================================================================
    # Fase A.2 — CPCA + Sincronizacion por repeticion
    # ==================================================================

    def _run_repetitions(
        self,
        rng:            np.random.Generator,
        tam:            int,
        top_matrix:     pd.DataFrame,
        expr_full_norm: pd.DataFrame,
        core_genes:     list[str],
        names_top:      list[str],
        n_top:          int,
        n_samp:         int,
        K:              int,
        theta_pre:      np.ndarray,
    ) -> tuple:
        """Por cada rep: sorteo de genes + CPCA(submatriz) + Synchronizer +
        ajuste FMM.

        Las muestras se preservan: el CPCA interno NO elimina outliers
        adicionales (``max_outlier_fraction=0``) — los outliers ya se quitaron
        en el CPCA principal del pipeline. Esto preserva ``n_samp`` constante
        en todas las reps.

        Legitimidad del orden circular candidato
        -----------------------------------------
        Tras sincronizar cada repeticion se obtiene ``esc_fin``, la escala
        circular de esa repeticion. Se exige que sea "legitima":

          (a) sus valores cubren mas de media circunferencia, es decir, el
              hueco angular maximo entre dos valores consecutivos es < pi;
          (b) ese hueco maximo no es mayor que el del orden preliminar de
              referencia (``theta_pre``, la escala de la Etapa 2).

        Si el subconjunto de genes sorteado para una repeticion produce una
        ``esc_fin`` que no cumple (a) y (b), se descarta y se vuelve a
        sortear un nuevo subconjunto de genes para esa misma repeticion,
        hasta ``_MAX_LEGITIMACY_ATTEMPTS`` intentos. Si se agota el limite sin
        exito, se registra un aviso y se sigue adelante con el ultimo intento
        en vez de interrumpir el pipeline completo por una unica repeticion
        problematica. Solo se aplica aqui, en la seleccion de las K
        ordenaciones de esta etapa, no al orden preliminar de la Etapa 2.
        """
        self._log("  --- Ejecutando repeticiones ---")

        # FMM "base" para core genes dentro de la sincronizacion previa
        # (en stab se reasigna por rep porque depende de la escala k)
        fmm = FMMModel(**self._fmm_kwargs) if self.fmm_method == "base" else None

        # Pares de consistencia a priori (Paso 2.2), resueltos una vez. Por
        # repeticion/intento se filtran a los que caigan dentro de ``genes_k``
        # (el subconjunto TOP sorteado para esa CPCA): un par cuyos dos genes
        # no estan en el muestreo de esa repeticion no aporta evidencia sobre
        # el orden que ese muestreo en concreto ha producido, y usarlo de
        # todas formas sesgaria por igual las K repeticiones en la misma
        # direccion, anulando la diversidad que el muestreo aleatorio busca.
        pairs_all = (
            list(self.consistency_pairs) if self.consistency_pairs is not None
            else list(DEFAULT_CONSISTENCY_PAIRS)
        )

        names_top_arr = np.asarray(names_top, dtype=object)
        theta_pre = np.asarray(theta_pre, dtype=np.float64)
        max_gap_pre = _max_circular_gap(theta_pre) if theta_pre.size == n_samp else None
        if max_gap_pre is None:
            self._log(
                "  [legitimidad] sin escala preliminar valida; solo se exige "
                "cobertura de mas de media circunferencia."
            )

        sample_orders = np.zeros((K, n_samp), dtype=int)
        circ_scales   = np.zeros((K, n_samp), dtype=np.float64)
        flipped_arr   = np.zeros(K, dtype=bool)
        stats_table   = np.zeros((K * n_top, 12), dtype=np.float64)
        per_rep_fmm   = []
        sel_names     = np.empty((K, tam), dtype=object)

        for k in range(K):
            self._log(f"    --- Repeticion {k+1}/{K} ---")

            attempt = 0
            while True:
                attempt += 1
                genes_k_arr = names_top_arr[rng.choice(n_top, size=tam, replace=False)]
                genes_k = [
                    g for g in genes_k_arr.tolist()
                    if g is not None and g in expr_full_norm.index
                ]
                if len(genes_k) < 2:
                    raise ValueError(
                        f"Repeticion {k}: menos de 2 genes disponibles para "
                        f"CPCA. genes_k = {genes_k}"
                    )

                genes_k_set = set(genes_k)
                pairs_k = [
                    (a, b) for (a, b) in pairs_all
                    if a in genes_k_set and b in genes_k_set
                ]
                prelim = CircularSynchronizer(
                    anchor_gene       = self.anchor_gene,
                    direction_gene    = self.direction_gene,
                    consistency_pairs = pairs_k,
                    verbose           = False,
                )

                # ── 1. CPCA sobre submatriz seleccionada ────────────────
                # max_outlier_fraction=0: NO elimina muestras adicionales
                # (los outliers ya se quitaron en el CPCA principal).
                cpca_stab_kwargs = dict(self.fmm_stab_kwargs)
                if (self.fmm_method == "stabilized"
                        and not self.recalibrate_lambda_per_rep
                        and self.lambda_star_inherited is not None
                        and "lambda_star" not in cpca_stab_kwargs):
                    cpca_stab_kwargs["lambda_star"] = float(self.lambda_star_inherited)

                sub_df = expr_full_norm.loc[genes_k]
                cpca_k = CPCA(
                    core_genes           = genes_k,
                    max_outlier_fraction = 0.0,
                    fmm_method           = self.fmm_method,
                    fmm_stab_kwargs      = cpca_stab_kwargs,
                    compute_fmm_fits     = False,  # solo se usan sample_order/circular_scale
                    verbose              = False,
                ).run(sub_df)
                order_k = cpca_k.sample_order
                esc_k   = cpca_k.circular_scale

                # ── 2. Reordenar la matriz TOP ──────────────────────────
                # ``iloc`` preserva los nombres originales de muestra; no
                # reseteamos columnas para que la salida final emparejada
                # por nombre sea posible aguas arriba.
                top_k = top_matrix.iloc[:, order_k].copy()

                # ── 3. Construir modelo FMM y ajustar core genes ────────
                if self.fmm_method == "stabilized":
                    stab_kwargs = dict(self.fmm_stab_kwargs)
                    if (not self.recalibrate_lambda_per_rep
                            and self.lambda_star_inherited is not None
                            and "lambda_star" not in stab_kwargs):
                        stab_kwargs["lambda_star"] = float(self.lambda_star_inherited)
                    prepared_k = prepare_dataset(esc_k, **stab_kwargs)
                    fmm = StabilizedFMMModel(prepared=prepared_k, refine=self.fmm_refine)

                core_present = [g for g in core_genes if g in top_k.index]
                fmm_fits = {}
                for g in core_present:
                    fmm_fits[g] = fmm.fit(
                        top_k.loc[g].values.astype(np.float64), esc_k,
                    )

                sync_k = prelim.run(
                    expr_ordered   = top_k,
                    circular_scale = esc_k,
                    core_genes     = core_present,
                    fmm_fits       = fmm_fits,
                )
                esc_fin = sync_k.circular_scale

                # ── Legitimidad del orden circular candidato ────────────
                max_gap_k = _max_circular_gap(esc_fin)
                legitimo = max_gap_k < pi and (
                    max_gap_pre is None or max_gap_k <= max_gap_pre
                )
                if legitimo or attempt >= _MAX_LEGITIMACY_ATTEMPTS:
                    if not legitimo:
                        ref_txt = (
                            f", referencia={max_gap_pre:.4f} rad"
                            if max_gap_pre is not None else ""
                        )
                        self._log(
                            f"  AVISO [legitimidad] repeticion {k+1}: ningun "
                            f"subconjunto cumplio las condiciones tras "
                            f"{attempt} intentos (hueco max={max_gap_k:.4f} "
                            f"rad{ref_txt}). Se acepta el ultimo intento."
                        )
                    break

            orientation_changed = bool(sync_k.direction_flipped)
            flipped_arr[k] = orientation_changed

            # Orden global (en indices originales 0..n-1)
            global_order = order_k[sync_k.sample_order]
            sample_orders[k] = global_order
            circ_scales[k]   = esc_fin
            sel_names[k]     = genes_k_arr

            # Pico ARNTL para rotacion de la escala
            arntl_peak_raw = (
                fmm_fits[self.anchor_gene].peak_time
                if self.anchor_gene in fmm_fits else 0.0
            )

            # ── 4. Re-ajuste FMM a cada gen del TOP ──────────
            fmm_mat = np.zeros((n_top, n_samp))

            esc_rot   = (esc_k - arntl_peak_raw + pi) % (2.0 * pi)
            order_rot = np.argsort(esc_rot)

            # Para stab reusamos el modelo ya preparado; para base, dict de kwargs
            fmm_obj = fmm if self.fmm_method == "stabilized" else self._fmm_kwargs

            gene_args = [
                (
                    i, gene,
                    top_k.loc[gene].values.astype(np.float64),
                    esc_k,
                    fmm_obj,
                    arntl_peak_raw,
                    orientation_changed,
                    order_rot,
                    n_samp,
                )
                for i, gene in enumerate(names_top)
            ]

            n_workers = os.cpu_count() if self.n_jobs == -1 else self.n_jobs
            if n_workers == 1:
                gene_results = [_fit_gene(a) for a in gene_args]
            elif self.fmm_method == "stabilized":
                # Threading: comparte prepared en memoria; BLAS/Numba liberan GIL
                from concurrent.futures import ThreadPoolExecutor
                with ThreadPoolExecutor(max_workers=n_workers) as pool:
                    gene_results = list(pool.map(_fit_gene, gene_args))
            else:
                # Multiproceso: argumentos pequenos por gen
                with ProcessPoolExecutor(max_workers=n_workers) as pool:
                    gene_results = list(pool.map(_fit_gene, gene_args))

            for res in gene_results:
                (i, gene, stat_fmm, fitted_fmm_reord) = res
                fmm_mat[i] = fitted_fmm_reord
                stats_table[k * n_top + i] = stat_fmm

            per_rep_fmm.append(pd.DataFrame(fmm_mat, index=names_top))

        return (sample_orders, circ_scales, flipped_arr,
                stats_table, per_rep_fmm, sel_names)

    # ==================================================================
    # Fase A.2b — Correccion de orientaciones
    # ==================================================================

    def _fix_orientations(
        self,
        sample_orders: np.ndarray,
        circ_scales:   np.ndarray,
        stats_table:   np.ndarray | None,
        per_rep_fmm:   list | None,
        theta_pre:     np.ndarray,
        n_samp:        int,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray | None, list | None]:
        """Corrige las repeticiones cuya orientacion salio invertida.

        La etapa de sincronizacion resuelve por separado el sentido de recorrido
        de cada uno de los K ordenes, y su regla puede fallar cuando el pico del
        gen de direccion queda cerca de la frontera entre semicirculos. Una
        repeticion asi invertida contamina todo lo que se construya despues, no
        solo la media circular de las fases sino tambien la matriz de costes de
        TSP3 y la hipermatriz de HODs, de modo que la correccion se aplica aqui,
        en un unico punto y antes de que las repeticiones se usen para nada.

        Criterio
        --------
        Una repeticion esta invertida si su correlacion circular de Fisher y
        Lee con el orden preliminar es negativa. Es el signo, no un
        umbral calibrable: equivale a escoger la orientacion que maximiza la
        concordancia con la referencia. El estadistico es invariante frente a
        rotaciones, asi que la decision no depende del origen de cada escala.

        Antes se usaba la correlacion de Jammalamadaka, que **tambien** es
        invariante a rotacion pero cuyo SIGNO no es fiable en este regimen:
        resta a cada escala su propia direccion media, y cuando las muestras
        cubren todo el ciclo la resultante es pequena (R ~ 0,13 con n=48) y
        esa direccion queda mal determinada. Medido sobre repeticiones
        genuinamente invertidas del pipeline, devolvia valores POSITIVOS de
        +0,11 a +0,36, con lo que la condicion nunca se cumplia y no se
        corregia ninguna repeticion. Medido sobre 720 escenarios de
        inversion conocida, Jammalamadaka acierta el signo el 95,3 % de las
        veces y Fisher-Lee el 100 %.

        Se corrige unicamente la orientacion, invirtiendo el orden y reflejando
        la escala, y nunca el origen. Rotar destruiria el anclaje biologico que
        la sincronizacion acaba de fijar situando el pico del gen ancla en pi.

        Que se refleja
        --------------
        La escala de muestras no es lo unico expresado en la orientacion
        antigua. Los parametros FMM por gen de esa repeticion, que son los que
        acaban en el atlas, tambien lo estan, y hay que reflejarlos aqui o el
        atlas sale especular respecto del orden que se publica a su lado.
        Bajo ``theta -> 2*pi - theta`` la onda FMM se transforma como::

            (M, A, alpha, beta, omega) -> (M, A, -alpha, -beta, omega)

        porque la tangente es impar, de modo que ``2*arctan(omega*tan(
        (-theta-alpha)/2)) = -2*arctan(omega*tan((theta+alpha)/2))`` y el
        signo pasa dentro del coseno, que es par. M y A no cambian, y R2
        tampoco: reflejar el eje temporal no altera el residuo. Pico y valle
        se recalculan de los parametros ya reflejados en vez de negarlos, que
        daria lo mismo pero por un camino que puede separarse del modelo.

        Sin esta parte el atlas quedaba invertido siempre que la correccion
        actuaba sobre la repeticion elegida. Medido en simulacion con n=48 y
        A/sigma=4, un caso en que el orden sale practicamente perfecto
        (MAE 0,25 h), el RMSE de la curva del atlas era 0,871 y al reflejarlo
        bajaba a 0,085.

        El orden preliminar sirve de referencia por estar anclado de forma
        independiente. Tomar como referencia la media de las propias K
        repeticiones da el mismo resultado mientras la mayoria este bien
        orientada, pero consolidaria la orientacion equivocada si no lo
        estuviera.

        En simulacion, corregir de este modo es inocuo cuando todas las
        repeticiones estan bien orientadas y recupera el consenso cuando alguna
        no lo esta. Solo perjudica cuando la dispersion entre repeticiones
        supera el radian, un regimen en el que el consenso ya carece de valor
        con correccion y sin ella.
        """
        K = sample_orders.shape[0]
        theta_pre = np.asarray(theta_pre, dtype=np.float64)
        if theta_pre.size != n_samp:
            self._log("  [orientacion] sin escala preliminar valida; no se corrige.")
            return sample_orders, circ_scales, stats_table, per_rep_fmm

        Theta = _build_theta(sample_orders, circ_scales, n_samp)
        corr  = np.array([_fisher_lee(Theta[k], theta_pre)
                          for k in range(K)])
        invertidas = np.flatnonzero(corr < 0.0)

        if invertidas.size == 0:
            self._log(
                f"  [orientacion] las {K} repeticiones concuerdan con el orden "
                f"preliminar (corr min = {corr.min():+.4f})"
            )
            return sample_orders, circ_scales, stats_table, per_rep_fmm

        so = sample_orders.copy()
        cs = circ_scales.copy()
        st = None if stats_table is None else np.asarray(stats_table).copy()
        pf = list(per_rep_fmm) if per_rep_fmm is not None else None
        n_top = (st.shape[0] // K) if (st is not None and K) else 0

        for k in invertidas:
            so[k] = sample_orders[k][::-1]
            cs[k] = (2.0 * pi - circ_scales[k][::-1]) % (2.0 * pi)

            # Parametros FMM de la rep, que van al atlas. Ver el docstring.
            bloque = st[k * n_top:(k + 1) * n_top] if n_top else None
            if bloque is not None:
                bloque[:, 2] = (-bloque[:, 2]) % (2.0 * pi)   # alpha
                bloque[:, 3] = (-bloque[:, 3]) % (2.0 * pi)   # beta
                for i in range(n_top):
                    al, be, om = bloque[i, 2], bloque[i, 3], bloque[i, 4]
                    bloque[i, 5] = FMMModel.peak_time(al, be, om)
                    bloque[i, 6] = FMMModel.trough_time(al, be, om)
                    bloque[i, 7] = bloque[i, 5] / (2.0 * pi) * 100.0
                    bloque[i, 8] = bloque[i, 6] / (2.0 * pi) * 100.0

            # Curvas ajustadas, indexadas por posicion en el orden de muestras
            if pf is not None and k < len(pf):
                pf[k] = pf[k].iloc[:, ::-1]
                pf[k].columns = range(pf[k].shape[1])

        self._log(
            f"  AVISO [orientacion]: {invertidas.size}/{K} repeticion(es) "
            f"invertida(s) respecto al orden preliminar, indices "
            f"{invertidas.tolist()} con correlaciones "
            f"{np.round(corr[invertidas], 4).tolist()}. Se corrigen. Revisa la "
            f"etapa de sincronizacion, en particular el gen de direccion "
            f"({self.direction_gene})."
        )
        return so, cs, st, pf

    # ==================================================================
    # Fase A.3 — Mediana R² (FMM) por repeticion
    # ==================================================================

    def _compute_median_r2(
        self,
        stats_table: np.ndarray,
        names_top:   list[str],
        eval_genes:  list[str],
        K:           int,
        n_top:       int,
    ) -> np.ndarray:
        """Mediana de R² (FMM) sobre ``eval_genes`` para cada repeticion.

        La columna r2_fmm en ``stats_table`` esta en posicion 11 (0-indexed)
        dentro de las 12 columnas FMM (5 params + 4 peaks + 3 metricas).
        """
        eval_indices = [i for i, name in enumerate(names_top) if name in eval_genes]

        median_r2 = np.zeros(K, dtype=np.float64)
        for k in range(K):
            start = k * n_top
            r2_vals = [stats_table[start + i, 11] for i in eval_indices]
            if r2_vals:
                median_r2[k] = float(np.median(r2_vals))
        return median_r2

    def _compute_n_well_fitted(
        self,
        stats_table: np.ndarray,
        names_top:   list[str],
        eval_genes:  list[str],
        K:           int,
        n_top:       int,
    ) -> np.ndarray:
        """Numero de genes de evaluacion con R² (FMM) >= ``selection_r2_threshold``.

        Mismo conjunto de genes y misma columna que ``_compute_median_r2``, de
        modo que ambos criterios sean comparables entre si.
        """
        eval_indices = [i for i, name in enumerate(names_top) if name in eval_genes]
        tau = self.selection_r2_threshold

        counts = np.zeros(K, dtype=np.float64)
        for k in range(K):
            start = k * n_top
            r2_vals = np.array(
                [stats_table[start + i, 11] for i in eval_indices], dtype=np.float64,
            )
            if r2_vals.size:
                counts[k] = float(np.count_nonzero(r2_vals >= tau))
        return counts

    def _rep_fmm_atlas(
        self,
        stats_table: np.ndarray,
        names_top:   list[str],
        rep:         int,
        n_top:       int,
    ) -> pd.DataFrame:
        """Atlas de best_k: parametros FMM por gen de UNA repeticion (la elegida
        ``rep``), sin agregar. Son los parametros del propio ajuste de esa rep,
        consistentes con su orden de muestras (a diferencia de promediar sobre
        las K, que mezclaria los ordenes ya descartados por peores).

        Columnas de ``stats_table`` (FMM): ``0=M, 1=A, 2=alpha, 3=beta,
        4=omega, 5=t_U, 6=t_L, 11=R2``.
        """
        block = stats_table[rep * n_top:(rep + 1) * n_top]
        return pd.DataFrame({
            "gene_name":        list(names_top),
            "M":                block[:, 0],
            "amplitude":        block[:, 1],   # A
            "alpha":            block[:, 2],
            "beta":             block[:, 3],
            "omega":            block[:, 4],
            "peak_phase_rad":   block[:, 5],   # t_U
            "trough_phase_rad": block[:, 6],   # t_L
            "r2_fmm":           block[:, 11],
        })

    def _median_fmm_atlas(
        self,
        stats_table: np.ndarray,
        names_top:   list[str],
        K:           int,
        n_top:       int,
    ) -> pd.DataFrame:
        """Atlas de CIRCUST original: mediana por gen sobre las K repeticiones.

        Es el paso 4 de Larriba et al. 2023: «[M_TOP] is the matrix that
        contains for the genes in the TOP the median of the FMM features: R2,
        tU and omega». Reproduce ``robustSincroDBP_v3_cores`` de la
        implementacion de referencia, que agrega los parametros por gen y
        descarta los K ordenes de muestras.

        Los parametros angulares (alpha, beta, t_U, t_L) se agregan con
        **mediana circular**, no aritmetica: son angulos en [0, 2*pi) y la
        mediana ordinaria se rompe en el corte 0/2*pi (la mediana de 0.1 y
        6.2 radianes es 3.15, diametralmente opuesta a las dos). El resto (M,
        A, omega, R2) son escalares y llevan mediana ordinaria.

        Columnas de ``stats_table`` (FMM): ``0=M, 1=A, 2=alpha, 3=beta,
        4=omega, 5=t_U, 6=t_L, 11=R2``.
        """
        # (K*n_top, 12) -> (K, n_top, 12) para agregar en el eje de las reps.
        # stats_table se rellena por bloques de rep: filas [k*n_top, (k+1)*n_top).
        cubo = stats_table[:K * n_top].reshape(K, n_top, -1)

        def med_lineal(col: int) -> np.ndarray:
            return np.median(cubo[:, :, col], axis=0)

        def med_circular(col: int) -> np.ndarray:
            return np.array([circular_median(cubo[:, i, col]) for i in range(n_top)])

        return pd.DataFrame({
            "gene_name":        list(names_top),
            "M":                med_lineal(0),
            "amplitude":        med_lineal(1),
            "alpha":            med_circular(2),
            "beta":             med_circular(3),
            "omega":            med_lineal(4),
            "peak_phase_rad":   med_circular(5),   # t_U
            "trough_phase_rad": med_circular(6),   # t_L
            "r2_fmm":           med_lineal(11),
        })

    # ==================================================================
    # Fase B — TAREA 1A: best_k por R² mediano
    # ==================================================================

    def _select_best_k(
        self,
        sample_orders: np.ndarray,
        circ_scales:   np.ndarray,
        median_r2:     np.ndarray,
        n_well_fitted: np.ndarray,
        n_samp:        int,
    ) -> tuple[np.ndarray, np.ndarray, int, float]:
        """Elige una de las K repeticiones segun ``selection_criterion``.

        Es la opcion mas simple y rapida: sin agregacion. Devuelve
        ``(sample_order_k*, circular_scale_k*, k*, score_k*)``, donde
        ``score_k*`` es el valor del criterio en la repeticion elegida.
        """
        crit = self.selection_criterion

        if crit == "r2_median":
            best_k = int(np.argmax(median_r2))
            score  = float(median_r2[best_k])
            self._log(
                f"  best_k [r2_median]: rep {best_k + 1} seleccionada "
                f"(mediana R² = {score:.4f})"
            )

        elif crit == "n_well_fitted":
            # Maximo numero de genes bien ajustados. Ante empate (frecuente con
            # pocos genes de evaluacion) se desempata por mediana de R², que
            # ordena dentro del grupo empatado sin volver a mirar el umbral.
            best_count = int(n_well_fitted.max())
            tied       = np.flatnonzero(n_well_fitted == best_count)
            best_k     = int(tied[np.argmax(median_r2[tied])])
            score      = float(n_well_fitted[best_k])
            msg_tie    = f", {tied.size} empatadas resueltas por mediana R²" if tied.size > 1 else ""
            self._log(
                f"  best_k [n_well_fitted]: rep {best_k + 1} seleccionada "
                f"({best_count} genes con R² >= {self.selection_r2_threshold:g}"
                f"{msg_tie})"
            )

        else:  # "medoid"
            mean_dist = self._compute_mean_cosine_distance(
                sample_orders, circ_scales, n_samp,
            )
            best_k = int(np.argmin(mean_dist))
            score  = float(mean_dist[best_k])
            self._log(
                f"  best_k [medoid]: rep {best_k + 1} seleccionada "
                f"(distancia coseno media = {score:.4f})"
            )

        return (
            sample_orders[best_k].copy(),
            circ_scales[best_k].copy(),
            best_k,
            score,
        )

    @staticmethod
    def _compute_mean_cosine_distance(
        sample_orders: np.ndarray,
        circ_scales:   np.ndarray,
        n_samp:        int,
    ) -> np.ndarray:
        """Distancia coseno media de cada repeticion al resto de las K.

        Para cada par de repeticiones (k, l) se promedia sobre las muestras la
        distancia coseno ``d(a, b) = 1 - cos(a - b)``, y despues se promedia
        sobre ``l``. La comparacion muestra a muestra es licita porque las K
        escalas ya han pasado por la etapa de sincronizacion, que fija un
        origen y un sentido de recorrido comunes.

        Devuelve un vector de longitud K (el termino l = k aporta 0).
        """
        Theta = _build_theta(sample_orders, circ_scales, n_samp)
        K     = Theta.shape[0]
        if K < 2:
            return np.zeros(K, dtype=np.float64)

        # D[k, l] = media_i (1 - cos(Theta[k,i] - Theta[l,i]))
        diff = Theta[:, None, :] - Theta[None, :, :]
        D    = 1.0 - np.cos(diff)
        return D.mean(axis=2).sum(axis=1) / (K - 1)

    # ==================================================================
    # Fase B — TAREA 1A (extra): consenso del orden por media circular
    #                            (opcion A) + contraste vs preliminar
    # ==================================================================

    def _circmean_consensus_and_contrast(
        self,
        sample_orders: np.ndarray,
        circ_scales:   np.ndarray,
        best_order:    np.ndarray,
        best_scale:    np.ndarray,
        theta_pre:     np.ndarray,
        n_samp:        int,
    ) -> dict:
        """Opcion A (solo best_k): agrega las K ordenaciones por media circular
        y contrasta el consenso frente al orden preliminar CPCA.

        Pasos:
          1. ``Theta`` (K x n_samp) con la fase de cada muestra por repeticion.
          2. Consenso por media circular alineada (Procrustes) al preliminar;
             ``order_A = argsort(theta_bar)`` y ``R_i`` como confianza por muestra.
          3. Contraste preliminar via correlacion circular (Jammalamadaka),
             tanto del consenso A como del orden best_k, para poder comparar
             ambas estrategias contra el mismo referente.

        Devuelve un dict con las claves que rellenan ``RobustOrderResult`` (o
        ``{}`` si no hay escala preliminar utilizable).
        """
        theta_pre = np.asarray(theta_pre, dtype=np.float64)
        if theta_pre.size != n_samp:
            self._log("  [circmean] sin escala preliminar valida; se omite consenso A.")
            return {}

        Theta = _build_theta(sample_orders, circ_scales, n_samp)
        theta_bar, resultant = _circular_mean_consensus(Theta)

        order_A = np.argsort(theta_bar)
        scale_A = theta_bar[order_A]
        resultant_A = resultant[order_A]

        # Fases por muestra del orden best_k, para un contraste comparable
        theta_bestk = np.zeros(n_samp, dtype=np.float64)
        theta_bestk[np.asarray(best_order, dtype=int)] = np.asarray(best_scale, dtype=np.float64)

        out = {
            "circmean_order":           order_A,
            "circmean_scale":           scale_A,
            "circmean_resultant":       resultant_A,
            "circmean_vs_prelim_corr":  _fisher_lee(theta_bar, theta_pre),
            "bestk_vs_prelim_corr":     _fisher_lee(theta_bestk, theta_pre),
        }
        self._log(
            f"  [circmean] consenso A vs prelim: corr="
            f"{out['circmean_vs_prelim_corr']:.4f}"
            f"  |  best_k vs prelim: corr={out['bestk_vs_prelim_corr']:.4f}"
        )
        self._log(
            f"  [circmean] R medio={float(resultant.mean()):.3f} "
            f"(min={float(resultant.min()):.3f})"
        )
        return out

    # ==================================================================
    # Fase B — TAREA 1B: aggregate via consenso circular
    # ==================================================================

    def _aggregate_orders_via_consensus(
        self,
        sample_orders: np.ndarray,
        circ_scales:   np.ndarray,
        n_samp:        int,
    ) -> tuple[np.ndarray, np.ndarray, float, str]:
        """Agrega los K ordenes de MUESTRAS en uno solo via TSP3/HODs.

        Aplicacion de la tecnica de Barragan et al. 2015 al problema de
        orden temporal de muestras (propuesta del TFG): los items que se
        agregan son las n muestras, y la matriz ``Theta`` tiene forma
        ``(K, n_samp)`` donde ``Theta[k, i]`` es la fase asignada a la
        muestra ``i`` en la repeticion ``k``.

        Despues deriva la escala circular asociada al orden consenso segun
        ``consensus_scale_mode``:

        - ``"circmean"`` (defecto): media circular de las K fases de cada
          muestra. Es el estimador de fase, no una relectura del orden, de modo
          que ``argsort(circular_scale)`` puede no coincidir con
          ``sample_order``. La discrepancia se mide y se registra.
        - ``"cire"``: proyeccion de esa media sobre el conjunto de vectores que
          verifican el orden consenso (CIRE de Rueda et al.). Garantiza
          ``phi_{o_1} <= ... <= phi_{o_n}`` como pide el anexo, a costa de
          agrupar muestras en bloques con la misma fase.

        Sobre la eleccion del modo: en simulacion, la proyeccion CIRE multiplica
        por entre 2 y 5 el error mediano de fase frente a la media circular
        (p.ej. 0,18 h -> 0,80 h con n=100, K=5), porque obliga a las fases a
        respetar un orden que es a su vez una estimacion ruidosa. Por eso el
        defecto es ``"circmean"`` y la consistencia estricta es opcional.
        """
        K = sample_orders.shape[0]

        # Construir Theta: para la rep k, el angulo asignado a la muestra i
        # es circ_scales[k, pos] donde pos es la posicion de i en
        # sample_orders[k] (Theta[k, sample_orders[k]] = circ_scales[k]).
        Theta = _build_theta(sample_orders, circ_scales, n_samp)

        self._log(f"    Construyendo matriz Theta ({K} x {n_samp} muestras) y "
                  f"agregando con method={self.aggregation_method}...")
        out = aggregate_circular_orders(
            Theta=Theta, weights=None, method=self.aggregation_method,
        )
        consensus_order = out["order"]
        msce_consensus  = out["msce"]
        selected_method = out["selected_method"]
        if self.aggregation_method == "auto":
            detalle = ", ".join(
                f"{m}={v:.4f}" for m, v in out["msce_by_method"].items()
            )
            self._log(f"    Auto-seleccion por MSCE [{detalle}] -> gana {selected_method}")
        self._log(f"    Consenso (muestras): MSCE = {msce_consensus:.4f} "
                  f"({selected_method})")

        # Escala circular del orden consenso: media circular, por muestra, de
        # las K fases asignadas. ``Theta[k, j]`` ya es la fase de la muestra j
        # en la repeticion k, asi que basta reordenar por el orden consenso.
        # (Antes se recorria un bucle K x n con una busqueda ``np.where``
        # dentro, O(K n^2), para calcular exactamente esto mismo.)
        Th = Theta[:, consensus_order]                      # (K, n_samp)
        final_scale = np.arctan2(
            np.sin(Th).mean(axis=0), np.cos(Th).mean(axis=0),
        ) % (2.0 * pi)

        if self.consensus_scale_mode == "cire":
            # Proyeccion sobre el conjunto de vectores que verifican el orden
            # consenso: garantiza phi_{o_1} <= ... <= phi_{o_n} (anexo CPCA).
            theta_by_item = np.empty(n_samp, dtype=np.float64)
            theta_by_item[consensus_order] = final_scale
            fitted, _ = circular_isotonic_regression(theta_by_item, consensus_order)
            final_scale = fitted[consensus_order]

        # Diagnostico: en modo "circmean" el orden que induce la escala no tiene
        # por que coincidir con el orden consenso. Se mide para poder reportarlo.
        n_desc = int(np.sum(np.diff(final_scale) < 0.0))
        self._log(
            f"    Escala consenso [{self.consensus_scale_mode}]: {n_desc} "
            f"descenso(s) a lo largo del orden consenso (1 = solo el cierre "
            f"del circulo)"
        )
        if self.consensus_scale_mode != "cire" and n_desc > 1:
            self._log(
                "    AVISO: sample_order (consenso TSP3/HODs) y "
                "argsort(circular_scale) no coinciden. Son dos estimadores "
                "distintos. Usa consensus_scale_mode='cire' si necesitas que "
                "la escala verifique el orden por construccion."
            )

        return consensus_order, final_scale, msce_consensus, selected_method

    # ==================================================================
    # Utilidades
    # ==================================================================

    def _log(self, message: str) -> None:
        if self.verbose:
            print(message, flush=True)
