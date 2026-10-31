"""
circust/top_genes.py
====================
Etapa 3: Determinacion de genes TOP ritmicos especificos del tejido.

Implementa la seleccion de genes TOP segun la Seccion 3 de CIRCUST:

  1. Ajuste ORI: para cada gen, ajustar regresion isotonica circular
     unimodal (PAVA) y calcular R2_ORI = 1 - MSE_unimodal / MSE_plano.

  2. Filtro ORI: descartar genes con R2_ORI < 0.5 (ritmicidad no
     parametrica insuficiente).

  3. Filtro FMM: ajustar el modelo FMM a los supervivientes y
     seleccionar aquellos que sean
     (i) no puntiagudos (omega > 0.1),
     (ii) altamente ritmicos (R2_FMM >= 0.5),
     (iii) picos (t_U) distribuidos en los 4 cuadrantes de [0, 2pi).

  4. Genes semilla forzados: los 12 genes reloj centrales se
     incluyen siempre; si alguno no paso los filtros, se fuerza su
     inclusion.

Posicion en el pipeline
-----------------------
    CircularSynchronizer  ->  **TopGeneSelector**  ->  RandomSelector
                                                   ->  RobustEstimator
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from dataclasses import dataclass, field
from typing import Optional

from math import pi

from circust.fitting.ori import ORIModel
from circust.fitting.fmm            import FMMModel
from circust.fitting.circafmm import CircaFMMModel, prepare_dataset


# ═══════════════════════════════════════════════════════════════════════════
# Funciones auxiliares de estadistica circular
# ═══════════════════════════════════════════════════════════════════════════

def circular_median(angles: np.ndarray) -> float:
    """
    Mediana circular: el angulo theta que minimiza  sum d(theta, theta_i)
    donde d es la distancia geodesica en el circulo unitario.

    Parameters
    ----------
    angles : array (n,) en [0, 2pi).

    Returns
    -------
    float — mediana circular en [0, 2pi).
    """
    angles = np.asarray(angles, dtype=np.float64) % (2.0 * pi)
    n = len(angles)
    if n == 0:
        return 0.0
    if n == 1:
        return float(angles[0])

    best_cost  = np.inf
    best_angle = angles[0]
    for candidate in angles:
        diffs = np.abs(angles - candidate)
        diffs = np.minimum(diffs, 2.0 * pi - diffs)
        cost  = diffs.sum()
        if cost < best_cost:
            best_cost  = cost
            best_angle = candidate

    return float(best_angle % (2.0 * pi))


def assign_sectors(
    peaks: np.ndarray,
    reference: float,
    n_sectors: int = 8,
) -> np.ndarray:
    """
    Asigna cada angulo de pico a uno de *n_sectors* sectores equiespaciados
    centrados en *reference*.

    Los sectores se numeran 1 ... n_sectors.

    Parameters
    ----------
    peaks     : array (n,) de angulos en [0, 2pi).
    reference : float — mediana circular usada como origen.
    n_sectors : int — por defecto 8.

    Returns
    -------
    array int (n,) — etiquetas de sector en {1, 2, ..., n_sectors}.
    """
    peaks = np.asarray(peaks, dtype=np.float64)
    rotated = (peaks - reference) % (2.0 * pi)
    sector_width = 2.0 * pi / n_sectors
    sectors = n_sectors - np.floor(rotated / sector_width).astype(int)
    sectors = np.clip(sectors, 1, n_sectors)
    return sectors


def sector_centres(reference: float, n_sectors: int = 8) -> np.ndarray:
    """
    Centros angulares de los sectores definidos por :func:`assign_sectors`.

    La indexacion es la que cabe esperar: ``sector_centres(ref, n)[L - 1]`` es
    el centro del sector cuya etiqueta es ``L`` en {1, ..., n}. La version
    anterior devolvia las fronteras inferiores y con las etiquetas desplazadas
    en una posicion, de modo que el centro asociado a la etiqueta ``L`` no
    coincidia con el elemento ``L - 1``.

    Deduccion: ``assign_sectors`` asigna la etiqueta ``L = n - floor(r / w)``
    con ``r = (pico - ref) mod 2pi`` y ``w = 2pi/n``, de manera que la etiqueta
    ``L`` cubre ``r`` en ``[(n - L) w, (n - L + 1) w)`` y su centro esta en
    ``ref + (n - L + 1/2) w``.

    Parameters
    ----------
    reference : float — mediana circular usada como origen.
    n_sectors : int — por defecto 8.

    Returns
    -------
    array (n_sectors,) — centro del sector ``L`` en la posicion ``L - 1``.
    """
    width = 2.0 * pi / n_sectors
    return np.array([
        (reference + (n_sectors - L + 0.5) * width) % (2.0 * pi)
        for L in range(1, n_sectors + 1)
    ])


def apply_sector_quota(
    sel:         np.ndarray,
    sectors:     np.ndarray,
    r2_fmm:      np.ndarray,
    coverage:    np.ndarray,
    n_sectors:   int,
    target_size: int,
) -> tuple[np.ndarray, dict[int, int]]:
    """Reduce ``sel`` a ``target_size`` genes repartiendo el cupo por sector.

    El cupo se reparte a partes iguales entre los sectores no vacios y dentro
    de cada uno se ordena por ``R2_FMM`` descendente. Las plazas que un sector
    no llegue a cubrir, por tener menos genes que su cupo, se redistribuyen
    entre los que aun tienen candidatos, de modo que se alcanza ``target_size``
    siempre que existan genes suficientes.

    Se reparte por sectores y no con un ranking global porque un ranking global
    tiende a concentrarse en la region del ciclo donde la senal es mas fuerte,
    deshaciendo la cobertura circular que exige el paper.

    Esta funcion es la unica implementacion del cupo: la usan tanto
    :class:`TopGeneSelector` como los scripts de calibracion de
    ``(n_sectors, target_size)``, para que el barrido de sensibilidad no pueda
    divergir del comportamiento real del pipeline.

    Parameters
    ----------
    sel : np.ndarray de int
        Indices de los genes candidatos (los que superan los filtros, mas los
        incorporados por la regla de cobertura).

    sectors : np.ndarray de int
        Etiqueta de sector de *cada* gen, indexable por los valores de ``sel``.

    r2_fmm : np.ndarray
        R2 del ajuste FMM de cada gen, indexable por los valores de ``sel``.

    coverage : np.ndarray de int
        Indices de los genes que entraron solo por la regla de cobertura. Se
        les da plaza garantizada: son el unico gen de su sector y descartarlos
        dejaria un hueco en la circunferencia.

    n_sectors, target_size : int

    Returns
    -------
    (np.ndarray de int, dict)
        Indices retenidos en orden creciente, y el cupo concedido a cada
        sector, para poder registrarlo en los logs.
    """
    # Candidatos de cada sector ordenados por R2_FMM descendente, con los
    # genes de cobertura delante para que su plaza quede garantizada.
    cov_set = set(coverage.tolist())
    by_sector: dict[int, np.ndarray] = {}
    for s in range(1, n_sectors + 1):
        idx = sel[sectors[sel] == s]
        if idx.size == 0:
            continue
        idx  = idx[np.argsort(-r2_fmm[idx], kind="stable")]
        prio = np.array([i in cov_set for i in idx])
        by_sector[s] = np.concatenate([idx[prio], idx[~prio]])

    quota  = {s: 0 for s in by_sector}
    libres = target_size
    while libres > 0:
        con_plazas = [s for s in by_sector if quota[s] < by_sector[s].size]
        if not con_plazas:
            break
        reparto = max(1, libres // len(con_plazas))
        for s in con_plazas:
            if libres == 0:
                break
            concedidas = min(reparto, by_sector[s].size - quota[s], libres)
            quota[s] += concedidas
            libres   -= concedidas

    chosen: list[int] = []
    for s, idx in by_sector.items():
        chosen.extend(idx[: quota[s]].tolist())

    return np.array(sorted(chosen), dtype=int), quota


# ═══════════════════════════════════════════════════════════════════════════
# Workers a nivel de modulo (serializables por joblib)
# ═══════════════════════════════════════════════════════════════════════════

def _fit_ori_one(
    row: np.ndarray,
    circ: np.ndarray,
) -> float:
    """
    Ajusta ORI a un unico gen. Devuelve R2_ORI.

    Funcion pura a nivel de modulo para que ``joblib`` pueda serializarla
    hacia procesos worker.
    """
    fr = ORIModel().fit(row, circ)
    return fr.r2


def _fit_fmm_one(
    row: np.ndarray,
    circ: np.ndarray,
    alpha_grid: int,
    omega_grid: int,
    num_reps: int,
) -> tuple[float, float, float, float, float, float, float]:
    """Ajusta FMM base a un unico gen. Devuelve ``(peak_time, omega, r2, A)``.

    Funcion pura a nivel de modulo para que ``joblib`` pueda serializarla
    hacia procesos worker (FMM base puede paralelizar con muchos genes).

    Devuelve ``(peak_time, omega, r2, A, M, alpha, beta)``. Los tres ultimos
    se anaden al final para que los consumidores que indexan por posicion no
    se rompan. Hacen falta para que el atlas pueda publicar los cinco
    parametros del modelo: con solo (A, t_U, omega) la onda no queda
    determinada, porque beta gobierna la asimetria entre subida y bajada.

    No devuelve el vector ``fitted``: no lo consume nadie y, con decenas de
    miles de genes, obliga a serializar de vuelta un vector de n muestras por
    gen desde cada proceso worker.
    """
    model = FMMModel(
        length_alpha_grid=alpha_grid,
        length_omega_grid=omega_grid,
        num_reps=num_reps,
    )
    fr = model.fit(row, circ)
    return (fr.peak_time, fr.params["omega"], fr.r2, fr.params["A"],
            fr.params["M"], fr.params["alpha"], fr.params["beta"])


def _fit_fmm_one_circafmm(
    row:      np.ndarray,
    circ:     np.ndarray,
    circafmm_model: CircaFMMModel,
) -> tuple[float, float, float, float, float, float, float]:
    """Ajusta CircaFMM reusando un ``CircaFMMModel`` ya preparado.

    Mismo formato de salida que :func:`_fit_fmm_one`.
    """
    fr = circafmm_model.fit(row, circ)
    return (fr.peak_time, fr.params["omega"], fr.r2, fr.params["A"],
            fr.params["M"], fr.params["alpha"], fr.params["beta"])


# ═══════════════════════════════════════════════════════════════════════════
# Dataclass de resultado
# ═══════════════════════════════════════════════════════════════════════════

@dataclass
class TopGeneResult:
    """
    Salida de :class:`TopGeneSelector`.

    Los campos principales son compatibles por duck-typing con
    ``ReferenceSetResult``, de modo que ``RandomSelector`` puede
    consumir este resultado directamente.

    Atributos — pipeline
    --------------------
    gene_names : list[str]
        Simbolos de los genes TOP seleccionados.

    fmm_peaks : np.ndarray (n_top,)
        Tiempos de pico FMM (t_U) por gen.

    r2_par : np.ndarray (n_top,)
        R2 parametrico (FMM) por gen.

    sector_labels : np.ndarray int (n_top,)
        Cuadrante asignado (1–4) por gen.

    sector_centres : np.ndarray (n_sectors,)
        Centros de los cuadrantes en [0, 2pi).

    candidate_matrix : pd.DataFrame (n_top, n_muestras)
        Matriz de expresion de los genes TOP (en el orden circular
        preliminar). Equivale a [X_TOP] del paper.

    circular_scale : np.ndarray (n_muestras,)
        Eje temporal circular del ordenamiento preliminar.

    added_genes : list[str]
        Genes semilla que fueron forzados (no pasaron los filtros
        por si solos). Se excluyen del pool de muestreo aleatorio
        en ``RandomSelector``.

    sample_order : np.ndarray int (n_muestras,)
        Orden de muestras del sincronizado preliminar.

    Atributos — diagnostico
    -----------------------
    omega : np.ndarray (n_top,)
        Parametro omega de FMM por gen.

    r2_ori : np.ndarray (n_top,)
        R2 ORI (no parametrico) de cada gen TOP.

    n_total_genes : int
        Numero total de genes antes de cualquier filtro.

    n_after_ori : int
        Genes que superaron el filtro R2_ORI >= umbral.

    n_after_fmm : int
        Genes que superaron los filtros FMM (omega + R2_FMM)
        *antes* de forzar genes semilla.
    """

    # Pipeline (compatibles con ReferenceSetResult por duck-typing)
    gene_names:       list          = field(default_factory=list)
    fmm_peaks:        np.ndarray    = field(default_factory=lambda: np.array([]))
    r2_par:           np.ndarray    = field(default_factory=lambda: np.array([]))
    sector_labels:    np.ndarray    = field(default_factory=lambda: np.array([], dtype=int))
    sector_centres:   np.ndarray    = field(default_factory=lambda: np.zeros(4))
    candidate_matrix: pd.DataFrame  = field(default_factory=pd.DataFrame)
    circular_scale:   np.ndarray    = field(default_factory=lambda: np.array([]))
    added_genes:      list          = field(default_factory=list)
    sample_order:     np.ndarray    = field(default_factory=lambda: np.array([], dtype=int))

    # Diagnostico
    omega:            np.ndarray    = field(default_factory=lambda: np.array([]))
    r2_ori:           np.ndarray    = field(default_factory=lambda: np.array([]))
    n_total_genes:    int           = 0
    n_after_ori:      int           = 0
    n_after_fmm:      int           = 0

    def summary(self) -> str:
        n = len(self.gene_names)
        n_sec = len(self.sector_centres)
        n_covered = len(np.unique(self.sector_labels)) if n > 0 else 0
        lines = [
            "=== Resumen de Genes TOP ===",
            f"  Genes totales            : {self.n_total_genes}",
            f"  Tras filtro R2_ORI>=0.5  : {self.n_after_ori}",
            f"  Tras filtros FMM         : {self.n_after_fmm}",
            f"  Genes semilla forzados   : {self.added_genes or 'ninguno'}",
            f"  Genes TOP finales        : {n}",
            f"  Cuadrantes cubiertos     : {n_covered} / {n_sec}",
        ]
        if n > 0:
            for s in sorted(np.unique(self.sector_labels)):
                cnt = int(np.sum(self.sector_labels == s))
                lines.append(f"    Cuadrante {s}: {cnt} gen(es)")
        return "\n".join(lines)


# ═══════════════════════════════════════════════════════════════════════════
# TopGeneSelector
# ═══════════════════════════════════════════════════════════════════════════

class TopGeneSelector:
    """
    Selecciona los genes TOP ritmicos especificos del tejido.

    Algoritmo (Seccion 3 del paper CIRCUST)
    ----------------------------------------
    1. Ajustar ORI a cada gen y calcular R2_ORI.
    2. Descartar genes con R2_ORI < ``r2_ori_threshold`` (0.5).
    3. Ajustar FMM a los supervivientes.
    4. Filtrar: omega > ``omega_min`` (0.1) Y R2_FMM >= ``r2_fmm_threshold`` (0.5).
    5. Repartir los picos FMM en ``n_sectors`` sectores y quedarse, dentro de
       cada uno, con los de mayor R2_FMM hasta completar el cupo, que es la
       fraccion ``top_fraction`` del total de genes de la matriz.
    6. Forzar inclusion de los genes core si no estan presentes.
    7. Construir la matriz TOP ([X_TOP] en el paper).

    Revision antes de publicar
    Los umbrales del paso 2 y 4 son los que fija el paper (0.5, 0.5, 0.1) y no
    controlan el tamano del conjunto: con pocas muestras un FMM de 5 parametros
    supera R2 > 0.5 en buena parte del transcriptoma. El tamano lo fija el cupo
    del paso 5, que es la lectura correcta del esquema del material suplementario
    (Fig B, pasos 3.2-3.4: dividir el circulo en intervalos, recorrer cada uno
    por R2 decreciente y quedarse con los de mayor R2_FMM del tejido). El
    articulo describe esa construccion pero no cuantifica ni cuantos intervalos
    ni cuantos genes, de ahi que ``n_sectors`` y ``top_fraction`` sean los dos
    parametros a calibrar.

    Parameters
    ----------
    r2_ori_threshold : float
        Umbral minimo de R2 ORI para el pre-filtro no parametrico.
        Paper: 0.5.

    r2_fmm_threshold : float
        Umbral minimo de R2 FMM para la ritmicidad parametrica.
        Paper: 0.5.

    omega_min : float
        omega minimo de FMM; rechaza ajustes puntiagudos.
        Paper: 0.1.

    n_sectors : int
        Numero de sectores circulares en los que se divide la circunferencia,
        con origen en la mediana circular de los picos. Cumplen dos funciones:
        detectar sectores sin ningun gen ritmico (se avisa por consola, no se
        rellena con un gen que no supere omega_min/r2_fmm_threshold, ver
        docstring de ``_select_top_indices``) y estratificar el cupo de
        el cupo.

        El texto principal del paper exige cubrir los cuatro cuadrantes, de
        modo que 4 es el minimo admisible. Por defecto 8.

    top_fraction : float
        Fraccion, en (0, 1], del numero total de genes de la matriz que se
        retiene en el conjunto TOP: el cupo es ``round(top_fraction * n_total)``
        y nunca baja de ``n_sectors``. Solo se retiene esa fraccion, no todos
        los genes que superan los filtros (si son menos que el cupo se
        conservan todos). El cupo se reparte entre los sectores no vacios y
        dentro de cada uno se ordena por R2_FMM descendente. El reparto por
        sectores, y no un ranking global, evita que el cupo deshaga la
        cobertura: un ranking global tiende a concentrarse en la region del
        ciclo donde la senal es mas fuerte.

        Los genes core forzados (paso 6) se anaden *despues*, igual que en la
        implementacion de referencia, de modo que el tamano final es el cupo
        mas los core que no hubiesen entrado por si solos.

    fmm_length_alpha_grid, fmm_length_omega_grid, fmm_num_reps : int
        Parametros de la rejilla de ajuste FMM.

    fmm_refine : bool
        Refinamiento local opcional de ``(alpha, omega)`` de CircaFMM
        (Seccion S3.4.4 del suplemento de CircaFMM). Por defecto False. Se
        ignora con ``fmm_method="base"``.

    n_jobs : int
        Numero de procesos para el ajuste en paralelo (ORI y FMM).
        ``-1`` usa todos los nucleos disponibles.

    verbose : bool
        Imprimir mensajes de progreso.
    """

    def __init__(
        self,
        r2_ori_threshold:      float = 0.5,
        r2_fmm_threshold:      float = 0.5,
        omega_min:             float = 0.1,
        n_sectors:             int   = 8,
        top_fraction:          float = 0.1,
        fmm_length_alpha_grid: int   = 48,
        fmm_length_omega_grid: int   = 24,
        fmm_num_reps:          int   = 3,
        fmm_method:            str   = "base",
        circafmm_kwargs:       dict|None = None,
        lambda_star:           Optional[float] = None,
        fmm_refine:            bool  = False,
        n_jobs:                int   = -1,
        verbose:               bool  = True,
    ) -> None:
        if fmm_method not in ("base", "circafmm"):
            raise ValueError(f"fmm_method debe ser 'base' o 'circafmm', recibido: {fmm_method!r}")
        self.r2_ori_threshold      = r2_ori_threshold
        self.r2_fmm_threshold      = r2_fmm_threshold
        self.omega_min             = omega_min
        if n_sectors < 4:
            raise ValueError(
                "n_sectors debe ser al menos 4: el paper exige que los picos "
                f"cubran los cuatro cuadrantes. Recibido: {n_sectors!r}"
            )
        self.n_sectors             = n_sectors
        if not 0.0 < top_fraction <= 1.0:
            raise ValueError(
                f"top_fraction debe estar en (0, 1]. Recibido: {top_fraction!r}"
            )
        self.top_fraction          = top_fraction
        self.fmm_length_alpha_grid = fmm_length_alpha_grid
        self.fmm_length_omega_grid = fmm_length_omega_grid
        self.fmm_num_reps          = fmm_num_reps
        self.fmm_method            = fmm_method
        # Si se proporciona lambda_star, lo inyectamos en circafmm_kwargs para
        # saltar la calibracion (es lo que hace el pipeline tras heredarlo de CPCA)
        self.circafmm_kwargs       = dict(circafmm_kwargs or {})
        if lambda_star is not None and "lambda_star" not in self.circafmm_kwargs:
            self.circafmm_kwargs["lambda_star"] = float(lambda_star)
        self.fmm_refine            = fmm_refine
        self.n_jobs                = n_jobs
        self.verbose               = verbose

    # ------------------------------------------------------------------
    # API publica
    # ------------------------------------------------------------------

    def run(
        self,
        expr_norm:       pd.DataFrame,
        circular_scale:  np.ndarray,
        seed_genes:      list[str],
        sample_order:    np.ndarray | None = None,
    ) -> TopGeneResult:
        """
        Selecciona los genes TOP.

        Internamente ajusta ORI (regresion isotonica circular unimodal)
        a cada gen para calcular R2_ORI, aplica los filtros y construye
        el conjunto TOP.

        Parameters
        ----------
        expr_norm : pd.DataFrame (n_genes, n_muestras)
            Matriz de expresion normalizada completa, ya ordenada
            por el tiempo circular preliminar.

        circular_scale : array (n_muestras,)
            Eje temporal circular del ordenamiento preliminar.

        seed_genes : list[str]
            Genes core del pipeline (los mismos seleccionados por
            ``CoreGeneSelector``). Se fuerzan en el conjunto TOP
            si no superan los filtros.

        sample_order : array int (n_muestras,) o None
            Indices de muestra del ordenamiento preliminar.
            Si es None se asume identidad.

        Returns
        -------
        TopGeneResult
        """
        seed_genes = list(dict.fromkeys(seed_genes))
        self._log("=== Etapa 3: Seleccion de Genes TOP ===")

        gene_names = list(expr_norm.index)
        circ       = np.asarray(circular_scale, dtype=np.float64)
        n_total    = len(gene_names)
        values     = expr_norm.values.astype(np.float64)

        if sample_order is None:
            sample_order = np.arange(expr_norm.shape[1], dtype=int)

        # ── 1. Ajuste ORI (regresion isotonica circular unimodal) ──
        self._log(f"  Genes totales: {n_total}")
        self._log("  Ajustando ORI (PAVA unimodal) ...")

        r2_ori = self._fit_ori_batch(values, circ)

        # ── 2. Filtro R2_ORI ───────────────────────────────────────
        self._log(
            f"  Filtro R2_ORI >= {self.r2_ori_threshold} ..."
        )

        mask_ori     = r2_ori >= self.r2_ori_threshold  # criterio: R2_ORI >= tau_ORI
        surv_genes   = [g for g, m in zip(gene_names, mask_ori) if m]
        surv_r2_ori  = r2_ori[mask_ori]
        n_after_ori  = len(surv_genes)

        self._log(f"  Supervivientes R2_ORI: {n_after_ori}")

        # ── 3. Ajustar FMM a los supervivientes (paralelo) ─────────
        self._log("  Ajustando FMM ...")

        surv_rows = np.array(
            [expr_norm.loc[g].values.astype(np.float64) for g in surv_genes]
        )

        fmm_results = self._fit_fmm_batch(surv_rows, circ)

        fmm_peaks_all = np.array([r[0] for r in fmm_results])
        fmm_omega_all = np.array([r[1] for r in fmm_results])
        fmm_r2_all    = np.array([r[2] for r in fmm_results])

        # ── 4. Filtro FMM: omega > min Y R2_FMM >= umbral ─────────────
        mask_fmm = (
            (fmm_omega_all > self.omega_min)
            & (fmm_r2_all >= self.r2_fmm_threshold)
        )
        n_after_fmm = int(mask_fmm.sum())

        self._log(
            f"  Filtro FMM (omega>{self.omega_min}, "
            f"R2>={self.r2_fmm_threshold}): {n_after_fmm} supervivientes"
        )

        # ── 4b. Cobertura circular y cupo por sector ───────────────────
        # Los tres criterios que definen el conjunto TOP son: R2_FMM alto,
        # perfil no spiked (omega > 0.1) y picos t_U que cubren de forma
        # aproximadamente uniforme la circunferencia. Los dos primeros son los
        # filtros anteriores; el tercero se comprueba aqui dividiendo el
        # circulo en ``n_sectors`` sectores con origen en la mediana circular
        # de los picos estimados. Si un sector queda vacio se avisa, pero no
        # se fuerza ningun gen que no supere ya los dos primeros criterios
        # (comprobado con datos reales: rellenar con R2<r2_fmm_threshold
        # mete en el TOP genes que el propio filtro ya habia rechazado). El
        # cupo (``top_fraction`` del total de genes) se reparte despues entre
        # los sectores que si tienen genes, tomando en cada uno los de mayor
        # R2_FMM.
        sel_idx, cov_idx = self._select_top_indices(
            fmm_peaks_all, fmm_omega_all, fmm_r2_all, mask_fmm, n_total,
        )
        coverage_names = [surv_genes[i] for i in cov_idx]

        # Acumular genes aceptados
        acc_names:   list[str]        = []
        acc_r2_ori:  list[float]      = []
        acc_r2_fmm:  list[float]      = []
        acc_omega:   list[float]      = []
        acc_peaks:   list[float]      = []
        acc_rows:    list[np.ndarray] = []

        for i in sel_idx:
            acc_names.append(surv_genes[i])
            acc_r2_ori.append(float(surv_r2_ori[i]))
            acc_r2_fmm.append(float(fmm_r2_all[i]))
            acc_omega.append(float(fmm_omega_all[i]))
            acc_peaks.append(float(fmm_peaks_all[i]))
            acc_rows.append(surv_rows[i])

        # ── 5. Forzar genes core (semilla) ────────────────────────────
        forced: list[str] = []
        accepted_set = set(acc_names)

        for seed in seed_genes:
            if seed in accepted_set:
                continue
            if seed not in expr_norm.index:
                self._log(
                    f"  AVISO: gen semilla {seed} ausente de la matriz."
                )
                continue

            row = expr_norm.loc[seed].values.astype(np.float64)
            if self.fmm_method == "circafmm":
                # Reusa el CircaFMMModel cacheado en _fit_fmm_batch
                pk, om, r2, *_ = _fit_fmm_one_circafmm(
                    row, circ, self._cached_circafmm_model,
                )
            else:
                pk, om, r2, *_ = _fit_fmm_one(
                    row, circ,
                    self.fmm_length_alpha_grid,
                    self.fmm_length_omega_grid,
                    self.fmm_num_reps,
                )
            idx_all = (
                gene_names.index(seed) if seed in gene_names else -1
            )
            seed_r2_ori = float(r2_ori[idx_all]) if idx_all >= 0 else 0.0

            acc_names.append(seed)
            acc_r2_ori.append(seed_r2_ori)
            acc_r2_fmm.append(float(r2))
            acc_omega.append(float(om))
            acc_peaks.append(float(pk))
            acc_rows.append(row)
            forced.append(seed)

        if forced:
            self._log(f"  Genes semilla forzados: {forced}")

        n_top = len(acc_names)
        self._log(f"  Genes TOP finales: {n_top}")

        # ── 6. Asignar cuadrantes ──────────────────────────────────────
        peaks_arr = np.array(acc_peaks)

        if n_top > 0:
            peak_ref = circular_median(peaks_arr)
            centres  = sector_centres(peak_ref, self.n_sectors)
            sectors  = assign_sectors(
                peaks_arr, peak_ref, self.n_sectors,
            )
        else:
            centres = np.zeros(self.n_sectors)
            sectors = np.array([], dtype=int)

        # Verificacion final. Tras _select_top_indices ningun sector deberia
        # estar vacio, salvo que no exista ningun gen ritmico con pico en el
        # (caso ya avisado alli). Los genes semilla forzados pueden desplazar
        # ligeramente la mediana circular de referencia, asi que se recomprueba.
        covered = set(sectors.tolist()) if n_top > 0 else set()
        expected = set(range(1, self.n_sectors + 1))
        missing = expected - covered
        if missing:
            self._log(
                f"  AVISO: cuadrantes sin cobertura: {sorted(missing)}. "
                "El conjunto TOP no cubre uniformemente la circunferencia, de "
                "modo que la estructura circular del CPCA posterior puede "
                "degradarse."
            )

        # ── 7. Construir matriz TOP ────────────────────────────────────
        if n_top > 0:
            top_matrix = pd.DataFrame(
                np.vstack(acc_rows),
                index=acc_names,
                columns=expr_norm.columns,
            )
        else:
            top_matrix = pd.DataFrame(
                dtype=float, columns=expr_norm.columns,
            )

        result = TopGeneResult(
            gene_names       = acc_names,
            fmm_peaks        = peaks_arr,
            r2_par           = np.array(acc_r2_fmm),
            sector_labels    = sectors,
            sector_centres   = centres,
            candidate_matrix = top_matrix,
            circular_scale   = circ,
            added_genes      = forced,
            sample_order     = np.asarray(sample_order, dtype=int),
            omega            = np.array(acc_omega),
            r2_ori           = np.array(acc_r2_ori),
            n_total_genes    = n_total,
            n_after_ori      = n_after_ori,
            n_after_fmm      = n_after_fmm,
        )

        self._log(result.summary())
        return result

    # ------------------------------------------------------------------
    # Cobertura circular de los picos y cuota por porcentaje
    # ------------------------------------------------------------------

    def _select_top_indices(
        self,
        peaks:    np.ndarray,
        omega:    np.ndarray,
        r2_fmm:   np.ndarray,
        mask_fmm: np.ndarray,
        n_total:  int,
    ) -> np.ndarray:
        """Indices del conjunto TOP tras imponer cobertura circular y cuota.

        Parte de los genes que superan los filtros FMM (``mask_fmm``) y aplica,
        en este orden:

        1. **Cobertura circular.** Se divide la circunferencia en
           ``n_sectors`` sectores tomando como origen la mediana circular de
           los picos ``t_U`` de los genes filtrados. Si algun sector queda
           vacio se avisa por consola, pero **no se rellena con ningun gen
           que no supere ya ``omega_min`` Y ``r2_fmm_threshold``** (es decir,
           que no pertenezca a ``mask_fmm``)

        2. **Cupo por sector.** Se retienen como mucho
           ``max(n_sectors, round(top_fraction * n_total))`` genes: el cupo se reparte a partes iguales
           entre los sectores no vacios y dentro de cada uno se ordena por
           ``R2_FMM`` descendente. Las plazas que un sector no llegue a cubrir
           (porque tiene menos genes que su cupo) se redistribuyen entre los
           que aun tienen candidatos, de forma que se alcance el cupo
           siempre que haya genes suficientes.

           El reparto por sectores, y no un ranking global, evita que el cupo
           deshaga la cobertura ya presente de forma natural: un ranking
           global por ``R2_FMM`` tiende a concentrarse en la region del ciclo
           donde la señal es mas fuerte.

        Returns
        -------
        np.ndarray de int : indices (sobre los supervivientes del filtro ORI)
        que forman el conjunto TOP, en orden creciente, y un array vacio
        (se mantiene por compatibilidad con el resto del pipeline, que ya
        no recibe genes de cobertura: nunca se fuerza ninguno).
        """
        sel = np.flatnonzero(mask_fmm)
        cov: np.ndarray = np.array([], dtype=int)
        if sel.size == 0:
            self._log("  AVISO: ningun gen supera los filtros FMM.")
            return sel, cov

        n_sec = self.n_sectors

        # Origen de los sectores: mediana circular de los picos filtrados
        peak_ref = circular_median(peaks[sel])
        sectors  = assign_sectors(peaks, peak_ref, n_sec)

        # ── 1. Cobertura: se registra si algun sector queda vacio, sin ────
        # forzar ningun gen que no supere ya omega Y R2_FMM (ver docstring).
        for s in range(1, n_sec + 1):
            if not np.any(sectors[sel] == s):
                self._log(
                    f"  AVISO: el sector {s} no contiene ningun gen con "
                    f"omega>{self.omega_min} y R2_FMM>={self.r2_fmm_threshold}; "
                    "queda sin cobertura."
                )

        # ── 2. Cupo absoluto, estratificado por sector ────────────────────
        cov_arr = cov
        target_size = max(n_sec, int(round(self.top_fraction * n_total)))
        if target_size >= sel.size:
            return sel, cov_arr

        chosen_arr, quota = apply_sector_quota(
            sel, sectors, r2_fmm, cov_arr, n_sec, target_size,
        )
        self._log(
            f"  Cupo {target_size} genes (top_fraction={self.top_fraction} "
            f"de {n_total}): "
            f"{chosen_arr.size}/{sel.size} genes retenidos "
            f"(por sector: { {s: int(quota[s]) for s in sorted(quota)} })"
        )
        return chosen_arr, cov_arr[np.isin(cov_arr, chosen_arr)]

    # ------------------------------------------------------------------
    # Ajuste ORI en lote (paralelo)
    # ------------------------------------------------------------------

    def _fit_ori_batch(
        self,
        rows: np.ndarray,
        circ: np.ndarray,
    ) -> np.ndarray:
        """Ajusta ORI a multiples genes en paralelo. Devuelve R2_ORI (n_genes,)."""
        n = rows.shape[0]

        if self.n_jobs == 1 or n <= 8:
            r2_vals = [
                _fit_ori_one(rows[i], circ) for i in range(n)
            ]
        else:
            from joblib import Parallel, delayed
            self._log(
                f"  Paralelizando {n} ajustes ORI en {self.n_jobs} jobs ..."
            )
            r2_vals = Parallel(n_jobs=self.n_jobs)(
                delayed(_fit_ori_one)(rows[i], circ) for i in range(n)
            )

        return np.array(r2_vals, dtype=np.float64)

    # ------------------------------------------------------------------
    # Ajuste FMM en lote (paralelo)
    # ------------------------------------------------------------------

    def _fit_fmm_batch(
        self,
        rows: np.ndarray,
        circ: np.ndarray,
    ) -> list[tuple[float, float, float, np.ndarray]]:
        """Ajusta FMM a multiples genes en paralelo.

        El backend de joblib se elige segun ``fmm_method``:

        - **base**: ``loky`` (multiproceso). Argumentos pequenos por gen
          (vector ``row`` + escalares); el coste de serializar es bajo y
          se gana paralelismo real saltando el GIL.

        - **circafmm**: ``threading``. Los workers comparten ``prepared``
          (~23 MB de matrices precomputadas) por memoria sin serializar.
          El fit interno usa BLAS/LAPACK y kernels ``@njit``, que liberan
          el GIL, asi que threading da paralelismo efectivo aqui.
        """
        n = rows.shape[0]
        from joblib import Parallel, delayed

        # ── Camino CircaFMM: prep + threading ────────────────────────
        if self.fmm_method == "circafmm":
            self._log(f"  Preparando CircaFMM (calibracion lambda)...")
            prepared = prepare_dataset(circ, **self.circafmm_kwargs)
            circafmm_model = CircaFMMModel(prepared=prepared, refine=self.fmm_refine)
            self._cached_circafmm_model = circafmm_model        # reutilizable para seeds forzados
            self._log(f"  lambda* = {prepared.lambda_star:.4f}")

            if self.n_jobs == 1 or n <= 8:
                return [_fit_fmm_one_circafmm(rows[i], circ, circafmm_model) for i in range(n)]

            self._log(f"  Paralelizando {n} ajustes FMM (CircaFMM, threads) en {self.n_jobs} jobs ...")
            return Parallel(n_jobs=self.n_jobs, backend="threading")(
                delayed(_fit_fmm_one_circafmm)(rows[i], circ, circafmm_model)
                for i in range(n)
            )

        # ── Camino base: loky (multiproceso, default) ────────────────────
        ag = self.fmm_length_alpha_grid
        og = self.fmm_length_omega_grid
        nr = self.fmm_num_reps

        if self.n_jobs == 1 or n <= 8:
            return [_fit_fmm_one(rows[i], circ, ag, og, nr) for i in range(n)]

        self._log(f"  Paralelizando {n} ajustes FMM en {self.n_jobs} jobs ...")
        return Parallel(n_jobs=self.n_jobs)(
            delayed(_fit_fmm_one)(rows[i], circ, ag, og, nr)
            for i in range(n)
        )

    # ------------------------------------------------------------------
    # Utilidad
    # ------------------------------------------------------------------

    def _log(self, message: str) -> None:
        if self.verbose:
            print(message, flush=True)
