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
     (ii) altamente ritmicos (R2_FMM > 0.5),
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
from circust.fitting.fmm_stabilized import StabilizedFMMModel, prepare_dataset


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
) -> tuple[float, float, float, float]:
    """Ajusta FMM base a un unico gen. Devuelve ``(peak_time, omega, r2, A)``.

    Funcion pura a nivel de modulo para que ``joblib`` pueda serializarla
    hacia procesos worker (FMM base puede paralelizar con muchos genes).

    Devuelve escalares y no el vector ``fitted``: este ultimo no lo consume
    nadie y, con decenas de miles de genes, obliga a serializar de vuelta un
    vector de n muestras por gen desde cada proceso worker. La amplitud, en
    cambio, es necesaria para construir el atlas.
    """
    model = FMMModel(
        length_alpha_grid=alpha_grid,
        length_omega_grid=omega_grid,
        num_reps=num_reps,
    )
    fr = model.fit(row, circ)
    return fr.peak_time, fr.params["omega"], fr.r2, fr.params["A"]


def _fit_fmm_one_stab(
    row:      np.ndarray,
    circ:     np.ndarray,
    stab_model: StabilizedFMMModel,
) -> tuple[float, float, float, float]:
    """Ajusta FMM estabilizado reusando un ``StabilizedFMMModel`` ya preparado.

    Mismo formato de salida que :func:`_fit_fmm_one`.
    """
    fr = stab_model.fit(row, circ)
    return fr.peak_time, fr.params["omega"], fr.r2, fr.params["A"]


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

    coverage_genes : list[str]
        Genes incorporados unicamente para que ningun sector circular quedase
        vacio, con el umbral de R2_FMM relajado. Se excluyen tambien del pool
        de muestreo aleatorio: al no alcanzar el R2 minimo harian fallar el
        criterio de calidad de la etapa 5 en cuanto entrasen en una extraccion.

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
    coverage_genes:   list          = field(default_factory=list)
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
    4. Filtrar: omega > ``omega_min`` (0.1) Y R2_FMM > ``r2_fmm_threshold`` (0.5).
    5. Asignar picos FMM a ``n_sectors`` cuadrantes y verificar cobertura.
    6. Forzar inclusion de los genes core si no estan presentes.
    7. Construir la matriz TOP ([X_TOP] en el paper).

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
        con origen en la mediana circular de los picos, para imponer la
        cobertura uniforme. Paper: 4 (cuadrantes). Ningun sector puede quedar
        vacio: si lo esta, se rellena con el gen de mayor R2_FMM de ese sector
        entre los que cumplen omega > omega_min.

    top_fraction : float | None
        Porcentaje prefijado de genes que forman el conjunto TOP, expresado
        como fraccion en (0, 1]. Se seleccionan los de mayor R2_FMM repartidos
        entre los sectores para preservar la cobertura circular. Si es None
        (defecto) se conservan todos los genes que superan los filtros, sin
        aplicar cuota.

    fmm_length_alpha_grid, fmm_length_omega_grid, fmm_num_reps : int
        Parametros de la rejilla de ajuste FMM.

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
        n_sectors:             int   = 4,
        top_fraction:          Optional[float] = None,
        fmm_length_alpha_grid: int   = 48,
        fmm_length_omega_grid: int   = 24,
        fmm_num_reps:          int   = 3,
        fmm_method:            str   = "base",
        fmm_stab_kwargs:       dict|None = None,
        lambda_star:           Optional[float] = None,
        n_jobs:                int   = -1,
        verbose:               bool  = True,
    ) -> None:
        if fmm_method not in ("base", "stabilized"):
            raise ValueError(f"fmm_method debe ser 'base' o 'stabilized', recibido: {fmm_method!r}")
        self.r2_ori_threshold      = r2_ori_threshold
        self.r2_fmm_threshold      = r2_fmm_threshold
        self.omega_min             = omega_min
        self.n_sectors             = n_sectors
        if top_fraction is not None and not (0.0 < top_fraction <= 1.0):
            raise ValueError(
                f"top_fraction debe estar en (0, 1], recibido: {top_fraction!r}"
            )
        self.top_fraction          = top_fraction
        self.fmm_length_alpha_grid = fmm_length_alpha_grid
        self.fmm_length_omega_grid = fmm_length_omega_grid
        self.fmm_num_reps          = fmm_num_reps
        self.fmm_method            = fmm_method
        # Si se proporciona lambda_star, lo inyectamos en stab_kwargs para
        # saltar la calibracion (es lo que hace el pipeline tras heredarlo de CPCA)
        self.fmm_stab_kwargs       = dict(fmm_stab_kwargs or {})
        if lambda_star is not None and "lambda_star" not in self.fmm_stab_kwargs:
            self.fmm_stab_kwargs["lambda_star"] = float(lambda_star)
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
        seed_genes = list(seed_genes)
        self._log("=== Etapa 3: Seleccion de Genes TOP ===")

        gene_names = list(expr_norm.index)
        circ       = np.asarray(circular_scale, dtype=np.float64)
        n_total    = len(gene_names)
        values     = expr_norm.values.astype(np.float64)

        if sample_order is None:
            sample_order = np.arange(expr_norm.shape[1], dtype=int)

        # ── 1. Ajuste ORI (regresion isotonica circular unimodal) ──────
        self._log(f"  Genes totales: {n_total}")
        self._log("  Ajustando ORI (PAVA unimodal) ...")

        r2_ori = self._fit_ori_batch(values, circ)

        # ── 2. Filtro R2_ORI ────────────────────────────────────────────
        self._log(
            f"  Filtro R2_ORI >= {self.r2_ori_threshold} ..."
        )

        mask_ori     = r2_ori >= self.r2_ori_threshold      # criterio: R2_ORI >= tau_ORI
        surv_genes   = [g for g, m in zip(gene_names, mask_ori) if m]
        surv_r2_ori  = r2_ori[mask_ori]
        n_after_ori  = len(surv_genes)

        self._log(f"  Supervivientes R2_ORI: {n_after_ori}")

        # ── 3. Ajustar FMM a los supervivientes (paralelo) ────────────
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

        # ── 4b. Cobertura circular y cuota por porcentaje ──────────────
        # Los tres criterios que definen el conjunto TOP son: R2_FMM alto,
        # perfil no spiked (omega > 0.1) y picos t_U que cubren de forma
        # aproximadamente uniforme la circunferencia. Los dos primeros son los
        # filtros anteriores; el tercero se impone aqui dividiendo el circulo en
        # ``n_sectors`` sectores con origen en la mediana circular de los picos
        # estimados y exigiendo que ninguno quede vacio. La cuota
        # ``top_fraction`` selecciona despues el porcentaje prefijado de genes
        # con mayor R2_FMM, repartido entre sectores para no romper la
        # cobertura recien impuesta.
        sel_idx, cov_idx = self._select_top_indices(
            fmm_peaks_all, fmm_omega_all, fmm_r2_all, mask_fmm,
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
            if self.fmm_method == "stabilized":
                # Reusa el StabilizedFMMModel cacheado en _fit_fmm_batch
                pk, om, r2, _ = _fit_fmm_one_stab(
                    row, circ, self._cached_stab_model,
                )
            else:
                pk, om, r2, _ = _fit_fmm_one(
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
            coverage_genes   = coverage_names,
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
    ) -> np.ndarray:
        """Indices del conjunto TOP tras imponer cobertura circular y cuota.

        Parte de los genes que superan los filtros FMM (``mask_fmm``) y aplica,
        en este orden:

        1. **Cobertura circular.** Se divide la circunferencia en
           ``n_sectors`` sectores tomando como origen la mediana circular de
           los picos ``t_U`` de los genes filtrados. Si algun sector queda
           vacio, se rellena con el gen de mayor ``R2_FMM`` cuyo pico caiga en
           el, entre los que cumplen ``omega > omega_min``. Es decir, se relaja
           unicamente el umbral de ``R2_FMM``, nunca el criterio de perfil no
           spiked. Un sector solo puede quedar vacio si ningun gen ritmico tiene
           su pico en el, situacion que se avisa por consola.

        2. **Cuota por porcentaje.** Si ``top_fraction`` no es None, se retiene
           esa fraccion de los genes, repartida proporcionalmente entre los
           sectores no vacios y ordenando por ``R2_FMM`` descendente dentro de
           cada uno. El reparto por sectores, y no global, evita que la cuota
           deshaga la cobertura que acaba de imponerse: un ranking global por
           ``R2_FMM`` tiende a concentrarse en la region del ciclo donde la
           señal es mas fuerte.

        Returns
        -------
        np.ndarray de int : indices (sobre los supervivientes del filtro ORI)
        que forman el conjunto TOP, en orden creciente, y los indices de los
        genes que han entrado unicamente por la regla de cobertura, que se
        excluyen despues del muestreo aleatorio de la etapa 5 por haber
        entrado con el umbral de R2 relajado.
        """
        sel = np.flatnonzero(mask_fmm)
        cov: list[int] = []
        if sel.size == 0:
            self._log("  AVISO: ningun gen supera los filtros FMM.")
            return sel, np.array([], dtype=int)

        n_sec = self.n_sectors

        # Origen de los sectores: mediana circular de los picos filtrados
        peak_ref = circular_median(peaks[sel])
        sectors  = assign_sectors(peaks, peak_ref, n_sec)

        # ── 1. Cobertura: ningun sector puede quedar vacio ────────────────
        keep = set(sel.tolist())
        pool = np.flatnonzero(omega > self.omega_min)     # se relaja solo R2_FMM
        for s in range(1, n_sec + 1):
            if np.any(sectors[sel] == s):
                continue
            cand = pool[sectors[pool] == s]
            if cand.size == 0:
                self._log(
                    f"  AVISO: el sector {s} no contiene ningun gen con "
                    f"omega>{self.omega_min}; no puede cubrirse."
                )
                continue
            best = int(cand[int(np.argmax(r2_fmm[cand]))])
            keep.add(best)
            cov.append(best)
            self._log(
                f"  Cobertura: sector {s} vacio, se incorpora el gen de mayor "
                f"R2_FMM del sector (R2={r2_fmm[best]:.3f}, "
                f"omega={omega[best]:.3f})"
            )
        sel = np.array(sorted(keep), dtype=int)

        # ── 2. Cuota por porcentaje, estratificada por sector ─────────────
        cov_arr = np.array(sorted(cov), dtype=int)
        if self.top_fraction is None:
            return sel, cov_arr

        n_target = int(np.ceil(self.top_fraction * sel.size))
        if n_target >= sel.size:
            return sel, cov_arr

        by_sector = {
            s: sel[sectors[sel] == s] for s in range(1, n_sec + 1)
        }
        by_sector = {s: idx for s, idx in by_sector.items() if idx.size > 0}

        # Reparto proporcional al tamano de cada sector, con al menos un gen
        # por sector no vacio para no perder la cobertura.
        n_sec_used = len(by_sector)
        quota = {s: max(1, int(np.floor(n_target * idx.size / sel.size)))
                 for s, idx in by_sector.items()}
        # Ajuste de redondeo: reparte o retira las plazas sobrantes por sector
        # de mayor a menor tamano, sin bajar de una plaza.
        orden = sorted(by_sector, key=lambda s: -by_sector[s].size)
        while sum(quota.values()) != n_target and n_target >= n_sec_used:
            delta = 1 if sum(quota.values()) < n_target else -1
            movido = False
            for s in orden:
                if delta > 0 and quota[s] < by_sector[s].size:
                    quota[s] += 1; movido = True
                elif delta < 0 and quota[s] > 1:
                    quota[s] -= 1; movido = True
                if movido and sum(quota.values()) == n_target:
                    break
            if not movido:
                break

        chosen: list[int] = []
        for s, idx in by_sector.items():
            top = idx[np.argsort(-r2_fmm[idx], kind="stable")][: quota[s]]
            chosen.extend(top.tolist())

        self._log(
            f"  Cuota top_fraction={self.top_fraction:g}: "
            f"{len(chosen)}/{sel.size} genes retenidos "
            f"(por sector: { {s: int(quota[s]) for s in sorted(quota)} })"
        )
        chosen_arr = np.array(sorted(chosen), dtype=int)
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

        - **stabilized**: ``threading``. Los workers comparten ``prepared``
          (~23 MB de matrices precomputadas) por memoria sin serializar.
          El fit interno usa BLAS/LAPACK y kernels ``@njit``, que liberan
          el GIL, asi que threading da paralelismo efectivo aqui.
        """
        n = rows.shape[0]
        from joblib import Parallel, delayed

        # ── Camino estabilizado: prep + threading ────────────────────────
        if self.fmm_method == "stabilized":
            self._log(f"  Preparando FMM estabilizado (calibracion lambda)...")
            prepared = prepare_dataset(circ, **self.fmm_stab_kwargs)
            stab_model = StabilizedFMMModel(prepared=prepared)
            self._cached_stab_model = stab_model        # reutilizable para seeds forzados
            self._log(f"  lambda* = {prepared.lambda_star:.4f}")

            if self.n_jobs == 1 or n <= 8:
                return [_fit_fmm_one_stab(rows[i], circ, stab_model) for i in range(n)]

            self._log(f"  Paralelizando {n} ajustes FMM (stab, threads) en {self.n_jobs} jobs ...")
            return Parallel(n_jobs=self.n_jobs, backend="threading")(
                delayed(_fit_fmm_one_stab)(rows[i], circ, stab_model)
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
