r"""
circust/cpca.py
===============
Etapas 1.1 + 1.2: Ordenamiento circular (CPCA) y deteccion de outliers.

Que es CPCA
------------
PCA sobre los 12(o mas) genes core del reloj. Los loadings de PC1 y PC2 de cada
muestra se tratan como coordenadas (x, y) en el plano, el angulo(phi) se obtiene:

    phi = atan2(PC2, PC1)   en [0, 2pi)

Y su norma(LEi):

    norm12 = sqrt(PC1² + PC2²)

da un *ordenamiento circular* que refleja el ritmo circadiano subyacente.
Muestras cercanas al origen tienen norma pequena y son potenciales outliers.

Deteccion de outliers
----------------------
Procedimiento estandar (spline + porcentaje de genes): por cada gen core
se ajusta un spline cubico periodico de nodos fijos (ver
fitting/cyclic_spline.py), independiente del modelo FMM usado despues en
el pipeline. Un gen es "candidata" si su residuo estandarizado supera
spline_gene_threshold; la muestra se declara atipica si es candidata en
al menos spline_outlier_fraction de los genes core (no basta con un solo
gen, evita que una desviacion propia de un gen tire la muestra entera).
Valido a partir de n>=24 aprox. (ver memoria de proyecto
"fmm-sobreajuste-oculta-outliers-n-bajo").

Criterio radial CPCA (opcional, DESACTIVADO por defecto): LEi < tight_radius
(0.10). Su umbral fijo no escala bien con n (no es data-driven todavia,
ver memoria de proyecto), asi que no participa en la decision salvo que
se active explicitamente (use_radial_criterion=True). Queda preparado
para combinarse con el spline via AND/OR mas adelante, tal como se
plantea en la correspondencia con la tutora.

Si se detectan outliers, se eliminan, se renormaliza la matriz core y se
re-ejecuta CPCA para obtener el ordenamiento final limpio.

Posicion en el pipeline
-----------------------
    Preprocessor  ->  CPCA  ->  CircularSynchronizer  ->  ...
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from dataclasses import dataclass, field
from typing import Optional

from circust.core_genes import SEED_GENES_DEFAULT
from circust.fitting.fmm            import FMMModel
from circust.fitting.circafmm import CircaFMMModel, prepare_dataset
from circust.fitting.cyclic_spline import fit_cyclic_spline_residuals
from circust.fitting.rhythm_model import FitResult
from circust.preprocessing import normalise_matrix


# ===========================================================================
# Dataclass de resultado
# ===========================================================================

@dataclass
class CPCAResult:
    r"""
    Salida de :class:`CPCA` — incluye el ordenamiento circular final y la
    deteccion de outliers (Etapas 1.1 + 1.2 del pipeline).

    Campos del pipeline (usados por etapas posteriores)
    ---------------------------------------------------
    sample_order : np.ndarray de int, forma (n_muestras,)
        Indices que ordenan las muestras por fase circular final.

    circular_scale : np.ndarray de float, forma (n_muestras,)
        Valores de phi ordenados en [0, 2pi). Eje temporal para ajustes.

    core_genes_found : list[str]
        Genes centrales presentes en la matriz y usados para CPCA.

    variance_explained : np.ndarray de float, forma (3,)
        Fraccion de varianza explicada por PC1, PC2, PC3 (CPCA final).

    expr_norm_final : pd.DataFrame, forma (n_genes, n_muestras_limpias)
        Matriz de expresion completa con muestras outlier eliminadas,
        reordenada por fase circular y renormalizada.

    core_norm_final : pd.DataFrame, forma (n_genes_core, n_muestras_limpias)
        Submatriz core normalizada, sin outliers, ordenada por fase.

    fmm_fits_final : dict[str, FitResult]
        Ajustes FMM sobre genes core usando el ordenamiento final.

    samples_dropped : list[int]
        Indices originales de columna de las muestras eliminadas.

    spline_outliers : list[int]
        Subconjunto de samples_dropped detectados exclusivamente por el
        criterio de residuos del spline periodico + porcentaje de genes.

    Campos de diagnostico (para visualizacion)
    ------------------------------------------
    pc1, pc2, pc3 : np.ndarray
        Loadings de PC1/PC2/PC3 del CPCA final (un valor por muestra).

    std_residuals_spline : pd.DataFrame o None
        Residuos estandarizados del spline periodico por gen core (orden
        inicial). Es lo que realmente usa la deteccion de outliers.
    """

    # ── Pipeline ─────────────────────────────────────────────────────────
    sample_order:          np.ndarray
    circular_scale:        np.ndarray
    core_genes_found:      list[str]
    variance_explained:    np.ndarray
    expr_norm_final:       pd.DataFrame
    core_norm_final:       pd.DataFrame
    fmm_fits_final:        dict
    samples_dropped:       list[int]   = field(default_factory=list)
    spline_outliers:       list[int]   = field(default_factory=list)

    # ── Diagnostico ───────────────────────────────────────────────────────
    pc1:                   np.ndarray                = field(default_factory=lambda: np.array([]))
    pc2:                   np.ndarray                = field(default_factory=lambda: np.array([]))
    pc3:                   np.ndarray                = field(default_factory=lambda: np.array([]))
    std_residuals_spline:  Optional[pd.DataFrame]    = None
    pc1_initial:           np.ndarray                = field(default_factory=lambda: np.array([]))
    pc2_initial:           np.ndarray                = field(default_factory=lambda: np.array([]))

    # ── CircaFMM: lambda calibrado, propagable a etapas siguientes ──
    lambda_star:           Optional[float]           = None

    def summary(self) -> str:
        n_cpca = len(self.samples_dropped) - len(self.spline_outliers)
        lines = [
            "=== Resumen CPCA + Deteccion de Outliers ===",
            f"  Genes core usados           : {len(self.core_genes_found)}  {self.core_genes_found}",
            f"  Varianza PC1 / PC2 / PC3    : "
            f"{self.variance_explained[0]:.1%} / "
            f"{self.variance_explained[1]:.1%} / "
            f"{self.variance_explained[2]:.1%}",
            f"  Muestras eliminadas         : {len(self.samples_dropped)}",
            f"    Criterio CPCA radial      : {n_cpca}",
            f"    Criterio spline residuos  : {len(self.spline_outliers)}",
            f"  Muestras finales            : {self.expr_norm_final.shape[1]}",
        ]
        return "\n".join(lines)


# ===========================================================================
# Clase CPCA
# ===========================================================================

class CPCA:
    r"""
    Etapas 1.1 y 1.2: ordenamiento circular de muestras y deteccion de
    outliers residuales, segun la Seccion 3.3 del suplementario de CIRCUST.

    Parametros — CPCA
    -----------------
    core_genes : list[str], opcional
        Genes ancla circadianos. Por defecto los 12 del articulo CIRCUST.

    tight_radius : float
        Umbral radial primario. Muestras con LEi <= tight_radius se confirman
        como outliers CPCA. Por defecto: 0.10.

    loose_radius : float
        Umbral de respaldo cuando ninguna muestra califica con tight_radius.
        Por defecto: 0.15.

    use_radial_criterion : bool
        Si True, el criterio radial participa en la deteccion de outliers
        (union OR con el spline). Por defecto: False (desactivado): el
        umbral fijo tight_radius/loose_radius no escala bien con n y aun
        no tiene sustituto data-driven validado (ver memoria de proyecto).
        El procedimiento estandar del TFG usa solo el criterio de spline.

    Parametros — deteccion de outliers (spline)
    ---------------------------------------------
    spline_n_knots : int
        Numero de nodos del spline periodico por gen core. Por defecto: 5.

    spline_gene_threshold : float
        Umbral de residuo estandarizado del spline para que un gen
        declare a una muestra "candidata". Por defecto: 3.5.

    spline_outlier_fraction : float
        Fraccion minima de genes core que deben declarar "candidata" a
        una muestra para que se elimine definitivamente. Por defecto:
        0.25 (25% de los genes).

    max_outlier_fraction : float
        Limite maximo de muestras eliminables: ceil(fraccion x n_muestras).
        Por defecto: 0.05.

    Parametros — ajuste FMM
    ------------------------
    fmm_length_alpha_grid : int
        Resolucion de rejilla para el parametro alfa de FMM. Por defecto: 48.

    fmm_length_omega_grid : int
        Resolucion de rejilla para el parametro omega de FMM. Por defecto: 24.

    fmm_num_reps : int
        Iteraciones de refinamiento FMM. Por defecto: 3.

    fmm_method : str
        ``"base"`` (default) usa ``FMMModel`` (RSS + Nelder-Mead).
        ``"circafmm"`` usa ``CircaFMMModel`` (criterio de potencia
        penalizado por rugosidad, pesos por densidad temporal). Cuando es
        ``"circafmm"``, los parametros ``fmm_length_*`` y ``fmm_num_reps``
        se ignoran y se usan los defaults de CircaFMM o los pasados via
        ``circafmm_kwargs``.

    circafmm_kwargs : dict, opcional
        Parametros adicionales para ``prepare_dataset`` cuando
        ``fmm_method="circafmm"`` (ej: ``n_simulations``, ``n_repeats``).

    fmm_refine : bool
        Activa el refinamiento local opcional de ``(alpha, omega)`` de
        CircaFMM (Seccion S3.4.4, Algoritmo S2 del suplemento de
        CircaFMM). Por defecto False, como en ``runCircaFMM.R``. Se ignora con
        ``fmm_method="base"``.

    verbose : bool
        Imprimir mensajes de progreso.

    Example
    -------
    >>> from circust.preprocessing import load_expression_matrix, Preprocessor
    >>> from circust.cpca import CPCA
    >>>
    >>> matrix = load_expression_matrix("data/raw/expression.csv")
    >>> prep   = Preprocessor().run(matrix)
    >>> result = CPCA().run(prep.expr_norm)
    >>> print(result.summary())
    """

    def __init__(
        self,
        core_genes:               list[str] = None,
        # CPCA
        tight_radius:             float     = 0.10,
        loose_radius:             float     = 0.15,
        use_radial_criterion:     bool      = False,
        # Outlier detection
        spline_n_knots:           int       = 5,
        spline_gene_threshold:    float     = 3.5,
        spline_outlier_fraction:  float     = 0.25,
        max_outlier_fraction:     float     = 0.05,
        # FMM fitting
        fmm_length_alpha_grid:    int       = 48,
        fmm_length_omega_grid:    int       = 24,
        fmm_num_reps:             int       = 3,
        fmm_method:               str       = "base",
        circafmm_kwargs:          dict|None = None,
        fmm_refine:               bool      = False,
        compute_fmm_fits:         bool      = True,
        verbose:                  bool      = True,
    ) -> None:
        if fmm_method not in ("base", "circafmm"):
            raise ValueError(f"fmm_method debe ser 'base' o 'circafmm', recibido: {fmm_method!r}")
        self.core_genes               = core_genes if core_genes is not None else SEED_GENES_DEFAULT
        # compute_fmm_fits=False salta el ajuste FMM final de los genes core y
        # la deteccion de outliers. Solo valido cuando al llamante le basta con
        # sample_order/circular_scale (vienen de la descomposicion CPCA, no del
        # FMM) y no elimina outliers (max_outlier_fraction=0). Lo usa
        # RobustOrder por replica para evitar ~el 77% del coste, que se tiraba.
        self.compute_fmm_fits         = compute_fmm_fits
        self.tight_radius             = tight_radius
        self.loose_radius             = loose_radius
        # El criterio radial (LEi) esta desactivado por defecto: su umbral
        # fijo (tight_radius/loose_radius) no escala bien con n (ver memoria
        # de proyecto), y esta pendiente de un umbral data-driven que aun no
        # se ha validado. Se calcula igualmente (barato) para dejarlo listo
        # para combinarse AND/OR con el spline mas adelante, pero no
        # participa en la decision final salvo que se active explicitamente.
        self.use_radial_criterion     = use_radial_criterion
        self.spline_n_knots           = spline_n_knots
        self.spline_gene_threshold    = spline_gene_threshold
        self.spline_outlier_fraction  = spline_outlier_fraction
        self.max_outlier_fraction     = max_outlier_fraction
        self.fmm_length_alpha_grid    = fmm_length_alpha_grid
        self.fmm_length_omega_grid    = fmm_length_omega_grid
        self.fmm_num_reps             = fmm_num_reps
        self.fmm_method               = fmm_method
        self.circafmm_kwargs          = circafmm_kwargs or {}
        self.fmm_refine               = fmm_refine
        self.verbose                  = verbose

        # Lambda calibrado en _fit_final (primera llamada a _make_fmm_model) y
        # expuesto en CPCAResult.lambda_star para etapas posteriores del pipeline.
        self._cached_lambda_star: Optional[float] = None

    def _make_fmm_model(self, times: np.ndarray):
        """Factory interno: construye el modelo FMM segun ``fmm_method``.

        Para ``"circafmm"``, calibra lambda solo en la primera llamada
        (en ``_fit_final``) y lo cachea para exponerlo en ``CPCAResult``.
        Si el usuario ya proporciono ``lambda_star`` en ``circafmm_kwargs``,
        salta la calibracion automatica desde el inicio.
        """
        if self.fmm_method == "circafmm":
            circafmm_kwargs = dict(self.circafmm_kwargs)

            # Si ya hay un lambda cacheado (segunda+ llamada), reusarlo
            if self._cached_lambda_star is not None and "lambda_star" not in circafmm_kwargs:
                circafmm_kwargs["lambda_star"] = self._cached_lambda_star

            prepared = prepare_dataset(times, **circafmm_kwargs)

            # Cachear la primera vez (o si el usuario forzo un lambda concreto)
            if self._cached_lambda_star is None:
                self._cached_lambda_star = float(prepared.lambda_star)

            return CircaFMMModel(prepared=prepared, refine=self.fmm_refine)
        return FMMModel(
            length_alpha_grid = self.fmm_length_alpha_grid,
            length_omega_grid = self.fmm_length_omega_grid,
            num_reps          = self.fmm_num_reps,
        )

    # -----------------------------------------------------------------------
    # API publica
    # -----------------------------------------------------------------------

    def run(self, expr_norm: pd.DataFrame) -> CPCAResult:
        """
        Ejecuta CPCA + deteccion de outliers sobre la matriz completa.

        Parameters
        ----------
        expr_norm : pd.DataFrame, forma (n_genes, n_muestras)
            Matriz de expresion normalizada completa (valores en [-1, 1]).
            Es ``PreprocessingResult.expr_norm``.

        Returns
        -------
        CPCAResult
        """
        self._log("=== Etapa 1.1: CPCA ===")

        # ── Paso 1: extraer submatriz de genes core ───────────────────────
        core_matrix, genes_found = self._extract_core_genes(expr_norm)

        # ── Paso 2: CPCA inicial ──────────────────────────────────────────
        (sample_order, circular_scale,
         pc1, pc2, pc3, var_exp,
         outlier_idx) = self._cpca_step(core_matrix)

        self._log(
            f"  Varianza PC1/PC2/PC3: {var_exp[0]:.1%} / "
            f"{var_exp[1]:.1%} / {var_exp[2]:.1%}"
        )
        self._log(f"  Outliers radiales CPCA (LEi < {self.tight_radius}): {len(outlier_idx)}")

        # ── Paso 3: spline en orden inicial + deteccion de outliers ──────
        # Atajo: si no se piden los ajustes FMM (compute_fmm_fits=False), se
        # salta la deteccion de outliers + _fit_final. Solo es
        # valido cuando al llamante le basta sample_order/circular_scale y no
        # se eliminan outliers; sample_order/circular_scale no dependen del FMM.
        if not self.compute_fmm_fits:
            cpca_outs, spline_outs, all_dropped = [], [], []
            pc1_initial, pc2_initial = pc1.copy(), pc2.copy()
            core_norm_clean = core_matrix
            expr_norm_final = self._order_full_matrix(expr_norm, all_dropped, sample_order)
            result = CPCAResult(
                sample_order         = sample_order,
                circular_scale       = circular_scale,
                core_genes_found     = genes_found,
                variance_explained   = var_exp,
                expr_norm_final      = expr_norm_final,
                core_norm_final      = core_norm_clean.iloc[:, sample_order],
                fmm_fits_final       = {},
                samples_dropped      = list(all_dropped),
                spline_outliers      = spline_outs,
                pc1                  = pc1,
                pc2                  = pc2,
                pc3                  = pc3,
                std_residuals_spline = None,
                pc1_initial          = pc1_initial,
                pc2_initial          = pc2_initial,
                lambda_star          = self._cached_lambda_star,
            )
            return result

        self._log("=== Etapa 1.2: Deteccion de Outliers ===")
        self._log("  Ajustando spline periodico en genes core (deteccion, independiente de FMM) ...")
        std_res_spline_df = self._fit_spline_initial(
            core_matrix, genes_found, sample_order, circular_scale,
        )

        # ── Paso 4: detectar outliers (CPCA radial OR spline residuos) ────
        cpca_outs, spline_outs, all_dropped = self._detect_outliers(
            outlier_idx, sample_order, std_res_spline_df, genes_found,
        )
        self._log(f"  Criterio CPCA radial     : {len(cpca_outs)} muestra(s)")
        self._log(f"  Criterio spline residuos : {len(spline_outs)} muestra(s)")
        self._log(f"  Total eliminados         : {len(all_dropped)}")

        # Guardar proyecciones iniciales (antes de posible re-ejecucion)
        # para visualizacion del scatter con outliers
        pc1_initial = pc1.copy()
        pc2_initial = pc2.copy()

        # ── Paso 5: limpiar + re-ejecutar CPCA si hay outliers ───────────
        if all_dropped:
            self._log("  Re-ejecutando CPCA sobre datos limpios ...")
            core_clean = self._drop_samples(core_matrix, all_dropped)
            core_norm_clean = normalise_matrix(core_clean)
            (sample_order, circular_scale,
             pc1, pc2, pc3, var_exp, _) = self._cpca_step(core_norm_clean)
        else:
            self._log("  Sin outliers — manteniendo CPCA inicial.")
            core_norm_clean = core_matrix

        # ── Paso 6: ordenar y renormalizar la matriz completa ─────────────
        self._log("  Ordenando y renormalizando la matriz completa ...")
        expr_norm_final = self._order_full_matrix(expr_norm, all_dropped, sample_order)

        # ── Paso 7: ajuste FMM final sobre genes core ─────────────────────
        self._log("  Ajustando FMM en genes core (orden final) ...")
        fmm_fits_fin = self._fit_final(
            expr_norm_final, genes_found, circular_scale,
        )

        result = CPCAResult(
            sample_order         = sample_order,
            circular_scale       = circular_scale,
            core_genes_found     = genes_found,
            variance_explained   = var_exp,
            expr_norm_final      = expr_norm_final,
            core_norm_final      = core_norm_clean.iloc[:, sample_order],
            fmm_fits_final       = fmm_fits_fin,
            samples_dropped      = list(all_dropped),
            spline_outliers      = spline_outs,
            # diagnostico
            pc1                  = pc1,
            pc2                  = pc2,
            pc3                  = pc3,
            std_residuals_spline = std_res_spline_df,
            pc1_initial          = pc1_initial,
            pc2_initial          = pc2_initial,
            lambda_star          = self._cached_lambda_star,
        )

        self._log(result.summary())
        return result

    # -----------------------------------------------------------------------
    # Pasos privados
    # -----------------------------------------------------------------------

    def _extract_core_genes(
        self, expr_norm: pd.DataFrame,
    ) -> tuple[pd.DataFrame, list[str]]:
        """Selecciona las filas de genes core de la matriz normalizada."""
        genes_found, missing = [], []
        for gene in self.core_genes:
            (genes_found if gene in expr_norm.index else missing).append(gene)

        if missing:
            self._log(
                f"  AVISO: {len(missing)} gen(es) core no encontrado(s): {missing}"
            )
        if len(genes_found) < 2:
            raise ValueError(
                f"CPCA requiere al menos 2 genes core. "
                f"Solo se encontraron: {genes_found}."
            )

        self._log(
            f"  Genes core extraidos: {len(genes_found)} / {len(self.core_genes)}"
        )
        return expr_norm.loc[genes_found], genes_found

    def _cpca_step(
        self, core_matrix: pd.DataFrame,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray,
               np.ndarray, np.ndarray]:
        """
        Ejecuta un paso de CPCA sobre la matriz core normalizada.

        Returns
        -------
        (sample_order, circular_scale, pc1, pc2, pc3, var_exp, outlier_idx)
        """
        values = core_matrix.values.astype(float)
        n_genes, n_samples = values.shape

        # Centrar por filas, escalar por columnas (R: prcomp scale.=TRUE center=FALSE)
        centred = values - np.nanmean(values, axis=1, keepdims=True)
        col_rms = np.sqrt(
            np.nansum(centred**2,axis=0) / (n_genes - 1)
        )
        col_rms[col_rms == 0] = 1.0
        scaled   = centred / col_rms

        n_comp = min(3, n_genes, n_samples)
        _, sigma_all, Vt = np.linalg.svd(scaled, full_matrices=False)

        pc1 = Vt[0]
        pc2 = Vt[1]
        pc3 = Vt[2] if n_comp >= 3 else np.zeros(n_samples)

        var_exp = np.zeros(3)
        var_exp[:n_comp] = (sigma_all[:n_comp]**2) / np.sum(sigma_all**2)

        norm12    = np.sqrt(pc1**2 + pc2**2)
        safe_norm = np.where(norm12 == 0, 1.0, norm12)
        phi       = np.arctan2(pc2 / safe_norm, pc1 / safe_norm) % (2 * np.pi)

        sample_order   = np.argsort(phi)
        circular_scale = phi[sample_order]

        # Outliers radiales: LEi <= tight_radius sobre TODAS las muestras
        tight_mask = norm12 <= self.tight_radius

        if tight_mask.any():
            outlier_idx = np.where(tight_mask)[0]
        else:
            loose_mask  = norm12 <= self.loose_radius
            outlier_idx = np.where(loose_mask)[0]
        # De menor a mayor norma: si hay mas candidatas que el tope de
        # _detect_outliers, se quedan las mas proximas al origen y no las de
        # menor indice de columna.
        outlier_idx = outlier_idx[np.argsort(norm12[outlier_idx], kind="stable")]

        return sample_order, circular_scale, pc1, pc2, pc3, var_exp, outlier_idx

    def _fit_spline_initial(
        self,
        core_matrix:    pd.DataFrame,
        genes_found:    list[str],
        sample_order:   np.ndarray,
        circular_scale: np.ndarray,
    ) -> pd.DataFrame:
        """
        Ajusta el spline periodico (independiente de FMM) en cada gen
        core usando el ordenamiento CPCA inicial. Solo genes core, no
        eigengenes: esto alimenta
        exclusivamente la deteccion de outliers.

        Returns
        -------
        std_res_df : DataFrame (genes x muestras) de residuos
        estandarizados del spline.
        """
        cm_vals = core_matrix.values
        std_res_rows: list[np.ndarray] = []

        for i, gene in enumerate(genes_found):
            data = cm_vals[i, sample_order]
            std_res_rows.append(
                fit_cyclic_spline_residuals(data, circular_scale, n_knots=self.spline_n_knots)
            )

        return pd.DataFrame(
            np.vstack(std_res_rows),
            index   = genes_found,
            columns = np.arange(len(circular_scale)),
        )

    def _detect_outliers(
        self,
        outlier_idx:        np.ndarray,
        sample_order:       np.ndarray,
        std_res_spline_df:  pd.DataFrame,
        genes_found:        list[str],
    ) -> tuple[list[int], list[int], list[int]]:
        r"""
        Aplica criterios de outlier en logica OR:

          1. CPCA radial (LEi < tight_radius): ya en outlier_idx. DESACTIVADO
             por defecto (use_radial_criterion=False); su umbral fijo no
             escala bien con n y aun no tiene sustituto data-driven
             validado. Se deja preparado para activarse/combinarse mas
             adelante sin cambiar la firma de este metodo.
          2. Spline periodico + porcentaje de genes (procedimiento estandar
             del TFG): un gen es "candidata" si su residuo estandarizado
             supera spline_gene_threshold; la muestra se declara atipica si
             es candidata en >= spline_outlier_fraction de los genes core.

        Si hay mas candidatas que el tope ``max_outlier_fraction``, se
        eliminan las mas atipicas: las radiales de menor norma y, despues, las
        del spline con mayor fraccion de genes que las marcan (a igualdad,
        mayor residuo medio absoluto).

        Returns
        -------
        cpca_outs   : eliminados por criterio radial (vacio si esta desactivado)
        spline_outs : eliminados solo por criterio spline
        all_dropped : union ordenada (con tope max_outlier_fraction)
        """
        n_samples = len(sample_order)
        cap       = int(np.ceil(self.max_outlier_fraction * n_samples))

        # Criterio 1: CPCA radial (opcional, desactivado por defecto)
        active_outlier_idx = outlier_idx if self.use_radial_criterion else np.array([], dtype=int)
        cpca_set = set(active_outlier_idx.tolist())

        # Criterio 2: spline periodico + porcentaje de genes. Valido a
        # partir de n>=24 aprox (ver memoria "fmm-sobreajuste-oculta-outliers-n-bajo",
        # el mismo problema de sobreajuste con pocos puntos afecta tambien
        # al spline, aunque con margen algo mayor que FMM).
        if n_samples < 24:
            self._log(
                f"  AVISO: n={n_samples} < 24, el criterio de residuos del "
                f"spline puede no ser fiable (ver limitacion conocida en "
                f"la memoria del TFG)."
            )

        core_res      = std_res_spline_df.loc[genes_found].values  # (n_genes, n_samples)
        per_gene_flag = np.abs(core_res) > self.spline_gene_threshold
        frac_flagged  = per_gene_flag.mean(axis=0)
        mean_abs_res  = np.abs(core_res).mean(axis=0)

        # Candidatas de mas a menos atipicas: primero la fraccion de genes que
        # las marcan y, a igualdad, el residuo medio absoluto en los genes core.
        # Antes se recorrian en el orden circular, de modo que, si habia mas
        # candidatas que el tope, se eliminaban las mas proximas al origen del
        # circulo y no las mas atipicas.
        pos_cands = np.where(frac_flagged >= self.spline_outlier_fraction)[0]
        pos_cands = pos_cands[np.lexsort((-mean_abs_res[pos_cands],
                                          -frac_flagged[pos_cands]))]
        spline_cands: list[int] = [int(sample_order[pos]) for pos in pos_cands]

        # Union OR con tope (chequea ANTES de agregar para que cap=0 elimine 0)
        seen:       set[int]  = set()
        all_dropped: list[int] = []
        for idx in list(active_outlier_idx) + spline_cands:
            if len(all_dropped) >= cap:
                break
            if idx not in seen:
                seen.add(idx)
                all_dropped.append(idx)

        cpca_outs   = [i for i in all_dropped if i in cpca_set]
        spline_outs = [i for i in all_dropped if i not in cpca_set]

        return cpca_outs, spline_outs, all_dropped

    def _drop_samples(
        self, matrix: pd.DataFrame, dropped: list[int],
    ) -> pd.DataFrame:
        """Elimina columnas por indice original de muestra."""
        all_cols  = np.arange(matrix.shape[1])
        keep_mask = ~np.isin(all_cols, dropped)
        return matrix.iloc[:, keep_mask]

    def _order_full_matrix(
        self,
        expr_norm:    pd.DataFrame,
        dropped:      list[int],
        sample_order: np.ndarray,
    ) -> pd.DataFrame:
        """
        Elimina outliers de la matriz completa, reordena por fase circular
        final y renormaliza gen a gen.
        """
        mat_clean   = self._drop_samples(expr_norm, dropped)
        mat_ordered = mat_clean.iloc[:, sample_order]
        return normalise_matrix(mat_ordered)

    def _fit_final(
        self,
        expr_norm_final: pd.DataFrame,
        genes_found:     list[str],
        circular_scale:  np.ndarray,
    ) -> dict[str, FitResult]:
        """
        Ajusta FMM en cada gen core con el ordenamiento final.

        Returns
        -------
        fmm_fits_final : {gen: FitResult}
        """
        fmm_model = self._make_fmm_model(circular_scale)

        fmm_fits: dict[str, FitResult] = {}

        n_total = len(genes_found)
        for idx, gene in enumerate(genes_found, 1):
            self._log(f"    [{idx}/{n_total}] {gene}", end="\r")
            if gene not in expr_norm_final.index:
                self._log(f"    AVISO: {gene} no en la matriz final, omitido.")
                continue
            data = expr_norm_final.loc[gene].values
            fr   = fmm_model.fit(data, circular_scale)
            fmm_fits[gene] = fr

        self._log(f"    Ajuste final completado: {len(fmm_fits)} genes core.  ")
        return fmm_fits

    # -----------------------------------------------------------------------
    # Utilidades
    # -----------------------------------------------------------------------

    def _log(self, message: str, end: str = "\n") -> None:
        if self.verbose:
            print(message, end=end, flush=True)
