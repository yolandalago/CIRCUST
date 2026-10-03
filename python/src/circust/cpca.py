r"""
circust/cpca.py
========================
Etapa 2.1 de CIRCUST: CPCA sobre los genes core y deteccion de muestras
atipicas, tal y como la describe el material suplementario (Seccion 3.3) y
la implementacion original (``obtainCPCA13`` y ``giveMatIniNP_v3_cores``).

CPCA
----
PCA sobre la submatriz de genes core. Los loadings de los dos primeros
eigengenes de cada muestra son coordenadas (x, y) en el plano:

    phi_i = atan2(y_i, x_i)  en [0, 2pi)        d_i = sqrt(x_i^2 + y_i^2)

El orden creciente de phi es el orden circular preliminar y d_i es la
distancia de la muestra al origen (``LEi`` en el suplementario).

Deteccion de outliers
---------------------
Dos criterios que se combinan:

1. Radial. Se ordenan las muestras por d_i y se toman las
   ``ceil(0.05 n)`` mas cercanas al origen. Entre ellas son candidatas CPCA
   las que cumplen d_i <= ``tight_radius`` (0.10). Si ninguna lo cumple, se
   usa ``loose_radius`` (0.15).

2. Residuos del FMM. Se ajusta el FMM a cada gen core en el orden CPCA y se
   estandarizan los residuos, Z = (r - media) / sd. Una muestra se elimina
   si en algun gen core |Z| > 4, o si |Z| > 3 y ademas es candidata CPCA.

Las muestras eliminadas se quitan de toda la matriz, se renormaliza gen a
gen y se repite el CPCA una vez.

Tope. El proceso recorre los genes core en orden y solo atiende a un gen
mientras el numero de eliminaciones acumuladas no supere ``ceil(0.05 n)``.
La comprobacion se hace antes de procesar cada gen, como en el original, de
modo que un gen puede aportar varias muestras y el total puede rebasar el
tope en esas ultimas. Cada muestra marcada por |Z| > 3 con candidatura suma
1 y cada gen con alguna muestra con |Z| > 4 suma 1.

El numero de muestras que se examinan por el criterio radial es
``ceil(0.05 n)``. El codigo R original fija ese numero en 8.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from circust.fitting.fmm import FMMModel
from circust.fitting.rhythm_model import FitResult
from circust.preprocessing import normalise_matrix


@dataclass
class CPCAResult:
    """Salida de la Etapa 2.1.

    Attributes
    ----------
    sample_order : indices de columna de ``expr_norm`` (sin outliers) en el
        orden CPCA final.
    circular_scale : fases phi ordenadas, en [0, 2pi).
    core_genes_found : genes core presentes en la matriz.
    expr_norm_final : matriz completa sin outliers, ordenada y renormalizada.
    fmm_fits_final : ajuste FMM de cada gen core en el orden final.
    samples_dropped : indices (de la matriz de entrada) de las muestras eliminadas.
    samples_dropped_names : sus nombres.
    radial_candidates : indices de las candidatas CPCA (criterio radial).
    outliers_multi, outliers_uni : eliminadas por |Z|>3 con candidatura y por |Z|>4.
    variance_explained : varianza de PC1, PC2 y PC3.
    """

    sample_order:          np.ndarray
    circular_scale:        np.ndarray
    core_genes_found:      list[str]
    expr_norm_final:       pd.DataFrame
    fmm_fits_final:        dict[str, FitResult]
    samples_dropped:       list[int]
    samples_dropped_names: list[str]
    radial_candidates:     list[int]
    outliers_multi:        list[int]
    outliers_uni:          list[int]
    variance_explained:    np.ndarray
    n_radial_examined:     int = 0
    pc1:                   np.ndarray = field(default_factory=lambda: np.array([]))
    pc2:                   np.ndarray = field(default_factory=lambda: np.array([]))
    radial_distance:       np.ndarray = field(default_factory=lambda: np.array([]))

    def summary(self) -> str:
        return "\n".join([
            "=== Resumen CPCA + deteccion de outliers ===",
            f"  Genes core usados          : {len(self.core_genes_found)}",
            f"  Muestras examinadas (radial): {self.n_radial_examined}",
            f"  Candidatas CPCA            : {len(self.radial_candidates)}",
            f"  Eliminadas |Z|>3 + candidata: {len(self.outliers_multi)}",
            f"  Eliminadas |Z|>4           : {len(self.outliers_uni)}",
            f"  Total eliminadas           : {len(self.samples_dropped)}",
            f"  Varianza PC1/PC2/PC3       : {self.variance_explained[0]:.1%} / "
            f"{self.variance_explained[1]:.1%} / {self.variance_explained[2]:.1%}",
        ])


def cpca_step(values: np.ndarray) -> tuple:
    """CPCA de una matriz genes x muestras ya normalizada.

    Equivale a ``prcomp(centrado(m), scale.=TRUE, center=FALSE)`` de R.

    Returns
    -------
    (sample_order, circular_scale, pc1, pc2, pc3, var_exp, d)
    """
    values = np.asarray(values, dtype=float)
    n_genes, n_samples = values.shape
    centred = values - np.nanmean(values, axis=1, keepdims=True)
    col_rms = np.sqrt(np.nansum(centred ** 2, axis=0) / (n_genes - 1))
    col_rms[col_rms == 0] = 1.0
    scaled = centred / col_rms

    n_comp = min(3, n_genes, n_samples)
    _, sigma, vt = np.linalg.svd(scaled, full_matrices=False)
    pc1, pc2 = vt[0], vt[1]
    pc3 = vt[2] if n_comp >= 3 else np.zeros(n_samples)

    var_exp = np.zeros(3)
    var_exp[:n_comp] = sigma[:n_comp] ** 2 / np.sum(sigma ** 2)

    d = np.sqrt(pc1 ** 2 + pc2 ** 2)
    safe = np.where(d == 0, 1.0, d)
    phi = np.arctan2(pc2 / safe, pc1 / safe) % (2 * np.pi)
    order = np.argsort(phi, kind="stable")
    return order, phi[order], pc1, pc2, pc3, var_exp, d


def standardized_residuals(data: np.ndarray, fitted: np.ndarray) -> np.ndarray:
    """(r - media(r)) / sd(r), con sd muestral como ``sd`` de R."""
    r = np.asarray(data) - np.asarray(fitted)
    sd = r.std(ddof=1)
    return (r - r.mean()) / sd if sd > 0 else np.zeros_like(r)


class CPCA:
    """Etapa 2.1: CPCA, deteccion de outliers y orden preliminar.

    Parameters
    ----------
    core_genes : lista de genes core.
    tight_radius, loose_radius : umbrales de d_i para las candidatas (0.10, 0.15).
    z_multi, z_uni : umbrales de |Z| (3 con candidatura, 4 sin ella).
    max_outlier_fraction : fraccion de muestras que se examinan por el criterio
        radial y tope de eliminaciones, ``ceil(fraccion * n)``. Por defecto 0.05.
    fmm_length_alpha_grid, fmm_length_omega_grid, fmm_num_reps : rejilla FMM.
    """

    def __init__(
        self,
        core_genes:            list[str],
        tight_radius:          float = 0.10,
        loose_radius:          float = 0.15,
        z_multi:               float = 3.0,
        z_uni:                 float = 4.0,
        max_outlier_fraction:  float = 0.05,
        fmm_length_alpha_grid: int   = 48,
        fmm_length_omega_grid: int   = 24,
        fmm_num_reps:          int   = 3,
        verbose:               bool  = True,
    ) -> None:
        self.core_genes            = list(core_genes)
        self.tight_radius          = tight_radius
        self.loose_radius          = loose_radius
        self.z_multi               = z_multi
        self.z_uni                 = z_uni
        self.max_outlier_fraction  = max_outlier_fraction
        self.fmm_alpha             = fmm_length_alpha_grid
        self.fmm_omega             = fmm_length_omega_grid
        self.fmm_reps              = fmm_num_reps
        self.verbose               = verbose

    # ------------------------------------------------------------------
    def run(self, expr_norm: pd.DataFrame) -> CPCAResult:
        self._log("=== Etapa 2.1: CPCA y deteccion de outliers ===")
        found   = [g for g in self.core_genes if g in expr_norm.index]
        missing = [g for g in self.core_genes if g not in expr_norm.index]
        if missing:
            self._log(f"  AVISO: genes core ausentes: {missing}")
        if len(found) < 2:
            raise ValueError(f"CPCA requiere al menos 2 genes core, hay {found}.")
        core = expr_norm.loc[found]
        n = core.shape[1]

        order, scale, pc1, pc2, pc3, var_exp, d = cpca_step(core.values)
        self._log(f"  Varianza PC1/PC2/PC3: {var_exp[0]:.1%} / {var_exp[1]:.1%} / {var_exp[2]:.1%}")

        # ── Criterio radial: candidatas CPCA ──────────────────────────
        cap = int(np.ceil(self.max_outlier_fraction * n))
        n_obs = max(1, min(cap, n))
        obs = np.argsort(d, kind="stable")[:n_obs]
        cand = [int(i) for i in obs if d[i] <= self.tight_radius]
        if not cand:
            cand = [int(i) for i in obs if d[i] <= self.loose_radius]
        self._log(f"  Criterio radial: {n_obs} muestras examinadas, {len(cand)} candidatas")

        # ── Criterio de residuos FMM sobre cada gen core ──────────────
        model = FMMModel(length_alpha_grid=self.fmm_alpha,
                         length_omega_grid=self.fmm_omega, num_reps=self.fmm_reps)
        z_rows = []
        for i in range(len(found)):
            data = core.values[i, order]
            fr = model.fit(data, scale)
            z_rows.append(standardized_residuals(data, fr.fitted))

        multi, uni = self._detect(np.asarray(z_rows), order, set(cand), cap)
        dropped = list(dict.fromkeys(uni + multi))   # union sin repetidos
        self._log(f"  |Z|>{self.z_uni:g}: {len(uni)} muestra(s)  |  "
                  f"|Z|>{self.z_multi:g} y candidata: {len(multi)} muestra(s)  |  "
                  f"total eliminadas: {len(dropped)}")

        # ── Eliminar, renormalizar y repetir CPCA ─────────────────────
        keep = np.setdiff1d(np.arange(n), dropped)
        if dropped:
            core_clean = normalise_matrix(core.iloc[:, keep])
            order, scale, pc1, pc2, pc3, var_exp, d = cpca_step(core_clean.values)
        full = expr_norm.iloc[:, keep].iloc[:, order]
        full = normalise_matrix(full)

        fits: dict[str, FitResult] = {}
        for g in found:
            fits[g] = model.fit(full.loc[g].values, scale)

        result = CPCAResult(
            sample_order=order, circular_scale=scale, core_genes_found=found,
            expr_norm_final=full, fmm_fits_final=fits,
            samples_dropped=[int(i) for i in dropped],
            samples_dropped_names=[str(expr_norm.columns[i]) for i in dropped],
            radial_candidates=cand, outliers_multi=multi, outliers_uni=uni,
            variance_explained=var_exp, n_radial_examined=n_obs,
            pc1=pc1, pc2=pc2, radial_distance=d,
        )
        self._log(result.summary())
        return result

    # ------------------------------------------------------------------
    def _detect(
        self, z: np.ndarray, order: np.ndarray, cand: set[int], cap: int,
    ) -> tuple[list[int], list[int]]:
        """Recorre los genes core con la regla del codigo original.

        ``z`` es (genes core x muestras) en el orden CPCA y ``order`` traduce
        posicion a indice de muestra. El tope se comprueba antes de atender a
        cada gen y no despues.
        """
        multi: list[int] = []
        uni:   list[int] = []
        total = 0
        for row in z:
            if np.any(np.abs(row) > self.z_multi) and total <= cap:
                for pos in np.where(np.abs(row) > self.z_multi)[0]:
                    s = int(order[pos])
                    if s in cand:
                        multi.append(s)
                        total += 1
            if np.any(np.abs(row) > self.z_uni) and total <= cap:
                uni.extend(int(order[p]) for p in np.where(np.abs(row) > self.z_uni)[0])
                total += 1
        return list(dict.fromkeys(multi)), list(dict.fromkeys(uni))

    def _log(self, msg: str) -> None:
        if self.verbose:
            print(msg, flush=True)
