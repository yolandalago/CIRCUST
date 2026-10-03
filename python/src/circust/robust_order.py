"""
circust/robust_order.py
================================
Etapas 3 (seleccion de los K subconjuntos) y 4 (estimacion robusta) de CIRCUST.

Seleccion aleatoria (paper, "3. TOP rhythmic orderings")
--------------------------------------------------------
Se extraen subconjuntos aleatorios de tamano 2/3 de los genes TOP y se
repite hasta reunir K que cumplan:

  (a) los valores angulares theta del CPCA del subconjunto se reparten por
      mas de media circunferencia (el arco que cubren supera pi),
  (b) la distancia maxima entre dos angulos consecutivos no supera la
      maxima observada en el orden preliminar.

El muestreo se hace sobre los genes TOP que superaron los filtros por si
mismos. Los genes core que se forzaron en el conjunto no entran en el sorteo,
aunque si se usan despues para sincronizar cada orden y se ajustan como el resto.
Si tras ``max_attempts`` sorteos no se reunen K validos, se completan con los
de menor distancia maxima entre angulos y se avisa.

Para cada subconjunto se recalcula el CPCA, se sincroniza el orden resultante
(ARNTL en pi y sentido, ver ``synchronizer``) y se reordena la matriz TOP.
Esos K ordenes sincronizados, o^(1), ..., o^(K), son el resultado.

Se ajusta el FMM a todos los genes TOP en cada orden. Cada gen tiene K
estimaciones de sus parametros y de R2. El atlas [M^TOP] contiene la mediana
de esas K estimaciones. Los parametros angulares (alpha, beta, t_U, t_L) se
resumen con la mediana circular y el resto con la mediana ordinaria. No se
construye ningun orden de muestras unico: CIRCUST entrega los K ordenes.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from math import pi

import numpy as np
import pandas as pd

from circust.cpca import cpca_step
from circust.fitting.fmm import FMMModel
from circust.synchronizer import CircularSynchronizer
from circust.top_genes import circular_median

TWO_PI = 2.0 * pi
PARAM_COLUMNS = ["M", "A", "alpha", "beta", "omega", "t_U", "t_L", "R2"]
ANGLE_COLUMNS = ["alpha", "beta", "t_U", "t_L"]


def max_circular_gap(theta: np.ndarray) -> float:
    """Mayor distancia entre dos angulos consecutivos sobre la circunferencia."""
    t = np.sort(np.asarray(theta, dtype=float) % TWO_PI)
    if t.size < 2:
        return TWO_PI
    return float(max(np.max(np.diff(t)), TWO_PI - t[-1] + t[0]))


def _fit_gene(args):
    """Ajusta FMM a un gen en un orden y devuelve sus parametros sincronizados."""
    (gene, data, scale, grid, anchor_peak_raw, net_flipped, order) = args
    alpha_grid, omega_grid, reps = grid
    fr = FMMModel(length_alpha_grid=alpha_grid, length_omega_grid=omega_grid,
                  num_reps=reps).fit(data, scale)
    p = fr.params
    al = (p["alpha"] - anchor_peak_raw + pi) % TWO_PI
    if not net_flipped:
        alpha, beta = al, p["beta"]
    else:
        alpha, beta = (TWO_PI - al) % TWO_PI, (TWO_PI - p["beta"]) % TWO_PI
    t_u = FMMModel.peak_time(alpha, beta, p["omega"])
    t_l = FMMModel.trough_time(alpha, beta, p["omega"])
    return (gene, p["M"], p["A"], alpha, beta, p["omega"], t_u, t_l, float(fr.r2),
            fr.fitted[order])


@dataclass
class RobustOrderResult:
    """Salida de las Etapas 3 (seleccion) y 4.

    sample_orders : (K, n) indices de columna de la matriz TOP preliminar que
        dan cada orden sincronizado. La muestra en la posicion j del orden k
        es ``columns[sample_orders[k, j]]``.
    circular_scales : (K, n) fases de cada orden, en [0, 2pi).
    top_matrices : K matrices TOP (genes TOP x muestras) reordenadas.
    fitted : K matrices con los valores ajustados del FMM de los genes TOP.
    parameters : tabla larga (gen, k, M, A, alpha, beta, omega, t_U, t_L, R2).
    atlas : [M^TOP], mediana de los K parametros por gen TOP.
    selections : genes del subconjunto aleatorio de cada repeticion.
    direction_flipped : si el sentido de cada orden quedo invertido respecto al CPCA.
    """

    sample_names:      list[str]
    sample_orders:     np.ndarray
    circular_scales:   np.ndarray
    top_matrices:      list[pd.DataFrame]
    fitted:            list[pd.DataFrame]
    parameters:        pd.DataFrame
    atlas:             pd.DataFrame
    selections:        list[list[str]]
    direction_flipped: np.ndarray
    max_gaps:          np.ndarray
    prelim_max_gap:    float
    n_attempts:        int
    n_relaxed:         int = 0

    @property
    def n_reps(self) -> int:
        return self.sample_orders.shape[0]

    def ordered_names(self, k: int) -> list[str]:
        return [self.sample_names[i] for i in self.sample_orders[k]]


class RobustOrderEstimator:
    """Etapas 3 (K subconjuntos) y 4 (estimacion robusta).

    Parameters
    ----------
    n_reps : numero K de ordenes. Por defecto 5.
    sample_size_fraction : fraccion de genes por subconjunto. Por defecto 2/3.
    max_attempts : sorteos maximos antes de relajar las condiciones. Original: 5000.
    anchor_gene, direction_gene, reference_gene : genes de la sincronizacion.
    seed : semilla del generador. El original no fija ninguna.
    """

    def __init__(
        self,
        n_reps:                int   = 5,
        sample_size_fraction:  float = 2.0 / 3.0,
        max_attempts:          int   = 5000,
        anchor_gene:           str   = "ARNTL",
        direction_gene:        str   = "DBP",
        reference_gene:        str   = "CRY1",
        fmm_length_alpha_grid: int   = 48,
        fmm_length_omega_grid: int   = 24,
        fmm_num_reps:          int   = 3,
        seed:                  int   = 42,
        n_jobs:                int   = -1,
        verbose:               bool  = True,
    ) -> None:
        self.n_reps               = n_reps
        self.sample_size_fraction = sample_size_fraction
        self.max_attempts         = max_attempts
        self.anchor_gene          = anchor_gene
        self.direction_gene       = direction_gene
        self.reference_gene       = reference_gene
        self.grid                 = (fmm_length_alpha_grid, fmm_length_omega_grid, fmm_num_reps)
        self.seed                 = seed
        self.n_jobs               = n_jobs
        self.verbose              = verbose

    # ------------------------------------------------------------------
    def run(
        self,
        top_matrix:     pd.DataFrame,
        prelim_scale:   np.ndarray,
        core_genes:     list[str],
        forced_genes:   list[str] | None = None,
    ) -> RobustOrderResult:
        """
        Parameters
        ----------
        top_matrix : [X_TOP], genes TOP x muestras en el orden preliminar.
        prelim_scale : fases del orden preliminar.
        core_genes : genes core presentes, para sincronizar cada orden.
        forced_genes : core forzados en la Etapa 3, que no entran en el sorteo.
        """
        self._log("=== Etapas 3 (subconjuntos) y 4 (estimacion robusta) ===")
        forced = set(forced_genes or [])
        names = list(top_matrix.index)
        pool = [i for i, g in enumerate(names) if g not in forced]
        if len(pool) < 3:
            raise ValueError("Quedan muy pocos genes TOP para sortear subconjuntos.")
        tam = max(2, int(np.ceil(self.sample_size_fraction * len(pool))))
        prelim_gap = max_circular_gap(prelim_scale)
        self._log(f"  Genes TOP: {len(names)}  | en el sorteo: {len(pool)}  | tam={tam}  | K={self.n_reps}")

        selections, gaps, n_attempts, n_relaxed = self._select(
            top_matrix.values, np.array(pool), tam, prelim_gap)

        cols = list(top_matrix.columns)
        n = len(cols)
        sample_orders = np.zeros((self.n_reps, n), dtype=int)
        circ_scales   = np.zeros((self.n_reps, n))
        flips = np.zeros(self.n_reps, dtype=bool)
        top_mats, fitted_mats, long_rows = [], [], []
        core_found = [g for g in core_genes if g in names]

        for k, sel in enumerate(selections):
            self._log(f"  --- Orden {k + 1}/{self.n_reps} ---")
            order_k, scale_k = cpca_step(top_matrix.values[sel])[:2]
            x_k = top_matrix.iloc[:, order_k]

            model = FMMModel(length_alpha_grid=self.grid[0],
                             length_omega_grid=self.grid[1], num_reps=self.grid[2])
            fits = {g: model.fit(x_k.loc[g].values, scale_k) for g in core_found}
            sync = CircularSynchronizer(
                self.anchor_gene, self.direction_gene, self.reference_gene,
                verbose=self.verbose,
            ).run(x_k, scale_k, core_found, fits)

            args = [(g, x_k.values[i], scale_k, self.grid, sync.anchor_peak_raw,
                     sync.direction_flipped, sync.sample_order)
                    for i, g in enumerate(names)]
            out = self._parallel(args)

            for (g, M, A, al, be, om, tu, tl, r2, _f) in out:
                long_rows.append((g, k + 1, M, A, al, be, om, tu, tl, r2))
            perm = order_k[sync.sample_order]
            sample_orders[k] = perm
            circ_scales[k]   = sync.circular_scale
            flips[k]         = sync.direction_flipped
            col_names = [cols[i] for i in perm]
            top_mats.append(pd.DataFrame(top_matrix.values[:, perm], index=names, columns=col_names))
            fitted_mats.append(pd.DataFrame(np.vstack([o[-1] for o in out]),
                                            index=names, columns=col_names))

        params = pd.DataFrame(long_rows, columns=["gene", "k"] + PARAM_COLUMNS)
        return RobustOrderResult(
            sample_names=[str(c) for c in cols], sample_orders=sample_orders,
            circular_scales=circ_scales, top_matrices=top_mats, fitted=fitted_mats,
            parameters=params, atlas=self.median_atlas(params, names),
            selections=[[names[i] for i in s] for s in selections],
            direction_flipped=flips, max_gaps=np.array(gaps),
            prelim_max_gap=prelim_gap, n_attempts=n_attempts, n_relaxed=n_relaxed,
        )

    # ------------------------------------------------------------------
    @staticmethod
    def median_atlas(params: pd.DataFrame, gene_order: list[str] | None = None) -> pd.DataFrame:
        """[M^TOP]: mediana de los K parametros de cada gen."""
        rows = []
        for gene, d in params.groupby("gene", sort=False):
            row = {"gene_name": gene}
            for c in PARAM_COLUMNS:
                row[c] = (circular_median(d[c].values) if c in ANGLE_COLUMNS
                          else float(np.median(d[c].values)))
            rows.append(row)
        atlas = pd.DataFrame(rows)
        if gene_order is not None:
            atlas = atlas.set_index("gene_name").loc[gene_order].reset_index()
        atlas = atlas.rename(columns={"R2": "r2_fmm"})
        atlas["classification"] = np.where(atlas["t_U"] < pi, "day", "night")
        return atlas

    # ------------------------------------------------------------------
    def _select(self, values: np.ndarray, pool: np.ndarray, tam: int, prelim_gap: float):
        """Sortea hasta reunir K subconjuntos que cumplan (a) y (b)."""
        rng = np.random.default_rng(self.seed)
        valid, valid_gaps = [], []
        stored, stored_gaps = [], []
        attempts = 0
        while len(valid) < self.n_reps and attempts < self.max_attempts:
            attempts += 1
            sel = np.sort(rng.choice(pool, size=tam, replace=False))
            theta = cpca_step(values[sel])[1]
            gap = max_circular_gap(theta)
            span = TWO_PI - gap
            if span > pi and gap <= prelim_gap:
                valid.append(sel); valid_gaps.append(gap)
            else:
                stored.append(sel); stored_gaps.append(gap)
        n_relaxed = 0
        if len(valid) < self.n_reps:
            need = self.n_reps - len(valid)
            best = np.argsort(stored_gaps)[:need]
            valid += [stored[i] for i in best]
            valid_gaps += [stored_gaps[i] for i in best]
            n_relaxed = need
            self._log(f"  AVISO: tras {attempts} sorteos solo {self.n_reps - need} cumplen (a) y (b). "
                      f"Se completan {need} con los de menor distancia maxima entre angulos.")
        self._log(f"  Subconjuntos validos: {self.n_reps - n_relaxed}/{self.n_reps} "
                  f"en {attempts} sorteos (hueco maximo preliminar {prelim_gap:.3f} rad)")
        return valid, valid_gaps, attempts, n_relaxed

    def _parallel(self, args: list[tuple]) -> list:
        if self.n_jobs == 1 or len(args) <= 8:
            return [_fit_gene(a) for a in args]
        from joblib import Parallel, delayed
        return Parallel(n_jobs=self.n_jobs)(delayed(_fit_gene)(a) for a in args)

    def _log(self, msg: str) -> None:
        if self.verbose:
            print(msg, flush=True)
