"""
circust/top_genes.py
=============================
Etapa 3 de CIRCUST: genes TOP ritmicos del tejido.

Sigue el texto principal del paper (seccion "3. TOP rhythmic orderings") y el
esquema de la Fig. B del suplementario:

  3.1  Filtro ORI. Se descartan los genes con R2_ORI < 0.5.
  3.2  Se divide la circunferencia en ``n_intervals`` intervalos de igual
       amplitud segun la fase del pico de cada gen.
  3.3  En cada intervalo se ajusta el FMM empezando por los genes de mayor
       R2_ORI. Un gen entra si es no puntiagudo (omega > 0.1) y esta bien
       ajustado (R2_FMM > 0.5). Se para al reunir ``genes_per_interval``.
  3.4  De los genes admitidos en cada intervalo se conservan los de mayor
       R2_FMM del tejido.
  Los genes core que no hayan entrado se fuerzan en el conjunto, como indica
  el paper. Se avisa si la mediana de R2_FMM de los genes core es < 0.3.

El pico que sitÃƒÆ’Ã‚Âºa a cada gen en un intervalo es el del Cosinor, igual que en
el original. Los picos del FMM se usan despues para resumir la cobertura.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from math import pi

import numpy as np
import pandas as pd

from circust.fitting.cosinor import CosinorModel
from circust.fitting.fmm import FMMModel
from circust.fitting.ori import ORIModel

TWO_PI = 2.0 * pi


def circular_median(angles: np.ndarray) -> float:
    """Mediana circular: el dato que minimiza la suma de distancias geodesicas."""
    a = np.asarray(angles, dtype=float) % TWO_PI
    if a.size == 0:
        return 0.0
    if a.size == 1:
        return float(a[0])
    diffs = np.abs(a[:, None] - a[None, :])
    diffs = np.minimum(diffs, TWO_PI - diffs)
    return float(a[np.argmin(diffs.sum(axis=1))] % TWO_PI)

def circular_median_r_compatible(angles: np.ndarray) -> float:
    """Mediana circular compatible con circular 0.5-2 de R."""
    a = np.asarray(angles, dtype=float) % TWO_PI
    n = len(a)

    if n == 0:
        return np.nan

    minimum = np.pi
    medians = []
    eps = np.finfo(float).eps

    for candidate in a:
        distances = np.pi - np.abs(np.pi - np.abs(a - candidate))
        deviation = np.sum(distances) / n

        if (deviation - minimum) / n < -eps:
            minimum = deviation
            medians = [candidate]
        elif abs(deviation - minimum) <= 1e-8:
            medians.append(candidate)

    sin_sum = np.sum(np.sin(medians))
    cos_sum = np.sum(np.cos(medians))

    if np.hypot(sin_sum, cos_sum) / len(medians) <= eps:
        return np.nan

    return float(np.arctan2(sin_sum, cos_sum) % TWO_PI)

def assign_intervals(peaks: np.ndarray, origin: float, n_intervals: int) -> np.ndarray:
    """Intervalo (1..n) de cada pico, con ``n`` arcos iguales desde ``origin``."""
    width = TWO_PI / n_intervals
    rel = (np.asarray(peaks, dtype=float) - origin) % TWO_PI
    return np.minimum((rel // width).astype(int), n_intervals - 1) + 1


def _fit_ori_one(row: np.ndarray, circ: np.ndarray) -> float:
    return ORIModel().fit(row, circ).r2


def _fit_fmm_one(row, circ, alpha_grid, omega_grid, num_reps):
    fr = FMMModel(length_alpha_grid=alpha_grid, length_omega_grid=omega_grid,
                  num_reps=num_reps).fit(row, circ)
    return (fr.peak_time, fr.params["omega"], fr.r2)


@dataclass
class TopGeneResult:
    """Salida de la Etapa 3.

    candidate_matrix : matriz de los genes TOP (genes x muestras), en el orden
        del orden preliminar. Es [X_TOP] del paper.
    gene_names : genes TOP, primero los de cada intervalo y despues los core
        forzados.
    forced_genes : genes core que no superaron los filtros y se incluyeron.
    fmm_peaks, omega, r2_fmm, r2_ori : valores de cada gen TOP en el orden preliminar.
    interval : intervalo al que pertenece cada gen TOP (0 si fue forzado).
    interval_counts : genes admitidos en 3.3 y conservados en 3.4, por intervalo.
    """

    candidate_matrix: pd.DataFrame
    circular_scale:   np.ndarray
    gene_names:       list[str]
    forced_genes:     list[str]
    fmm_peaks:        np.ndarray
    omega:            np.ndarray
    r2_fmm:           np.ndarray
    r2_ori:           np.ndarray
    interval:         np.ndarray
    n_total_genes:    int
    n_after_ori:      int
    n_admitted:       int
    interval_counts:  list[tuple[int, int, int]] = field(default_factory=list)
    interval_origin:  float = 0.0
    cosinor_peaks: np.ndarray | None = None
    def summary(self) -> str:
        lines = [
            "=== Resumen de la Etapa 3 ===",
            f"  Genes totales              : {self.n_total_genes}",
            f"  Tras filtro ORI            : {self.n_after_ori}",
            f"  Admitidos en los intervalos: {self.n_admitted}",
            f"  Genes TOP finales          : {len(self.gene_names)}"
            f"  (core forzados: {len(self.forced_genes)})",
            "  Intervalo | candidatos | admitidos | conservados",
        ]
        for k, cand, adm, kept in self.interval_counts:
            lines.append(f"      {k:>2}    | {cand:>9}  | {adm:>8}  | {kept:>10}")
        return "\n".join(lines)


class TopGeneSelector:
    """Selecciona los genes TOP del tejido (ver el docstring del modulo).

    Parameters
    ----------
    r2_ori_threshold : R2_ORI minimo del filtro 3.1. Paper: 0.5.
    r2_fmm_threshold : R2_FMM minimo para entrar. Paper: 0.5.
    omega_min : omega minimo, para excluir perfiles puntiagudos. Paper: 0.1.
    n_intervals : numero de intervalos de la circunferencia. Original en R: 8.
    genes_per_interval : genes admitidos por intervalo en 3.3. Original: 50.
    reduce_above, top_quantile : paso 3.4. Si un intervalo admite mas de
        ``reduce_above`` genes se conservan los de R2_FMM >= cuantil
        ``top_quantile``. Original: 25 y 0.75.
    n_jobs : procesos para los ajustes. -1 usa todos los nucleos.
    """

    def __init__(
        self,
        r2_ori_threshold:      float = 0.5,
        r2_fmm_threshold:      float = 0.5,
        omega_min:             float = 0.1,
        n_intervals:           int   = 8,
        genes_per_interval:    int   = 50,
        reduce_above:          int   = 25,
        top_quantile:          float = 0.75,
        fmm_length_alpha_grid: int   = 48,
        fmm_length_omega_grid: int   = 24,
        fmm_num_reps:          int   = 3,
        n_jobs:                int   = -1,
        verbose:               bool  = True,
        ori_filter_method:     str   = "fixed",
        fmm_filter_method:     str   = "fixed",
        circular_median_method: str  = "observed",
        sector_numbering: str = "original",
        admission_method: str = "fixed",
    ) -> None:
        if n_intervals < 4:
            raise ValueError("n_intervals debe ser al menos 4: el paper exige "
                             "cubrir los cuatro cuadrantes de la circunferencia.")
        self.r2_ori_threshold   = r2_ori_threshold
        self.r2_fmm_threshold   = r2_fmm_threshold
        self.omega_min          = omega_min
        self.n_intervals        = n_intervals
        self.genes_per_interval = genes_per_interval
        self.reduce_above       = reduce_above
        self.top_quantile       = top_quantile
        self.fmm_alpha          = fmm_length_alpha_grid
        self.fmm_omega          = fmm_length_omega_grid
        self.fmm_reps           = fmm_num_reps
        self.n_jobs             = n_jobs
        self.verbose            = verbose

        if ori_filter_method not in ("fixed", "median"):
            raise ValueError(
                "ori_filter_method debe ser 'fixed' o 'median'."
            )
        self.ori_filter_method = ori_filter_method

        if fmm_filter_method not in ("fixed", "core_quantile"):
            raise ValueError(
                "fmm_filter_method debe ser 'fixed' o 'core_quantile'."
            )
        self.fmm_filter_method = fmm_filter_method
        if circular_median_method not in ("observed", "circular_mean"):
            raise ValueError(
                "circular_median_method debe ser 'observed' o 'circular_mean'."
            )
        self.circular_median_method = circular_median_method
        if sector_numbering not in ("original", "r_compatible"):
            raise ValueError(
                "sector_numbering debe ser 'original' o 'r_compatible'."
            )
        self.sector_numbering = sector_numbering
        if admission_method not in ("fixed", "r_compatible"):
            raise ValueError(
                "admission_method debe ser 'fixed' o 'r_compatible'."
            )
        self.admission_method = admission_method
    # ------------------------------------------------------------------
    def run(
        self,
        expr_norm:      pd.DataFrame,
        circular_scale: np.ndarray,
        seed_genes:     list[str],
        r2_ori_precomputed: pd.Series | None = None,
	r2_fmm_core_precomputed: pd.Series | None = None,
    ) -> TopGeneResult:
        """
        Parameters
        ----------
        expr_norm : matriz normalizada (genes x muestras) en el orden preliminar.
        circular_scale : fases del orden preliminar, alineadas con las columnas.
        seed_genes : genes core presentes en la matriz.
        """
        self._log("=== Etapa 3: genes TOP ===")
        seed_genes = list(dict.fromkeys(seed_genes))
        names = list(expr_norm.index)
        circ  = np.asarray(circular_scale, dtype=float)
        vals  = expr_norm.values.astype(float)

        
        if r2_ori_precomputed is None:
            self._log(f"  Ajustando ORI a {len(names)} genes ...")
            r2_ori = np.array(self._parallel(
                _fit_ori_one,
                [(vals[i], circ) for i in range(len(names))]
            ))
        else:
            if not isinstance(r2_ori_precomputed, pd.Series):
                raise TypeError("r2_ori_precomputed debe ser una pd.Series.")

            if not r2_ori_precomputed.index.is_unique:
                raise ValueError("Hay nombres de genes duplicados en los resultados ORI.")

            if set(r2_ori_precomputed.index) != set(names):
                raise ValueError("Los genes de los resultados ORI no coinciden con la matriz.")

            r2_ori = r2_ori_precomputed.reindex(names).to_numpy(dtype=float)

            if not np.all(np.isfinite(r2_ori)):
                raise ValueError("Los resultados ORI contienen valores no finitos.")

            if np.any((r2_ori < 0) | (r2_ori > 1)):
                raise ValueError("Los valores R2_ORI deben estar entre 0 y 1.")

            self._log(f"  Reutilizando ORI precalculado para {len(names)} genes.")

        if self.ori_filter_method == "median":
            ori_threshold = float(np.median(r2_ori))
            keep = np.where(r2_ori > ori_threshold)[0]
        else:
            ori_threshold = self.r2_ori_threshold
            keep = np.where(r2_ori >= ori_threshold)[0]
        operador_ori = ">" if self.ori_filter_method == "median" else ">="
        self._log(f"  R2_ORI {operador_ori} {ori_threshold:.6f}: {keep.size} genes")
        if keep.size == 0:
            raise ValueError("Ningun gen supera el filtro R2_ORI.")


        cos = CosinorModel()
        cos_peaks = np.array([cos.fit(vals[i], circ).peak_time for i in keep])
        cos_peaks_by_idx = dict(zip(keep, cos_peaks))
        if self.circular_median_method == "circular_mean":
            origin = circular_median_r_compatible(cos_peaks)
        else:
            origin = circular_median(cos_peaks)
        interval = assign_intervals(cos_peaks, origin, self.n_intervals)
        if self.sector_numbering == "r_compatible":
            interval = self.n_intervals + 1 - interval


        if self.fmm_filter_method == "core_quantile":
            if r2_fmm_core_precomputed is None:
                raise ValueError(
                    "El metodo 'core_quantile' requiere r2_fmm_core_precomputed."
                )

            if not isinstance(r2_fmm_core_precomputed, pd.Series):
                raise TypeError(
                    "r2_fmm_core_precomputed debe ser una pd.Series."
                )

            if not r2_fmm_core_precomputed.index.is_unique:
                raise ValueError("Hay genes core duplicados en los resultados FMM.")

            if set(r2_fmm_core_precomputed.index) != set(seed_genes):
                raise ValueError(
                    "Los genes de los resultados FMM no coinciden con los genes core."
                )

            core_r2_values = r2_fmm_core_precomputed.reindex(
                seed_genes
            ).to_numpy(dtype=float)

            if not np.all(np.isfinite(core_r2_values)):
                raise ValueError("Los resultados FMM core contienen valores no finitos.")

            if np.any((core_r2_values < 0) | (core_r2_values > 1)):
                raise ValueError("Los valores R2_FMM core deben estar entre 0 y 1.")

            effective_fmm_threshold = min(
                float(np.quantile(core_r2_values, 0.25)), 0.4
            )
        else:
            effective_fmm_threshold = self.r2_fmm_threshold

        self._log(f"  Umbral R2_FMM aplicado: {effective_fmm_threshold:.6f}")

        admitted: list[dict] = []
        counts: list[tuple[int, int, int]] = []
        for k in range(1, self.n_intervals + 1):
            members = keep[interval == k]
            members = members[np.argsort(-r2_ori[members], kind="stable")]
            got = self._admit(
                vals,
                circ,
                members,
                r2_fmm_threshold=effective_fmm_threshold,
                r2_ori=r2_ori,
                q75_global=float(np.quantile(r2_ori, 0.75)),
            )
            for g in got:
                g["interval"] = k
            admitted.extend(got)
            counts.append((k, int(members.size), len(got), 0))
            if not got:
                self._log(f"  AVISO: el intervalo {k} no admite ningun gen.")

        n_admitted = len(admitted)


        final: list[dict] = []
        for k in range(1, self.n_intervals + 1):
            grp = [g for g in admitted if g["interval"] == k]
            if len(grp) > self.reduce_above:
                cut = np.quantile([g["r2_fmm"] for g in grp], self.top_quantile)
                grp = [g for g in grp if g["r2_fmm"] >= cut]
            final.extend(grp)
            c = counts[k - 1]
            counts[k - 1] = (c[0], c[1], c[2], len(grp))


        in_top = {names[g["idx"]] for g in final}
        forced = [s for s in seed_genes if s in names and s not in in_top]
        if forced:
            self._log(f"  Genes core forzados (no superaron los filtros): {forced}")
            fz = [names.index(s) for s in forced]
            fits = self._parallel(_fit_fmm_one, [
                (vals[i], circ, self.fmm_alpha, self.fmm_omega, self.fmm_reps) for i in fz])
            for i, (pk, om, r2) in zip(fz, fits):
                final.append({"idx": i, "peak": pk, "omega": om, "r2_fmm": r2, "interval": 0})

        idx = [g["idx"] for g in final]
        gene_names = [names[i] for i in idx]
        cosinor_peaks_final = np.array([
            cos_peaks_by_idx[i]
            if i in cos_peaks_by_idx
            else cos.fit(vals[i], circ).peak_time
            for i in idx
        ])
        core_r2 = [g["r2_fmm"] for g in final if names[g["idx"]] in seed_genes]
        if core_r2 and np.median(core_r2) < 0.3:
            self._log("  AVISO: la mediana de R2_FMM de los genes core es < 0.3 tras el "
                      "orden preliminar. Segun el paper, el analisis posterior puede "
                      "ser poco fiable.")

        result = TopGeneResult(
            candidate_matrix=expr_norm.iloc[idx],
            circular_scale=circ, gene_names=gene_names, forced_genes=forced,
            fmm_peaks=np.array([g["peak"] for g in final]),
            omega=np.array([g["omega"] for g in final]),
            r2_fmm=np.array([g["r2_fmm"] for g in final]),
            r2_ori=r2_ori[idx],
            interval=np.array([g["interval"] for g in final], dtype=int),
            cosinor_peaks=cosinor_peaks_final,
            n_total_genes=len(names), n_after_ori=int(keep.size),
            n_admitted=n_admitted, interval_counts=counts, interval_origin=origin,
        )
        self._log(result.summary())
        return result

    # ------------------------------------------------------------------
    def _admit(
        self,
        vals: np.ndarray,
        circ: np.ndarray,
        members: np.ndarray,
        r2_fmm_threshold: float | None = None,
        r2_ori: np.ndarray | None = None,
        q75_global: float | None = None,
    ) -> list[dict]:
        """Recorre ``members`` (ya ordenados) y admite hasta la cuota."""
        got: list[dict] = []
        if self.admission_method == "r_compatible":
            if r2_ori is None or q75_global is None:
                raise ValueError(
                    "El modo r_compatible requiere r2_ori y q75_global."
                )

            if members.size == 0:
                return got

            q25_sector = float(np.quantile(r2_ori[members], 0.25))
            umbral_ori = max(q25_sector, 0.4)
            umbral_fmm = (
                self.r2_fmm_threshold
                if r2_fmm_threshold is None
                else r2_fmm_threshold
            )

            for i in members:
                valor_ori = float(r2_ori[i])

                if valor_ori <= umbral_ori:
                    break

                if (
                    len(got) >= self.genes_per_interval
                    and valor_ori >= q75_global
                ):
                    break

                pk, om, r2 = _fit_fmm_one(
                    vals[i],
                    circ,
                    self.fmm_alpha,
                    self.fmm_omega,
                    self.fmm_reps,
                )

                if om > self.omega_min and r2 > umbral_fmm:
                    got.append({
                        "idx": int(i),
                        "peak": pk,
                        "omega": om,
                        "r2_fmm": r2,
                    })

            return got
        chunk = max(2 * self.genes_per_interval, 32)
        pos = 0
        while pos < members.size and len(got) < self.genes_per_interval:
            block = members[pos:pos + chunk]
            fits = self._parallel(_fit_fmm_one, [
                (vals[i], circ, self.fmm_alpha, self.fmm_omega, self.fmm_reps) for i in block])
            for i, (pk, om, r2) in zip(block, fits):
                if om > self.omega_min and r2 > (
                  self.r2_fmm_threshold if r2_fmm_threshold is None else r2_fmm_threshold
                ):
                    got.append({"idx": int(i), "peak": pk, "omega": om, "r2_fmm": r2})
                    if len(got) == self.genes_per_interval:
                        break
            pos += chunk
        return got

    def _parallel(self, fn, arglist: list[tuple]) -> list:
        if self.n_jobs == 1 or len(arglist) <= 8:
            return [fn(*a) for a in arglist]
        from joblib import Parallel, delayed
        return Parallel(n_jobs=self.n_jobs)(delayed(fn)(*a) for a in arglist)

    def _log(self, msg: str) -> None:
        if self.verbose:
            print(msg, flush=True)
