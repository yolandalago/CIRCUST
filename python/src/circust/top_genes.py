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

El pico que sitúa a cada gen en un intervalo es el del Cosinor, igual que en
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

    # ------------------------------------------------------------------
    def run(
        self,
        expr_norm:      pd.DataFrame,
        circular_scale: np.ndarray,
        seed_genes:     list[str],
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

        # ── 3.1 filtro ORI ────────────────────────────────────────────
        self._log(f"  Ajustando ORI a {len(names)} genes ...")
        r2_ori = np.array(self._parallel(_fit_ori_one, [(vals[i], circ) for i in range(len(names))]))
        keep = np.where(r2_ori >= self.r2_ori_threshold)[0]
        self._log(f"  R2_ORI >= {self.r2_ori_threshold}: {keep.size} genes")
        if keep.size == 0:
            raise ValueError("Ningun gen supera el filtro R2_ORI.")

        # ── 3.2 intervalos segun el pico de cada gen ──────────────────
        cos = CosinorModel()
        cos_peaks = np.array([cos.fit(vals[i], circ).peak_time for i in keep])
        origin = circular_median(cos_peaks)
        interval = assign_intervals(cos_peaks, origin, self.n_intervals)

        # ── 3.3 ajuste FMM por intervalo, de mayor a menor R2_ORI ─────
        admitted: list[dict] = []
        counts: list[tuple[int, int, int]] = []
        for k in range(1, self.n_intervals + 1):
            members = keep[interval == k]
            members = members[np.argsort(-r2_ori[members], kind="stable")]
            got = self._admit(vals, circ, members)
            for g in got:
                g["interval"] = k
            admitted.extend(got)
            counts.append((k, int(members.size), len(got), 0))
            if not got:
                self._log(f"  AVISO: el intervalo {k} no admite ningun gen.")

        n_admitted = len(admitted)

        # ── 3.4 los de mayor R2_FMM de cada intervalo ─────────────────
        final: list[dict] = []
        for k in range(1, self.n_intervals + 1):
            grp = [g for g in admitted if g["interval"] == k]
            if len(grp) > self.reduce_above:
                cut = np.quantile([g["r2_fmm"] for g in grp], self.top_quantile)
                grp = [g for g in grp if g["r2_fmm"] >= cut]
            final.extend(grp)
            c = counts[k - 1]
            counts[k - 1] = (c[0], c[1], c[2], len(grp))

        # ── Genes core forzados ───────────────────────────────────────
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
            n_total_genes=len(names), n_after_ori=int(keep.size),
            n_admitted=n_admitted, interval_counts=counts, interval_origin=origin,
        )
        self._log(result.summary())
        return result

    # ------------------------------------------------------------------
    def _admit(self, vals: np.ndarray, circ: np.ndarray, members: np.ndarray) -> list[dict]:
        """Recorre ``members`` (ya ordenados) y admite hasta la cuota."""
        got: list[dict] = []
        chunk = max(2 * self.genes_per_interval, 32)
        pos = 0
        while pos < members.size and len(got) < self.genes_per_interval:
            block = members[pos:pos + chunk]
            fits = self._parallel(_fit_fmm_one, [
                (vals[i], circ, self.fmm_alpha, self.fmm_omega, self.fmm_reps) for i in block])
            for i, (pk, om, r2) in zip(block, fits):
                if om > self.omega_min and r2 > self.r2_fmm_threshold:
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
