"""
circust/synchronizer.py
================================
Sincronizacion del orden circular (punto de partida y sentido).

El orden de CPCA es invariante a rotacion y su sentido es arbitrario: cada
orden circular representa 2m configuraciones temporales segun el punto de
arranque y la direccion (suplementario, Seccion 3.2). CIRCUST fija ambos con
tres supuestos biologicos:

  (1) el pico de ARNTL se sitúa en pi, lo que fija el punto de partida,
  (2) DBP pica despues de los genes de fase ROR (ARNTL, NPAS, CLOCK),
  (3) la mayoria de los genes core pica en el periodo activo [0, pi).

Los supuestos (2) y (3) determinan el sentido. La implementacion sigue las
funciones ``basicPreOder_cores`` y ``basicOder_cores`` del codigo original.

Paso 1. Rotacion y direccion.
    Se rota la escala de forma que el pico de ARNTL quede en pi. Si
    ``(pico_DBP - pico_ARNTL + pi) mod 2pi < pi`` (DBP pica despues de ARNTL)
    se mantiene el sentido y si no se invierte.

Paso 2. Verificacion.
    Con los picos ya rotados se evalua un aviso::

        aviso = (1 - cos(pico_DBP) <= 0.1) or (1 - cos(pico_DBP - pi) <= 0.1)
                or (n_genes_core_en_[0,pi) < (n_core - 2) / 2)

    donde el recuento excluye ARNTL y DBP. Si hay aviso y ademas NO se
    cumple la cadena ``pico_DBP < pico_CRY1 < pico_ARNTL``, se invierte de
    nuevo el sentido. Es una comparacion lineal de los picos, igual que en el
    original. Si CRY1 no esta entre los genes core la cadena no se puede
    evaluar y no se invierte.
"""
from __future__ import annotations

from dataclasses import dataclass
from math import pi

import numpy as np
import pandas as pd

from circust.fitting.fmm import FMMModel
from circust.fitting.rhythm_model import FitResult

TWO_PI = 2.0 * pi


@dataclass
class SynchronizationResult:
    """Salida de la sincronizacion.

    sample_order : indices de columna de la matriz de entrada, en el orden final.
    circular_scale : fases en [0, 2pi) tras la correccion.
    expr_ordered : matriz de entrada reordenada.
    peak_times, r2_fmm, fmm_params : picos, R2 y parametros [M, A, alpha, beta,
        omega] de los genes core en el marco final.
    direction_flipped : True si el sentido neto quedo invertido respecto al CPCA.
    """

    sample_order:      np.ndarray
    circular_scale:    np.ndarray
    expr_ordered:      pd.DataFrame
    peak_times:        np.ndarray
    r2_fmm:            np.ndarray
    fmm_params:        np.ndarray
    direction_flipped: bool
    core_genes:        list[str]
    anchor_peak_raw:   float = 0.0
    reversed_step1:    bool = False
    flipped_step2:     bool = False


class CircularSynchronizer:
    """Rotacion al ancla y eleccion del sentido (regla del codigo original)."""

    def __init__(
        self,
        anchor_gene:     str  = "ARNTL",
        direction_gene:  str  = "DBP",
        reference_gene:  str  = "CRY1",
        verbose:         bool = True,
    ) -> None:
        self.anchor_gene    = anchor_gene
        self.direction_gene = direction_gene
        self.reference_gene = reference_gene
        self.verbose        = verbose

    # ------------------------------------------------------------------
    def run(
        self,
        expr_ordered:   pd.DataFrame,
        circular_scale: np.ndarray,
        core_genes:     list[str],
        fmm_fits:       dict[str, FitResult],
    ) -> SynchronizationResult:
        core_genes = list(core_genes)
        for g in (self.anchor_gene, self.direction_gene):
            if g not in fmm_fits:
                raise ValueError(f"La sincronizacion necesita el gen {g} entre los genes core.")

        scale = np.asarray(circular_scale, dtype=float)
        n_core = len(core_genes)
        r2 = np.array([fmm_fits[g].r2 for g in core_genes])
        anchor_peak = fmm_fits[self.anchor_gene].peak_time
        dir_peak    = fmm_fits[self.direction_gene].peak_time

        # ── Paso 1: rotar al ancla y fijar el sentido con DBP ──────────
        rotated = (scale - anchor_peak + pi) % TWO_PI
        esc_t   = np.argsort(rotated, kind="stable")
        esc2    = rotated[esc_t]
        forward = ((dir_peak - anchor_peak + pi) % TWO_PI) < pi

        peaks = np.zeros(n_core)
        pars  = np.zeros((n_core, 5))
        for i, g in enumerate(core_genes):
            p = fmm_fits[g].params
            al = (p["alpha"] - anchor_peak + pi) % TWO_PI
            if forward:
                peaks[i] = FMMModel.peak_time(al, p["beta"], p["omega"])
                pars[i]  = [p["M"], p["A"], al, p["beta"], p["omega"]]
            else:
                a2, b2 = (TWO_PI - al) % TWO_PI, (TWO_PI - p["beta"]) % TWO_PI
                peaks[i] = FMMModel.peak_time(a2, b2, p["omega"])
                pars[i]  = [p["M"], p["A"], a2, b2, p["omega"]]

        if forward:
            order, esc = esc_t, esc2
        else:
            order, esc = esc_t[::-1].copy(), (TWO_PI - esc2[::-1]) % TWO_PI

        # ── Paso 2: verificacion y posible segunda inversion ───────────
        flipped2 = self._needs_flip(core_genes, peaks)
        if flipped2:
            peaks = (TWO_PI - peaks) % TWO_PI
            order = order[::-1].copy()
            esc   = (TWO_PI - esc[::-1]) % TWO_PI
            pars[:, 2] = (TWO_PI - pars[:, 2]) % TWO_PI
            pars[:, 3] = (TWO_PI - pars[:, 3]) % TWO_PI

        self._log(f"    {self.anchor_gene} pico {anchor_peak:.3f} rad | sentido paso 1: "
                  f"{'directo' if forward else 'invertido'} | paso 2: "
                  f"{'invertido' if flipped2 else 'sin cambio'}")

        return SynchronizationResult(
            sample_order=order, circular_scale=esc,
            expr_ordered=expr_ordered.iloc[:, order].copy(),
            peak_times=peaks, r2_fmm=r2, fmm_params=pars,
            direction_flipped=bool((not forward) ^ flipped2),
            core_genes=core_genes, anchor_peak_raw=float(anchor_peak),
            reversed_step1=not forward, flipped_step2=flipped2,
        )

    # ------------------------------------------------------------------
    def _needs_flip(self, core_genes: list[str], peaks: np.ndarray) -> bool:
        """Paso 2 de ``basicOder_cores``."""
        idx = {g: i for i, g in enumerate(core_genes)}
        dir_peak = peaks[idx[self.direction_gene]]
        n_day = sum(
            1 for g, i in idx.items()
            if g not in (self.anchor_gene, self.direction_gene) and 0.0 < peaks[i] < pi
        )
        half = (len(core_genes) - 2) / 2.0
        p6am = 1.0 - np.cos(dir_peak)
        p6pm = 1.0 - np.cos(dir_peak - pi)
        warning = (p6am <= 0.1) or (p6pm <= 0.1) or (n_day < half)
        if not warning:
            return False
        if self.reference_gene not in idx:
            self._log(f"    AVISO: {self.reference_gene} no esta entre los genes core, "
                      "no se puede evaluar la cadena DBP < CRY1 < ARNTL. Sin inversion.")
            return False
        chain = (peaks[idx[self.direction_gene]] < peaks[idx[self.reference_gene]]
                 < peaks[idx[self.anchor_gene]])
        return not chain

    def _log(self, msg: str) -> None:
        if self.verbose:
            print(msg, flush=True)
