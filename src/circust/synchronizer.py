"""
circust/synchronizer.py
=======================
Etapa 2: Sincronizacion del orden circular.

Objetivo
--------
El ordenamiento circular derivado de CPCA es invariante a rotacion y tiene
una direccion arbitraria (sentido horario o antihorario). Esta etapa corrige
ambos problemas estableciendo un marco de referencia biologico.

Sincronizacion manual basada en genes ancla conocidos
------------------------------------------------------
* Un gen de referencia (``anchor_gene``, por defecto ARNTL) se coloca en
  π en la escala circular — su pico define el punto de arranque.
* Un gen de direccion (``direction_gene``, por defecto DBP) determina la
  orientacion (CW vs CCW): si su pico cae en la mitad incorrecta se invierte.
* Un gen de consistencia (``consistency_gene``, por defecto CRY1) participa
  en una verificacion de consistencia biologica refinada (Paso 2.2).

Paso 2.1 — Rotacion y determinacion de direccion:
    1. Rotar escala circular: peak(anchor) → π.
    2. Si (direction_peak − anchor_peak + π) mod 2π ≥ π → invertir.
    3. Reparametrizar ajustes FMM de genes core en el nuevo marco.
    4. Clasificar genes core como *dia* (pico ∈ [0, π)) o *noche* (pico ∈ [π, 2π)).

Paso 2.2 — Verificacion de consistencia biologica:
    Si direction_gene esta demasiado cerca de 0 o π, o menos de la mitad
    de los genes sin anchor/direction tienen pico en [0, π), y el orden
    direction < consistency < anchor NO se cumple, entonces invertir.

Posicion en el pipeline
-----------------------
    OutlierRefiner  →  CircularSynchronizer  →  (Etapa 3)
"""
import numpy as np
import pandas as pd
from dataclasses import dataclass, field
from math import pi

from circust.cpca import CPCAResult
from circust.fitting.fmm import FMMModel
from circust.fitting.rhythm_model import FitResult


# ---------------------------------------------------------------------------
# Dataclass de resultado
# ---------------------------------------------------------------------------

@dataclass
class SynchronizationResult:
    """
    Salida combinada de la Etapa 2.

    Attributes
    ----------
    sample_order : np.ndarray de int, forma (n_muestras,)
        Indices finales de muestra en la matriz de expresion limpia.
        Equivalente en R: ``basicOrdRefG2[[1]]``  (= ``indSin``).

    circular_scale : np.ndarray de float, forma (n_muestras,)
        Eje temporal circular en [0, 2π) tras correccion de direccion.
        Equivalente en R: ``basicOrdRefG2[[2]]``  (= ``escSincroRefG``).

    expr_ordered : pd.DataFrame, forma (n_genes, n_muestras)
        Matriz de expresion normalizada completa en el orden circular final.
        Equivalente en R: ``basicOrdRefG2[[3]]``  (= ``matNewNew``).

    peak_times : np.ndarray de float, forma (n_genes_core,)
        Tiempo de pico FMM para cada gen central en el marco de coordenadas final.
        Equivalente en R: ``basicOrdRefG2[[4]]``  (= ``peaksPreNew``).

    r2_fmm : np.ndarray de float, forma (n_genes_core,)
        R² del ajuste FMM para cada gen central.
        Equivalente en R: ``basicOrdRefG2[[5]]``.

    fmm_params : np.ndarray de float, forma (n_genes_core, 5)
        Parametros FMM [M, A, α, β, ω] por gen central en el nuevo marco.
        Equivalente en R: ``basicOrdRefG2[[6]]``  (= ``pars``).

    direction_flipped : bool
        True si la direccion fue invertida (por el Paso 2.1 o 2.2).
        Equivalente en R: ``basicOrdRefG2[[7]]``  (= ``cambiaOri``).

    within_group_indices : np.ndarray de int, forma (n_muestras,)
        Indices posicionales base 0 dentro del grupo tras posible inversion.
        Equivalente en R: ``basicOrdRefG2[[8]]``  (= ``indNewNew``).

    core_genes : list[str]
        Simbolos de genes centrales en el orden utilizado.

    day_genes : list[str]
        Genes centrales (exc. anchor y direction) con pico en [0, π).

    night_genes : list[str]
        Genes centrales (exc. anchor y direction) con pico en [π, 2π).

    pre_sample_order : np.ndarray de int
        Orden de muestras del Paso 2.1 (antes de la verificacion del Paso 2.2).
        Equivalente en R: ``preOrdRefG2[[1]]``.

    pre_circular_scale : np.ndarray de float
        Escala circular del Paso 2.1.
        Equivalente en R: ``preOrdRefG2[[2]]``.

    pre_expr_ordered : pd.DataFrame
        Matriz de expresion ordenada por el Paso 2.1.
        Equivalente en R: ``preOrdRefG2[[3]]``.

    pre_peak_times : np.ndarray de float
        Tiempos de pico FMM del Paso 2.1.
        Equivalente en R: ``preOrdRefG2[[4]]``.

    pre_fmm_params : np.ndarray de float, forma (n_genes_core, 5)
        Parametros FMM del Paso 2.1.
        Equivalente en R: ``preOrdRefG2[[6]]``  (= ``parCore``).

    pre_direction_reversed : bool
        True si el Paso 2.1 invirtio la direccion (direction_gene en la mitad incorrecta).
        Equivalente en R: ``preOrdRefG2[[13]]``  (= ``reverse``).
    """

    sample_order:           np.ndarray
    circular_scale:         np.ndarray
    expr_ordered:           pd.DataFrame
    peak_times:             np.ndarray
    r2_fmm:                 np.ndarray
    fmm_params:             np.ndarray
    direction_flipped:      bool
    within_group_indices:   np.ndarray
    core_genes:             list

    day_genes:              list       = field(default_factory=list)
    night_genes:            list       = field(default_factory=list)

    pre_sample_order:       np.ndarray = field(default_factory=lambda: np.array([], dtype=int))
    pre_circular_scale:     np.ndarray = field(default_factory=lambda: np.array([]))
    pre_expr_ordered:       pd.DataFrame = field(default_factory=pd.DataFrame)
    pre_peak_times:         np.ndarray = field(default_factory=lambda: np.array([]))
    pre_fmm_params:         np.ndarray = field(default_factory=lambda: np.zeros((0, 5)))
    pre_direction_reversed: bool       = False

    def summary(self) -> str:
        lines = [
            "=== Resumen de Sincronizacion Circular ===",
            f"  Direccion invertida  : {self.direction_flipped}",
            f"  Paso 2.1 invertido  : {self.pre_direction_reversed}",
            f"  Genes dia  (0..π)   : {self.day_genes}",
            f"  Genes noche (π..2π) : {self.night_genes}",
            f"  Picos genes core    : {np.round(self.peak_times, 3).tolist()}",
        ]
        return "\n".join(lines)


# ---------------------------------------------------------------------------
# CircularSynchronizer
# ---------------------------------------------------------------------------

class CircularSynchronizer:
    """
    Etapa 2 de CIRCUST: establece el marco de referencia biologico del
    ordenamiento circular.

    Parameters
    ----------
    anchor_gene : str
        Gen cuyo pico se coloca en π para definir el punto de arranque del
        ordenamiento. En el reloj mamifero, ARNTL tiene su pico cerca de la
        medianoche subjetiva (Zeitgeber 0).
        Por defecto: ``"ARNTL"``.

    direction_gene : str
        Gen que determina la orientacion CW/CCW del ordenamiento. Se espera
        que su pico caiga en [0, π) (fase activa) tras la rotacion.
        Por defecto: ``"DBP"``.

    consistency_gene : str
        Gen auxiliar usado en la verificacion de consistencia biologica
        del Paso 2.2. Se comprueba el orden direction < consistency < anchor.
        Por defecto: ``"CRY1"``.

    verbose : bool
        Imprimir mensajes de progreso.

    Example
    -------
    >>> syncer = CircularSynchronizer()
    >>> result = syncer.run(refined)
    >>> print(result.summary())
    """

    def __init__(
        self,
        anchor_gene:      str  = "ARNTL",
        direction_gene:   str  = "DBP",
        consistency_gene: str  = "CRY1",
        verbose:          bool = True,
    ) -> None:
        self.anchor_gene      = anchor_gene
        self.direction_gene   = direction_gene
        self.consistency_gene = consistency_gene
        self.verbose          = verbose

    # ------------------------------------------------------------------
    # API publica
    # ------------------------------------------------------------------

    def run(
        self,
        refined: CPCAResult,
    ) -> SynchronizationResult:
        """
        Ejecuta la sincronizacion circular.

        Parameters
        ----------
        refined : CPCAResult
            Salida de ``CPCA.run()``. Los genes core se toman de
            ``refined.core_genes_found``.

        Returns
        -------
        SynchronizationResult
        """
        self._log("=== Etapa 2: Sincronizacion Circular ===")
        core_genes = list(refined.core_genes_found)

        # ── Paso 2.1: rotacion + determinacion de direccion ───────────────
        self._log("  Paso 2.1 — rotando al marco de referencia ...")
        (o_pre, esc_pre, mat_pre, peaks_pre, r2_fmm, par_pre,
         names_day, names_night, reversed_21) = self._pre_order(
            refined, core_genes
        )

        # ── Paso 2.2: verificacion de consistencia biologica ──────────────
        self._log("  Paso 2.2 — verificando consistencia biologica ...")
        (o_fin, esc_fin, mat_fin, peaks_fin,
         pars_fin, flipped_22, ind_new) = self._basic_order(
            o_pre, esc_pre, mat_pre, peaks_pre, par_pre, core_genes
        )

        direction_flipped = reversed_21 ^ flipped_22

        result = SynchronizationResult(
            sample_order           = o_fin,
            circular_scale         = esc_fin,
            expr_ordered           = mat_fin,
            peak_times             = peaks_fin,
            r2_fmm                 = r2_fmm,
            fmm_params             = pars_fin,
            direction_flipped      = direction_flipped,
            within_group_indices   = ind_new,
            core_genes             = core_genes,
            day_genes              = names_day,
            night_genes            = names_night,
            pre_sample_order       = o_pre,
            pre_circular_scale     = esc_pre,
            pre_expr_ordered       = mat_pre,
            pre_peak_times         = peaks_pre,
            pre_fmm_params         = par_pre,
            pre_direction_reversed = reversed_21,
        )

        self._log(result.summary())
        return result

    # ------------------------------------------------------------------
    # Pasos internos
    # ------------------------------------------------------------------

    def _pre_order(
        self,
        refined:    CPCAResult,
        core_genes: list[str],
    ) -> tuple:
        """
        Returns
        -------
        o_new       : np.ndarray  — indices de ordenamiento de muestras
        esc_new     : np.ndarray  — escala circular en [0, 2π)
        mat_new     : pd.DataFrame — matriz completa en nuevo orden
        peaks       : np.ndarray  — tiempos de pico FMM por gen central
        r2_fmm      : np.ndarray  — R² por gen central
        par_core    : np.ndarray  forma (n_core, 5) — [M, A, α, β, ω]
        names_day   : list[str]   — genes con pico en [0, π)
        names_night : list[str]   — genes con pico en [π, 2π)
        reversed    : bool        — True si la direccion se invirtio
        """
        circular_scale = refined.circular_scale
        expr_ordered   = refined.expr_norm_final
        fmm_fits       = refined.fmm_fits_final
        peak_times_fin = refined.fmm_peak_times_final

        n_core = len(core_genes)

        # ── R² de ajustes FMM con ordenamiento final ──────────────────────
        r2_fmm = np.array([
            fmm_fits[g].r2 if g in fmm_fits else float("nan")
            for g in core_genes
        ])

        # ── Tiempos de pico anchor y direction en el marco CPCA original ──
        anchor_peak    = peak_times_fin[self.anchor_gene]
        direction_peak = peak_times_fin[self.direction_gene]

        # ── Rotar: anchor → π ─────────────────────────────────────────────
        rotated = (circular_scale - anchor_peak + pi) % (2.0 * pi)
        esc_T   = np.argsort(rotated)
        esc2    = rotated[esc_T]

        # ── Determinar orientacion ─────────────────────────────────────────
        direction_rotated = (direction_peak - anchor_peak + pi) % (2.0 * pi)
        forward           = direction_rotated < pi

        # ── Reparametrizar FMM de cada gen central en el nuevo marco ──────
        peaks    = np.zeros(n_core)
        par_core = np.zeros((n_core, 5))

        for i, gene in enumerate(core_genes):
            if gene not in fmm_fits:
                continue
            fr = fmm_fits[gene]
            alpha = fr.params["alpha"]
            beta  = fr.params["beta"]
            omega = fr.params["omega"]
            M     = fr.params["M"]
            A     = fr.params["A"]

            new_alpha = (alpha - anchor_peak + pi) % (2.0 * pi)

            if forward:
                peaks[i]    = FMMModel.peak_time(new_alpha, beta, omega)
                par_core[i] = [M, A, new_alpha, beta, omega]
            else:
                na2 = (2.0 * pi - new_alpha) % (2.0 * pi)
                nb2 = (2.0 * pi - beta)      % (2.0 * pi)
                peaks[i]    = FMMModel.peak_time(na2, nb2, omega)
                par_core[i] = [M, A, na2, nb2, omega]

        # ── Reordenar columnas de la matriz ───────────────────────────────
        # Nota: ``iloc`` preserva los nombres originales de muestra; **no**
        # reseteamos a ``range(n)`` para que los consumidores (R, scripts
        # de validacion, etc.) puedan emparejar por nombre directamente.
        if forward:
            o_new   = esc_T
            esc_new = esc2
            mat_new = expr_ordered.iloc[:, esc_T].copy()
        else:
            o_new   = esc_T[::-1].copy()
            esc_new = (2.0 * pi - esc2[::-1]) % (2.0 * pi)
            mat_new = expr_ordered.iloc[:, esc_T[::-1]].copy()

        # ── Clasificar como dia / noche ────────────────────────────────────
        names_day   = []
        names_night = []
        for i, gene in enumerate(core_genes):
            if gene in (self.anchor_gene, self.direction_gene):
                continue
            if 0 <= peaks[i] < pi:
                names_day.append(gene)
            else:
                names_night.append(gene)

        self._log(
            f"    {self.anchor_gene} peak: {anchor_peak:.3f} rad  |  "
            f"{self.direction_gene} rotado: {direction_rotated:.3f} rad  |  "
            f"Direccion: {'directa' if forward else 'INVERTIDA'}"
        )

        return (o_new, esc_new, mat_new, peaks, r2_fmm, par_core,
                names_day, names_night, not forward)

    def _basic_order(
        self,
        o:          np.ndarray,
        esc:        np.ndarray,
        mat:        pd.DataFrame,
        peaks:      np.ndarray,
        par:        np.ndarray,
        core_genes: list[str],
    ) -> tuple:
        """
        Paso 2.2 — Equivalente en R: ``basicOder_cores``.

        Returns
        -------
        o_new, esc_new, mat_new, peaks_new, pars_new, flipped, ind_new
        """
        anchor_i      = core_genes.index(self.anchor_gene)      if self.anchor_gene      in core_genes else None
        direction_i   = core_genes.index(self.direction_gene)   if self.direction_gene   in core_genes else None
        consistency_i = core_genes.index(self.consistency_gene) if self.consistency_gene in core_genes else None

        if direction_i is None:
            self._log(
                f"  AVISO: {self.direction_gene} no encontrado en core_genes; "
                "omitiendo verificacion de direccion."
            )
            n = len(o)
            return o, esc, mat, peaks, par, False, np.arange(n)

        direction_peak = peaks[direction_i]
        n_core         = len(core_genes)

        # Contar genes con pico en el semicirculo [0, pi), excluyendo anchor y
        # direction. El intervalo es el mismo que usa la clasificacion dia/noche
        # de ``_pre_order``: antes uno abria en 0 y el otro lo cerraba, de modo
        # que un pico exactamente en 0 contaba como "dia" pero no computaba
        # aqui.
        peak_d = sum(
            1 for i, g in enumerate(core_genes)
            if g not in (self.anchor_gene, self.direction_gene)
            and 0 <= peaks[i] < pi
        )
        mitad = int(np.floor((n_core - 2) / 2))

        p6am  = 1.0 - np.cos(direction_peak)
        p6pm  = 1.0 - np.cos(direction_peak - pi)
        aviso = (p6am <= 0.1) or (p6pm <= 0.1) or (peak_d < mitad)

        if aviso and consistency_i is not None and anchor_i is not None:
            excl = (peaks[direction_i] < peaks[consistency_i] < peaks[anchor_i])
            flip = not excl
        elif aviso:
            # Sin gen de consistencia la orientacion se invierte por defecto.
            # Se avisa porque es una decision tomada sin evidencia: ocurre, por
            # ejemplo, con el conjunto core de Ruben et al. (2018), que incluye
            # CRY2 pero no CRY1.
            self._log(
                f"  AVISO: se ha producido una situacion de aviso sobre la "
                f"orientacion y el gen de consistencia ({self.consistency_gene}) "
                f"no esta entre los genes core disponibles. La orientacion se "
                f"invierte por defecto, sin verificar la relacion de fases."
            )
            flip = True
        else:
            flip = False

        n = len(o)

        if flip:
            self._log(
                f"    Invirtiendo direccion  "
                f"(aviso=True, {self.direction_gene}_peak={direction_peak:.3f}, "
                f"peak_d={peak_d}/{mitad})"
            )
            peaks_new       = (2.0 * pi - peaks) % (2.0 * pi)
            o_new           = o[::-1].copy()
            esc_new         = (2.0 * pi - esc[::-1]) % (2.0 * pi)
            mat_new         = mat.iloc[:, ::-1].copy()
            # ``iloc`` preserva los nombres de muestra; no reseteamos columnas
            pars_new        = par.copy()
            pars_new[:, 2]  = (2.0 * pi - par[:, 2]) % (2.0 * pi)  # α
            pars_new[:, 3]  = (2.0 * pi - par[:, 3]) % (2.0 * pi)  # β
            ind_new         = np.arange(n - 1, -1, -1)
        else:
            self._log("    La direccion es consistente — no se necesita inversion.")
            peaks_new = peaks.copy()
            o_new     = o.copy()
            esc_new   = esc.copy()
            mat_new   = mat.copy()
            pars_new  = par.copy()
            ind_new   = np.arange(n)

        return o_new, esc_new, mat_new, peaks_new, pars_new, flip, ind_new

    # ------------------------------------------------------------------
    # Utilidad
    # ------------------------------------------------------------------

    def _log(self, message: str) -> None:
        if self.verbose:
            print(message)
