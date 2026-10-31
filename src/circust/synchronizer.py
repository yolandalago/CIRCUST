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
* Una lista de pares de genes con orden de pico conocido a priori
  (``consistency_pairs``, por defecto solo ``DBP ≺ CRY1 ≺ ARNTL``) participa
  en una verificacion de consistencia biologica refinada (Paso 2.2).

Paso 2.1 — Rotacion y determinacion de direccion:
    1. Rotar escala circular: peak(anchor) → π.
    2. Si (direction_peak − anchor_peak + π) mod 2π ≥ π → invertir.
    3. Reparametrizar ajustes FMM de genes core en el nuevo marco.
    4. Clasificar genes core como *dia* (pico ∈ [0, π)) o *noche* (pico ∈ [π, 2π)).

Paso 2.2 — Verificacion de consistencia biologica:
    Si direction_gene esta demasiado cerca de 0 o π, o menos de la mitad
    de los genes sin anchor/direction tienen pico en [0, π), se consulta
    ``consistency_pairs``: por cada par ``(g1, g2)`` presente en los genes
    core y con R² >= ``consistency_r2_min`` para ambos, se acumula
    ``sin(peak(g2) − peak(g1))``. Si la suma sale negativa se invierte; si
    no hay ningun par valido, se mantiene la orientacion del Paso 2.1.

    El seno es la version circular de "g2 viene despues de g1", y sustituye
    a la comparacion ``peak(g1) < peak(g2)`` que se usaba antes. Aquella era
    lineal sobre datos circulares y por tanto dependia de donde estuviera el
    corte: con g1 a las 23 h y g2 a la 1 h, g2 pica dos horas despues pero
    ``23 < 1`` es falso. Ademas ``|sin|`` pondera solo cuanta informacion
    aporta cada par, maxima con 6 h de separacion y **cero con 12 h**,
    porque en un ciclo de 24 h "12 horas despues" y "12 horas antes" son la
    misma posicion. Con la regla anterior un par antipodal emitia un voto de
    pleno derecho, y como la reflexion θ → 2π − θ deja quietos exactamente a
    0 y a π, un par anclado ahi seguia "cumpliendose" tras invertir el
    orden. Con dos pares eso empataba la votacion y dejaba pasar la
    inversion, que es la situacion de ARNTL y DBP, casi antipodas.

    Esto generaliza la tripleta fija original (direction < consistency <
    anchor, un unico par encadenado) a una lista configurable de relaciones
    de orden a priori: el propio script R de CIRCUST evaluaba muchas
    tripletas candidatas por tejido (ver ``rscripts/functionGTEX_cores.R``,
    lineas 1907-1990) en vez de confiar en una unica combinacion fija; aqui
    se deja preparado para anadir mas conocimiento biologico segun se valide.

Posicion en el pipeline
-----------------------
    CPCA  →  CircularSynchronizer  →  (Etapa 3)
"""
import numpy as np
import pandas as pd
from dataclasses import dataclass, field
from math import pi

from circust.fitting.fmm import FMMModel
from circust.fitting.rhythm_model import FitResult


# Unica relacion de orden a priori con respaldo directo: es la cadena que ya
# usaba el metodo original (direction < consistency < anchor, con DBP, CRY1 y
# ARNTL). No se anaden mas genes por defecto (RORA, STAT3, CRY2, PER1-3) sin
# verificar antes su relacion de fase en la bibliografia: una relacion a
# priori incorrecta es peor que no tener ninguna.
DEFAULT_CONSISTENCY_PAIRS: list[tuple[str, str]] = [
    ("DBP",  "CRY1"),
    ("CRY1", "ARNTL"),
]


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
        Indices finales de muestra en la matriz de expresion sincronizada.
        
    circular_scale : np.ndarray de float, forma (n_muestras,)
        Eje temporal circular en [0, 2π) tras correccion de direccion.

    expr_ordered : pd.DataFrame, forma (n_genes, n_muestras)
        Matriz de expresion normalizada completa en el orden circular final.

    peak_times : np.ndarray de float, forma (n_genes_core,)
        Tiempo de pico FMM para cada gen central en el marco de coordenadas final.

    r2_fmm : np.ndarray de float, forma (n_genes_core,)
        R² del ajuste FMM para cada gen central.

    fmm_params : np.ndarray de float, forma (n_genes_core, 5)
        Parametros FMM [M, A, α, β, ω] por gen central en el nuevo marco.

    direction_flipped : bool
        True si la direccion fue invertida (por el Paso 2.1 o 2.2).

    day_genes : list[str]
        Genes centrales (exc. anchor y direction) con pico en [0, π).

    night_genes : list[str]
        Genes centrales (exc. anchor y direction) con pico en [π, 2π).

    pre_sample_order : np.ndarray de int
        Orden de muestras del Paso 2.1 (antes de la verificacion del Paso 2.2).

    pre_circular_scale : np.ndarray de float
        Escala circular del Paso 2.1.

    pre_expr_ordered : pd.DataFrame
        Matriz de expresion ordenada por el Paso 2.1.

    pre_peak_times : np.ndarray de float
        Tiempos de pico FMM del Paso 2.1.

    pre_fmm_params : np.ndarray de float, forma (n_genes_core, 5)
        Parametros FMM del Paso 2.1.

    pre_direction_reversed : bool
        True si el Paso 2.1 invirtio la direccion (direction_gene en la mitad incorrecta).
    """

    sample_order:           np.ndarray
    circular_scale:         np.ndarray
    expr_ordered:           pd.DataFrame
    peak_times:             np.ndarray
    r2_fmm:                 np.ndarray
    fmm_params:             np.ndarray
    direction_flipped:      bool

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

    consistency_pairs : list[tuple[str, str]], opcional
        Lista de pares ``(g1, g2)`` con orden de pico conocido a priori
        (se espera que g2 pique despues de g1), usada en la verificacion
        de consistencia biologica del Paso 2.2. Cada par solo cuenta si
        ambos genes estan en los genes core y su R² de ajuste FMM es
        >= ``consistency_r2_min``.

        Cada par aporta ``sin(peak(g2) - peak(g1))`` y se invierte la
        orientacion si la suma sale negativa. El seno es la version
        circular de "viene despues": es positivo cuando g2 esta entre 0 y
        media vuelta por delante de g1, y negativo si esta por detras. Se
        usa en vez de la comparacion ``peak(g1) < peak(g2)`` porque esa es
        lineal sobre datos circulares y depende de donde este el corte:
        con g1 estimado a las 23 h y g2 a la 1 h, g2 pica dos horas
        despues pero ``23 < 1`` es falso.

        La magnitud ``|sin|`` pondera sola cuanta informacion aporta el
        par. Es maxima con 6 h de separacion y **cero con 12 h**, porque
        en un ciclo de 24 h "12 horas despues" y "12 horas antes" son la
        misma posicion y ese par no puede distinguir la direccion. Con la
        regla anterior un par asi emitia un voto de pleno derecho que
        podia empatar la votacion y dejar pasar una inversion, que es lo
        que ocurre en cuanto los dos genes de anclaje son casi antipodas.
        Por defecto: ``DEFAULT_CONSISTENCY_PAIRS`` (``DBP < CRY1 < ARNTL``,
        la unica relacion con respaldo directo del metodo original).

    consistency_r2_min : float
        R² minimo (ambos genes del par) para que un par de
        ``consistency_pairs`` cuente como evidencia valida en el Paso 2.2.
        Por defecto: 0.5.

    verbose : bool
        Imprimir mensajes de progreso.

    Example
    -------
    >>> syncer = CircularSynchronizer()
    >>> result = syncer.run(
    ...     expr_ordered   = cpca.expr_norm_final,
    ...     circular_scale = cpca.circular_scale,
    ...     core_genes     = cpca.core_genes_found,
    ...     fmm_fits       = cpca.fmm_fits_final,
    ... )
    >>> print(result.summary())
    """

    def __init__(
        self,
        anchor_gene:        str  = "ARNTL",
        direction_gene:     str  = "DBP",
        consistency_pairs:  list[tuple[str, str]] | None = None,
        consistency_r2_min: float = 0.5,
        verbose:            bool = True,
    ) -> None:
        self.anchor_gene        = anchor_gene
        self.direction_gene     = direction_gene
        self.consistency_pairs  = (
            list(consistency_pairs) if consistency_pairs is not None
            else list(DEFAULT_CONSISTENCY_PAIRS)
        )
        self.consistency_r2_min = consistency_r2_min
        self.verbose            = verbose

    # ------------------------------------------------------------------
    # API publica
    # ------------------------------------------------------------------

    def run(
        self,
        expr_ordered:   pd.DataFrame,
        circular_scale: np.ndarray,
        core_genes:     list[str],
        fmm_fits:       dict,
    ) -> SynchronizationResult:
        """
        Ejecuta la sincronizacion circular.

        Recibe datos planos y no el resultado de otro paso, de modo que
        cualquier llamante que disponga de una matriz ordenada, su escala y
        los ajustes FMM de unos genes core pueda sincronizar, venga esa
        entrada del CPCA principal o de una submatriz remuestreada como en la
        Etapa 4.

        Parameters
        ----------
        expr_ordered : pd.DataFrame (n_genes, n_muestras)
            Matriz de expresion en el orden circular del paso anterior.
        circular_scale : np.ndarray (n_muestras,)
            Fases del orden, alineadas con las columnas de ``expr_ordered``.
        core_genes : list[str]
            Genes core **presentes** en ``expr_ordered``, no los solicitados.
        fmm_fits : dict[str, FitResult]
            Ajuste FMM de cada gen de ``core_genes``. Los tiempos de pico se
            derivan de aqui (``fit.peak_time``) en vez de recibirse aparte,
            que permitiria que ambos quedasen desincronizados.

        Returns
        -------
        SynchronizationResult
        """
        self._log("=== Etapa 2: Sincronizacion Circular ===")
        core_genes = list(core_genes)
        peak_times = {g: fr.peak_time for g, fr in fmm_fits.items()}

        # ── Paso 2.1: rotacion + determinacion de direccion ───────────────
        self._log("  Paso 2.1 — rotando al marco de referencia ...")
        (o_pre, esc_pre, mat_pre, peaks_pre, r2_fmm, par_pre,
         names_day, names_night, reversed_21) = self._pre_order(
            expr_ordered, circular_scale, fmm_fits, peak_times, core_genes
        )

        # ── Paso 2.2: verificacion de consistencia biologica ──────────────
        self._log("  Paso 2.2 — verificando consistencia biologica ...")
        (o_fin, esc_fin, mat_fin, peaks_fin,
         pars_fin, flipped_22) = self._basic_order(
            o_pre, esc_pre, mat_pre, peaks_pre, par_pre, r2_fmm, core_genes
        )

        direction_flipped = reversed_21 ^ flipped_22 #Se aplica xor

        result = SynchronizationResult(
            sample_order           = o_fin,
            circular_scale         = esc_fin,
            expr_ordered           = mat_fin,
            peak_times             = peaks_fin,
            r2_fmm                 = r2_fmm,
            fmm_params             = pars_fin,
            direction_flipped      = direction_flipped,
            day_genes              = names_day,
            night_genes            = names_night,
            pre_sample_order       = o_pre,
            pre_circular_scale     = esc_pre,
            pre_expr_ordered       = mat_pre,
            pre_peak_times         = peaks_pre,
            pre_fmm_params         = par_pre,
            pre_direction_reversed = reversed_21,
        )

        # El resumen se construye solo si se va a imprimir: la Etapa 4 llama
        # aqui una vez por repeticion y por reintento.
        if self.verbose:
            self._log(result.summary())
        return result

    # ------------------------------------------------------------------
    # Pasos internos
    # ------------------------------------------------------------------

    def _pre_order(
        self,
        expr_ordered:   pd.DataFrame,
        circular_scale: np.ndarray,
        fmm_fits:       dict,
        peak_times_fin: dict,
        core_genes:     list[str],
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
        r2_fmm:     np.ndarray,
        core_genes: list[str],
    ) -> tuple:
        """
        Paso 2.2 .

        Returns
        -------
        o_new, esc_new, mat_new, peaks_new, pars_new, flipped
        """
        direction_i = core_genes.index(self.direction_gene) if self.direction_gene in core_genes else None

        if direction_i is None:
            self._log(
                f"  AVISO: {self.direction_gene} no encontrado en core_genes; "
                "omitiendo verificacion de direccion."
            )
            return o, esc, mat, peaks, par, False

        direction_peak = peaks[direction_i]
        n_core         = len(core_genes)

        # Contar genes con pico en el semicirculo [0, pi), excluyendo anchor y
        # direction. El intervalo es el mismo que usa la clasificacion dia/noche
        # antes uno abria en 0 y el otro lo cerraba, de modo
        # que un pico exactamente en 0 contaba como "dia" pero no computaba
        # aqui.
        peak_d = sum(
            1 for i, g in enumerate(core_genes)
            if g not in (self.anchor_gene, self.direction_gene)
            and 0 <= peaks[i] < pi
        )
        mitad = int(np.floor((n_core - 2) / 2))

        # Diagnostico de si el Paso 2.1 tenia base para decidir. Ya no decide
        # si se consulta el conocimiento a priori, solo se registra.
        p6am  = 1.0 - np.cos(direction_peak)
        p6pm  = 1.0 - np.cos(direction_peak - pi)
        aviso = (p6am <= 0.1) or (p6pm <= 0.1) or (peak_d < mitad)

        gene_idx = {g: i for i, g in enumerate(core_genes)}
        score    = 0.0     # suma de sin(delta), su signo decide
        evidence = 0.0     # suma de |sin(delta)|, cuanta informacion habia
        n_valid  = 0
        for g_before, g_after in self.consistency_pairs:
            i_before = gene_idx.get(g_before)
            i_after  = gene_idx.get(g_after)
            if i_before is None or i_after is None:
                continue
            if (r2_fmm[i_before] < self.consistency_r2_min
                    or r2_fmm[i_after] < self.consistency_r2_min):
                continue
            n_valid  += 1
            aporte    = np.sin(peaks[i_after] - peaks[i_before])
            score    += aporte
            evidence += abs(aporte)

        if n_valid > 0:
            # bool de Python y no np.bool_: ``direction_flipped`` esta
            # documentado como bool y aguas abajo se combina con XOR.
            flip = bool(score < 0.0)
            self._log(
                f"    Consistencia a priori: score={score:+.3f} sobre "
                f"{n_valid} pares (evidencia {evidence:.3f} de {n_valid:.1f} "
                f"posible, aviso del Paso 2.1={aviso})."
            )
        else:
            # Sin ningun par valido (genes ausentes o mal ajustados) no hay
            # evidencia a priori con la que decidir, asi que se mantiene la
            # orientacion que ya establecio el Paso 2.1 en vez de invertir a
            # ciegas. Ocurre, por ejemplo, con el conjunto core de Ruben et
            # al. (2018), que incluye CRY2 pero no CRY1, o cuando ninguno de
            # los dos genes de un par cayo en el subconjunto muestreado de
            # una repeticion de la Etapa 4.
            self._log(
                f"  AVISO: ningun par de consistency_pairs esta "
                f"disponible/bien ajustado (R2 >= {self.consistency_r2_min}) "
                f"entre los genes core. Se mantiene la orientacion del "
                f"Paso 2.1, decidida solo con {self.direction_gene}."
            )
            flip = False

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
            pars_new        = par.copy()
            pars_new[:, 2]  = (2.0 * pi - par[:, 2]) % (2.0 * pi)  # α
            pars_new[:, 3]  = (2.0 * pi - par[:, 3]) % (2.0 * pi)  # β
        else:
            self._log("    La direccion es consistente — no se necesita inversion.")
            peaks_new = peaks.copy()
            o_new     = o.copy()
            esc_new   = esc.copy()
            mat_new   = mat.copy()
            pars_new  = par.copy()

        return o_new, esc_new, mat_new, peaks_new, pars_new, flip

    # ------------------------------------------------------------------
    # Utilidad
    # ------------------------------------------------------------------

    def _log(self, message: str) -> None:
        if self.verbose:
            print(message)
