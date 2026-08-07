#!/usr/bin/env python3
"""
experiments.py
==============
Definicion **declarativa** del estudio de simulacion del Capitulo 3.

Este modulo no ejecuta nada: solo describe *que* se contrasta. El runner
(``run_experiment.py``) lo consume y se encarga de generar, ejecutar y
validar. Separarlo permite cambiar el diseno del estudio sin tocar la
maquinaria, y deja el diseno en un unico sitio citable desde la memoria.

Dos planos de parametros
------------------------
El estudio cruza dos cosas de naturaleza distinta que conviene no mezclar:

1. **Factores del dato** (:mod:`generate_synthetic`): definen el escenario
   simulado (``sampling``, ``omega``, ``sigma``, ``prop_rhythmic``,
   ``n_samples``, ``n_total``). No dependen del algoritmo.

2. **Opciones del algoritmo** (flags de ``run_pipeline.py``): son las
   variantes que se contrastan, es decir las aportaciones del TFG frente
   al procedimiento original. Ver :data:`VARIANTS`.

Diseno del estudio
------------------
En lugar de un unico factorial completo (384 celdas x 6 variantes, que no
es viable en tiempo y produce tablas ilegibles), el estudio se organiza en
sub-experimentos dirigidos. Cada aportacion se interroga con los factores
del dato que realmente la ponen a prueba, y cada uno lleva su **control
negativo**, esto es, la celda donde por construccion no deberia haber
diferencia entre variantes.

=====  ==========================================  ==========================
Exp.   Pregunta                                    Subgrid
=====  ==========================================  ==========================
E1     La regularizacion del FMM (Sec. 2.2.1),     ``sampling`` x ``omega``
       concebida contra el muestreo irregular y
       los perfiles apuntados, aporta algo?
E2     La agregacion de ordenes (Sec. 2.2.2)       ``sigma`` x ``n_samples``
       mejora sobre elegir una de las K reps?      x ``prop_rhythmic``
E3     Coste computacional Python vs el CIRCUST    ``n_samples`` x ``n_total``
       original en R, y precision de ambos.
E4     Interaccion entre las dos aportaciones.     solo escenario base
T      Calibracion del umbral tau del criterio     celdas de E2
       ``n_well_fitted`` (hoy provisional 0,65).
=====  ==========================================  ==========================

Uso
---
    from experiments import EXPERIMENTS, VARIANTS

    exp = EXPERIMENTS["E1"]
    exp.cells()          # lista de dicts de factores del dato
    exp.variant_objs()   # lista de Variant
    exp.n_runs(reps=5)   # celdas x variantes x reps
"""
from __future__ import annotations

import itertools
from dataclasses import dataclass, field
from typing import Any, Callable, Iterable


# ===========================================================================
# Escenario base
# ===========================================================================

#: Celda de referencia del estudio. Cada sub-experimento mueve solo los
#: factores que le competen y mantiene el resto en estos valores, de modo
#: que todos los sub-experimentos comparten un punto comun y sus resultados
#: son comparables entre si.
#:
#: Los niveles son los estandar en cronobiologia: muestreo equiespaciado,
#: onda moderadamente asimetrica (omega = 0,5), ruido medio (sigma = 1, es
#: decir ruido de la misma magnitud que la amplitud) y un 5 % de genes
#: ritmicos, que es la proporcion habitual en tejido humano.
BASE_SCENARIO: dict[str, Any] = {
    "sampling":      "equispaced",
    "omega":         0.5,
    "sigma":         1.0,
    "prop_rhythmic": 0.05,
    "n_samples":     48,
    "n_total":       1200,
}


# ===========================================================================
# Variantes del algoritmo
# ===========================================================================

@dataclass(frozen=True)
class Variant:
    """
    Una configuracion del algoritmo a contrastar.

    Attributes
    ----------
    id : str
        Identificador corto, se usa como nombre de carpeta y como clave
        en ``results.csv``. Sin espacios.
    label : str
        Etiqueta legible para las tablas y figuras de la memoria.
    flags : dict
        Argumentos que se pasan a ``run_pipeline.py``. Las claves son los
        nombres largos sin ``--`` (``fmm_method`` -> ``--fmm-method``).
    engine : str
        ``"python"`` para la reimplementacion de este TFG, ``"R"`` para el
        CIRCUST original. El runner despacha a un ejecutable distinto
        segun este campo.
    family : str
        Agrupacion a efectos de analisis: ``"fmm"`` (aportacion 2.2.1),
        ``"order"`` (aportacion 2.2.2), ``"external"`` (referencia).
    """
    id:     str
    label:  str
    flags:  dict[str, Any] = field(default_factory=dict)
    engine: str            = "python"
    family: str            = "order"

    def cli_args(self) -> list[str]:
        """Traduce ``flags`` a la lista de argumentos de ``run_pipeline.py``."""
        args: list[str] = []
        for key, value in self.flags.items():
            args += [f"--{key.replace('_', '-')}", str(value)]
        return args


# --- Aportacion 2.2.1: regularizacion del modelo FMM ----------------------
# Al contrastar el modelo FMM se mantiene FIJA la eleccion del orden final
# (best_k con mediana de R2, que es el criterio por defecto y el mas proximo
# al procedimiento original), de modo que la unica diferencia entre las dos
# variantes sea el procedimiento de ajuste.

_ORDER_HELD_FIXED = {
    "robust_method":       "best_k",
    "selection_criterion": "r2_median",
}

V_FMM_BASE = Variant(
    id     = "fmm_base",
    label  = "FMM base",
    flags  = {"fmm_method": "base", **_ORDER_HELD_FIXED},
    family = "fmm",
)

V_FMM_STAB = Variant(
    id     = "fmm_stabilized",
    label  = "FMM regularizado",
    flags  = {"fmm_method": "stabilized", **_ORDER_HELD_FIXED},
    family = "fmm",
)


# --- Aportacion 2.2.2: eleccion del orden circular optimo ------------------
# Aqui se mantiene FIJO el modelo FMM en 'base', de modo que la unica
# diferencia sea como se resume la coleccion de K ordenaciones. Las tres
# primeras son criterios de *seleccion* (el orden final es una de las K) y
# las tres ultimas construyen un orden *nuevo* por agregacion.

_FMM_HELD_FIXED = {"fmm_method": "base"}

#: Umbral tau por defecto del criterio ``n_well_fitted``. Valor provisional
#: a la espera de la calibracion del experimento :data:`EXP_TAU`.
TAU_DEFAULT: float = 0.65

V_BESTK_R2 = Variant(
    id     = "bestk_r2median",
    label  = "best-K (mediana $R^2$)",
    flags  = {**_FMM_HELD_FIXED,
              "robust_method": "best_k", "selection_criterion": "r2_median"},
    family = "order",
)

V_BESTK_NWF = Variant(
    id     = "bestk_nwellfitted",
    label  = "best-K (genes bien ajustados)",
    flags  = {**_FMM_HELD_FIXED,
              "robust_method": "best_k", "selection_criterion": "n_well_fitted",
              "selection_r2_threshold": TAU_DEFAULT},
    family = "order",
)

V_BESTK_MEDOID = Variant(
    id     = "bestk_medoid",
    label  = "best-K (medoide)",
    flags  = {**_FMM_HELD_FIXED,
              "robust_method": "best_k", "selection_criterion": "medoid"},
    family = "order",
)

V_AGG_TSP3 = Variant(
    id     = "agg_tsp3",
    label  = "Agregacion TSP3",
    flags  = {**_FMM_HELD_FIXED,
              "robust_method": "aggregate", "aggregation": "tsp3"},
    family = "order",
)

V_AGG_HODS = Variant(
    id     = "agg_hods",
    label  = "Agregacion HODs",
    flags  = {**_FMM_HELD_FIXED,
              "robust_method": "aggregate", "aggregation": "hods"},
    family = "order",
)

V_AGG_AUTO = Variant(
    id     = "agg_auto",
    label  = "Agregacion automatica (MSCE)",
    flags  = {**_FMM_HELD_FIXED,
              "robust_method": "aggregate", "aggregation": "auto"},
    family = "order",
)


# --- Referencia externa: el CIRCUST original en R --------------------------
# El procedimiento original no resume las K ordenaciones en un orden de
# consenso: la estimacion de fase por muestra que publica es la del orden
# preliminar (ver la nota de runCIRCUST.R sobre basicOrdRefG2[[2]]). Esa es
# por tanto la salida que se compara, y es exactamente el punto de partida
# que la aportacion 2.2.2 pretende mejorar.

V_R_ORIGINAL = Variant(
    id     = "R_original",
    label  = "CIRCUST original (R)",
    flags  = {},
    engine = "R",
    family = "external",
)


#: Registro de todas las variantes, indexado por ``id``.
VARIANTS: dict[str, Variant] = {
    v.id: v for v in (
        V_FMM_BASE, V_FMM_STAB,
        V_BESTK_R2, V_BESTK_NWF, V_BESTK_MEDOID,
        V_AGG_TSP3, V_AGG_HODS, V_AGG_AUTO,
        V_R_ORIGINAL,
    )
}


# ===========================================================================
# Experimentos
# ===========================================================================

@dataclass(frozen=True)
class Experiment:
    """
    Un sub-experimento: un subgrid de factores del dato x un conjunto de
    variantes del algoritmo.

    Attributes
    ----------
    id, title, question
        Identificacion y pregunta que responde, para los encabezados de
        las tablas de la memoria.
    grid : dict
        Factores que **varian**, con sus niveles. El resto se toma de
        :data:`BASE_SCENARIO`.
    variants : list[str]
        Ids de :data:`VARIANTS` que se evaluan en cada celda.
    negative_control : dict | None
        Celda en la que, por construccion, no se espera diferencia entre
        variantes. Sirve de comprobacion de que el montaje experimental no
        esta inventando efectos. ``None`` si no aplica.
    primary_metrics : tuple[str, ...]
        Metricas de ``metrics.py`` sobre las que se lee el resultado. Se
        listan en orden de importancia: la primera es la que lidera las
        tablas.
    """
    id:               str
    title:            str
    question:         str
    grid:             dict[str, tuple]
    variants:         list[str]
    negative_control: dict[str, Any] | None    = None
    primary_metrics:  tuple[str, ...]          = ("medae_h", "sdae_h", "pct_within_2h")
    notes:            str                      = ""

    def cells(self) -> list[dict[str, Any]]:
        """
        Producto cartesiano de :attr:`grid`, completado con los valores de
        :data:`BASE_SCENARIO` para los factores que no varian. Devuelve
        una lista de dicts, cada uno describiendo por completo un dataset.
        """
        keys = list(self.grid.keys())
        cells = []
        for combo in itertools.product(*(self.grid[k] for k in keys)):
            cell = {**BASE_SCENARIO, **dict(zip(keys, combo))}
            cell["cell_id"] = self.cell_id(cell)
            cells.append(cell)
        return cells

    def cell_id(self, cell: dict[str, Any]) -> str:
        """
        Identificador de carpeta. Solo incluye los factores que varian en
        este experimento, de forma que las rutas queden cortas y legibles.
        """
        parts = []
        for key in self.grid:
            value = cell[key]
            token = value if isinstance(value, str) else f"{value:g}"
            parts.append(f"{_ABBREV.get(key, key)}-{token}")
        return "__".join(parts) if parts else "base"

    def variant_objs(self) -> list[Variant]:
        """Las variantes de este experimento, resueltas a objetos."""
        return [VARIANTS[v] for v in self.variants]

    def is_negative_control(self, cell: dict[str, Any]) -> bool:
        """True si ``cell`` es la celda de control negativo del experimento."""
        if self.negative_control is None:
            return False
        return all(cell.get(k) == v for k, v in self.negative_control.items())

    def n_runs(self, reps: int) -> int:
        """Numero total de ejecuciones del pipeline: celdas x variantes x reps."""
        return len(self.cells()) * len(self.variants) * reps


#: Abreviaturas de los factores para construir nombres de carpeta cortos.
_ABBREV: dict[str, str] = {
    "sampling":      "samp",
    "omega":         "om",
    "sigma":         "sig",
    "prop_rhythmic": "prop",
    "n_samples":     "n",
    "n_total":       "G",
}


# --- E1: regularizacion del modelo FMM ------------------------------------

EXP_FMM = Experiment(
    id       = "E1",
    title    = "Regularizacion del modelo FMM",
    question = ("Aporta la regularizacion (ponderacion por densidad temporal "
                "+ criterio de potencia penalizado) precision frente al ajuste "
                "base cuando el muestreo es irregular o el perfil es apuntado?"),
    grid     = {
        # Irregularidad creciente de la cobertura del ciclo. Es el factor
        # que motiva la ponderacion por densidad temporal von Mises.
        "sampling": ("equispaced", "jitter", "clustered", "clustered_gap"),
        # Apuntamiento de la onda. omega pequeno es donde el criterio de
        # minimos cuadrados puro se vuelve inestable y donde entra la
        # penalizacion por rugosidad. omega = 1 es la sinusoide pura.
        "omega":    (0.1, 0.5, 0.9, 1.0),
    },
    variants = ["fmm_base", "fmm_stabilized"],
    # Muestreo uniforme y onda sinusoidal: ni la ponderacion ni la
    # penalizacion tienen nada que corregir, las dos variantes deberian
    # empatar. Si aqui aparece diferencia, el montaje esta sesgado.
    negative_control = {"sampling": "equispaced", "omega": 1.0},
    primary_metrics  = ("medae_h", "sdae_h", "pct_within_2h", "ccc"),
    notes = ("El orden final se mantiene fijo en best_k/r2_median para que la "
             "unica diferencia entre variantes sea el ajuste del FMM."),
)


# --- E2: eleccion del orden circular optimo -------------------------------

EXP_ORDER = Experiment(
    id       = "E2",
    title    = "Eleccion del orden circular optimo",
    question = ("Construir un orden de consenso por agregacion mejora sobre "
                "seleccionar una de las K repeticiones, y bajo que condiciones "
                "de senal?"),
    grid     = {
        # La ventaja de agregar deberia crecer al bajar la senal, porque es
        # cuando las K repeticiones discrepan mas entre si.
        "sigma":         (0.5, 1.0, 2.0),
        "n_samples":     (24, 48, 75, 96),
        "prop_rhythmic": (0.01, 0.05),
    },
    variants = ["bestk_r2median", "bestk_nwellfitted", "bestk_medoid",
                "agg_tsp3", "agg_hods", "agg_auto"],
    # Senal fuerte, muestra grande y muchos genes ritmicos: las K
    # repeticiones deberian coincidir casi exactamente, de modo que da
    # igual como se resuman. Ahi las 6 variantes deben empatar.
    negative_control = {"sigma": 2.0, "n_samples": 96, "prop_rhythmic": 0.05},
    primary_metrics  = ("medae_h", "sdae_h", "pct_within_2h", "circ_corr", "cre"),
    notes = ("El modelo FMM se mantiene fijo en 'base' para aislar el efecto "
             "de la etapa de resumen de las K ordenaciones."),
)


# --- E3: coste computacional y precision frente a R ------------------------

EXP_COST = Experiment(
    id       = "E3",
    title    = "Coste computacional y precision frente al CIRCUST original",
    question = ("Como escalan en tiempo la reimplementacion en Python y el "
                "codigo original en R con el numero de muestras y de genes, y "
                "que precision alcanza cada uno?"),
    grid     = {
        "n_samples": (24, 48, 75, 96),
        "n_total":   (1200, 5000, 15000),
    },
    variants = ["fmm_base", "fmm_stabilized", "R_original"],
    negative_control = None,
    # Aqui lo que lidera es el tiempo, pero se recogen igualmente las
    # metricas de precision porque el interes es que la aceleracion no
    # venga a costa de peores estimaciones.
    primary_metrics = ("t_pipeline_s", "medae_h", "sdae_h"),
    notes = ("La salida por muestra del CIRCUST original es la del orden "
             "preliminar, que es el punto de partida que mejora la Sec. 2.2.2. "
             "Requiere R con el codigo de rscripts/ operativo."),
)


# --- E4: interaccion entre las dos aportaciones ----------------------------

EXP_INTERACTION = Experiment(
    id       = "E4",
    title    = "Interaccion entre regularizacion FMM y agregacion de ordenes",
    question = ("Las dos aportaciones son aditivas, o el beneficio de una "
                "absorbe al de la otra cuando se aplican juntas?"),
    # Un unico escenario, pero uno exigente: muestreo agrupado con hueco y
    # onda apuntada, que es donde ambas aportaciones deberian tener margen.
    grid     = {
        "sampling": ("clustered_gap",),
        "omega":    (0.1,),
        "sigma":    (0.5,),
    },
    # El cruce 2 x 6 se resuelve en el runner expandiendo estas variantes
    # con los dos modelos FMM (ver _expand_interaction en run_experiment).
    variants = ["bestk_r2median", "bestk_medoid",
                "agg_tsp3", "agg_hods", "agg_auto"],
    negative_control = None,
    primary_metrics  = ("medae_h", "sdae_h", "pct_within_2h"),
    notes = ("El runner cruza estas variantes de orden con fmm-method en "
             "{base, stabilized}, generando el diseno 2 x 5."),
)


# --- T: calibracion del umbral tau ----------------------------------------

#: Rejilla de umbrales tau a barrer. Se excluye 0,5 del rango util porque
#: coincide con el umbral de construccion del propio conjunto TOP, de modo
#: que el criterio satura por construccion (ver Sec. 2.2.2). Se incluye
#: igualmente como extremo para poder mostrar esa saturacion en la memoria.
TAU_GRID: tuple[float, ...] = (0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90)


def tau_variants(taus: Iterable[float] = TAU_GRID) -> list[Variant]:
    """
    Construye una variante ``bestk_nwellfitted`` por cada umbral de
    ``taus``. Se usa en el experimento de calibracion :data:`EXP_TAU`,
    donde el umbral pasa a ser el factor bajo estudio.
    """
    out = []
    for tau in taus:
        out.append(Variant(
            id     = f"bestk_nwf_tau{tau:g}".replace(".", "p"),
            label  = rf"best-K ($\tau$ = {tau:g})",
            flags  = {**_FMM_HELD_FIXED,
                      "robust_method": "best_k",
                      "selection_criterion": "n_well_fitted",
                      "selection_r2_threshold": tau},
            family = "tau",
        ))
    return out


# Las variantes de tau se registran en VARIANTS para que el runner las
# resuelva igual que las demas.
for _v in tau_variants():
    VARIANTS[_v.id] = _v

EXP_TAU = Experiment(
    id       = "T",
    title    = "Calibracion del umbral tau del criterio de genes bien ajustados",
    question = ("Que valor de tau maximiza la precision del criterio "
                "'mayor numero de genes bien ajustados', y como de sensible "
                "es el resultado a esa eleccion?"),
    # Un subgrid reducido de E2: basta con recorrer el eje de senal, que es
    # donde el criterio discrimina. Fijar n y prop reduce el coste sin
    # perder la lectura, porque tau no interactua con el tamano muestral.
    grid     = {
        "sigma":     (0.5, 1.0, 2.0),
        "n_samples": (48, 96),
    },
    variants = [v.id for v in tau_variants()],
    # Se incluye la mediana de R2 como referencia: si ningun tau la bate,
    # el criterio n_well_fitted no aporta y eso tambien es un resultado.
    negative_control = None,
    primary_metrics  = ("medae_h", "sdae_h", "pct_within_2h"),
    notes = ("Cierra el valor provisional 0,65 usado en la Sec. 2.2.2. "
             "tau = 0,5 se incluye solo para evidenciar la saturacion."),
)
# Referencia de comparacion dentro del experimento de calibracion.
EXP_TAU.variants.append("bestk_r2median")


#: Registro de experimentos, indexado por ``id``.
EXPERIMENTS: dict[str, Experiment] = {
    e.id: e for e in (EXP_FMM, EXP_ORDER, EXP_COST, EXP_INTERACTION, EXP_TAU)
}


# ===========================================================================
# Resumen legible (para --list del runner)
# ===========================================================================

def describe(reps: int = 5) -> str:
    """Tabla resumen del estudio completo, con el coste en ejecuciones."""
    lines = [
        f"{'Exp':<5} {'Celdas':>7} {'Variantes':>10} {'Reps':>5} {'Runs':>7}  Titulo",
        "-" * 78,
    ]
    total = 0
    for exp in EXPERIMENTS.values():
        n_cells = len(exp.cells())
        n_var   = len(exp.variants)
        runs    = exp.n_runs(reps)
        total  += runs
        lines.append(
            f"{exp.id:<5} {n_cells:>7} {n_var:>10} {reps:>5} {runs:>7}  {exp.title}"
        )
    lines.append("-" * 78)
    lines.append(f"{'TOTAL':<5} {'':>7} {'':>10} {'':>5} {total:>7}")
    return "\n".join(lines)


if __name__ == "__main__":
    print(describe())
