"""
circust/core_genes.py
=====================
Modulo de seleccion de genes core (reloj circadiano) para el pipeline CIRCUST.

Responsabilidad
---------------
Este modulo resuelve que genes core se usan en la ejecucion y valida
que esos genes esten presentes en la matriz de expresion. Es el unico
punto del pipeline donde se toma esa decision.

Conjuntos predefinidos
----------------------
- "circust"  Larriba et al. 2023  — 12 genes  (defecto de CIRCUST)
- "zhang"    Zhang et al. 2014    — 10 genes
- "ruben"    Ruben et al. 2018    — 54 genes de ciclado ubicuo en humano

Uso tipico
----------
>>> from circust.core_genes import CoreGeneSelector
>>> from circust.preprocessing import Preprocessor, load_expression_matrix
>>>
>>> matrix = load_expression_matrix("data/raw/gtex_brain.csv")
>>> prep   = Preprocessor().run(matrix)
>>>
>>> # Preset por nombre
>>> sel    = CoreGeneSelector(preset="circust")
>>> result = sel.select(prep.expr_norm)
>>> result.genes          # lista validada
>>> result.missing        # genes no encontrados en la matriz
>>>
>>> # Lista personalizada
>>> sel    = CoreGeneSelector(custom_genes=["ARNTL", "DBP", "PER1", "PER2"])
>>> result = sel.select(prep.expr_norm)


Posicion en el pipeline
-----------------------
    Preprocessor  ->  **CoreGeneSelector**  ->  CPCA  ->  ...
"""
from __future__ import annotations

from dataclasses import dataclass, field

import pandas as pd

# ---------------------------------------------------------------------------
# Listas de genes core predefinidas
# ---------------------------------------------------------------------------

SEED_GENES_CIRCUST: list[str] = [
    "PER1", "PER2", "PER3",
    "CRY1", "CRY2",
    "ARNTL", "CLOCK",
    "NR1D1", "RORA",
    "DBP", "TEF", "STAT3",
]

SEED_GENES_ZHANG: list[str] = [
    "ARNTL", "DBP", "NR1D1", "NR1D2", "PER1", "PER2", "PER3",
    "USP2", "TSC22D3", "TSPAN4",
]

# Ruben et al. (2018), los 54 genes de ciclado ubicuo en humano.
#
# Fuente: material suplementario de la publicacion, hoja
# ``DatafileS2_HumanUbiquitousCycle`` de
# ``data/NIHMS1725947-supplement-File_S2 54 cores.xlsx``. Son los genes que
# ciclan en al menos 8 de los 13 tejidos del atlas. El orden es el del propio
# fichero, por numero de tejidos en los que ciclan.
#
# La lista se transcribe aqui y no se lee del Excel en tiempo de ejecucion para
# que importar este modulo no dependa de un fichero de datos ni de openpyxl, y
# para que quede versionada y a la vista.
#
# AVISO: este conjunto no contiene CRY1, que es el gen de consistencia por
# defecto de la etapa de sincronizacion. Con el conjunto de Ruben esa
# verificacion se queda sin tercer gen y, segun la regla, invierte la
# orientacion por defecto cuando se produce un aviso. ``CircularSynchronizer``
# lo notifica por consola cuando ocurre.
SEED_GENES_RUBEN: list[str] = [
    "DBP", "NR1D2", "CIART", "CRY2", "GSTM2", "PER2", "PER3", "TEF", "ARNTL",
    "FAM76A", "HLF", "BHLHE41", "IL11RA", "LDB1", "TCEA3", "TRIM45",
    "ACOX2", "CHD6", "FAM117B", "GLT8D2", "MARVELD1", "NMNAT3", "OLFML1",
    "PER1", "RTTN", "ANKRD6", "C1QTNF2", "CCNH", "CCNL1", "CHIC1", "DHRS12",
    "EDA", "ELL", "ERF", "EVC", "FAM117A", "GSTM5", "IGSF8", "KDM6B",
    "KIAA1328", "NFIL3", "NR3C2", "PARP2", "PPFIBP2", "PRDM5", "PRTFDC1",
    "RBM43", "SEC16B", "SUN1", "TCEA2", "TTLL3", "UBA7", "ZBTB44", "ZSWIM6",
]

SEED_GENES_DEFAULT = SEED_GENES_CIRCUST


# ---------------------------------------------------------------------------
# Conjuntos predefinidos
# ---------------------------------------------------------------------------

PRESETS: dict[str, list[str]] = {
    "circust": SEED_GENES_CIRCUST,
    "zhang":   SEED_GENES_ZHANG,
    "ruben":   SEED_GENES_RUBEN,
}


# ---------------------------------------------------------------------------
# Alias de genes (nombre_canonico → alias_alternativo)
# ---------------------------------------------------------------------------
# Si un dataset usa el alias en vez del canonico, ``apply_gene_aliases`` lo
# renombra. Para anadir un nuevo alias en el futuro, basta con incluir una
# linea aqui (ej: ``"PER1": "RIGUI"``).
GENE_ALIASES: dict[str, str] = {
    "ARNTL": "BMAL1",      # ARNTL = BMAL1 (Aryl hydrocarbon Receptor Nuclear Translocator-Like)
}


def apply_gene_aliases(matrix: pd.DataFrame) -> pd.DataFrame:
    """Renombra alias al nombre canonico en el indice de ``matrix``.

    Por cada par ``canonico ↔ alias``:
      - Si el canonico ya esta presente, no hace nada (los alias se ignoran).
      - Si solo el alias esta, lo renombra al canonico.
    """
    rename = {
        alias: canonical
        for canonical, alias in GENE_ALIASES.items()
        if alias in matrix.index and canonical not in matrix.index
    }
    return matrix.rename(index=rename) if rename else matrix


# ---------------------------------------------------------------------------
# Contenedor de resultados
# ---------------------------------------------------------------------------

@dataclass
class CoreGeneResult:
    """
    Resultado de la seleccion de genes core.

    Attributes
    ----------
    genes : list[str]
        Genes seleccionados y validados (presentes en la matriz de expresion).
        Es el valor que consume CPCA y el resto del pipeline.

    genes_requested : list[str]
        Genes solicitados originalmente, antes de validar contra la matriz.

    missing : list[str]
        Genes solicitados que no se encontraron en la matriz de expresion.
        Puede ser una lista vacia si todos los genes estaban presentes.

    preset : str
        Nombre del preset utilizado (``"circust"``, ``"zhang"``, ``"ruben"``)
        o ``"custom"`` si el usuario paso una lista explicita.
    """

    genes:            list[str]
    genes_requested:  list[str]
    missing:          list[str]   = field(default_factory=list)
    preset:           str         = "circust"

    def summary(self) -> str:
        lines = [
            "=== Seleccion de Genes Core ===",
            f"  Preset           : {self.preset}",
            f"  Solicitados      : {len(self.genes_requested)}",
            f"  Encontrados      : {len(self.genes)}  {self.genes}",
        ]
        if self.missing:
            lines.append(f"  No encontrados   : {self.missing}")
        return "\n".join(lines)


# ---------------------------------------------------------------------------
# Selector
# ---------------------------------------------------------------------------

class CoreGeneSelector:
    """
    Selecciona y valida los genes core del reloj circadiano.

    Parameters
    ----------
    preset : str o None
        Nombre del conjunto predefinido: "circust", "zhang", "ruben".
        Se ignora si se pasa custom_genes.
        Si ambos son None se usa **"circust"** por defecto.

    custom_genes : list[str] o None
        Lista explicita de simbolos genicos. Tiene prioridad sobre "preset".

    verbose : bool
        Si True, imprime el resumen de seleccion al llamar a :meth:`select`.

    Raises
    ------
    ValueError
        Si preset no es un nombre reconocido y no se paso custom_genes.
    """

    PRESETS = PRESETS

    def __init__(
        self,
        preset:       str | None       = "circust",
        custom_genes: list[str] | None = None,
        verbose:      bool             = True,
    ) -> None:
        if custom_genes is not None:
            self._genes_requested = list(custom_genes)
            self._preset_name     = "custom"
        elif preset is not None:
            preset_lower = preset.lower()
            if preset_lower not in self.PRESETS:
                raise ValueError(
                    f"Preset '{preset}' desconocido. "
                    f"Opciones disponibles: {list(self.PRESETS)}"
                )
            # Un preset registrado pero sin genes falla mas tarde, en ``select``,
            # con el mensaje "ningun gen core encontrado en la matriz", que
            # apunta a la matriz de entrada y no a la causa real. Se detecta
            # aqui para que el error senale al preset.
            if not self.PRESETS[preset_lower]:
                raise ValueError(
                    f"El preset '{preset_lower}' esta registrado pero no "
                    f"contiene ningun gen. Completa su lista en core_genes.py "
                    f"o utiliza otro preset: "
                    f"{[k for k, v in self.PRESETS.items() if v]}"
                )
            self._genes_requested = list(self.PRESETS[preset_lower])
            self._preset_name     = preset_lower
        else:
            self._genes_requested = list(SEED_GENES_DEFAULT)
            self._preset_name     = "circust"

        self.verbose = verbose

    # ------------------------------------------------------------------
    # API publica
    # ------------------------------------------------------------------

    def select(self, expr_norm: pd.DataFrame) -> CoreGeneResult:
        """
        Valida los genes solicitados contra la matriz de expresion.

        Parameters
        ----------
        expr_norm : pd.DataFrame
            Matriz de expresion normalizada, genes x muestras.
            Tipicamente PreprocessingResult.expr_norm.

        Returns
        -------
        CoreGeneResult
            result.genes contiene la lista validada lista para CPCA.

        Raises
        ------
        ValueError
            Si ningun gen solicitado esta presente en la matriz.
        """
        genes_in_matrix = set(expr_norm.index)
        found   = [g for g in self._genes_requested if g in genes_in_matrix]
        missing = [g for g in self._genes_requested if g not in genes_in_matrix]

        if not found:
            raise ValueError(
                f"Ningun gen core encontrado en la matriz de expresion.\n"
                f"Solicitados : {self._genes_requested}\n"
                f"Comprueba que los nombres coinciden con el indice de la "
                f"matriz (mayusculas/minusculas incluidas)."
            )

        result = CoreGeneResult(
            genes           = found,
            genes_requested = self._genes_requested,
            missing         = missing,
            preset          = self._preset_name,
        )

        if self.verbose:
            print(result.summary())

        return result

    # ------------------------------------------------------------------
    # Metodo de clase para construir desde un string de CLI
    # ------------------------------------------------------------------

    @classmethod
    def from_string(
        cls,
        raw:     str | None,
        verbose: bool = True,
    ) -> "CoreGeneSelector":
        """
        Construye un :class:`CoreGeneSelector` desde un argumento de linea de comandos.

        Acepta tres formas:
        - ``None``           → preset por defecto (``"circust"``)
        - ``"circust"``      → preset con ese nombre
        - ``"PER1,PER2,..."`` → lista separada por comas

        Parameters
        ----------
        raw : str o None
            Valor del argumento ``--core-genes`` de la CLI.

        verbose : bool
            Propagado al constructor.

        Returns
        -------
        CoreGeneSelector

        Raises
        ------
        ValueError
            Si la lista personalizada tiene menos de 2 genes.
        """
        if raw is None:
            return cls(preset="circust", verbose=verbose)

        if raw.lower() in PRESETS:
            return cls(preset=raw.lower(), verbose=verbose)

        # Lista separada por comas
        genes = [g.strip() for g in raw.split(",") if g.strip()]
        if len(genes) < 2:
            raise ValueError(
                f"--core-genes: se necesitan al menos 2 genes, se recibio: {genes}"
            )
        return cls(custom_genes=genes, verbose=verbose)