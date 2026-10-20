# preprocessing.py  ── seccion de carga y preprocesamiento
import pathlib
from typing import Optional
from dataclasses import dataclass, field
from scipy.stats import median_abs_deviation

import pandas as pd
import numpy as np


# ── formatos soportados ───────────────────────────────────────────────────
_LOADERS = {
    ".csv":     "csv",
    ".tsv":     "tsv",
    ".txt":     "tsv",       # se asume separado por tabuladores cuando .txt
    ".xlsx":    "excel",
    ".xls":     "excel",
    ".parquet": "parquet",
    ".gct":     "gct",       # Gene Cluster Text (especificacion GenePattern/GSEA)
    ".gct.gz":  "gct",       # GCT comprimido con gzip
}


def _detect_format(file_path: pathlib.Path) -> tuple[str | None, str | None]:
    """Detecta el formato y la extension efectiva.

    Maneja extensiones compuestas como ``.gct.gz`` que ``pathlib.suffix`` no
    captura por si solo (devolveria ``.gz``).

    Returns
    -------
    (extension, fmt) : (str, str) o (None, None) si no se reconoce.
    """
    name_lower = file_path.name.lower()
    # Comprueba primero las extensiones compuestas (mas largas)
    for ext in sorted(_LOADERS.keys(), key=len, reverse=True):
        if name_lower.endswith(ext):
            return ext, _LOADERS[ext]
    return None, None


def load_expression_matrix(
    path: str,
    gene_column: Optional[str] = None,
    chunk_size: Optional[int] = None,
) -> pd.DataFrame:
    """
    Carga una matriz de expresion genica desde disco en un DataFrame de pandas.

    Forma esperada tras la carga
    ----------------------------
    - Filas    = genes  (indice de fila = simbolos genicos, ej. "ARNTL", "PER1")
    - Columnas = muestras (cada columna = un individuo/muestra)
    - Valores  = conteos de expresion brutos (float); ceros y NaN permitidos

    Formatos soportados
    -------------------
    .csv            Valores separados por comas
    .tsv / .txt     Valores separados por tabuladores
    .xlsx / .xls    Libro Excel (se lee la primera hoja)
    .parquet        Apache Parquet — recomendado para conjuntos grandes
    .gct            Gene Cluster Text (formato GenePattern/GSEA): cabecera de
                    2 lineas (#1.2 + dimensiones), columnas Name+Description
                    seguidas de las muestras. Description pasa a indice de fila;
                    Name se descarta.
    .gct.gz         Igual que .gct pero comprimido con gzip.

    Parameters
    ----------
    path : str
        Ruta al archivo. El formato se detecta por la extension.

    gene_column : str, opcional
        Nombre de la columna que contiene los simbolos genicos.
        Especificar explicitamente si el archivo tiene un layout no estandar,
        ej. gene_column="gene_id".

    chunk_size : int, opcional
        Solo para archivos .parquet. Numero de row-groups a leer por vez.
        Util cuando el archivo es mas grande que la RAM disponible.
        Si None (defecto), se carga el archivo completo de una vez.
        Ignorado para CSV/TSV/Excel.

    Returns
    -------
    pd.DataFrame
        Matriz con simbolos genicos como indice de fila (index.name = "gene_symbol")
        e IDs de muestra como nombres de columna. Todos los valores son float64.
        Celdas que no se pudieron parsear como numeros se establecen a NaN.

    Raises
    ------
    FileNotFoundError
        Si no existe archivo en ``path``.
    ValueError
        Si la extension del archivo no esta en la lista soportada.
    ImportError
        Si se lee un .parquet y pyarrow no esta instalado.

    Example
    -------
    Cargar un CSV estandar:

    >>> df = load_expression_matrix("data/raw/gtex_brain.csv")
    >>> df.shape
    (56200, 479)
    >>> df.index[:3].tolist()
    ['ARNTL', 'PER1', 'CRY1']

    Cargar un parquet grande con lectura por bloques:

    >>> df = load_expression_matrix(
    ...     "data/raw/gtex_brain.parquet",
    ...     chunk_size=10_000,
    ... )
    """
    file_path = pathlib.Path(path)

    # ── verificacion de existencia ────────────────────────────────────────
    if not file_path.exists():
        raise FileNotFoundError(
            f"Archivo de matriz de expresion no encontrado: '{path}'\n"
            f"Verifica la ruta y asegurate de que el archivo esta en la carpeta data/."
        )

    # ── deteccion de formato (incluye extensiones compuestas .gct.gz) ─────
    extension, fmt = _detect_format(file_path)

    if fmt is None:
        supported = ", ".join(_LOADERS.keys())
        raise ValueError(
            f"Formato de archivo no soportado: '{file_path.suffix}'\n"
            f"Formatos soportados: {supported}"
        )

    # ── que columna se convierte en el indice de fila? ────────────────────
    # index_col=0 significa "usar la primera columna como indice de fila"
    # Layout estandar: primera columna = nombres de genes,
    # columnas restantes = muestras.
    # Si el usuario conoce su columna de genes por nombre, usar esa.
    if fmt in ("csv", "tsv", "excel"):
        index_col = gene_column if gene_column is not None else 0
    elif fmt == "parquet":
        if gene_column is None:
            raise ValueError(
                "Para archivos Parquet debes especificar el nombre de la columna "
                "que contiene los identificadores de genes mediante el parametro 'gene_column'."
            )
        index_col = gene_column
    elif fmt == "gct":
        # GCT tiene formato fijo: ``Description`` se usa como indice; ``gene_column``
        # se ignora porque la especificacion del formato lo determina.
        index_col = "Description"
    else:
        # no deberia ocurrir porque ya validamos fmt
        index_col = None

    # ── carga ─────────────────────────────────────────────────────────────
    if fmt == "csv":
        matrix = pd.read_csv(file_path, index_col=index_col)

    elif fmt == "tsv":
        matrix = pd.read_csv(file_path, sep="\t", index_col=index_col)

    elif fmt == "excel":
        matrix = pd.read_excel(file_path, index_col=index_col)

    elif fmt == "parquet":
        matrix = _load_parquet(file_path, index_col=index_col,
                               chunk_size=chunk_size)

    elif fmt == "gct":
        matrix = _load_gct(file_path, gzipped=extension.endswith(".gz"))

    # ── forzar valores numericos ──────────────────────────────────────────
    # Algunos Excel o CSV mal formateados contienen cadenas sueltas.
    # errors="coerce" las convierte en NaN en lugar de fallar.
    # Los archivos Parquet ya estan tipados asi que esto es un no-op para ellos.
    matrix = matrix.apply(pd.to_numeric, errors="coerce")

    # ── nombrar el indice ─────────────────────────────────────────────────
    matrix.index.name = "gene_symbol"

    return matrix


def _load_parquet(
    file_path: pathlib.Path,
    index_col,
    chunk_size: Optional[int],
) -> pd.DataFrame:
    try:
        import pyarrow.parquet as pq
    except ImportError:
        raise ImportError(
            "Leer archivos parquet requiere pyarrow.\n"
            "Instalalo con:  pip install pyarrow"
        )

    parquet_file = pq.ParquetFile(file_path)

    if chunk_size is None:
        # Usar pandas directamente — restaura correctamente el indice guardado
        matrix = pd.read_parquet(file_path)
    else:
        chunks = []
        for batch in parquet_file.iter_batches(batch_size=chunk_size):
            chunks.append(batch.to_pandas())
        matrix = pd.concat(chunks, ignore_index=False)

    # Solo establecer indice manualmente si no se restauro automaticamente
    # (es decir, el archivo se escribio SIN un indice con nombre)
    if matrix.index.name != "gene_symbol":
        if isinstance(index_col, str) and index_col in matrix.columns:
            matrix = matrix.set_index(index_col)
        elif isinstance(index_col, int) and index_col < len(matrix.columns):
            matrix = matrix.set_index(matrix.columns[index_col])

    return matrix


def _load_gct(file_path: pathlib.Path, gzipped: bool) -> pd.DataFrame:
    """Carga un archivo GCT (Gene Cluster Text) con o sin compresion gzip.

    Estructura del formato (especificacion GenePattern/GSEA):
      Linea 1:  ``#1.2`` (version)
      Linea 2:  ``<n_genes>\\t<n_samples>`` (dimensiones)
      Linea 3:  cabecera ``Name\\tDescription\\t<sample_1>\\t...\\t<sample_N>``
      Resto:    una fila por gen con valores numericos.

    El indice resultante es la columna ``Description`` (simbolo del gen). La columna
    ``Name`` se descarta porque es el ensemble id y no es una muestra.
    """
    compression = "gzip" if gzipped else None
    # skiprows=2: saltar la linea de version y la de dimensiones
    matrix = pd.read_csv(
        file_path,
        sep         = "\t",
        skiprows    = 2,
        index_col   = "Description",
        compression = compression,
    )
    # Descartar la columna Name si esta presente
    if "Name" in matrix.columns:
        matrix = matrix.drop(columns=["Name"])
    return matrix


@dataclass
class PreprocessingResult:
    """
    Todos los resultados producidos por :class:`Preprocessor`.

    Attributes
    ----------
    expr_norm : pd.DataFrame
        Matriz de expresion normalizada, genes x muestras, valores en [-1, 1].
        Es la matriz que consume cada paso (CPCA, FMM).

    expr_raw : pd.DataFrame
        Matriz bruta filtrada despues de eliminar genes malos e imputar los
        valores ausentes, pero antes de normalizar. Util para depuracion y
        validacion.

    dropped_sparse : List[str]
        Simbolos genicos eliminados porque > zero_threshold de muestras eran cero o NaN.

    dropped_duplicates : List[str]
        Simbolos genicos eliminados durante la resolucion de duplicados (filas con menor MAD).

    dropped_constant : List[str]
        Simbolos genicos eliminados por ser constantes (un unico valor distinto
        entre todas las muestras) tras la imputacion del paso 4.

    imputed_genes : List[str]
        Simbolos genicos que tenian algun valor ausente y han sido imputados.
        Estos genes superaron el filtro de dispersion, es decir, su fraccion de
        NaN no llega al umbral, de modo que se conservan.

    n_imputed_cells : int
        Numero total de celdas imputadas en toda la matriz.

    dropped_samples : list
        Columnas descartadas por no contener ningun valor numerico. Suelen ser
        columnas de metadatos interpretadas como muestras al cargar.

    n_genes_in : int
        Numero de genes que entraron al preprocesador.

    n_genes_out : int
        Numero de genes que sobrevivieron a todos los pasos de filtrado.

    n_samples : int
        Numero de muestras. No cambia nunca durante el preprocesamiento.
    """

    expr_norm:           pd.DataFrame
    expr_raw:            pd.DataFrame
    dropped_sparse:      list[str] = field(default_factory=list)
    dropped_duplicates:  list[str] = field(default_factory=list)
    dropped_constant:    list[str] = field(default_factory=list)
    imputed_genes:       list[str] = field(default_factory=list)
    n_imputed_cells:     int = 0
    dropped_samples:     list      = field(default_factory=list)
    n_genes_in:          int = 0
    n_genes_out:         int = 0
    n_samples:           int = 0

    def summary(self) -> str:
        """Devuelve un resumen legible de lo ocurrido durante el preprocesamiento."""
        lines = [
            "=== Resumen del Preprocesamiento ===",
            f"  Genes entrada        : {self.n_genes_in}",
            f"  Genes salida         : {self.n_genes_out}",
            f"  Muestras             : {self.n_samples}",
            f"  Eliminados (sparse)  : {len(self.dropped_sparse)}",
            f"  Eliminados (duplicados): {len(self.dropped_duplicates)}",
            f"  Eliminados (constantes): {len(self.dropped_constant)}",
            f"  Genes imputados      : {len(self.imputed_genes)} "
            f"({self.n_imputed_cells} celda(s))",
        ]
        if self.dropped_samples:
            lines.append(
                f"  Columnas descartadas : {len(self.dropped_samples)} "
                f"{self.dropped_samples[:5]}"
            )
        if self.dropped_sparse:
            shown = self.dropped_sparse[:5]
            tail  = f" ... (+{len(self.dropped_sparse)-5} mas)" if len(self.dropped_sparse) > 5 else ""
            lines.append(f"    genes sparse       : {shown}{tail}")
        if self.dropped_duplicates:
            unique_dupes = sorted(set(self.dropped_duplicates))[:5]
            tail = f" ... (+{len(set(self.dropped_duplicates))-5} mas)" if len(set(self.dropped_duplicates)) > 5 else ""
            lines.append(f"    genes duplicados   : {unique_dupes}{tail}")
        if self.dropped_constant:
            shown = self.dropped_constant[:5]
            tail  = (f" ... (+{len(self.dropped_constant)-5} mas)"
                     if len(self.dropped_constant) > 5 else "")
            lines.append(f"    genes constantes   : {shown}{tail}")
        if self.imputed_genes:
            shown = self.imputed_genes[:5]
            tail  = (f" ... (+{len(self.imputed_genes)-5} mas)"
                     if len(self.imputed_genes) > 5 else "")
            lines.append(f"    genes imputados    : {shown}{tail}")
        return "\n".join(lines)


def normalise_matrix(
    mat: pd.DataFrame,
    norm_min: float = -1.0,
    norm_max: float =  1.0,
) -> pd.DataFrame:
    """
    Normaliza min-max cada fila de ``mat`` al intervalo [norm_min, norm_max].

    Filas constantes (min == max) se mapean a todos ceros.

    Tratamiento de valores ausentes
    -------------------------------
    El minimo y el maximo se calculan con ``nanmin`` y ``nanmax``, es decir,
    sobre los valores observados de cada gen. Con ``min`` y ``max`` propagantes
    un unico NaN en la fila convertia ``row_min`` y ``row_max`` en NaN y, a
    traves de la division por ``span``, arrasaba las demas celdas de ese gen,
    que eran validas. El gen sobrevivia entonces al filtro del 30 % de
    ``Preprocessor`` y moria en este paso: con el modelo FMM base su R2 salia
    NaN y quedaba descartado en silencio por ``r2 >= umbral``, y si era un
    gen core el SVD del CPCA no convergia. La implementacion de referencia en R arrastra el
    mismo comportamiento en ``normalice2``.

    Las celdas que ya eran NaN siguen siendolo a la salida, porque normalizar
    no puede inventar el valor que falta. Quien decide que hacer con ellas es
    ``Preprocessor``, que las imputa antes de llegar aqui.

    Parameters
    ----------
    mat : pd.DataFrame
        Matriz genes × muestras. Valores pueden ser float o int.
    norm_min, norm_max : float
        Extremos del rango de salida. Por defecto [-1, 1].

    Returns
    -------
    pd.DataFrame con la misma forma e indices que ``mat``.
    """
    values = mat.values.astype(float)

    # Filas enteramente NaN: nanmin/nanmax emitirian un RuntimeWarning y
    # devolverian NaN. Se tratan como filas constantes y se mapean a ceros.
    all_nan = np.all(np.isnan(values), axis=1, keepdims=True)

    with np.errstate(invalid="ignore", divide="ignore"):
        safe    = np.where(all_nan, 0.0, values)
        row_min = np.nanmin(safe, axis=1, keepdims=True)
        row_max = np.nanmax(safe, axis=1, keepdims=True)
        span    = row_max - row_min
        center = (norm_min + norm_max) / 2.0
        normed = np.where(
            (span == 0) | all_nan,
            center,
            (norm_max - norm_min) * ((values - row_min) / span) + norm_min,
        )
    return pd.DataFrame(normed, index=mat.index, columns=mat.columns)


class Preprocessor:
    """
    Limpia y normaliza una matriz de expresion genica bruta.

    1. Eliminar genes sin nombre (NA rownames en R)
    2. Eliminar genes sparse (> ``sparse_threshold`` de muestras son cero O NaN)
    3. Resolver nombres de genes duplicados — conservar la fila con mayor MAD
    4. Imputar los valores ausentes que quedan
    4.5. Eliminar genes constantes (varianza cero tras la imputacion)
    5. Normalizar cada gen a [-1, 1]

    El paso 4 existe porque el filtro del paso 2 tolera hasta un
    ``sparse_threshold`` de NaN por gen, pero ni el SVD del CPCA ni el ajuste
    del modelo FMM admiten valores ausentes. Sin el, un gen que supera el
    filtro con un solo NaN acaba descartado en silencio o hace fallar el CPCA
    si pertenece al conjunto \textit{core}.

    Parameters
    ----------
    sparse_threshold : float
        Fraccion de muestras que pueden ser cero o NaN antes de que un gen
        se elimine. La comparacion es estricta (>), asi que un gen
        exactamente en el umbral se conserva.

    impute_strategy : str
        Que hacer con los valores ausentes de los genes que superan el filtro
        de dispersion.

        ``"median"`` (defecto)
            Sustituirlos por la mediana de los valores observados del gen. Se
            usa la mediana y no la media porque es robusta a los valores
            extremos, frecuentes en expresion genica.
        ``"drop"``
            Eliminar el gen entero. Equivale a exigir un 0 % de valores
            ausentes y hace inoperante la tolerancia de ``sparse_threshold``
            sobre los NaN.

    zero_threshold : float o None
        Override para el umbral solo de ceros. Si None (defecto), usa
        ``sparse_threshold`` para ceros.

    nan_threshold : float o None
        Override para el umbral solo de NaN. Si None (defecto), usa
        ``sparse_threshold`` para NaN.

    verbose : bool
        Si True, imprime mensajes de progreso tras cada paso.

    Example
    -------
    >>> import pandas as pd
    >>> import numpy as np
    >>>
    >>> rng = np.random.default_rng(0)
    >>> mat = pd.DataFrame(
    ...     rng.poisson(10, size=(200, 60)).astype(float),
    ...     index=[f"GENE_{i}" for i in range(200)],
    ...     columns=[f"sample_{j}" for j in range(60)],
    ... )
    >>>
    >>> result = Preprocessor().run(mat)
    >>> result.expr_norm.values.min() >= -1.0
    True
    >>> result.expr_norm.values.max() <= 1.0
    True
    """
    def __init__(
            self,
            sparse_threshold: float         = 0.3,
            zero_threshold:   Optional[float] = None,
            nan_threshold:    Optional[float] = None,
            norm_min:         float         = -1.0,
            norm_max:         float         =  1.0,
            verbose:          bool          = True,
    ) -> None:

        if not 0.0 < sparse_threshold < 1.0:
            raise ValueError(
                f"sparse_threshold debe estar entre 0 y 1, se recibio {sparse_threshold}"
            )
        if zero_threshold is not None:
            if not 0.0 < zero_threshold < 1.0:
                raise ValueError(
                    f"zero_threshold debe estar entre 0 y 1, se recibio {zero_threshold}"
                )
        else:
            zero_threshold = sparse_threshold

        if nan_threshold is not None:
            if not 0.0 < nan_threshold < 1.0:
                raise ValueError(
                   f"nan_threshold debe estar entre 0 y 1, se recibio {nan_threshold}"
                )
        else:
            nan_threshold = sparse_threshold
        
        self.sparse_threshold = sparse_threshold
        self.zero_threshold   = zero_threshold 
        self.nan_threshold    = nan_threshold
        self.norm_min         = norm_min
        self.norm_max         = norm_max
        self.verbose          = verbose

    def run(self, matrix: pd.DataFrame) -> PreprocessingResult:
        """
        Ejecuta los cuatro pasos de preprocesamiento en orden.

        Parameters
        ----------
        matrix : pd.DataFrame
            Matriz de expresion bruta cargada por ``load_expression_matrix()``.
            Filas = genes, columnas = muestras. Puede contener ceros y NaN.

        Returns
        -------
        PreprocessingResult
        """

        if matrix.empty:
            raise ValueError("La matriz de entrada esta vacia.")

        n_genes_in, n_samples = matrix.shape
        self._log(f"Matriz de entrada: {n_genes_in} genes x {n_samples} muestras")

        # paso 0 — descartar columnas enteramente ausentes
        matrix, dropped_cols = self._drop_empty_samples(matrix)
        n_samples = matrix.shape[1]

        # paso 1 — eliminar genes sin nombre
        mat = self._drop_unnamed(matrix)

        # paso 2 — eliminar genes sparse
        mat, dropped_sparse = self._drop_sparse(mat)

        # paso 3 — resolver nombres de genes duplicados
        mat, dropped_duplicates = self._resolve_duplicates(mat)

        # paso 4 — imputar los valores ausentes que han superado el filtro
        mat, imputed_genes, n_imputed = self._impute_missing(mat)

        # paso 4.5 — descartar genes constantes (varianza cero tras imputar)
        mat, dropped_constant = self._delete_constant(mat)

        # paso 5 — normalizar cada gen a [-1, 1]
        expr_norm = normalise_matrix(mat, self.norm_min, self.norm_max)

        result = PreprocessingResult(
            expr_norm          = expr_norm,
            expr_raw           = mat,
            dropped_sparse     = dropped_sparse,
            dropped_duplicates = dropped_duplicates,
            dropped_constant   = dropped_constant,
            imputed_genes      = imputed_genes,
            n_imputed_cells    = n_imputed,
            dropped_samples    = dropped_cols,
            n_genes_in         = n_genes_in,
            n_genes_out        = len(expr_norm),
            n_samples          = n_samples,
        )

        self._log(result.summary())
        return result

    def _drop_empty_samples(
        self, mat: pd.DataFrame,
    ) -> tuple[pd.DataFrame, list]:
        """
        Descarta las columnas cuyos valores son todos ausentes.
        """
        empty = mat.columns[mat.isna().all(axis=0)]
        if len(empty) == 0:
            return mat, []

        self._log(
            f"  Paso 0 — descartadas {len(empty)} columna(s) sin ningun valor "
            f"numerico: {list(empty)[:5]}"
            + (" ..." if len(empty) > 5 else "")
            + ". Revisa si son columnas de metadatos mal interpretadas como "
              "muestras."
        )
        return mat.loc[:, ~mat.columns.isin(empty)], list(empty)

    def _drop_unnamed(self, mat: pd.DataFrame) -> pd.DataFrame:
        """
        Elimina filas cuyo gen falta o esta efectivamente vacio.
        - NaN real en el indice
        - Cadena vacia ``""``
        - Cadena de solo espacios ``"   "``
        - La cadena literal ``"nan"`` (ocurre cuando NaN se convierte a str)
        """
        index_str = mat.index.astype(str)
        is_invalid = (
            mat.index.isna()                          # NaN real
            | (index_str.str.strip() == "")           # vacio / solo espacios
            | (index_str.str.lower() == "nan")        # NaN convertido a cadena
        )

        n_bad = int(is_invalid.sum())
        if n_bad > 0:
            self._log(f"  Paso 1 — eliminados {n_bad} gen(es) sin nombre")
        else:
            self._log("  Paso 1 — no se encontraron genes sin nombre")

        return mat.loc[~is_invalid]

    def _drop_sparse(self, mat: pd.DataFrame) -> tuple[pd.DataFrame, list[str]]:
        """
        Elimina genes donde mas de zero_threshold de muestras son cero O NaN.

        Nota: la condicion es estrictamente > asi que un gen exactamente en el
        30% se conserva. Preservado aqui: (frac > threshold), no (frac >= threshold).
        """
        # fraccion de ceros y NaN por gen (por fila)
        zero_frac = (mat == 0).mean(axis=1)
        nan_frac  = mat.isna().mean(axis=1)

        is_sparse = (zero_frac > self.zero_threshold) | (nan_frac > self.nan_threshold)

        dropped = mat.index[is_sparse].tolist()

        # replica el print exacto de R (linea 3925)
        self._log(
            f"  Paso 2 — eliminados {len(dropped)} gen(es) con mas del "
            f"{int(self.zero_threshold * 100)}% de ceros o NaN"
        )

        return mat.loc[~is_sparse].copy(), dropped

    def _impute_missing(
        self, mat: pd.DataFrame,
    ) -> tuple[pd.DataFrame, list[str], int]:
        """
        Resuelve los valores ausentes de los genes que superaron el filtro
        de dispersion mediante el valor mediano.

        Devuelve la matriz, los genes afectados y el numero de celdas
        imputadas.
        """
        na_per_gene = mat.isna().sum(axis=1)
        affected    = na_per_gene[na_per_gene > 0]

        if affected.empty:
            self._log("  Paso 4 — sin valores ausentes que imputar")
            return mat, [], 0

        genes     = affected.index.tolist()
        n_imputed = int(affected.sum())

        # Mediana por gen sobre los valores observados. Se opera sobre el array
        # y no con ``fillna``, que exigiria transponer y realinear por etiqueta
        # con decenas de miles de genes.
        values = mat.values.astype(float)
        mask   = np.isnan(values)
        with np.errstate(invalid="ignore"):
            medians = np.nanmedian(np.where(mask.all(axis=1, keepdims=True), 0.0, values),
                                   axis=1, keepdims=True)
        values[mask] = np.broadcast_to(medians, values.shape)[mask]
        filled = pd.DataFrame(values, index=mat.index, columns=mat.columns)

        self._log(
            f"  Paso 4 — imputadas {n_imputed} celda(s) en {len(genes)} gen(es) "
            f"con la mediana del gen"
        )
        return filled, genes, n_imputed

    def _resolve_duplicates(self, mat: pd.DataFrame) -> tuple[pd.DataFrame, list[str]]:
        """
        Para cada nombre de gen que aparece mas de una vez, conserva la fila
        con mayor MAD y descarta el resto.

        El MAD se calcula con ``nan_policy="omit"``. Con la politica por
        defecto, ``propagate``, una fila con algun valor ausente daba MAD igual
        a NaN y ``np.argmax`` devuelve la posicion del primer NaN, de modo que
        entre dos filas homonimas se conservaba precisamente la incompleta.
        """
        # encontrar nombres de genes que aparecen mas de una vez
        duplicated_mask  = mat.index.duplicated(keep=False)
        duplicated_names = mat.index[duplicated_mask].unique().tolist()

        if not duplicated_names:
            self._log("  Paso 3 — no se encontraron nombres de genes duplicados")
            return mat, []

        self._log(
            f"  Paso 3 — resolviendo {len(duplicated_names)} "
            f"nombre(s) de gen duplicado(s)"
        )

        # recopilar posiciones enteras de filas a eliminar
        rows_to_drop:  list[int] = []
        dropped_names: list[str] = []

        for gene in duplicated_names:
            # np.where da posiciones enteras necesario porque dos filas
            # comparten la misma etiqueta, asi que .loc devolveria ambas
            positions = np.where(mat.index == gene)[0]

            mads = np.array([
                float(median_abs_deviation(mat.iloc[pos].values, nan_policy="omit"))
                for pos in positions
            ])
            mads = np.nan_to_num(mads, nan=-np.inf)   # fila totalmente NaN: la ultima opcion

            # conservar la fila con mayor MAD, eliminar el resto
            best_position   = positions[int(np.argmax(mads))]
            loser_positions = [p for p in positions if p != best_position]

            rows_to_drop.extend(loser_positions)
            dropped_names.extend([gene] * len(loser_positions))

        # construir mascara booleana de conservacion y aplicar
        keep_mask = np.ones(len(mat), dtype=bool)
        keep_mask[rows_to_drop] = False

        return mat.iloc[keep_mask].copy(), dropped_names

    def _delete_constant(self, mat: pd.DataFrame) -> tuple[pd.DataFrame, list[str]]:
        """
        Elimina genes constantes, es decir, con un unico valor distinto entre
        todas las muestras tras la imputacion del paso 4.
        """
        constant_mask = mat.nunique(axis=1) <= 1
        if not constant_mask.any():
            self._log("  Paso 4.5 — no se encontraron genes constantes")
            return mat, []

        dropped = mat.index[constant_mask].tolist()
        self._log(f"  Paso 4.5 — eliminados {len(dropped)} gen(es) constantes")

        return mat.loc[~constant_mask].copy(), dropped

    # -----------------------------------------------------------------------
    # Utilidades
    # -----------------------------------------------------------------------

    def _log(self, message: str) -> None:
        """Imprime solo cuando verbose=True."""
        if self.verbose:
            print(message)
