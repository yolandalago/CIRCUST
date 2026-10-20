# CIRCUST — implementación en Python

Implementación de CIRCUST para reconstruir órdenes circadianos a partir de matrices de expresión génica. **Todos los cálculos se realizan desde los datos de entrada**: no se requieren archivos de R, checkpoints ni ajustes precalculados.

## 1. Estructura del repositorio

```text
README.md
requirements.txt
pyproject.toml
LICENSE
.gitignore
data/.gitkeep
scripts/run_circust.py
src/circust/
    preprocessing.py
    core_genes.py
    cpca.py
    synchronizer.py
    top_genes.py
    robust_order.py
    outputs.py
    fitting/
```

**Importante:** no existe una carpeta `python/`. Al descargar el ZIP de GitHub, entra en la carpeta descomprimida que contiene `pyproject.toml`, `requirements.txt`, `scripts` y `src`. Todos los comandos siguientes se ejecutan **desde esa carpeta**. No hay que copiar ni mover el contenido de ninguna subcarpeta `python/`.

## 2. Instalación en Windows (PowerShell)

Instala Python 3.10 o posterior. La versión utilizada en las pruebas de este proyecto fue Python 3.14.4 en Windows; no se han comprobado todas las versiones intermedias.

Abre PowerShell en la carpeta raíz del proyecto. Por ejemplo, para la carpeta utilizada en la prueba:

```powershell
cd "C:\Users\yolanda\Downloads\CIRCUST_Python_GitHub_v1"
```

**Sustituye esta ruta por la carpeta donde hayas descomprimido el repositorio.** Verifica que estás en la carpeta correcta:

```powershell
Test-Path .\pyproject.toml
Test-Path .\scripts\run_circust.py
```

**Ambos comandos deben devolver `True`.** Si alguno devuelve `False`, revisa la carpeta antes de continuar.

Crea el entorno virtual e instala el paquete y sus dependencias:

```powershell
py -m venv .venv
.\.venv\Scripts\python.exe -m pip install --upgrade pip
.\.venv\Scripts\python.exe -m pip install -e .
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
```

`pip install -e .` es imprescindible para que `from circust...` funcione al ejecutar `scripts/run_circust.py`. `requirements.txt` instala además Numba para acelerar los cálculos. No es necesario activar el entorno virtual.

Comprueba la instalación:

```powershell
.\.venv\Scripts\python.exe -c "import circust; print('CIRCUST importado correctamente')"
.\.venv\Scripts\python.exe .\scripts\run_circust.py --help
```

Si PowerShell indica que `py` no existe, prueba `python -m venv .venv` y sigue con los comandos que empiezan por `.\.venv\Scripts\python.exe`. Si aparece `ModuleNotFoundError: No module named 'circust'`, comprueba que estás en la raíz (donde está `pyproject.toml`) y repite `.\.venv\Scripts\python.exe -m pip install -e .`.

## 3. Preparación de datos

La matriz debe tener **genes en filas y muestras en columnas**. Los valores de expresión deben ser numéricos y la primera columna del CSV debe identificar los genes (o se indica otra con `--gene-column`). No incluyas columnas de metadatos de muestras. CIRCUST aplica su propio preprocesamiento y normalización; no normalices manualmente a [-1, 1].


El programa admite CSV, TSV, Parquet y Excel (para Excel puede ser necesario instalar `openpyxl`: `.\.venv\Scripts\python.exe -m pip install openpyxl`). Para Parquet se espera por defecto una columna `gene_id`. Si la columna de genes tiene otro nombre, añade `--gene-column NOMBRE`.

### CSV con punto y coma y coma decimal

Algunos CSV utilizan `;` como separador de campos y `,` como separador decimal. Si se leen con la configuración estándar de `pandas.read_csv`, todos los campos pueden quedar agrupados en una sola columna. **Convierte estos archivos a un CSV estándar antes de ejecutar CIRCUST.**

Si la primera columna ya se llama `gene_id`, no es necesario renombrarla. Si contiene otra variable en lugar de identificadores de genes, revisa el archivo antes de convertirlo.

### Datos RDS

El script no lee `.rds` directamente. Si el objeto de R es una matriz numérica con genes en los nombres de fila:

## 4. Ejecutar CIRCUST desde cero

Coloca tu matriz en `data/`. Por ejemplo:

```powershell
.\.venv\Scripts\python.exe .\scripts\run_circust.py --data data/matrixIn.csv --r-compatible --n-reps 5 --jobs -1 -o resultados_circust
```

## 5. Resultados

Con K=5, **una ejecución completada correctamente genera 20 archivos de resultados**:

- `sample_order_k1.csv` a `sample_order_k5.csv`: nombre, posición y fase circular de cada muestra.
- `ordered_matrix_k1.csv` a `ordered_matrix_k5.csv`: expresión observada, reordenada, de los genes preprocesados (CSV o Parquet según la opción elegida).
- `fitted_top_k1.csv` a `fitted_top_k5.csv`: expresión FMM ajustada de los genes TOP.
- `parameters_by_k.csv`: parámetros de cada gen TOP y cada repetición.
- `circust_atlas.csv`: resumen mediante medianas de las K estimaciones.
- `top_genes.csv`: genes TOP y estadísticas preliminares.
- `outliers.csv`: muestras excluidas.
- `manifest.json`: configuración, información de selección y tiempos.

Para ver el manifiesto:

```powershell
Get-Content .\resultados_circust\manifest.json
```

## 6. Tiempo de ejecución

En una ejecución previa de MatrixIn (Windows, Python 3.14.4) se procesaron desde cero 56.200 genes y 479 muestras, obteniendo 21.149 genes tras preprocesamiento, 22 muestras excluidas, 65 genes TOP y cinco órdenes. El tiempo comunicado fue 48,4 s hasta la escritura. Es una referencia de una sola ejecución, no una comparación formal de rendimiento con R.


## 7. Publicación en GitHub

Sube **directamente el contenido del ZIP a la raíz del repositorio** (`README.md`, `requirements.txt`, `pyproject.toml`, `scripts/`, `src/`, etc.). **No existe una carpeta `python/` y no hay que crearla.** No subas `.venv/`, datos de expresión, resultados, archivos temporales ni material de terceros sin autorización.
