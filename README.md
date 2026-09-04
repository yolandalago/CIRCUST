# CIRCUST - Python Implementation

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![DOI](https://img.shields.io/badge/DOI-10.1371%2Fjournal.pcbi.1011510-blue)](https://doi.org/10.1371/journal.pcbi.1011510)
[![Docs](https://img.shields.io/badge/docs-GitHub%20Pages-blue?logo=github)](https://samuel-ortega-bustos.github.io/CIRCUST/)

**Autores de la metodología original (R):** Yolanda Larriba, Ivy C. Mason, Richa Saxena, Frank A.J.L. Scheer, Cristina Rueda  
**Implementación en Python:** Samuel Ortega Bustos (TFG supervisado por Yolanda Larriba)

Este repositorio contiene una implementación en Python de la metodología **CIRCUST** (CIRCular-robUST), originalmente desarrollada en R por Yolanda Larriba y colaboradores (Larriba et al., 2023) para la reconstrucción del orden temporal de ritmos moleculares y la construcción de atlas de expresión génica circadiana en humanos.

## 🌐 Contexto

Este proyecto es un *fork* del [repositorio original en R](https://github.com/yolandalago/CIRCUST/) y forma parte de un Trabajo de Fin de Grado (TFG) supervisado por la autora principal. El objetivo es traducir y adaptar el flujo de trabajo a Python, aprovechando las ventajas de este ecosistema para el procesamiento de datos genómicos a gran escala, manteniendo la fidelidad al método original.

## 📄 Descripción

CIRCUST aborda dos problemas clave en cronobiología molecular a partir de datos de expresión génica con tiempos de muestreo desconocidos o imprecisos (por ejemplo, muestras post-mortem):

1. **Reconstrucción del orden temporal** de las muestras utilizando Análisis de Componentes Principales Circular (CPCA).
2. **Estimación robusta de ritmicidad** mediante el modelo paramétrico FMM (Frequency Modulated Möbius), que permite capturar patrones asimétricos no sinusoidales.

La versión original en R fue validada en múltiples conjuntos de datos (epidermis humana, biopsias de músculo, cerebro post-mortem, tejidos de babuino) y aplicada a 34 tejidos del proyecto GTEx, generando el atlas de ritmos diarios más completo hasta la fecha.

Esta implementación en Python sigue los pasos descritos en el artículo original:

- Preprocesamiento (filtrado de genes con baja expresión, normalización min-max).
- Detección de outliers mediante CPCA y residuos FMM.
- Ordenación preliminar usando genes semilla (12 genes centrales del reloj circadiano).
- Selección de genes TOP (altamente rítmicos y con fases distribuidas) y remuestreo para obtener órdenes robustos.
- Ajuste del modelo FMM y cálculo de parámetros (amplitud, fase pico, nitidez).
- Generación de atlas por tejido.

## 🔁 Relación con el repositorio original

- **Repositorio original (R)**: [yolandalago/CIRCUST](https://github.com/yolandalago/CIRCUST/)
- **Artículo de referencia**: [PLOS Comput Biol 19(9): e1011510](https://doi.org/10.1371/journal.pcbi.1011510)
- **Este repositorio**: Traducción a Python, manteniendo la lógica y los pasos del método original, pero con adaptaciones propias del entorno Python (uso de `pandas`, `numpy`, `scipy`, `matplotlib`, y `scikit-learn` para PCA).

Se ha partido de los archivos de datos y resultados generados por la autora (en formato RData) para verificar la corrección de la implementación.

## 🚀 Características

- Implementación modular (cada paso del método en un módulo separado).
- Soporte para leer datos desde archivos Parquet (recomendado). La lectura directa de RData vía `pyreadr` solo se usa en desarrollo, para validar contra los resultados del repositorio en R.
- Visualización de resultados (fases pico, patrones de expresión, comparación con tiempos reales).
- Generación de tablas resumen similares a las del artículo (genes TOP por tejido con sus parámetros FMM).
- Validación reproducible mediante notebooks que replican las figuras del paper.

## ⚙️ Instalación

Se necesita Python 3.10 o superior. La fuente de verdad de las versiones es `pyproject.toml`; los ficheros `requirements*.txt` son atajos equivalentes.

**Uso normal (algoritmo + aceleración en CPU con Numba), recomendado:**

```bash
git clone https://github.com/Samuel-Ortega-Bustos/CIRCUST.git
cd CIRCUST
pip install -e ".[accel]"          # o: pip install -r requirements.txt
python scripts/check_env.py        # comprueba que no falta nada
```

**Solo el algoritmo, sin Numba:**

```bash
pip install -e .
```

Funciona igual y da los mismos resultados. La diferencia es velocidad: el ajuste FMM base tarda 2 o 3 veces más porque la rejilla y el objetivo de Nelder-Mead se evalúan en NumPy puro en lugar de compilados.

**Ruta CUDA opcional del FMM base:**

```bash
pip install -e ".[accel,cuda]"     # o: pip install -r requirements.txt -r requirements-cuda.txt
```

Requiere una GPU NVIDIA con driver compatible con CUDA 12. El toolkit llega en ruedas pip, no hace falta instalarlo en el sistema. La ruta está desactivada por defecto y se activa con `FMMModel(use_cuda=True)`. En GPU de consumo la CPU compilada sigue siendo más rápida, porque la rejilla es pequeña y la doble precisión va a 1/64 de la simple, así que solo tiene interés para verificar ese camino.

Para usar CIRCUST no hace falta nada más. La suite de tests (`pytest`) y la lectura de los RData del repositorio original (`pyreadr`) solo sirven para desarrollar y validar la traducción, y se instalan con `pip install -e ".[dev]"`.

`python scripts/check_env.py` imprime qué hay instalado, qué falta y el comando para arreglarlo, además de si Numba ve la GPU. Con `--strict` devuelve código de salida 1 si falta alguna dependencia obligatoria, útil en CI.

## 📦 Estructura del repositorio
```bash
circust/
├── README.md
├── LICENSE
├── .gitignore
├── pyproject.toml # Metadatos del paquete y dependencias (fuente de verdad)
├── requirements.txt # Atajo: algoritmo + aceleración en CPU
├── requirements-cuda.txt # Atajo: extras para la ruta CUDA opcional
├── circust/ # Paquete principal
│ ├── __init__.py
│ ├── preprocessing.py
│ ├── cpca.py
│ ├── outlier_detection.py
│ ├── order_estimation.py
│ ├── fmm.py
│ ├── top_genes.py
│ ├── utils.py
│ └── constants.py
├── tests/ # Pruebas unitarias
│ └── test_*.py
├── notebooks/ # Jupyter notebooks con ejemplos y validación
│ ├── 
├── scripts/ # Scripts para ejecutar análisis completos
│ ├── run_pipeline.py
│ └── check_env.py # Comprueba dependencias, Numba y GPU
├── rscripts/ #Scripts de codigo original de R
│ ├──
└── data/ # Carpeta para datos (ignorada por git)
```