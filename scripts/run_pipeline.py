#!/usr/bin/env python3
"""
run_pipeline.py
===============
Ejecucion del algoritmo CIRCUST (Etapas 1-4). **Solo** ejecuta el pipeline
y vuelca los resultados a CSV/TXT/JSON para que los consuman otros
programas (scripts de validacion sintetica, generadores de graficos en
Python o R, analisis aguas abajo, etc.).

No produce graficos: el plotting vive en scripts independientes.

Etapas ejecutadas
-----------------
    1.0  Preprocesamiento (min-max, filtrado de genes sparse)
    1.1  CPCA + deteccion de outliers
    2    Sincronizador (ancla ARNTL en pi, direccion via DBP)
    3    Seleccion de genes TOP (R2_ORI, R2_FMM, omega, cobertura)
    4    Estimacion robusta del orden (K repeticiones + agregacion)

Ejemplos
--------
    # Por defecto: data/matrixIn.parquet con genes core de Larriba
    python scripts/run_pipeline.py

    # Entrada CSV
    python scripts/run_pipeline.py --data data/BA46_glut_sample_no_minmax.csv

    # Preset alternativo de genes core
    python scripts/run_pipeline.py --core-genes zhang

    # Lista personalizada
    python scripts/run_pipeline.py --core-genes ARNTL,DBP,NR1D1,PER1,PER2,PER3

Salida
------
    <output>/
    ├── manifest.json         config + semilla (reproducibilidad)
    │
    ├── sample_order.csv      orden de muestras — UNA FILA POR MUESTRA.
    │                          Es la salida principal del TFG (mejora
    │                          en la inferencia del orden temporal de
    │                          muestras via agregacion de Barragan 2015).
    │                          Se valida contra ground truth (TOD o
    │                          tiempos reales conocidos):
    │                            sample_name     : nombre real de la muestra
    │                            circular_phase  : fase en RADIANES, [0, 2*pi) — canonico
    │                            phase_h         : fase en HORAS,   [0, 24)  — interpretable
    │
    └── circust_atlas.csv     [M^TOP] — atlas circadiano del tejido,
                              UNA FILA POR GEN TOP. Salida cientifica
                              canonica de CIRCUST (Larriba et al. 2023):
                                gene_name      : simbolo del gen TOP
                                r2_fmm         : R^2 del ajuste FMM
                                amplitude      : amplitud A
                                peak_phase_rad : pico t_U en [0, 2*pi)
                                peak_phase_h   : pico en [0, 24) horas
                                omega          : asimetria / agudeza
                                classification : 'day' (t_U in [0, pi))
                                                 o 'night' (t_U in [pi, 2pi))

Solo se generan estos tres archivos. Las metricas y resumenes intermedios
se imprimen por stdout durante la ejecucion pero no se guardan en disco.

Notas sobre la agregacion
-------------------------
La **agregacion opera a nivel de MUESTRAS** (propuesta del TFG, seccion 2:
aplicar la tecnica de Barragan 2015 al orden de muestras de CIRCUST):

- ``--robust-method best_k``: se selecciona la repeticion con mayor
  mediana R^2 y se usan su orden de muestras + sus parametros FMM por gen.
- ``--robust-method aggregate``: TSP3/HODs se aplica a la matriz
  ``(K reps x n_muestras)`` de fases asignadas por muestra para producir
  UN orden de muestras consenso. Despues se re-ajusta FMM por gen sobre
  ese orden consenso para obtener los parametros del atlas.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

# ── sys.path para importar el paquete circust ──────────────────────────────
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from circust.preprocessing import load_expression_matrix, Preprocessor
from circust.cpca          import CPCA
from circust.synchronizer  import CircularSynchronizer
from circust.top_genes     import TopGeneSelector
from circust.robust_order  import RobustOrderEstimator
from circust.core_genes    import CoreGeneSelector


# ═══════════════════════════════════════════════════════════════════════════
# Configuracion por defecto — editar aqui o sobreescribir via CLI
# ═══════════════════════════════════════════════════════════════════════════

DEFAULT_CONFIG = {
    # Preprocesamiento
    "sparse_threshold":       0.3,      # fraccion max de ceros/NaN por gen
    # CPCA
    "tight_radius":           0.10,     # umbral radial primario
    "loose_radius":           0.15,     # umbral radial alternativo
    # Deteccion de outliers
    "spline_gene_threshold":   3.5,     # umbral residuo estandarizado del spline, por gen
    "spline_outlier_fraction": 0.25,    # fraccion de genes core que deben coincidir
    "max_outlier_fraction":    0.05,    # limite max de muestras eliminables
    # Ajuste FMM
    "fmm_alpha_grid":         48,       # resolucion de la rejilla alpha
    "fmm_omega_grid":         24,       # resolucion de la rejilla omega
    "fmm_reps":               3,        # iteraciones de refinamiento
    # Anclas biologicas
    "anchor_gene":            "ARNTL",  # gen situado en pi (amanecer)
    "direction_gene":         "DBP",    # gen para verificacion de direccion
}

# Nombre de la columna de genes por extension. Parquet suele tener una
# columna explicita; CSV/TSV usan la primera columna como indice.
_GENE_COL_DEFAULTS = {
    ".parquet": "gene_id",
    ".csv":     None,
    ".tsv":     None,
    ".txt":     None,
    ".xlsx":    None,
    ".xls":     None,
}


def _resolve_gene_column(data_path: Path, override: str | None) -> str | None:
    if override is not None:
        return override
    return _GENE_COL_DEFAULTS.get(data_path.suffix.lower(), None)


# ═══════════════════════════════════════════════════════════════════════════
# CLI
# ═══════════════════════════════════════════════════════════════════════════

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Ejecutar el algoritmo CIRCUST (Etapas 1-4) y volcar CSVs.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=f"""
ejemplos:
  python scripts/run_pipeline.py
  python scripts/run_pipeline.py --data data/BA46_glut_sample_no_minmax.csv
  python scripts/run_pipeline.py --core-genes zhang
  python scripts/run_pipeline.py --core-genes ARNTL,DBP,NR1D1,PER1,PER2,PER3

conjuntos de genes disponibles: {list(CoreGeneSelector.PRESETS.keys())}
""",
    )
    parser.add_argument(
        "--data", type=str, default=None,
        help="Ruta a la matriz de expresion (CSV, TSV, Parquet, Excel). "
             "Por defecto: data/matrixIn.parquet",
    )
    parser.add_argument(
        "--gene-column", type=str, default=None,
        help="Nombre de la columna identificadora de genes. Se detecta "
             "automaticamente por extension si se omite (Parquet -> 'gene_id', "
             "CSV -> primera col).",
    )
    parser.add_argument(
        "--core-genes", type=str, default=None,
        help=f"Genes core: preset {list(CoreGeneSelector.PRESETS.keys())} o "
             "lista separada por comas. Por defecto: circust.",
    )
    parser.add_argument(
        "--label", type=str, default=None,
        help="Etiqueta del dataset (queda en manifest.json). Por defecto: "
             "nombre del archivo.",
    )
    parser.add_argument(
        "-o", "--output", type=str, default="output",
        help="Directorio de salida. Por defecto: output/",
    )

    # ── Parametros del algoritmo (5 esenciales) ──────────────────────────
    parser.add_argument(
        "--fmm-method", choices=["base", "stabilized"], default="base",
        help="Modelo FMM: 'base' (Larriba/Rueda 2019, MLE puro) o "
             "'stabilized' (power-roughness + KDE von Mises). "
             "Cuando se usa 'stabilized', CPCA calibra lambda y se hereda "
             "a TopGenes y RobustOrder. Por defecto: base.",
    )
    parser.add_argument(
        "--robust-method", choices=["median", "best_k", "aggregate"],
        default="best_k",
        help="Como se resumen las K repeticiones: 'median' (CIRCUST original, "
             "mediana de los parametros FMM por gen; NO produce orden de "
             "muestras), 'best_k' (selecciona una repeticion) o 'aggregate' "
             "(consenso TSP3/HODs de los K ordenes). Por defecto: best_k.",
    )
    parser.add_argument(
        "--selection-criterion",
        choices=["r2_median", "n_well_fitted", "medoid"], default="r2_median",
        help="Criterio con el que se elige la repeticion. Solo aplica si "
             "--robust-method=best_k. 'r2_median' (mediana de R2 maxima), "
             "'n_well_fitted' (mas genes con R2 >= umbral) o 'medoid' (orden "
             "mas representativo). Por defecto: r2_median.",
    )
    parser.add_argument(
        "--selection-r2-threshold", type=float, default=0.65,
        help="Umbral tau del criterio 'n_well_fitted'. Valor provisional 0.65, "
             "pendiente de calibrar en el estudio de simulacion.",
    )
    parser.add_argument(
        "--aggregation", choices=["auto", "tsp3", "hods"], default="auto",
        help="Metodo de agregacion de ordenes. Solo aplica si "
             "--robust-method=aggregate. Por defecto: auto (corre TSP3 y "
             "HODs y elige el de menor MSCE).",
    )
    parser.add_argument(
        "--consensus-scale", choices=["circmean", "cire"], default="circmean",
        help="Como se construye la escala circular del orden consenso. Solo "
             "aplica si --robust-method=aggregate. 'circmean' (media circular "
             "de las K fases, mejor estimador de fase) o 'cire' (proyeccion "
             "isotonica circular sobre el orden consenso, que garantiza "
             "phi_o1 <= ... <= phi_on a costa de precision). Por defecto: "
             "circmean.",
    )
    parser.add_argument(
        "--target-size", type=int, default=100,
        help="Numero objetivo de genes del conjunto TOP. El cupo se reparte "
             "entre los sectores y en cada uno se toman los de mayor R2_FMM, "
             "para preservar la cobertura circular. Los genes core forzados se "
             "anaden despues. Por defecto: 100.",
    )
    parser.add_argument(
        "--n-sectors", type=int, default=8,
        help="Numero de sectores circulares para la cobertura y el reparto del "
             "cupo. El paper exige cubrir los 4 cuadrantes (minimo admisible); "
             "la implementacion de referencia usa 8. Por defecto: 8.",
    )
    parser.add_argument(
        "--n-reps", type=int, default=5,
        help="Numero K de repeticiones aleatorias en RobustOrder. "
             "Por defecto: 5.",
    )
    parser.add_argument(
        "--seed", type=int, default=42,
        help="Semilla del RNG (para las selecciones aleatorias de TOP "
             "en RobustOrder). Por defecto: 42.",
    )
    return parser.parse_args()


# ═══════════════════════════════════════════════════════════════════════════
# Helpers
# ═══════════════════════════════════════════════════════════════════════════

def _banner(msg: str) -> None:
    width = 60
    print()
    print("=" * width)
    print(f"  {msg}")
    print("=" * width)


def _build_atlas(
    robust_result,
    top_result,
    cpca_result,
    fmm_method:  str,
    cfg:         dict,
) -> pd.DataFrame:
    """
    Construye la tabla [M^TOP] del atlas circadiano del tejido
    (Larriba et al. 2023, Step 4): una fila por gen TOP con sus
    parametros FMM finales (R^2, A, t_U, omega) + clasificacion
    day/night, ordenada por R^2 descendente.

    La **agregacion opera a nivel de MUESTRAS** (propuesta del TFG,
    seccion 2: aplicar la tecnica de Barragan 2015 al orden de muestras):
    los K reps producen K orderings de muestras, y se consolida en UN
    orden consenso. Luego se ajusta FMM por gen sobre ese consenso.

    - ``method="best_k"``: [M^TOP] = parametros FMM **de la rep seleccionada**
      (``robust_result.selected_fmm``), consistentes con su orden de muestras.
      No se promedia sobre las K (eso mezclaria ordenes ya descartados por
      peores segun el R² mediano).
    - ``method="aggregate"``: el orden consenso (TSP3/HODs sobre los K
      orderings de muestras) no coincide con ninguna rep individual, asi
      que **se re-fittea FMM por gen** sobre el orden consenso para
      obtener los parametros canonicos.
    """
    top_names = list(robust_result.top_gene_names)

    if robust_result.method_used in ("best_k", "median"):
        # best_k : parametros de la rep elegida.
        # median : mediana por gen sobre las K (CIRCUST original). En los dos
        #          casos el atlas ya viene calculado en selected_fmm y no hay
        #          nada que re-ajustar.
        raw = robust_result.selected_fmm[
            ["gene_name", "r2_fmm", "amplitude", "peak_phase_rad", "omega"]
        ].copy()
    else:
        # Re-fittear FMM sobre el orden consenso de MUESTRAS
        top_matrix_final = top_result.candidate_matrix.iloc[
            :, robust_result.sample_order
        ]
        final_scale = np.asarray(robust_result.circular_scale, dtype=np.float64)

        # El ajuste de cada gen es independiente, asi que se reparte entre
        # procesos igual que en las etapas 3 y 4. Antes era un bucle secuencial
        # sobre los miles de genes TOP, y con la matriz completa suponia un
        # tercio del tiempo total del pipeline corriendo en un unico nucleo.
        from joblib import Parallel, delayed
        from circust.top_genes import _fit_fmm_one, _fit_fmm_one_stab

        rows_mat = top_matrix_final.loc[top_names].values.astype(np.float64)
        n_jobs   = os.cpu_count() or 1

        if fmm_method == "base":
            ag, og, nr = (cfg["fmm_alpha_grid"], cfg["fmm_omega_grid"], cfg["fmm_reps"])
            # loky (multiproceso): los argumentos por gen son pequenos y el
            # ajuste es intensivo en CPU sin liberar el GIL.
            fits = Parallel(n_jobs=n_jobs)(
                delayed(_fit_fmm_one)(rows_mat[i], final_scale, ag, og, nr)
                for i in range(len(top_names))
            )
        else:
            from circust.fitting.fmm_stabilized import (
                StabilizedFMMModel, prepare_dataset,
            )
            stab_kwargs: dict = {}
            if cpca_result.lambda_star is not None:
                stab_kwargs["lambda_star"] = float(cpca_result.lambda_star)
            prepared   = prepare_dataset(final_scale, **stab_kwargs)
            stab_model = StabilizedFMMModel(prepared=prepared)
            # threading: los workers comparten ``prepared`` en memoria en vez
            # de serializar sus matrices precomputadas hacia cada proceso.
            fits = Parallel(n_jobs=n_jobs, backend="threading")(
                delayed(_fit_fmm_one_stab)(rows_mat[i], final_scale, stab_model)
                for i in range(len(top_names))
            )

        raw = pd.DataFrame({
            "gene_name":      list(top_names),
            "r2_fmm":         [float(f[2]) for f in fits],
            "amplitude":      [float(f[3]) for f in fits],
            "peak_phase_rad": [float(f[0]) for f in fits],
            "omega":          [float(f[1]) for f in fits],
        })

    # Postprocesar: pico en horas + clasificacion day/night
    raw["peak_phase_h"]   = raw["peak_phase_rad"] * 24.0 / (2.0 * np.pi)
    raw["classification"] = np.where(
        raw["peak_phase_rad"] < np.pi, "day", "night",
    )
    atlas = raw[[
        "gene_name", "r2_fmm", "amplitude",
        "peak_phase_rad", "peak_phase_h",
        "omega", "classification",
    ]].sort_values("r2_fmm", ascending=False).reset_index(drop=True)
    return atlas


# ═══════════════════════════════════════════════════════════════════════════
# Pipeline principal
# ═══════════════════════════════════════════════════════════════════════════

def main() -> None:
    args = parse_args()
    cfg = DEFAULT_CONFIG.copy()
    t_start = time.time()

    # Cronometro por etapas (perf_counter): se vuelca a manifest.json para
    # los benchmarks de tiempo (Exp 1 R vs Python y Exp 2/3 del grid).
    timing: dict[str, float] = {}

    # ── Resolver entradas ────────────────────────────────────────────────
    data_path = Path(args.data) if args.data else (PROJECT_ROOT / "data" / "matrixIn.parquet")
    gene_column   = _resolve_gene_column(data_path, args.gene_column)
    core_selector = CoreGeneSelector.from_string(args.core_genes)
    label         = args.label if args.label else data_path.stem

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    # ── Mostrar configuracion ────────────────────────────────────────────
    _banner("Configuracion")
    print(f"  Archivo de datos : {data_path}")
    print(f"  Columna genes    : {gene_column or '(primera columna)'}")
    print(f"  Preset genes core: {core_selector._preset_name}")
    print(f"  Gen ancla        : {cfg['anchor_gene']}")
    print(f"  Gen direccion    : {cfg['direction_gene']}")
    print(f"  FMM method       : {args.fmm_method}")
    print(f"  Robust method    : {args.robust_method}"
          + (f"  (agregacion: {args.aggregation})"
             if args.robust_method == "aggregate" else ""))
    print(f"  K repeticiones   : {args.n_reps}")
    print(f"  Semilla          : {args.seed}")
    print(f"  Etiqueta         : {label}")
    print(f"  Salida           : {output_dir}")

    # ─────────────────────────────────────────────────────────────────────
    # Etapa 0 — Cargar datos
    # ─────────────────────────────────────────────────────────────────────
    _banner("Etapa 0: Cargando matriz de expresion")
    _t = time.perf_counter()
    raw_matrix = load_expression_matrix(str(data_path), gene_column=gene_column)
    timing["load"] = time.perf_counter() - _t
    print(f"  Cargado: {raw_matrix.shape[0]} genes x {raw_matrix.shape[1]} muestras")

    # ─────────────────────────────────────────────────────────────────────
    # Etapa 1.0 — Preprocesamiento
    # ─────────────────────────────────────────────────────────────────────
    _banner("Etapa 1.0: Preprocesamiento")
    _t = time.perf_counter()
    prep = Preprocessor(
        sparse_threshold=cfg["sparse_threshold"],
        verbose=True,
    ).run(raw_matrix)
    timing["preprocess"] = time.perf_counter() - _t

    # ─────────────────────────────────────────────────────────────────────
    # Etapa 1.1a — Seleccion de genes core
    # ─────────────────────────────────────────────────────────────────────
    _banner("Etapa 1.1a: Seleccion de Genes Core")
    _t = time.perf_counter()
    core_result = core_selector.select(prep.expr_norm)
    core_genes  = core_result.genes
    timing["core_genes"] = time.perf_counter() - _t

    # ─────────────────────────────────────────────────────────────────────
    # Etapa 1.1 — CPCA + deteccion de outliers
    # ─────────────────────────────────────────────────────────────────────
    _banner("Etapa 1.1: PCA Circular")
    _t = time.perf_counter()
    cpca_result = CPCA(
        core_genes            = core_genes,
        tight_radius          = cfg["tight_radius"],
        loose_radius          = cfg["loose_radius"],
        spline_gene_threshold   = cfg["spline_gene_threshold"],
        spline_outlier_fraction = cfg["spline_outlier_fraction"],
        max_outlier_fraction  = cfg["max_outlier_fraction"],
        fmm_length_alpha_grid = cfg["fmm_alpha_grid"],
        fmm_length_omega_grid = cfg["fmm_omega_grid"],
        fmm_num_reps          = cfg["fmm_reps"],
        fmm_method            = args.fmm_method,
        verbose               = True,
    ).run(prep.expr_norm)
    timing["cpca"] = time.perf_counter() - _t

    # ─────────────────────────────────────────────────────────────────────
    # Etapa 2 — Sincronizador (orden preliminar)
    # ─────────────────────────────────────────────────────────────────────
    _banner("Etapa 2: Ordenamiento Circular Preliminar")
    _t = time.perf_counter()
    core_genes_found = list(cpca_result.core_genes_found)
    order_result = CircularSynchronizer(
        anchor_gene    = cfg["anchor_gene"],
        direction_gene = cfg["direction_gene"],
        verbose        = True,
    ).run(cpca_result)
    timing["preliminary_order"] = time.perf_counter() - _t

    # ─────────────────────────────────────────────────────────────────────
    # Etapa 3 — Seleccion de genes TOP
    # ─────────────────────────────────────────────────────────────────────
    _banner("Etapa 3: Seleccion de Genes TOP")
    _t = time.perf_counter()
    top_result = TopGeneSelector(
        r2_ori_threshold      = 0.5,
        r2_fmm_threshold      = 0.5,
        omega_min             = 0.1,
        n_sectors             = args.n_sectors,
        target_size           = args.target_size,
        fmm_length_alpha_grid = cfg["fmm_alpha_grid"],
        fmm_length_omega_grid = cfg["fmm_omega_grid"],
        fmm_num_reps          = cfg["fmm_reps"],
        fmm_method            = args.fmm_method,
        lambda_star           = cpca_result.lambda_star,   # heredado si stabilized
        n_jobs                = -1,
        verbose               = True,
    ).run(
        expr_norm      = order_result.expr_ordered,
        circular_scale = order_result.circular_scale,
        seed_genes     = core_genes_found,
        sample_order   = order_result.sample_order,
    )
    timing["top_genes"] = time.perf_counter() - _t
    print(f"  Matriz TOP: {top_result.candidate_matrix.shape[0]} genes x "
          f"{top_result.candidate_matrix.shape[1]} muestras")

    # ─────────────────────────────────────────────────────────────────────
    # Etapa 4 — Estimacion robusta del orden
    # ─────────────────────────────────────────────────────────────────────
    _banner("Etapa 4: Estimacion Robusta del Orden")
    _t = time.perf_counter()
    robust_result = RobustOrderEstimator(
        n_reps               = args.n_reps,
        sample_size_fraction = 2.0 / 3.0,
        method               = args.robust_method,
        selection_criterion  = args.selection_criterion,
        selection_r2_threshold = args.selection_r2_threshold,
        aggregation_method   = args.aggregation,
        consensus_scale_mode = args.consensus_scale,
        r2_min               = 0.5,
        anchor_gene          = cfg["anchor_gene"],
        direction_gene       = cfg["direction_gene"],
        fmm_length_alpha_grid= cfg["fmm_alpha_grid"],
        fmm_length_omega_grid= cfg["fmm_omega_grid"],
        fmm_num_reps         = cfg["fmm_reps"],
        fmm_method           = args.fmm_method,
        lambda_star          = cpca_result.lambda_star,   # heredado si stabilized
        n_jobs               = -1,   # paralelo: el ajuste FMM por gen es el cuello de botella
        seed                 = args.seed,
        verbose              = True,
    ).run(
        top_result     = top_result,
        expr_full_norm = order_result.expr_ordered,
        core_genes     = core_genes_found,
    )
    timing["robust_order"] = time.perf_counter() - _t

    # ─────────────────────────────────────────────────────────────────────
    # Salida 1: orden de MUESTRAS — UNA FILA POR MUESTRA
    # ─────────────────────────────────────────────────────────────────────
    # Es la salida principal del TFG (mejora del orden de muestras via
    # aggregation de Barragan 2015 aplicada a las K iteraciones de CIRCUST,
    # seccion 2 de la propuesta). Es lo que se valida contra ground truth
    # (tiempos reales conocidos o TOD), como hace el paper original en
    # Fig 5 (skeletal muscle) y Fig 7A (baboons).
    if robust_result.sample_order.size == 0:
        # method="median" (CIRCUST original) no produce orden de muestras: las
        # K ordenaciones se usan para ajustar y se descartan. No se escribe un
        # fichero vacio que luego parezca un fallo.
        print("    sample_order.csv NO se escribe: el metodo 'median' no "
              "define un orden de muestras consenso.")
    else:
        robust_names = [order_result.expr_ordered.columns[i]
                        for i in robust_result.sample_order]
        phases_rad   = np.asarray(robust_result.circular_scale, dtype=np.float64)
        phases_h     = phases_rad * 24.0 / (2.0 * np.pi)
        pd.DataFrame({
            "sample_name":    robust_names,
            "circular_phase": phases_rad,   # radianes en [0, 2*pi) — canonico
            "phase_h":        phases_h,     # horas    en [0, 24)  — interpretable
        }).to_csv(output_dir / "sample_order.csv", index=False)
        print(f"    Saved: sample_order.csv  ({len(phases_rad)} muestras)")

    # ─────────────────────────────────────────────────────────────────────
    # Salida 2: [M^TOP] — atlas circadiano del tejido (UNA FILA POR GEN)
    # ─────────────────────────────────────────────────────────────────────
    # Salida cientifica canonica de CIRCUST (Larriba et al. 2023, Step 4):
    # para cada gen TOP, sus parametros FMM finales (R^2, amplitud, peak
    # phase, omega) y clasificacion day/night. Construido re-ajustando
    # FMM sobre el orden de muestras consenso (en aggregate) o tomando
    # los parametros de la rep seleccionada (en best_k).
    _t = time.perf_counter()
    atlas_df = _build_atlas(
        robust_result = robust_result,
        top_result    = top_result,
        cpca_result   = cpca_result,
        fmm_method    = args.fmm_method,
        cfg           = cfg,
    )
    timing["build_atlas"] = time.perf_counter() - _t
    atlas_df.to_csv(output_dir / "circust_atlas.csv", index=False)
    print(f"    Saved: circust_atlas.csv  ({len(atlas_df)} genes TOP)")

    # ─────────────────────────────────────────────────────────────────────
    # Manifest (config + semilla) para reproducibilidad
    # ─────────────────────────────────────────────────────────────────────
    manifest = {
        "input": {
            "data_path":   str(data_path),
            "label":       label,
            "gene_column": gene_column,
            "core_preset": core_selector._preset_name,
            "core_genes":  list(core_genes),
        },
        "config":    cfg,
        "algorithm": {
            "fmm_method":    args.fmm_method,
            "robust_method": args.robust_method,
            "aggregation":   (args.aggregation
                              if args.robust_method == "aggregate" else None),
            "n_reps":        args.n_reps,
            "n_sectors":     args.n_sectors,
            "target_size":   args.target_size,
            "lambda_star":   cpca_result.lambda_star,
        },
        "seed":      args.seed,
        "timing": {
            **{k: round(v, 4) for k, v in timing.items()},
            "elapsed_total_s": round(time.time() - t_start, 4),
        },
    }
    (output_dir / "manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )
    print(f"  Saved: manifest.json")

    # ─────────────────────────────────────────────────────────────────────
    # Informe final
    # ─────────────────────────────────────────────────────────────────────
    elapsed = time.time() - t_start

    _banner("Pipeline completado")
    print(f"  Dataset            : {data_path.name}")
    print(f"  Genes (entrada)    : {prep.n_genes_in}")
    print(f"  Genes (filtrados)  : {prep.n_genes_out}")
    print(f"  Muestras (limpias) : {cpca_result.expr_norm_final.shape[1]}")
    print(f"  Outliers eliminados: {len(cpca_result.samples_dropped)}")
    print(f"  Inversion direcc.  : {order_result.direction_flipped}")
    print(f"  Genes tras ORI     : {top_result.n_after_ori}")
    print(f"  Genes tras FMM     : {top_result.n_after_fmm}")
    print(f"  Genes TOP finales  : {len(top_result.gene_names)}")
    print(f"  Metodo orden rob.  : {robust_result.method_used}")
    print(f"  Repeticiones K     : {robust_result.sample_orders.shape[0]}")
    print(f"  Flips de orient.   : {int(robust_result.direction_flipped.sum())}")
    print(f"  Salida             : {output_dir.resolve()}")
    print(f"  Tiempo transcurrido: {elapsed:.1f}s")


if __name__ == "__main__":
    main()
