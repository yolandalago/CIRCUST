#!/usr/bin/env python3
"""Compara las aportaciones del TFG sobre datos con tiempos conocidos.

CIRCUST resuelve la incertidumbre del orden calculando K ordenaciones y
quedandose con la **mediana de los parametros FMM por gen** (Larriba et al.
2023, paso 4). Los K ordenes se usan para ajustar y se descartan: el metodo no
define ningun orden de muestras consenso.

Este trabajo propone agregar los K ordenes entre si (Barragan et al. 2015) para
obtener un unico orden consenso y calcular los parametros una sola vez sobre
el. Este script mide si eso mejora, en que, y a que coste.

Ejes comparados
---------------
1. Como se resumen las K repeticiones:

       median                CIRCUST original. NO produce orden de muestras.
       best_k:r2_median      la repeticion de mayor R2 mediano
       best_k:n_well_fitted  la de mas genes bien ajustados
       best_k:medoid         la mas central entre las K
       aggregate:tsp3        consenso de los K ordenes por TSP3
       aggregate:hods        consenso de los K ordenes por HODs

2. Modelo FMM: ``base`` (Rueda et al. 2019) o ``stabilized`` (CircaFMM).

3. Tamano del conjunto TOP. No se barre para elegirlo (el barrido de
   ``sweep_top.py`` mostro que la superficie es plana), sino para comprobar que
   el orden entre metodos no depende de esa eleccion.

Metricas
--------
**Atlas** — ``picos_medae_h``: error mediano de las fases de pico ``t_U``, en
horas. La referencia se obtiene ajustando FMM a cada gen usando los TIEMPOS
REALES como escala circular, que es lo que hace el articulo al validar sobre
GSE71620 ("the difference between TOD and CIRCUST for the seed genes is, on
median, lower than 1.671 hours"). Como las fases solo estan definidas salvo
rotacion y reflexion, se alinean antes de comparar. Es la unica metrica que
admite los seis metodos, y mide la salida canonica del algoritmo.

La referencia se ajusta **con el mismo modelo FMM** que la configuracion que se
esta evaluando. Asi el sesgo del estimador se cancela entre las dos partes y lo
que queda medido es el error del ORDEN, que es lo que se compara. El precio es
que ``picos_medae_h`` no sirve para decidir que modelo FMM es mejor en terminos
absolutos, solo cual reconstruye mejor el orden.

**Orden de muestras** — ``rho`` (regresion circular-circular) y ``mad_h``
(MAD frente a los tiempos reales). Son las metricas de las validaciones del
articulo sobre babuino y musculo esqueletico. ``median`` queda a NaN porque no
produce orden: no es un fallo, es la carencia que motiva el trabajo.

**Coste** — ``segundos`` de la etapa robusta. Mide el sobrecoste marginal de
agregar frente a resumir por medianas.

Todos los metodos de una misma unidad comparten preprocesado, CPCA, orden
preliminar, conjunto TOP, semilla y por tanto las MISMAS K repeticiones. Lo
unico que cambia es como se resumen, asi que cualquier diferencia es
atribuible a eso.

Uso
---
    python validation/compare_methods.py --datos data/babuinos
    python validation/compare_methods.py --datos data/GSE71620
    python validation/compare_methods.py --datos data/babuinos --analizar
    python validation/compare_methods.py --datos data/babuinos \
        --target-sizes 100 --fmm-methods base        # version corta
"""
from __future__ import annotations

import argparse
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))

from circust.preprocessing import load_expression_matrix, Preprocessor
from circust.cpca import CPCA
from circust.synchronizer import CircularSynchronizer
from circust.top_genes import TopGeneSelector, _fit_fmm_one, _fit_fmm_one_stab
from circust.robust_order import RobustOrderEstimator
from circust.core_genes import CoreGeneSelector
from validation.metrics import (
    align_phases_by_mad, circular_abs_diff, circular_regression_rho,
)
from validation.sweep_top import CFG, R2_ORI, R2_FMM, OMEGA_MIN, descubrir

warnings.filterwarnings("ignore")

K, SEED, N_SECTORS = 5, 42, 8

METODOS = [
    ("median",               dict(method="median")),
    ("best_k:r2_median",     dict(method="best_k", selection_criterion="r2_median")),
    ("best_k:n_well_fitted", dict(method="best_k", selection_criterion="n_well_fitted")),
    ("best_k:medoid",        dict(method="best_k", selection_criterion="medoid")),
    ("aggregate:tsp3",       dict(method="aggregate", aggregation_method="tsp3")),
    ("aggregate:hods",       dict(method="aggregate", aggregation_method="hods")),
]


def ruta_salida(datos: Path) -> Path:
    return ROOT / "validation" / "resultados" / f"compare_methods_{datos.name}.csv"


def ajustar_picos(filas: np.ndarray, escala: np.ndarray, fmm_method: str,
                  lambda_star: float | None) -> list[float]:
    """Fase de pico t_U de cada fila de expresion sobre la escala circular dada."""
    if fmm_method == "base":
        from joblib import Parallel, delayed
        ajustes = Parallel(n_jobs=-1)(
            delayed(_fit_fmm_one)(filas[i], escala, CFG["fmm_alpha_grid"],
                                  CFG["fmm_omega_grid"], CFG["fmm_reps"])
            for i in range(len(filas))
        )
    else:
        from circust.fitting.fmm_stabilized import StabilizedFMMModel, prepare_dataset
        kw = {} if lambda_star is None else {"lambda_star": float(lambda_star)}
        modelo = StabilizedFMMModel(prepared=prepare_dataset(escala, **kw))
        ajustes = [_fit_fmm_one_stab(filas[i], escala, modelo)
                   for i in range(len(filas))]
    return [float(a[0]) for a in ajustes]


def picos_de_referencia(expr: pd.DataFrame, genes: list[str], reales_h: np.ndarray,
                        fmm_method: str, lambda_star) -> dict[str, float]:
    """Fases de pico ajustando FMM sobre los TIEMPOS REALES.

    Es el "orden por TOD" contra el que el articulo contrasta el orden
    estimado. Las muestras se reordenan por tiempo real creciente, que es lo
    que espera el ajuste.
    """
    presentes = [g for g in genes if g in expr.index]
    if not presentes:
        return {}
    escala = np.asarray(reales_h, dtype=np.float64) * 2 * np.pi / 24.0
    orden  = np.argsort(escala)
    filas  = expr.loc[presentes].to_numpy(dtype=np.float64)[:, orden]
    picos  = ajustar_picos(filas, escala[orden], fmm_method, lambda_star)
    return dict(zip(presentes, picos))


def error_de_picos(atlas: pd.DataFrame, referencia: dict[str, float]) -> dict:
    """Error mediano de las fases de pico, en horas, tras alinear.

    Las fases solo estan definidas salvo rotacion y reflexion (el CPCA fija el
    origen por convencion, anclando ARNTL en pi), asi que compararlas en bruto
    daria un error arbitrario. Se alinea minimizando la mediana del error
    absoluto circular, probando ambos sentidos de giro.
    """
    comunes = [g for g in atlas["gene_name"] if g in referencia]
    if len(comunes) < 3:
        return {"picos_medae_h": np.nan, "picos_p75_h": np.nan,
                "picos_n_genes": len(comunes)}

    idx = atlas.set_index("gene_name")
    est = np.asarray([idx.loc[g, "peak_phase_rad"] for g in comunes]) * 24 / (2 * np.pi)
    ref = np.asarray([referencia[g] for g in comunes]) * 24 / (2 * np.pi)

    err = circular_abs_diff(align_phases_by_mad(est, ref).aligned_pred, ref)
    return {"picos_medae_h": float(np.median(err)),
            "picos_p75_h":   float(np.percentile(err, 75)),
            "picos_n_genes": len(comunes)}


def atlas_de(rob, top, fmm_method: str, lambda_star) -> pd.DataFrame:
    """Atlas del resultado robusto, reajustando solo si el metodo lo exige.

    ``median`` y ``best_k`` lo traen calculado en ``selected_fmm``.
    ``aggregate`` no, porque el orden consenso no coincide con ninguna
    repeticion y hay que reajustar FMM sobre el.
    """
    if rob.method_used in ("median", "best_k"):
        return rob.selected_fmm

    nombres = list(rob.top_gene_names)
    filas = (top.candidate_matrix.iloc[:, rob.sample_order]
             .loc[nombres].to_numpy(dtype=np.float64))
    escala = np.asarray(rob.circular_scale, dtype=np.float64)
    return pd.DataFrame({
        "gene_name":      nombres,
        "peak_phase_rad": ajustar_picos(filas, escala, fmm_method, lambda_star),
    })


def comparar_unidad(nombre: str, matriz: Path, tiempos_h: dict[str, float],
                    target_sizes: list[int], fmm_methods: list[str]) -> list[dict]:
    mat  = load_expression_matrix(str(matriz), gene_column=None)
    prep = Preprocessor(sparse_threshold=0.3, verbose=False).run(mat)
    core = CoreGeneSelector.from_string(None).select(prep.expr_norm).genes

    filas: list[dict] = []
    for fmm_method in fmm_methods:
        cpca = CPCA(core_genes=core, n_outlier_candidates=8, tight_radius=0.10,
                    loose_radius=0.15, fmm_threshold=3.0, max_outlier_fraction=0.05,
                    fmm_length_alpha_grid=CFG["fmm_alpha_grid"],
                    fmm_length_omega_grid=CFG["fmm_omega_grid"],
                    fmm_num_reps=CFG["fmm_reps"], fmm_method=fmm_method,
                    verbose=False).run(prep.expr_norm)
        encontrados = list(cpca.core_genes_found)
        orden = CircularSynchronizer(anchor_gene=CFG["anchor_gene"],
                                     direction_gene=CFG["direction_gene"],
                                     verbose=False).run(cpca, encontrados)
        columnas = list(orden.expr_ordered.columns)
        reales   = np.array([tiempos_h[c] for c in columnas], dtype=np.float64)

        # TODO(calibracion): fit_cache es andamiaje temporal de TopGeneSelector.
        cache: dict = {}
        for objetivo in target_sizes:
            top = TopGeneSelector(
                r2_ori_threshold=R2_ORI, r2_fmm_threshold=R2_FMM,
                omega_min=OMEGA_MIN, n_sectors=N_SECTORS, target_size=objetivo,
                fmm_length_alpha_grid=CFG["fmm_alpha_grid"],
                fmm_length_omega_grid=CFG["fmm_omega_grid"],
                fmm_num_reps=CFG["fmm_reps"], fmm_method=fmm_method,
                lambda_star=cpca.lambda_star, n_jobs=-1, verbose=False,
            ).run(expr_norm=orden.expr_ordered,
                  circular_scale=orden.circular_scale,
                  seed_genes=encontrados, sample_order=orden.sample_order,
                  fit_cache=cache)

            referencia = picos_de_referencia(orden.expr_ordered,
                                             list(top.gene_names), reales,
                                             fmm_method, cpca.lambda_star)

            for etiqueta, kwargs in METODOS:
                t0 = time.perf_counter()
                rob = RobustOrderEstimator(
                    n_reps=K, sample_size_fraction=2.0 / 3.0,
                    selection_r2_threshold=0.65, r2_min=0.5,
                    anchor_gene=CFG["anchor_gene"],
                    direction_gene=CFG["direction_gene"],
                    fmm_length_alpha_grid=CFG["fmm_alpha_grid"],
                    fmm_length_omega_grid=CFG["fmm_omega_grid"],
                    fmm_num_reps=CFG["fmm_reps"], fmm_method=fmm_method,
                    lambda_star=cpca.lambda_star, n_jobs=-1, seed=SEED,
                    verbose=False, **kwargs,
                ).run(top_result=top, expr_full_norm=orden.expr_ordered,
                      core_genes=encontrados)
                segundos = time.perf_counter() - t0

                fila = {"unidad": nombre, "n_muestras": len(columnas),
                        "fmm_method": fmm_method, "target_size": objetivo,
                        "n_top": len(top.gene_names), "metodo": etiqueta,
                        "segundos": segundos}
                fila.update(error_de_picos(
                    atlas_de(rob, top, fmm_method, cpca.lambda_star), referencia))

                if rob.sample_order.size == 0:
                    # 'median' no define orden: NaN a proposito, no por fallo.
                    fila.update({"rho": np.nan, "mad_h": np.nan})
                else:
                    est_h = np.asarray(rob.circular_scale) * 24 / (2 * np.pi)
                    ver_h = np.array([tiempos_h[columnas[i]]
                                      for i in rob.sample_order])
                    fila["rho"]   = float(circular_regression_rho(est_h, ver_h))
                    fila["mad_h"] = float(align_phases_by_mad(est_h, ver_h).mad)
                filas.append(fila)
    return filas


def analizar(datos: Path) -> None:
    salida = ruta_salida(datos)
    if not salida.exists():
        raise SystemExit(f"No existe {salida}. Lanza antes la comparacion.")
    d = pd.read_csv(salida)
    pd.set_option("display.width", 240)
    orden_m = [m for m, _ in METODOS]

    print(f"=== {datos.name}: {d['unidad'].nunique()} unidades, "
          f"{int(d['n_muestras'].median())} muestras (mediana) ===")

    for col, titulo, sentido in [
        ("picos_medae_h", "ATLAS - error mediano de las fases de pico (h)", "menor mejor"),
        ("rho",           "ORDEN - rho de regresion circular-circular",     "mayor mejor"),
        ("mad_h",         "ORDEN - MAD frente a los tiempos reales (h)",    "menor mejor"),
        ("segundos",      "COSTE - segundos de la etapa robusta",           "menor mejor"),
    ]:
        print(f"\n{titulo}   ({sentido})")
        t = d.pivot_table(index="metodo", columns=["fmm_method", "target_size"],
                          values=col, aggfunc="median").reindex(orden_m)
        print(t.round(3).to_string())

    # ── Comparacion pareada contra el metodo original ─────────────────
    # Pareada por unidad: cada metodo se compara con 'median' sobre la MISMA
    # unidad y configuracion, que es donde comparten las K repeticiones.
    print("\n" + "=" * 72)
    print("Frente al metodo original ('median'), pareado por unidad")
    for fmm in sorted(d["fmm_method"].unique()):
        for tgt in sorted(d["target_size"].unique()):
            sub  = d[(d.fmm_method == fmm) & (d.target_size == tgt)]
            base = sub[sub.metodo == "median"].set_index("unidad")
            if base.empty:
                continue
            print(f"\n  fmm={fmm}  target_size={tgt}")
            for etiqueta in orden_m[1:]:
                otro  = sub[sub.metodo == etiqueta].set_index("unidad")
                comun = base.index.intersection(otro.index)
                if not len(comun):
                    continue
                dif = (otro.loc[comun, "picos_medae_h"]
                       - base.loc[comun, "picos_medae_h"])
                print(f"    {etiqueta:22s} picos: mejora en "
                      f"{int((dif < 0).sum()):3d}/{len(comun)} unidades, "
                      f"diferencia mediana {dif.median():+.3f} h")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--datos", type=Path, required=True,
                   help="Directorio del conjunto etiquetado (ver sweep_top.py).")
    p.add_argument("--solo", type=str, default=None,
                   help="Lista separada por comas para restringir las unidades.")
    p.add_argument("--target-sizes", type=str, default="50,100,200")
    p.add_argument("--fmm-methods", type=str, default="base,stabilized")
    p.add_argument("--analizar", action="store_true")
    args = p.parse_args()

    datos = args.datos.resolve()
    if not datos.is_dir():
        raise SystemExit(f"No existe el directorio {datos}")
    if args.analizar:
        analizar(datos)
        return

    target_sizes = [int(s) for s in args.target_sizes.split(",")]
    fmm_methods  = [s.strip() for s in args.fmm_methods.split(",")]

    unidades = descubrir(datos)
    if args.solo:
        permitidas = {s.strip() for s in args.solo.split(",")}
        unidades = [u for u in unidades if u[0] in permitidas]
    if not unidades:
        raise SystemExit(f"No se ha encontrado ninguna unidad en {datos}")
    print(f"{len(unidades)} unidades x {len(fmm_methods)} modelos FMM x "
          f"{len(target_sizes)} tamanos x {len(METODOS)} metodos")

    salida = ruta_salida(datos)
    salida.parent.mkdir(parents=True, exist_ok=True)
    filas, hechos = [], set()
    if salida.exists():
        previo = pd.read_csv(salida)
        filas  = previo.to_dict("records")
        hechos = set(previo["unidad"])
        print(f"Reanudando: {len(hechos)} unidades ya hechas", flush=True)

    for nombre, matriz, ruta_tiempos in unidades:
        if nombre in hechos:
            continue
        tiempos = pd.read_csv(ruta_tiempos)
        tiempos_h = dict(zip(tiempos["sample_name"], tiempos["true_time_h"]))
        try:
            nuevas = comparar_unidad(nombre, matriz, tiempos_h,
                                     target_sizes, fmm_methods)
        except Exception as exc:                          # noqa: BLE001
            print(f"[{nombre}] ERROR: {type(exc).__name__}: {exc}", flush=True)
            continue
        filas.extend(nuevas)
        ref = [f for f in nuevas
               if f["target_size"] == target_sizes[0]
               and f["fmm_method"] == fmm_methods[0]]
        resumen = "  ".join(
            f"{f['metodo'].split(':')[-1][:7]}={f['picos_medae_h']:.2f}" for f in ref)
        print(f"[{nombre}] picos MedAE (h) | {resumen}", flush=True)
        pd.DataFrame(filas).to_csv(salida, index=False)

    print(f"\nGuardado: {salida}")
    print(f"Resumen con: python validation/compare_methods.py "
          f"--datos {args.datos} --analizar")


if __name__ == "__main__":
    main()
