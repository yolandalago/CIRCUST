#!/usr/bin/env python3
"""Barrido de sensibilidad de (n_sectors, target_size) sobre datos etiquetados.

El articulo CIRCUST describe la construccion del conjunto TOP (material
suplementario, Fig B, pasos 3.2-3.4: dividir el circulo en intervalos,
recorrer cada uno por R2 decreciente y quedarse con los de mayor R2_FMM del
tejido) pero no cuantifica ni cuantos intervalos ni cuantos genes. Los
umbrales si estan fijados en el articulo (R2_ORI >= 0.5, R2_FMM >= 0.5,
omega > 0.1) y aqui no se tocan.

Este script recorre la rejilla de los dos parametros que el articulo deja
abiertos y mide rho (regresion circular-circular frente a los tiempos reales)
en cada unidad de un conjunto etiquetado, para poder elegirlos con un criterio
empirico en vez de por copia de la implementacion de referencia.

Formato de entrada
------------------
``--datos`` apunta a un directorio. El script descubre como unidad cada
subdirectorio que contenga ``matriz.csv``, y toma los tiempos reales de su
propio ``tiempos.csv`` o, si no lo tiene, del ``tiempos.csv`` del directorio
padre (util cuando todas las unidades comparten calendario de muestreo). No
hay ninguna lista de conjuntos cableada: para anadir datos nuevos basta con
dejarlos en ese formato.

    <datos>/
        tiempos.csv          opcional, compartido por todas las unidades
        <unidad>/matriz.csv  genes x muestras, primera columna = gen
        <unidad>/tiempos.csv opcional, propio de la unidad

``tiempos.csv`` tiene columnas ``sample_name`` y ``true_time_h``.

La eleccion NO debe hacerse con el maximo de la mediana sobre las unidades de
un conjunto: eso es ajustar y reportar sobre el mismo dato. ``--analizar``
aplica leave-one-out sobre las unidades, y la configuracion elegida debe
validarse despues en datos que no hayan intervenido en el barrido.

Uso
---
    python validation/sweep_top.py --datos data/babuinos
    python validation/sweep_top.py --datos data/GSE71620
    python validation/sweep_top.py --datos data/babuinos --analizar
"""
from __future__ import annotations

import argparse
import sys
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
from circust.top_genes import TopGeneSelector
from circust.robust_order import RobustOrderEstimator
from circust.core_genes import CoreGeneSelector
from validation.metrics import circular_regression_rho

warnings.filterwarnings("ignore")

# Umbrales del articulo. No son parametros del barrido.
R2_ORI, R2_FMM, OMEGA_MIN = 0.5, 0.5, 0.1

# Rejilla. n_sectors >= 4 porque el articulo exige cubrir los cuatro
# cuadrantes; el techo lo impone el numero de muestras, que en babuino es 12.
N_SECTORS   = [4, 6, 8, 12]
TARGET_SIZE = [50, 100, 150, 200, 300, 400]

CFG = dict(fmm_alpha_grid=48, fmm_omega_grid=24, fmm_reps=3,
           anchor_gene="ARNTL", direction_gene="DBP")
K, SEED = 5, 42


def descubrir(datos: Path) -> list[tuple[str, Path, Path]]:
    """Unidades de un directorio de datos: (nombre, matriz.csv, tiempos.csv).

    Es unidad todo subdirectorio con ``matriz.csv``. Los tiempos salen de su
    propio ``tiempos.csv`` o, en su defecto, del compartido del directorio
    padre.
    """
    compartido = datos / "tiempos.csv"
    unidades = []
    for sub in sorted(p for p in datos.iterdir() if p.is_dir()):
        matriz = sub / "matriz.csv"
        if not matriz.exists():
            continue
        tiempos = sub / "tiempos.csv"
        if not tiempos.exists():
            tiempos = compartido
        if not tiempos.exists():
            print(f"  AVISO: {sub.name} no tiene tiempos.csv ni lo hereda; se omite.")
            continue
        unidades.append((sub.name, matriz, tiempos))
    return unidades


def ruta_salida(datos: Path) -> Path:
    return ROOT / "validation" / "resultados" / f"sweep_top_{datos.name}.csv"


def rho_de(nombres, fases_rad, tiempos_h) -> float:
    """rho entre las fases estimadas y los ZT reales de esas muestras.

    ``circular_regression_rho`` trabaja en horas (period=24), de ahi la
    conversion desde radianes.
    """
    reales = np.array([tiempos_h[n] for n in nombres], dtype=float)
    estimadas = np.asarray(fases_rad, float) * 24.0 / (2 * np.pi)
    return float(circular_regression_rho(estimadas, reales))


def etapas_previas(matriz: Path):
    """Etapas 1-2: preprocesado, CPCA y orden preliminar. Una vez por unidad."""
    mat = load_expression_matrix(str(matriz), gene_column=None)
    prep = Preprocessor(sparse_threshold=0.3, verbose=False).run(mat)
    core = CoreGeneSelector.from_string(None).select(prep.expr_norm).genes
    cpca = CPCA(core_genes=core, n_outlier_candidates=8, tight_radius=0.10,
                loose_radius=0.15, fmm_threshold=3.0, max_outlier_fraction=0.05,
                fmm_length_alpha_grid=CFG["fmm_alpha_grid"],
                fmm_length_omega_grid=CFG["fmm_omega_grid"],
                fmm_num_reps=CFG["fmm_reps"], fmm_method="base",
                verbose=False).run(prep.expr_norm)
    encontrados = list(cpca.core_genes_found)
    orden = CircularSynchronizer(anchor_gene=CFG["anchor_gene"],
                                 direction_gene=CFG["direction_gene"],
                                 verbose=False).run(cpca, encontrados)
    return cpca, encontrados, orden


def barrer_unidad(nombre: str, matriz: Path,
                  tiempos_h: dict[str, float]) -> list[dict]:
    cpca, encontrados, orden = etapas_previas(matriz)
    columnas = list(orden.expr_ordered.columns)
    rho_prelim = rho_de([columnas[i] for i in orden.sample_order],
                        orden.circular_scale, tiempos_h)

    # TODO(calibracion): `fit_cache` es andamiaje temporal de TopGeneSelector.
    # Cuando (n_sectors, target_size) queden fijados hay que eliminarlo alli y
    # aqui. Ver el bloque TODO(calibracion) en src/circust/top_genes.py.
    cache: dict = {}
    filas = []
    for n_sec in N_SECTORS:
        for objetivo in TARGET_SIZE:
            top = TopGeneSelector(
                r2_ori_threshold=R2_ORI, r2_fmm_threshold=R2_FMM,
                omega_min=OMEGA_MIN, n_sectors=n_sec, target_size=objetivo,
                fmm_length_alpha_grid=CFG["fmm_alpha_grid"],
                fmm_length_omega_grid=CFG["fmm_omega_grid"],
                fmm_num_reps=CFG["fmm_reps"], fmm_method="base",
                lambda_star=cpca.lambda_star, n_jobs=-1, verbose=False,
            ).run(expr_norm=orden.expr_ordered,
                  circular_scale=orden.circular_scale,
                  seed_genes=encontrados, sample_order=orden.sample_order,
                  fit_cache=cache)

            rob = RobustOrderEstimator(
                n_reps=K, sample_size_fraction=2.0 / 3.0, method="best_k",
                selection_criterion="n_well_fitted", selection_r2_threshold=0.65,
                aggregation_method="auto", r2_min=0.5,
                anchor_gene=CFG["anchor_gene"], direction_gene=CFG["direction_gene"],
                fmm_length_alpha_grid=CFG["fmm_alpha_grid"],
                fmm_length_omega_grid=CFG["fmm_omega_grid"],
                fmm_num_reps=CFG["fmm_reps"], fmm_method="base",
                lambda_star=cpca.lambda_star, n_jobs=-1, seed=SEED, verbose=False,
            ).run(top_result=top, expr_full_norm=orden.expr_ordered,
                  core_genes=encontrados)

            rho_reps = [rho_de([columnas[i] for i in rob.sample_orders[k]],
                               rob.circular_scales[k], tiempos_h)
                        for k in range(rob.sample_orders.shape[0])]
            filas.append({
                "unidad": nombre, "n_sectors": n_sec, "target_size": objetivo,
                "n_muestras": len(columnas), "n_top": len(top.gene_names),
                "rho": rho_de([columnas[i] for i in rob.sample_order],
                              rob.circular_scale, tiempos_h),
                "rho_prelim": rho_prelim,
                "sd_reps": float(np.std(rho_reps)),
            })
    return filas


def analizar(datos: Path) -> None:
    """Resume el barrido con leave-one-out sobre las unidades."""
    salida = ruta_salida(datos)
    if not salida.exists():
        raise SystemExit(f"No existe {salida}. Lanza antes el barrido.")
    d = pd.read_csv(salida)
    unidades = sorted(d["unidad"].unique())
    pd.set_option("display.width", 220)

    print(f"=== {datos.name}: {len(unidades)} unidades, "
          f"{int(d['n_muestras'].median())} muestras (mediana) ===\n")
    tabla = d.pivot_table(index="target_size", columns="n_sectors",
                          values="rho", aggfunc="median")
    print("Mediana de rho por configuracion:")
    print(tabla.round(3).to_string())

    print("\nGenes TOP resultantes (mediana):")
    print(d.pivot_table(index="target_size", columns="n_sectors",
                        values="n_top", aggfunc="median").round(0).to_string())

    print("\nsd de rho entre las K repeticiones (mediana):")
    print(d.pivot_table(index="target_size", columns="n_sectors",
                        values="sd_reps", aggfunc="median").round(4).to_string())

    # ── Leave-one-out sobre las unidades ──────────────────────────────
    # Para cada unidad se elige la configuracion con mayor mediana de rho
    # sobre las OTRAS, y se evalua en la que se dejo fuera. Asi el rho
    # reportado no ha intervenido en la eleccion de su configuracion.
    if len(unidades) < 3:
        print("\nMenos de 3 unidades: no se aplica leave-one-out (no seria "
              "informativo). Mira la tabla de medianas.")
        return

    elegidas, rho_loo = [], []
    for u in unidades:
        resto = d[d["unidad"] != u]
        mejor = resto.groupby(["n_sectors", "target_size"])["rho"].median().idxmax()
        fila  = d[(d["unidad"] == u) & (d["n_sectors"] == mejor[0])
                  & (d["target_size"] == mejor[1])]
        elegidas.append(mejor)
        rho_loo.append(float(fila["rho"].iloc[0]))

    conteo = pd.Series(elegidas).value_counts()
    print("\n=== Leave-one-out ===")
    print(f"Configuraciones elegidas (n_sectors, target_size): {dict(conteo)}")
    print(f"rho LOO       : mediana={np.median(rho_loo):.3f}")
    print(f"rho preliminar: mediana="
          f"{d.groupby('unidad')['rho_prelim'].first().median():.3f}")

    if len(conteo) == 1:
        print("\nLa configuracion elegida es la misma dejando fuera cualquier "
              "unidad: la eleccion no depende de ninguna en concreto.")
    else:
        print("\nAVISO: la configuracion elegida cambia segun la unidad que se "
              "deje fuera. Es senal de que la superficie es plana o ruidosa; "
              "mira la tabla de medianas antes de fijar un valor.")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--datos", type=Path, required=True,
                   help="Directorio del conjunto etiquetado. Cada subdirectorio "
                        "con matriz.csv es una unidad del barrido.")
    p.add_argument("--solo", type=str, default=None,
                   help="Lista separada por comas para restringir las unidades. "
                        "Por ejemplo, en el babuino, los 47 tejidos de la Tabla A.")
    p.add_argument("--analizar", action="store_true",
                   help="No barre: resume el CSV ya calculado.")
    args = p.parse_args()

    datos = args.datos.resolve()
    if not datos.is_dir():
        raise SystemExit(f"No existe el directorio {datos}")

    if args.analizar:
        analizar(datos)
        return

    unidades = descubrir(datos)
    if args.solo:
        permitidas = {s.strip() for s in args.solo.split(",")}
        unidades = [u for u in unidades if u[0] in permitidas]
    if not unidades:
        raise SystemExit(f"No se ha encontrado ninguna unidad en {datos}")
    print(f"{len(unidades)} unidades en {datos.name}")

    salida = ruta_salida(datos)
    salida.parent.mkdir(parents=True, exist_ok=True)
    filas, hechos = [], set()
    if salida.exists():                                   # reanudable
        previo = pd.read_csv(salida)
        filas  = previo.to_dict("records")
        hechos = set(previo["unidad"])
        print(f"Reanudando: {len(hechos)} unidades ya barridas", flush=True)

    n_cfg = len(N_SECTORS) * len(TARGET_SIZE)
    for nombre, matriz, ruta_tiempos in unidades:
        if nombre in hechos:
            continue
        tiempos = pd.read_csv(ruta_tiempos)
        tiempos_h = dict(zip(tiempos["sample_name"], tiempos["true_time_h"]))
        try:
            nuevas = barrer_unidad(nombre, matriz, tiempos_h)
        except Exception as exc:                          # noqa: BLE001
            print(f"[{nombre}] ERROR: {type(exc).__name__}: {exc}", flush=True)
            continue
        filas.extend(nuevas)
        mejor = max(nuevas, key=lambda f: f["rho"])
        print(f"[{nombre}] {n_cfg} cfgs | prelim={nuevas[0]['rho_prelim']:.3f} "
              f"mejor={mejor['rho']:.3f} "
              f"en (n_sec={mejor['n_sectors']}, tgt={mejor['target_size']})",
              flush=True)
        pd.DataFrame(filas).to_csv(salida, index=False)

    print(f"\nGuardado: {salida}")
    print(f"Resumen con: python validation/sweep_top.py --datos {args.datos} --analizar")


if __name__ == "__main__":
    main()
