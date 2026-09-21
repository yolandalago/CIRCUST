"""
Bloque computacional: coste de la implementacion en Python frente a la de R.

La metrica es el tiempo, no la calidad. Va aparte de la rejilla unificada
porque la implementacion original en R es unas 20 veces mas lenta y necesita
un proceso por celda, asi que no se puede cruzar con los demas factores. El
eje es el tamano de la matriz, que es lo que gobierna el coste.

Las dos implementaciones corren sobre la MISMA matriz sintetica, con la
configuracion comparable: FMM base y el resumen de las K repeticiones que hace
el CIRCUST original. El lado de Python usa ``runner.py`` como el resto del
estudio. El lado de R usa ``rscripts/benchmark_timing.R``, que cronometra las
mismas cuatro etapas.

Por que hace falta renombrar
----------------------------
La implementacion de R no recibe los genes core como dato, los tiene escritos
en el codigo, y ``robustSincroDBP_v3_cores`` depende de que exista
literalmente ``DBP``. Antes que operar el codigo de R, se exporta una segunda
copia de la misma matriz con las filas core renombradas a los nombres que R
espera: ``CORE_A`` pasa a ``ARNTL`` porque es el ancla, ``CORE_B`` a ``DBP``
porque es el gen de direccion, y el gen intermedio de la cadena de
consistencia a ``CRY1``, con lo que se preserva ``DBP < CRY1 < ARNTL``, que
es justo nuestra cadena. Los valores no cambian, solo las etiquetas de fila.

Uso
---
    python -m simulation.rbench --reps 3 --out output/simulacion
    python -m simulation.rbench --reps 3 --out output/simulacion --run-r
    python -m simulation.rbench --consolidar-r --out output/simulacion

El primero exporta las celdas y cronometra Python. El segundo hace ademas la
pasada de R celda a celda. El tercero solo recoge los ``timing_r.csv`` que R
haya dejado y los vuelca en ``tiempos.csv``.
"""
from __future__ import annotations

import argparse
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from simulation.generate import generate_dataset
from simulation.grid import build_grid
from simulation.runner import RobustConfig, run_stage_4, run_stages_0_3
from simulation.truth import ANCHOR_NAME, OPPOSITE_NAME, SimConfig, SimulationTruth

# Los 12 nombres que benchmark_timing.R tiene escritos a mano en su ``coreG``.
# El orden importa: los tres primeros de la asignacion son fijos.
R_CORE_NAMES: tuple[str, ...] = (
    "ARNTL", "DBP", "CRY1",
    "PER1", "PER2", "PER3", "CRY2", "CLOCK", "NR1D1", "RORA", "TEF", "STAT3",
)

SCRIPT_R = PROJECT_ROOT / "rscripts" / "benchmark_timing.R"

# Etapas que cronometra el lado de R, en el orden del pipeline. El lado de
# Python guarda las mismas con estos nombres en ``estado.timing``.
ETAPAS = ("preprocess", "preliminary_order", "top_genes", "robust_order")


def rename_map_for_r(truth: SimulationTruth) -> dict[str, str]:
    """
    Correspondencia entre los genes core sinteticos y los que espera R.

    Las tres asignaciones que no son arbitrarias: ``CORE_A`` a ``ARNTL``
    (ancla, pico en pi), ``CORE_B`` a ``DBP`` (direccion, pico en 0) y el gen
    intermedio de la cadena a ``CRY1``. Con eso se preserva
    ``DBP < CRY1 < ARNTL``. Los demas core se reparten por los nueve nombres
    restantes en orden.
    """
    if len(truth.core_names) < len(R_CORE_NAMES):
        raise ValueError(
            f"El renombrado para R necesita al menos 12 genes core, la celda "
            f"tiene {len(truth.core_names)}"
        )
    intermedio = truth.chain[1]
    mapa = {
        ANCHOR_NAME:   "ARNTL",
        OPPOSITE_NAME: "DBP",
        intermedio:    "CRY1",
    }
    restantes = [g for g in truth.core_names if g not in mapa]
    for sintetico, real in zip(restantes, R_CORE_NAMES[3:]):
        mapa[sintetico] = real
    return mapa


def export_for_r(
    matrix:  pd.DataFrame,
    truth:   SimulationTruth,
    out_dir: Path,
) -> dict[str, Path]:
    """
    Escribe la celda para los dos lados del benchmark.

    Cuatro ficheros: la matriz con los nombres sinteticos (lado Python), la
    lista de genes core declarados, los pares de consistencia, y una segunda
    copia de la **misma** matriz con las filas core renombradas a lo que
    espera R. Los tiempos verdaderos **no** se exportan.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    rutas = {
        "matriz":     out_dir / "matriz.csv",
        "core_genes": out_dir / "core_genes.txt",
        "pares":      out_dir / "pares.csv",
        "matriz_r":   out_dir / "matriz_r.csv",
    }
    matrix.to_csv(rutas["matriz"])
    tamano = max(truth.core_sets)
    rutas["core_genes"].write_text(
        "\n".join(truth.core_genes(tamano)) + "\n", encoding="utf-8",
    )
    pd.DataFrame(truth.consistency_pairs(),
                 columns=["antes", "despues"]).to_csv(rutas["pares"], index=False)
    matrix.rename(index=rename_map_for_r(truth)).to_csv(rutas["matriz_r"])
    return rutas


def time_python(matrix: pd.DataFrame, truth: SimulationTruth,
                n_core: int, celda_dir: Path | None = None) -> dict[str, float]:
    """
    Tiempos por etapa de una ejecucion completa en Python, FMM base y el
    metodo de resumen de referencia, sin construir el atlas: el R original
    tampoco lo construye dentro de las etapas que se cronometran.

    Con ``celda_dir`` vuelca ademas las dos salidas que se comparan con R:

    ``py_fases_preliminar.csv``
        La fase por muestra del ORDEN PRELIMINAR, que es la salida homologa a
        la que reporta R. El CIRCUST original resume las K repeticiones con
        una mediana de parametros FMM y no las colapsa en un orden de
        muestras, asi que la unica fase por muestra que publica es esta.
    ``py_fases_robusto.csv``
        La fase por muestra tras la Etapa 4. No tiene equivalente en R, y esa
        ausencia es justamente el hueco que llena este trabajo.

    El volcado va fuera de las dos regiones cronometradas, asi que no
    contamina los tiempos.
    """
    t = time.perf_counter()
    estado = run_stages_0_3(matrix, truth, n_core, "base")
    t_0_3 = time.perf_counter() - t

    t = time.perf_counter()
    etapa4 = run_stage_4(estado, truth, RobustConfig(), "base", con_atlas=False)
    t_4 = time.perf_counter() - t

    if celda_dir is not None:
        celda_dir = Path(celda_dir)
        celda_dir.mkdir(parents=True, exist_ok=True)
        pd.DataFrame({"sample": estado.sample_names,
                      "phase": np.asarray(estado.theta_pre, dtype=float)}) \
          .to_csv(celda_dir / "py_fases_preliminar.csv", index=False)
        etapa4.theta_hat.rename_axis("sample").reset_index(name="phase") \
              .to_csv(celda_dir / "py_fases_robusto.csv", index=False)
        pd.DataFrame({"sample": truth.sample_names,
                      "t_true": truth.t_true}) \
          .to_csv(celda_dir / "verdad.csv", index=False)

    tiempos = {k: float(v) for k, v in estado.timing.items()}
    tiempos["robust_order"] = t_4
    tiempos["stages_0_3"] = t_0_3
    tiempos["total"] = t_0_3 + t_4
    return tiempos


def run_r(celda_dir: Path, etiqueta: str) -> Path:
    """
    Lanza el benchmark de R sobre la matriz renombrada de una celda.

    Con ``--dump`` R escribe ademas ``r_fases.csv`` y ``r_atlas.csv`` en el
    directorio de la celda, que es lo que permite comparar su precision con la
    de Python y no solo su velocidad.

    Devuelve la ruta del ``timing_r.csv`` que R escribe. Un fallo de R lanza
    ``CalledProcessError`` con la salida capturada, para que quede en el log.
    """
    celda_dir = Path(celda_dir)
    destino = celda_dir / "timing_r.csv"
    subprocess.run(
        ["Rscript", str(SCRIPT_R),
         "--matrix", str(celda_dir / "matriz_r.csv"),
         "--reps", "1",
         "--label", etiqueta,
         "--out", str(destino),
         "--dump", str(celda_dir)],
        cwd=PROJECT_ROOT, check=True, capture_output=True, text=True,
    )
    return destino


def consolidar_r(raiz: Path) -> pd.DataFrame:
    """
    Recoge los ``timing_r.csv`` de cada celda y los pasa al formato ancho de
    ``tiempos.csv``, una fila por celda y repeticion con ``t_<etapa>_s``.
    """
    trozos = []
    for f in sorted(Path(raiz).glob("*/timing_r.csv")):
        t = pd.read_csv(f)
        ancho = t.pivot_table(index=["celda", "rep"], columns="stage",
                              values="seconds").reset_index()
        ancho.columns = [c if c in ("celda", "rep") else f"t_{c}_s"
                         for c in ancho.columns]
        ancho = ancho.rename(columns={"t_elapsed_total_s": "t_total_s"})
        ancho["implementacion"] = "r"
        partes = f.parent.name.split("_")
        ancho["n_genes"]   = int(partes[0][1:])
        ancho["n_samples"] = int(partes[1][1:])
        # la replica del CSV de R es la pasada interna (siempre 1); la del
        # estudio va en el nombre de la celda
        ancho["rep"] = int(partes[2][1:])
        trozos.append(ancho)
    if not trozos:
        return pd.DataFrame()
    return pd.concat(trozos, ignore_index=True)


def _nombre(cfg: SimConfig, rep: int) -> str:
    return f"g{cfg.n_genes}_n{cfg.n_samples}_r{rep}"


def _anadir(destino: Path, filas: pd.DataFrame) -> None:
    """
    Anade filas al CSV alineando las columnas con las que ya tiene.

    ``to_csv`` en modo append no reordena nada, escribe los valores en el
    orden en que vienen. Las filas de R y las de Python se construyen por
    caminos distintos y no salen en el mismo orden, asi que sin este
    realineado los tiempos de R aterrizan bajo las cabeceras de Python y el
    fichero queda inservible. Paso ya una vez.
    """
    if filas.empty:
        return
    if destino.exists():
        columnas = list(pd.read_csv(destino, nrows=0).columns)
        for falta in [c for c in filas.columns if c not in columnas]:
            columnas.append(falta)
        filas = filas.reindex(columns=columnas)
    filas.to_csv(destino, mode="a", header=not destino.exists(), index=False)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Bloque computacional: exporta las celdas para R, "
                    "cronometra Python y opcionalmente R.",
    )
    parser.add_argument("--reps", type=int, default=3)
    parser.add_argument("--out", type=str, default="output/simulacion")
    parser.add_argument("--export-only", action="store_true",
                        help="Solo escribe las celdas para R, sin cronometrar.")
    parser.add_argument("--run-r", action="store_true",
                        help="Lanza tambien el benchmark de R en cada celda.")
    parser.add_argument("--consolidar-r", action="store_true",
                        help="Solo recoge los timing_r.csv existentes.")
    args = parser.parse_args()

    salida = Path(args.out)
    raiz = salida / "rbench"
    raiz.mkdir(parents=True, exist_ok=True)
    destino = salida / "tiempos.csv"

    if args.consolidar_r:
        r = consolidar_r(raiz)
        _anadir(destino, r)
        print(f"Anadidas {len(r)} filas de R a {destino}")
        return

    filas_py = []
    filas_r  = []
    for celda in build_grid("computacional", args.reps):
        nombre = _nombre(celda.cfg, celda.rep)
        matriz, truth = generate_dataset(celda.cfg)
        export_for_r(matriz, truth, raiz / nombre)
        if args.export_only:
            continue

        tiempos = time_python(matriz, truth, celda.n_core_declared,
                              celda_dir=raiz / nombre)
        filas_py.append(dict(
            implementacion = "python",
            celda          = nombre,
            n_genes        = celda.cfg.n_genes,
            n_samples      = celda.cfg.n_samples,
            rep            = celda.rep,
            **{f"t_{k}_s": v for k, v in tiempos.items()},
        ))
        print(f"{nombre}: python {tiempos['total']:.2f} s", flush=True)

        if args.run_r:
            try:
                t0 = time.perf_counter()
                run_r(raiz / nombre, nombre)
                print(f"{nombre}: r {time.perf_counter() - t0:.1f} s", flush=True)
            except subprocess.CalledProcessError as exc:
                print(f"{nombre}: R FALLO\n{exc.stderr[-2000:]}", flush=True)

    _anadir(destino, pd.DataFrame(filas_py))
    if args.run_r:
        filas_r = consolidar_r(raiz)
        _anadir(destino, filas_r)
    print(f"Escritas {len(filas_py)} filas de Python y {len(filas_r)} de R "
          f"en {destino}")
    print(f"Celdas exportadas para R en {raiz}")


if __name__ == "__main__":
    main()
