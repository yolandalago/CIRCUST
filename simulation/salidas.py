"""
Salidas del pipeline en cada celda del estudio, guardadas para reevaluarlas.

El CSV de ``run.py`` lleva las medidas ya calculadas. Si se cambia una medida o
se anade otra, rehacer la rejilla cuesta dias. Por eso, junto al CSV, se
guarda lo que el pipeline produce en cada grupo y en cada metodo de resumen, y
``reevaluar.py`` recalcula cualquier medida a partir de ello sin volver a
correr el pipeline.

Lo que no se guarda es la matriz de expresion. Se regenera exactamente a
partir de la celda de la rejilla (la semilla es determinista) con
:func:`matriz_evaluacion`, que repite el preprocesado y la eliminacion de
atipicos con las muestras que quedaron.

Tablas, cada una en ``<salida>/salidas/<tabla>/`` como partes parquet:

``grupos``        una fila por grupo: lo comun a sus metodos.
``preliminar``    fase del orden preliminar (Etapa 3) por muestra.
``repeticiones``  fase de cada una de las K repeticiones por muestra, ya
                  sincronizadas y orientadas. Es la salida del CIRCUST
                  original en la Etapa 4.
``ajustes_rep``   parametros FMM de cada gen TOP en cada repeticion, las mismas
                  columnas que el atlas que publica la implementacion de R.
``fases``         fase final por muestra, una por metodo de resumen.
``atlas``         atlas final por gen TOP, uno por metodo de resumen.

Las de grupo se indexan por ``grupo_id`` y las de metodo por ``row_id``, el
mismo del CSV. Si una ejecucion se interrumpe entre escribir las partes y el
CSV, al reanudar el grupo se repite y sus filas aparecen dos veces, asi que
:func:`leer` se queda con la ultima.
"""
from __future__ import annotations

import time
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

TABLAS = ("grupos", "preliminar", "repeticiones", "ajustes_rep", "fases", "atlas")

_CLAVE = {
    "grupos":       ["grupo_id"],
    "preliminar":   ["grupo_id", "sample"],
    "repeticiones": ["grupo_id", "k", "sample"],
    "ajustes_rep":  ["grupo_id", "k", "gene"],
    "fases":        ["row_id", "sample"],
    "atlas":        ["row_id", "gene_name"],
}

# Columnas de ``RobustRepetitions.stats_table`` que se guardan, con su nombre.
# Son las del atlas por repeticion de la implementacion de R.
COLUMNAS_AJUSTE = {0: "M", 1: "A", 2: "alpha", 3: "beta", 4: "omega",
                   5: "t_U", 6: "t_L", 11: "R2"}


def grupo_id(bloque: str, cfg, n_core: int, fmm_method: str, rep: int) -> str:
    """Identidad de un grupo: la de sus filas sin el metodo de resumen."""
    return (
        f"{bloque}|n{cfg.n_samples}|{cfg.sampling}"
        f"|w{cfg.omega}|A{cfg.amplitude}|s{cfg.sigma}"
        f"|g{cfg.n_genes}|p{cfg.prop_rhythmic}"
        f"|c{n_core}|{fmm_method}|r{rep}"
    )


def fases_por_muestra(sample_orders: np.ndarray, circ_scales: np.ndarray) -> np.ndarray:
    """``Theta[k, i]``: fase de la muestra ``i`` en la repeticion ``k``.

    ``circ_scales[k]`` viene en el orden de ``sample_orders[k]``, no por
    muestra, asi que comparar columnas de ``circ_scales`` entre repeticiones
    compara muestras distintas.
    """
    K, n = circ_scales.shape
    Theta = np.empty((K, n))
    for k in range(K):
        Theta[k, sample_orders[k]] = circ_scales[k]
    return Theta


@dataclass
class Salidas:
    """Acumulador de las tablas de uno o varios grupos."""

    partes: dict[str, list[pd.DataFrame]] = field(
        default_factory=lambda: {t: [] for t in TABLAS})

    def anadir(self, tabla: str, df: pd.DataFrame) -> None:
        self.partes[tabla].append(df)

    def unir(self, otra: "Salidas") -> None:
        for t in TABLAS:
            self.partes[t].extend(otra.partes[t])

    def vacia(self) -> bool:
        return not any(self.partes.values())

    def escribir(self, raiz: Path) -> None:
        """Escribe una parte parquet por tabla y vacia el acumulador."""
        sello = time.time_ns()
        for t in TABLAS:
            if not self.partes[t]:
                continue
            carpeta = Path(raiz) / "salidas" / t
            carpeta.mkdir(parents=True, exist_ok=True)
            df = pd.concat(self.partes[t], ignore_index=True)
            for c in df.select_dtypes("float64").columns:
                df[c] = df[c].astype("float32")
            df.to_parquet(carpeta / f"parte-{sello}.parquet", index=False)
            self.partes[t] = []


def salidas_grupo(gid: str, estado, reps, t_reps: float) -> Salidas:
    """Las tablas comunes a todos los metodos de un grupo."""
    s = Salidas()
    muestras = list(estado.sample_names)
    s.anadir("grupos", pd.DataFrame([{
        "grupo_id":    gid,
        "n_muestras":  len(muestras),
        "n_top":       estado.n_top,
        "n_outliers":  estado.n_outliers,
        "lambda_star": estado.lambda_star,
        "top_genes":   ",".join(reps.names_top),
        "core_genes":  ",".join(estado.core_genes),
        "flipped":     ",".join(str(int(f)) for f in reps.flipped),
        "t_reps_s":    t_reps,
    }]))
    s.anadir("preliminar", pd.DataFrame({
        "grupo_id": gid, "sample": muestras,
        "theta_pre": np.asarray(estado.theta_pre, dtype=float),
    }))
    Theta = fases_por_muestra(reps.sample_orders, reps.circ_scales)
    K, n = Theta.shape
    s.anadir("repeticiones", pd.DataFrame({
        "grupo_id": gid,
        "k":        np.repeat(np.arange(K), n),
        "sample":   muestras * K,
        "phase":    Theta.ravel(),
    }))
    n_top = len(reps.names_top)
    ajustes = pd.DataFrame(reps.stats_table[:, list(COLUMNAS_AJUSTE)],
                           columns=list(COLUMNAS_AJUSTE.values()))
    ajustes.insert(0, "gene", list(reps.names_top) * K)
    ajustes.insert(0, "k", np.repeat(np.arange(K), n_top))
    ajustes.insert(0, "grupo_id", gid)
    s.anadir("ajustes_rep", ajustes)
    return s


def salidas_metodo(row_id: str, gid: str, etapa4) -> Salidas:
    """Las tablas propias de un metodo de resumen."""
    s = Salidas()
    s.anadir("fases", pd.DataFrame({
        "row_id": row_id, "grupo_id": gid,
        "sample": list(etapa4.theta_hat.index),
        "theta_hat": etapa4.theta_hat.to_numpy(dtype=float),
    }))
    if etapa4.atlas is not None and not etapa4.atlas.empty:
        atlas = etapa4.atlas.copy()
        atlas.insert(0, "grupo_id", gid)
        atlas.insert(0, "row_id", row_id)
        s.anadir("atlas", atlas)
    return s


def leer(raiz: Path, tabla: str) -> pd.DataFrame:
    """Una tabla completa, sin las filas repetidas por reanudaciones."""
    carpeta = Path(raiz) / "salidas" / tabla
    partes = sorted(carpeta.glob("parte-*.parquet"))
    if not partes:
        raise FileNotFoundError(f"No hay partes de {tabla!r} en {carpeta}")
    df = pd.concat([pd.read_parquet(p) for p in partes], ignore_index=True)
    return df.drop_duplicates(_CLAVE[tabla], keep="last").reset_index(drop=True)


def matriz_evaluacion(matriz: pd.DataFrame, muestras: list[str]) -> pd.DataFrame:
    """
    La matriz normalizada sobre la que se evaluo el grupo, ``expr_full_norm``.

    Repite lo que hace el pipeline: preprocesado de la matriz bruta, quitar
    las muestras atipicas (las que no estan en ``muestras``) y volver a
    normalizar cada gen a [-1, 1]. ``matriz`` es la que devuelve
    ``generate_dataset`` con la configuracion de la celda.
    """
    from circust.preprocessing import Preprocessor, normalise_matrix
    from run_pipeline import DEFAULT_CONFIG

    prep = Preprocessor(sparse_threshold=DEFAULT_CONFIG["sparse_threshold"],
                        verbose=False).run(matriz)
    return normalise_matrix(prep.expr_norm[list(muestras)])
