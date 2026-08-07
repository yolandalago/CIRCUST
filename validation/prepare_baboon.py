#!/usr/bin/env python3
"""
prepare_baboon.py
=================
Extrae del atlas de babuino (Mure et al. 2018, GSE98965) un dataset por
**tejido** en el formato que consumen ``run_pipeline.py`` y ``validate.py``,
exactamente igual que ``generate_synthetic.write_dataset`` hace para el grid
sintetico. Esta es la capa de adaptacion que permite validar el CIRCUST de
Python sobre datos reales con tiempo conocido (el ZT del muestreo), como hizo
el paper original (Larriba et al. 2023, Fig. 7A).

Estructura del CSV de entrada
-----------------------------
    EnsemblID,Symbol,ADC.ZT00,ADC.ZT02,...,ADM.ZT00,...

- Dos columnas de ID: ``EnsemblID`` (se descarta) y ``Symbol`` (pasa a indice
  de genes; es donde estan ARNTL, DBP, PER1... que necesita CIRCUST).
- Columnas de muestra ``TEJIDO.ZTxx``. Peculiaridades a parsear:
    * algunos tejidos usan ``_ZT`` en vez de ``.ZT`` (COR, IRI, RPE, SKI),
    * algunas columnas llevan sufijo ``.approx`` (p.ej. ``ANT.ZT04.approx``).
- 64 tejidos, cada uno con 11-12 puntos ZT cada 2 h (ZT00..ZT22 -> 0..22 h).
  **El ZT es el ground truth temporal directo.**

Salida (una carpeta por tejido bajo ``<out>/_tissues/<TEJIDO>/``)
-----------------------------------------------------------------
    matrix.csv       indice = Symbol, columnas = las 11-12 muestras del tejido
                     con su nombre original (``LIV.ZT06``...). Consumible por
                     run_pipeline.py (gen column por defecto = primera columna).
    true_times.csv   sample_name, true_time_h, true_phase_h, true_phase_rad.
                     Mismo esquema que generate_synthetic; lo cruza validate.py.

Ejemplos
--------
    # Extrae los 64 tejidos
    python validation/prepare_baboon.py -o output/baboon

    # Solo un subconjunto (para smoke test)
    python validation/prepare_baboon.py -o output/baboon --tissues LIV,KIM
"""
from __future__ import annotations

import argparse
import json
import re
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_CSV  = PROJECT_ROOT / "data" / "GSE98965_baboon_tissue_expression_FPKM.csv"
PERIOD_HOURS = 24.0

# Separador ``.`` o ``_`` antes de ``ZT``; sufijo ``.approx`` opcional.
#   ADC.ZT00          -> ('ADC', 0)
#   COR_ZT22          -> ('COR', 22)
#   ANT.ZT04.approx   -> ('ANT', 4)
_SAMPLE_RE = re.compile(r"^(?P<tissue>.+?)[._]ZT(?P<zt>\d+)(?:\.approx)?$")


def parse_sample(col: str) -> tuple[str, int] | None:
    """``'LIV.ZT06'`` -> ``('LIV', 6)``. Devuelve ``None`` si la columna no es
    una muestra (p.ej. las columnas de ID ``EnsemblID`` / ``Symbol``)."""
    m = _SAMPLE_RE.match(str(col))
    if m is None:
        return None
    return m.group("tissue"), int(m.group("zt"))


def load_baboon(csv_path: Path) -> pd.DataFrame:
    """Carga el CSV grande UNA vez: descarta ``EnsemblID`` y fija ``Symbol``
    como indice de genes. Los Symbols duplicados/sin nombre se dejan tal cual;
    el Preprocessor del pipeline ya los resuelve (mayor MAD) y filtra."""
    df = pd.read_csv(csv_path)
    if "Symbol" not in df.columns:
        raise SystemExit(
            f"'{csv_path}' no tiene columna 'Symbol'. "
            f"Columnas iniciales: {list(df.columns[:4])}"
        )
    df = df.drop(columns=[c for c in ("EnsemblID",) if c in df.columns])
    df = df.set_index("Symbol")
    df.index.name = "gene_symbol"
    return df


def group_by_tissue(columns) -> dict[str, list[str]]:
    """Agrupa las columnas de muestra por tejido, preservando el orden ZT."""
    tissues: dict[str, list[str]] = {}
    for col in columns:
        parsed = parse_sample(col)
        if parsed is None:
            continue
        tissue, _zt = parsed
        tissues.setdefault(tissue, []).append(col)
    # Ordenar las muestras de cada tejido por su ZT (orden temporal real)
    for tissue, cols in tissues.items():
        cols.sort(key=lambda c: parse_sample(c)[1])
    return tissues


def write_tissue(
    df:      pd.DataFrame,
    tissue:  str,
    cols:    list[str],
    out_dir: Path,
) -> dict:
    """Escribe ``matrix.csv`` + ``true_times.csv`` de un tejido. Devuelve un
    dict resumen (tejido, n muestras, ZTs)."""
    out_dir.mkdir(parents=True, exist_ok=True)

    # 1) Matriz genes x muestras (nombres originales -> casan con el GT)
    sub = df[cols].copy()
    sub.index.name = "gene_symbol"
    sub.to_csv(out_dir / "matrix.csv")

    # 2) Ground truth temporal a partir del ZT
    zt_h      = np.array([parse_sample(c)[1] for c in cols], dtype=np.float64)
    phase_h   = zt_h % PERIOD_HOURS
    phase_rad = phase_h * (2.0 * np.pi / PERIOD_HOURS)
    pd.DataFrame({
        "sample_name":    cols,
        "true_time_h":    zt_h,
        "true_phase_h":   phase_h,
        "true_phase_rad": phase_rad,
    }).to_csv(out_dir / "true_times.csv", index=False)

    return {"tissue": tissue, "n_samples": len(cols),
            "zt": [int(z) for z in zt_h]}


# ===========================================================================
# API de alto nivel (reutilizable por el orquestador)
# ===========================================================================

def extract_all(
    csv_path:    Path,
    out_root:    Path,
    only:        list[str] | None = None,
    verbose:     bool             = True,
) -> Path:
    """Extrae los datasets por tejido a ``<out_root>/_tissues/<TEJIDO>/`` y
    devuelve la ruta ``_tissues``. ``only`` restringe a un subconjunto."""
    tissues_dir = out_root / "_tissues"
    tissues_dir.mkdir(parents=True, exist_ok=True)

    if verbose:
        print(f"  Cargando atlas: {csv_path}  ...", flush=True)
    t0 = time.time()
    df = load_baboon(csv_path)
    if verbose:
        print(f"  Cargado: {df.shape[0]} genes x {df.shape[1]} columnas "
              f"de muestra  ({time.time() - t0:.1f}s)")

    groups = group_by_tissue(df.columns)
    if only:
        groups = {t: c for t, c in groups.items() if t in set(only)}
        missing = set(only) - set(groups)
        if missing:
            raise SystemExit(f"Tejidos no encontrados en el atlas: {sorted(missing)}")

    summaries = []
    for tissue, cols in sorted(groups.items()):
        summary = write_tissue(df, tissue, cols, tissues_dir / tissue)
        summaries.append(summary)
        if verbose:
            print(f"    {tissue:6s}  n={summary['n_samples']:2d}  "
                  f"ZT={summary['zt'][0]:02d}..{summary['zt'][-1]:02d}")

    # Indice maestro para el orquestador
    manifest = {
        "source_csv": str(csv_path),
        "period_h":   PERIOD_HOURS,
        "n_tissues":  len(summaries),
        "tissues":    summaries,
    }
    (tissues_dir / "tissues_manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8")

    if verbose:
        print(f"  ✓ {len(summaries)} tejidos extraidos en {tissues_dir.resolve()}")
    return tissues_dir


# ===========================================================================
# CLI
# ===========================================================================

def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=("Extrae un dataset por tejido del atlas de babuino "
                     "(Mure 2018) para validar CIRCUST contra el ZT real."),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("-o", "--output", required=True, type=str,
                   help="Directorio raiz (los tejidos van a <output>/_tissues/).")
    p.add_argument("--data", default=str(DEFAULT_CSV), type=str,
                   help=f"CSV del atlas. Defecto: {DEFAULT_CSV.name}")
    p.add_argument("--tissues", default=None, type=str,
                   help="Subconjunto de tejidos separado por comas "
                        "(p.ej. LIV,KIM). Sin flag: los 64.")
    return p.parse_args()


def main() -> None:
    args = _parse_args()
    csv_path = Path(args.data)
    if not csv_path.exists():
        raise SystemExit(f"No existe el CSV del atlas: {csv_path}")

    only = ([t.strip() for t in args.tissues.split(",") if t.strip()]
            if args.tissues else None)

    print("╔══════════════════════════════════════════════════════════════════╗")
    print("  Preparacion del atlas de babuino (Mure 2018) para CIRCUST")
    print("╚══════════════════════════════════════════════════════════════════╝")
    extract_all(csv_path, Path(args.output), only=only, verbose=True)


if __name__ == "__main__":
    sys.exit(main())
