#!/usr/bin/env python
"""
scripts/check_env.py
====================
Comprueba que el entorno tiene todo lo necesario para ejecutar CIRCUST y
dice, para cada pieza, si esta, que version hay y como instalarla si falta.

Uso
---
    python scripts/check_env.py            # informe
    python scripts/check_env.py --strict   # ademas sale con codigo 1 si falta
                                           # alguna dependencia obligatoria

Tres niveles:

  obligatorias   sin ellas el algoritmo no arranca
  aceleracion    Numba: mismos resultados, FMM base 2-3 veces mas rapido
  cuda           ruta GPU opcional del FMM base (desactivada por defecto)
"""
from __future__ import annotations

import argparse
import importlib
import platform
import sys
from importlib import metadata

# (modulo importable, nombre del paquete pip, version minima)
OBLIGATORIAS = [
    ("numpy",      "numpy",        "1.24"),
    ("scipy",      "scipy",        "1.10"),
    ("pandas",     "pandas",       "2.0"),
    ("sklearn",    "scikit-learn", "1.2"),
    ("matplotlib", "matplotlib",   "3.7"),
    ("joblib",     "joblib",       "1.2"),
    ("pyarrow",    "pyarrow",      "12.0"),
]
ACELERACION = [
    ("numba", "numba", "0.59"),
]
CUDA = [
    ("numba_cuda",       "numba-cuda",           "0.10"),
    (None,               "nvidia-cuda-nvcc-cu12", "12.8"),
]
PY_MIN = (3, 10)


def _version(pip_name: str, module: str | None) -> str | None:
    """Version instalada del paquete, o None si no esta."""
    try:
        return metadata.version(pip_name)
    except metadata.PackageNotFoundError:
        pass
    if module is None:
        return None
    try:
        mod = importlib.import_module(module)
    except Exception:
        return None
    return getattr(mod, "__version__", "?")


def _older(v: str | None, vmin: str) -> bool:
    if v is None or v == "?":
        return False
    def key(s: str):
        out = []
        for p in s.split("."):
            digits = "".join(ch for ch in p if ch.isdigit())
            out.append(int(digits) if digits else 0)
        return out
    return key(v) < key(vmin)


def _report(bloque: str, items, pip_extra: str) -> tuple[list[str], list[str]]:
    """Imprime el bloque; devuelve (faltan, anticuadas) como nombres pip."""
    faltan, viejas = [], []
    print(f"\n{bloque}")
    for module, pip_name, vmin in items:
        v = _version(pip_name, module)
        if v is None:
            print(f"  [FALTA]  {pip_name:24s} (>= {vmin})")
            faltan.append(pip_name)
        elif _older(v, vmin):
            print(f"  [VIEJA]  {pip_name:24s} {v}  (se necesita >= {vmin})")
            viejas.append(pip_name)
        else:
            print(f"  [OK]     {pip_name:24s} {v}")
    if faltan or viejas:
        print(f"  -> pip install -U {' '.join(faltan + viejas)}"
              f"      (o: pip install -e \".[{pip_extra}]\")" if pip_extra else
              f"  -> pip install -U {' '.join(faltan + viejas)}")
    return faltan, viejas


def _numba_cuda_status() -> str:
    """Que ve Numba respecto a la GPU, sin lanzar excepciones al usuario."""
    try:
        from numba import cuda
    except Exception as e:                       # numba no instalado
        return f"no evaluable ({type(e).__name__})"
    try:
        if not cuda.is_available():
            # Distinguir "no hay GPU" de "hay GPU pero falta el toolkit"
            try:
                from numba.cuda.cudadrv import driver as drv
                drv.driver  # fuerza la carga de libcuda
                n = len(cuda.list_devices())
            except Exception:
                n = 0
            if n == 0:
                return "sin GPU CUDA visible (o sin driver): la ruta CUDA no se puede usar"
            return ("GPU visible pero Numba no encuentra el toolkit CUDA "
                    "(libcudart / NVVM): pip install -r requirements-cuda.txt")
        dev = cuda.get_current_device()
        name = dev.name.decode() if isinstance(dev.name, bytes) else str(dev.name)
        cc = dev.compute_capability
        return f"disponible: {name}, capacidad {cc[0]}.{cc[1]}"
    except Exception as e:
        return f"error al inicializar CUDA: {type(e).__name__}: {e}"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[3])
    ap.add_argument("--strict", action="store_true",
                    help="codigo de salida 1 si falta alguna dependencia obligatoria")
    args = ap.parse_args()

    print(f"Python {platform.python_version()} ({sys.executable})")
    py_ok = sys.version_info[:2] >= PY_MIN
    if not py_ok:
        print(f"  [ERROR] se necesita Python >= {PY_MIN[0]}.{PY_MIN[1]}")

    faltan_obl, viejas_obl = _report("Obligatorias (algoritmo)", OBLIGATORIAS, "")
    faltan_acc, _          = _report("Aceleracion en CPU (opcional, recomendada)", ACELERACION, "accel")
    faltan_cuda, _         = _report("CUDA (opcional, desactivada por defecto)", CUDA, "cuda")

    print("\nEstado efectivo de CIRCUST en este entorno")
    if faltan_obl or viejas_obl or not py_ok:
        print("  algoritmo:      NO ejecutable hasta instalar lo marcado como FALTA/VIEJA")
    else:
        print("  algoritmo:      ejecutable")
    if faltan_acc:
        print("  aceleracion:    sin Numba -> camino NumPy puro (mismos resultados, "
              "FMM base 2-3x mas lento)")
    else:
        print("  aceleracion:    Numba activo (la primera ejecucion compila y cachea "
              "en __pycache__)")
    if faltan_acc:
        print("  cuda:           no evaluable sin Numba")
    else:
        print(f"  cuda:           {_numba_cuda_status()}")
        if not faltan_cuda:
            print("                  (recordatorio: solo se usa con FMMModel(use_cuda=True); "
                  "en GPU de consumo la CPU es mas rapida)")

    if args.strict and (faltan_obl or viejas_obl or not py_ok):
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
