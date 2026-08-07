#!/usr/bin/env python3
"""
generate_synthetic.py
=====================
Genera un dataset sintetico de expresion genica circadiana y lo guarda
**en disco** como artefactos que pueden alimentar el resto del flujo.

Diseno factorial
----------------
Este modulo esta pensado para un **estudio de simulacion disenado**: en
lugar de muestrear los parametros de cada gen de rangos aleatorios, se
fijan ``M = 0`` (nivel basal) y ``A = 1`` (amplitud), y se cruzan como
**factores** de escenario:

- **muestreo temporal** (``sampling``): equiespaciado, equiespaciado con
  jitter, agrupado (clustered) sin huecos, agrupado con huecos. Representan
  grados crecientes de irregularidad en la cobertura del ciclo circadiano.
- **omega** (``omega``): grado de asimetria de la senal FMM. ``omega = 1``
  reduce el FMM a una sinusoide pura (referencia cosinor).
- **sigma** (``sigma``): intensidad de la senal. El ruido es
  ``(A/sigma)*N(0,1)``, de modo que sigma alto = menos ruido.
- **proporcion ritmica** (``prop_rhythmic``) y **tamano muestral**
  (``n_samples``).

Forma de la senal (FMM, Rueda et al. 2019)
------------------------------------------
``f(t) = M + A * cos(beta + 2*arctan(omega * tan((t_n - alpha)/2))) + ruido``
permite picos asimetricos y estrechos. Es el modelo que CIRCUST ajusta
internamente, asi que validar con datos generados de esta forma evalua la
capacidad de CIRCUST en su escenario de diseno. Los angulos ``alpha`` y
``beta`` se sortean ``U(0, 2*pi)`` por gen (no se favorece ninguna
localizacion de pico). Con ``omega = 1`` el FMM coincide con
``A*cos(2*pi*t/24 - alpha + beta)``, una sinusoide pura.

Salidas a disco
---------------
    ``matrix.csv``        Matriz genes x muestras. Formato compatible
                          con ``run_pipeline.py`` (primera columna =
                          gene_symbol, resto = nombres de muestras).

    ``true_times.csv``    Ground truth temporal:
                              sample_name   : nombre de la muestra
                              true_time_h   : tiempo real en horas
                                              (sin plegar — puede exceder
                                              24 h si time_span > 24)
                              true_phase_h  : tiempo real plegado a [0, 24)
                              true_phase_rad: misma fase en [0, 2*pi)
                          Es el archivo que consume ``validate.py``.

    ``gene_params.csv``   Parametros usados por gen al generar (model,
                          amp, snr=sigma, alpha/beta/omega). Util para
                          diagnostico post-hoc.

    ``manifest.json``     Configuracion completa del generador
                          (reproducibilidad).

Workflow tipico
---------------
    # 1) Generar
    python scripts/generate_synthetic.py --dataset SynDST4 --seed 42 \\
            -o output/synth/SynDST4

    # 2) Correr CIRCUST sobre la matriz generada
    python scripts/run_pipeline.py \\
            --data output/synth/SynDST4/matrix.csv \\
            -o output/run/SynDST4

    # 3) Validar contra el ground truth
    python scripts/validate.py \\
            --predicted output/run/SynDST4/sample_order.csv \\
            --true      output/synth/SynDST4/true_times.csv \\
            -o output/validate/SynDST4

Notas
-----
- La **periodo del ritmo** es siempre 24 h, aunque la ventana de muestreo
  (``time_span``) sea 14, 18, 24 o 48 h.
- Para datasets de 48 h los tiempos cubren dos ciclos; ``true_times.csv``
  guarda tanto el tiempo bruto (``true_time_h``) como la fase plegada
  modulo 24 (``true_phase_h``, ``true_phase_rad``).
"""
from __future__ import annotations

import argparse
import json
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd


# ===========================================================================
# Constantes
# ===========================================================================

PERIOD_HOURS: float = 24.0

#: Genes semilla CIRCUST que el pipeline buscara siempre. El generador
#: los marca como ritmicos y los coloca primero en la matriz.
SEED_GENES: list[str] = [
    "PER1", "PER2", "PER3",
    "CRY1", "CRY2",
    "ARNTL", "CLOCK",
    "NR1D1", "RORA",
    "DBP", "TEF", "STAT3",
]

#: Esquemas de muestreo temporal soportados (grados crecientes de
#: irregularidad en la cobertura del ciclo de 24 h).
SAMPLING_SCHEMES: tuple[str, ...] = (
    "equispaced",     # uniforme sobre [0, 24)
    "jitter",         # equiespaciado + perturbacion uniforme pequena
    "clustered",      # concentracion local, sin huecos
    "clustered_gap",  # concentracion local + ventana(s) infrarepresentada(s)
)

#: Parametros por defecto de cada esquema de muestreo. Se pueden
#: sobreescribir parcialmente via ``sampling_params``.
SAMPLING_DEFAULTS: dict[str, dict] = {
    "equispaced":    dict(),
    "jitter":        dict(jitter_frac=0.4),
    "clustered":     dict(frac=0.4, center=12.0, sd=2.0),
    "clustered_gap": dict(frac=0.4, center=12.0, sd=2.0, gaps=((0.0, 6.0),)),
}

#: Niveles de los factores del estudio factorial (los importa el runner
#: del grid). ``sigma`` se confirma en el exploratorio de tiempos.
FACTOR_LEVELS: dict[str, tuple] = {
    "sampling":      SAMPLING_SCHEMES,
    "omega":         (0.1, 0.5, 0.9, 1.0),
    "sigma":         (0.5, 1.0, 2.0),
    "prop_rhythmic": (0.01, 0.05),
    "n_samples":     (24, 48, 75, 96),
}


# ===========================================================================
# Generador
# ===========================================================================

@dataclass
class SyntheticDataset:
    """Contenedor del resultado del generador (in-memory)."""
    name:        str
    expr:        pd.DataFrame
    true_times:  pd.Series
    is_rhythmic: pd.Series
    gene_params: pd.DataFrame
    config:      dict


def _sampling_times(
    n_samples: int,
    sampling:  str,
    params:    dict,
    rng:       np.random.Generator,
    span:      float = PERIOD_HOURS,
) -> np.ndarray:
    """
    Genera ``n_samples`` tiempos de muestreo en ``[0, span)`` segun el
    esquema ``sampling`` (ver :data:`SAMPLING_SCHEMES`). Devuelve los
    tiempos ordenados.

    Esquemas
    --------
    ``equispaced``
        ``linspace(0, span, n, endpoint=False)``. Cobertura uniforme.
    ``jitter``
        Equiespaciado + ``U(-f*delta, f*delta)`` por punto, con
        ``delta = span/n`` y ``f = jitter_frac``. Horarios planificados
        con pequenas desviaciones operativas.
    ``clustered``
        Una fraccion ``frac`` de los puntos se sortea de
        ``N(center, sd)`` (concentracion local); el resto se reparte
        ``U(0, span)`` para mantener cobertura de todo el ciclo.
    ``clustered_gap``
        Como ``clustered`` pero los puntos que caen dentro de alguna
        ventana de ``gaps`` se reubican fuera, dejando esas franjas
        practicamente vacias.
    """
    if sampling == "equispaced":
        T = np.linspace(0.0, span, n_samples, endpoint=False)

    elif sampling == "jitter":
        f     = params.get("jitter_frac", 0.4)
        delta = span / n_samples
        base  = np.linspace(0.0, span, n_samples, endpoint=False)
        T     = np.mod(base + rng.uniform(-f * delta, f * delta, size=n_samples), span)

    elif sampling in ("clustered", "clustered_gap"):
        frac   = params.get("frac",   0.4)
        center = params.get("center", 12.0)
        sd     = params.get("sd",     2.0)
        n_clu  = int(round(frac * n_samples))
        n_uni  = n_samples - n_clu
        uni    = rng.uniform(0.0, span, size=n_uni)
        clu    = np.mod(rng.normal(center, sd, size=n_clu), span)
        T      = np.concatenate([uni, clu])

        if sampling == "clustered_gap":
            gaps = params.get("gaps", ((0.0, 6.0),))
            in_gap = lambda x: np.any(
                [(x >= g0) & (x < g1) for g0, g1 in gaps], axis=0
            )
            mask = in_gap(T)
            # Reubica iterativamente los puntos que caigan en un hueco.
            while mask.any():
                T[mask] = rng.uniform(0.0, span, size=int(mask.sum()))
                mask    = in_gap(T)

    else:
        raise ValueError(
            f"sampling debe ser uno de {SAMPLING_SCHEMES}, recibido: {sampling!r}"
        )

    return np.sort(T)


def _fmm_signal(
    t_hours: np.ndarray,
    M:       float,
    A:       float,
    alpha:   float,
    beta:    float,
    omega:   float,
    period:  float = PERIOD_HOURS,
) -> np.ndarray:
    r"""
    Senal del modelo FMM (Rueda et al. 2019)::

        mu(t) = M + A * cos(beta + 2*arctan(omega*tan((t_n - alpha)/2)))

    con :math:`t_n = 2\pi t / period`. Implementado via forma racional
    estable (mismo enfoque que ``mobius_cos_sin`` del FMM estabilizado),
    evitando la singularidad de ``tan()`` cuando :math:`(t_n-\alpha)/2 = \pm\pi/2`.
    """
    t_rad     = 2.0 * np.pi * t_hours / period
    u         = (t_rad - alpha) / 2.0
    cu, su    = np.cos(u), np.sin(u)
    cu2, su2  = cu * cu, su * su
    om2_su2   = (omega * omega) * su2
    den       = cu2 + om2_su2
    den       = np.where(den < 1e-12, 1e-12, den)
    C         = (cu2 - om2_su2) / den
    S         = (2.0 * omega * su * cu) / den
    return M + A * (np.cos(beta) * C - np.sin(beta) * S)


def generate(
    name:            str,
    n_samples:       int,
    *,
    sampling:        str             = "equispaced",
    sampling_params: Optional[dict]  = None,
    omega:           float           = 0.5,
    sigma:           float           = 1.0,
    prop_rhythmic:   float           = 0.05,
    n_total:         int             = 1200,
    amplitude:       float           = 1.0,
    baseline:        float           = 0.0,
    time_span:       float           = PERIOD_HOURS,
    seed:            Optional[int]   = None,
) -> SyntheticDataset:
    """
    Genera un dataset sintetico segun el diseno factorial. Ver docstring
    del modulo para el detalle de los factores.

    Genes ritmicos (FMM, A=``amplitude``, M=``baseline`` fijos)::

        f(t) = M + A*cos(beta + 2*atan(omega*tan((t_n - alpha)/2))) + (A/sigma)*S

    con ``alpha, beta ~ U(0, 2*pi)`` por gen y ``S ~ N(0, 1)``. Genes no
    ritmicos: solo ruido al mismo nivel, ``f(t) = (A/sigma)*S``.

    El numero de genes ritmicos es ``round(prop_rhythmic * n_total)`` y
    debe ser ``>= len(SEED_GENES)`` para que el pipeline encuentre sus
    genes ancla.

    Raises
    ------
    ValueError
        Si ``sampling`` no es valido, si ``n_rhythmic < len(SEED_GENES)``
        o si ``n_total < n_rhythmic``.
    """
    if sampling not in SAMPLING_SCHEMES:
        raise ValueError(
            f"sampling debe ser uno de {SAMPLING_SCHEMES}, recibido: {sampling!r}"
        )

    n_rhythmic = int(round(prop_rhythmic * n_total))
    if n_rhythmic < len(SEED_GENES):
        raise ValueError(
            f"n_rhythmic={n_rhythmic} (= round({prop_rhythmic}*{n_total})) "
            f"debe ser >= {len(SEED_GENES)} (genes semilla CIRCUST). "
            f"Sube n_total o prop_rhythmic."
        )
    if n_total < n_rhythmic:
        raise ValueError(
            f"n_total={n_total} debe ser >= n_rhythmic={n_rhythmic}"
        )

    rng = np.random.default_rng(seed)

    # 1) Tiempos de muestreo segun el esquema elegido
    params = {**SAMPLING_DEFAULTS[sampling], **(sampling_params or {})}
    T = _sampling_times(n_samples, sampling, params, rng, span=time_span)
    n_samples = len(T)

    # 2) Nombres de genes y muestras
    n_other_rhy = n_rhythmic - len(SEED_GENES)
    n_noise     = n_total - n_rhythmic

    other_rhythmic = [f"rhythm_{i:04d}" for i in range(n_other_rhy)]
    noise_genes    = [f"noise_{i:04d}"  for i in range(n_noise)]
    gene_names     = SEED_GENES + other_rhythmic + noise_genes

    is_rhy = np.concatenate([
        np.ones(n_rhythmic, dtype=bool),
        np.zeros(n_noise,   dtype=bool),
    ])

    sample_names = [f"s{i:04d}" for i in range(n_samples)]

    # 3) Pre-allocar matriz y registro de parametros por gen
    expr = np.zeros((n_total, n_samples), dtype=np.float64)
    params_rows: list[dict] = []

    A          = float(amplitude)
    noise_std  = A / sigma

    # 4) Genes ritmicos: senal FMM + ruido (alpha, beta aleatorios por gen)
    for g in range(n_rhythmic):
        alpha  = float(rng.uniform(0.0, 2.0 * np.pi))
        beta   = float(rng.uniform(0.0, 2.0 * np.pi))
        signal = _fmm_signal(T, M=baseline, A=A, alpha=alpha, beta=beta, omega=omega)
        noise  = noise_std * rng.standard_normal(n_samples)
        expr[g] = signal + noise
        params_rows.append(dict(
            gene  = gene_names[g], model = "fmm",
            amp   = A,             snr   = sigma,
            alpha = alpha,         beta  = beta,    omega = omega,
        ))

    # 5) Genes no ritmicos: solo ruido al mismo nivel
    for g in range(n_rhythmic, n_total):
        expr[g] = noise_std * rng.standard_normal(n_samples)
        params_rows.append(dict(
            gene  = gene_names[g], model = "noise",
            amp   = A,             snr   = sigma,
            alpha = np.nan,        beta  = np.nan,  omega = np.nan,
        ))

    # 6) Empaquetar
    expr_df     = pd.DataFrame(expr, index=gene_names, columns=sample_names)
    true_t      = pd.Series(T, index=sample_names, name="true_time_h")
    is_rhy_s    = pd.Series(is_rhy, index=gene_names, name="is_rhythmic")
    gene_params = pd.DataFrame(params_rows).set_index("gene")

    return SyntheticDataset(
        name        = name,
        expr        = expr_df,
        true_times  = true_t,
        is_rhythmic = is_rhy_s,
        gene_params = gene_params,
        config      = dict(
            n_samples       = n_samples,
            time_span       = time_span,
            sampling        = sampling,
            sampling_params = params,
            omega           = omega,
            sigma           = sigma,
            prop_rhythmic   = prop_rhythmic,
            n_rhythmic      = n_rhythmic,
            n_total         = n_total,
            amplitude       = A,
            baseline        = baseline,
            seed            = seed,
        ),
    )


# ===========================================================================
# CLI
# ===========================================================================

def write_dataset(ds: SyntheticDataset, out: Path) -> dict:
    """
    Escribe en ``out`` los 4 artefactos del dataset (``matrix.csv``,
    ``true_times.csv``, ``gene_params.csv``, ``manifest.json``) en el
    formato que consumen ``run_pipeline.py`` y ``validate.py``. Devuelve
    el manifest como dict. Reutilizable por el runner del grid.
    """
    out.mkdir(parents=True, exist_ok=True)

    # 1) Matriz en formato consumible por run_pipeline.py
    matrix_path = out / "matrix.csv"
    ds.expr.index.name = "gene_symbol"
    ds.expr.to_csv(matrix_path)

    # 2) Ground truth temporal
    times_rad = (ds.true_times.values % PERIOD_HOURS) * (2.0 * np.pi / PERIOD_HOURS)
    times_df = pd.DataFrame({
        "sample_name":    ds.true_times.index,
        "true_time_h":    ds.true_times.values,
        "true_phase_h":   ds.true_times.values % PERIOD_HOURS,
        "true_phase_rad": times_rad,
    })
    times_path = out / "true_times.csv"
    times_df.to_csv(times_path, index=False)

    # 3) Parametros por gen
    params_path = out / "gene_params.csv"
    ds.gene_params.to_csv(params_path)

    # 4) Manifest (config completa + period)
    manifest = {
        "name":     ds.name,
        "config":   ds.config,
        "period_h": PERIOD_HOURS,
        "outputs": {
            "matrix":      str(matrix_path.name),
            "true_times":  str(times_path.name),
            "gene_params": str(params_path.name),
        },
    }
    with open(out / "manifest.json", "w") as f:
        json.dump(manifest, f, indent=2, default=str)

    return manifest


def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Genera UNA celda del estudio factorial (un escenario) y la "
            "guarda en disco como matrix.csv + true_times.csv (consumibles "
            "por run_pipeline.py y validate.py). Para el grid completo usar "
            "run_simulation_grid.py."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--name",      default="cell",
                   help="Nombre del escenario. Defecto: 'cell'.")
    p.add_argument("--n-samples", type=int,   default=48,
                   help="Numero de muestras. Defecto: 48.")
    p.add_argument("--sampling",  default="equispaced",
                   choices=list(SAMPLING_SCHEMES),
                   help="Esquema de muestreo temporal. Defecto: equispaced.")
    p.add_argument("--omega",     type=float, default=0.5,
                   help="Asimetria FMM. omega=1 => sinusoide. Defecto: 0.5.")
    p.add_argument("--sigma",     type=float, default=1.0,
                   help="Intensidad de senal (ruido=A/sigma). Defecto: 1.0.")
    p.add_argument("--prop-rhythmic", type=float, default=0.05,
                   help="Proporcion de genes ritmicos. Defecto: 0.05.")
    p.add_argument("--n-total",   type=int,   default=1200,
                   help="Genes totales (ritmicos + ruido). Defecto: 1200.")
    p.add_argument("--seed",      type=int,   default=42,
                   help="Semilla del RNG. Defecto: 42.")
    p.add_argument("-o", "--output", type=str, required=True,
                   help="Directorio de salida (se crea si no existe).")
    return p.parse_args()


def main() -> None:
    args = _parse_args()
    out  = Path(args.output)

    n_rhythmic = int(round(args.prop_rhythmic * args.n_total))

    print(f"╔══════════════════════════════════════════════════════════════════╗")
    print(f"  Escenario  : {args.name}")
    print(f"  n_samples  : {args.n_samples}")
    print(f"  sampling   : {args.sampling}")
    print(f"  omega      : {args.omega}")
    print(f"  sigma      : {args.sigma}  (ruido = A/sigma = {1.0/args.sigma:.3f})")
    print(f"  prop rhy.  : {args.prop_rhythmic}  -> {n_rhythmic} ritmicos / {args.n_total} totales")
    print(f"  seed       : {args.seed}")
    print(f"  output     : {out.resolve()}")
    print(f"╚══════════════════════════════════════════════════════════════════╝")

    t_start = time.time()

    ds = generate(
        name          = args.name,
        n_samples     = args.n_samples,
        sampling      = args.sampling,
        omega         = args.omega,
        sigma         = args.sigma,
        prop_rhythmic = args.prop_rhythmic,
        n_total       = args.n_total,
        seed          = args.seed,
    )
    write_dataset(ds, out)

    elapsed = time.time() - t_start
    n_fmm   = int((ds.gene_params["model"] == "fmm").sum())

    print()
    print(f"  ✓ Generado en {elapsed:.1f}s")
    print(f"    matrix      : {ds.expr.shape[0]} genes x {ds.expr.shape[1]} muestras "
          f"({n_fmm} fmm + ruido)")
    print(f"    Archivos en {out.resolve()}:")
    print(f"      - matrix.csv")
    print(f"      - true_times.csv")
    print(f"      - gene_params.csv")
    print(f"      - manifest.json")


if __name__ == "__main__":
    main()
