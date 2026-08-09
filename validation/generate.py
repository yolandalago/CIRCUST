#!/usr/bin/env python3
"""
validation/generate.py
======================
Generador de conjuntos sinteticos de expresion circadiana para el estudio de
simulacion.

Idea general
------------
Se sigue la estructura del generador de DCPR, que fija la biologia y mueve el
diseno de muestreo, pero la senal ritmica se genera con el modelo FMM en lugar
de con un coseno. El motivo es que un coseno equivale a fijar el parametro de
apuntamiento en ``omega = 1``, de modo que no habria asimetria de forma de onda
que capturar y no se podria evaluar la estimacion estabilizada del modelo.

Que produce
-----------
Una matriz de expresion ``(n_genes, n_muestras)`` junto con los tiempos reales
de cada muestra y los parametros con los que se genero cada gen. Los tiempos
reales son el ground truth del orden, y los parametros por gen son el ground
truth de la caracterizacion circadiana.

Composicion de la matriz
------------------------
1. Los genes \\textit{core} del reloj, con fases relativas tomadas de un atlas
   real (ver ``CORE_PHASE_H``). Sin ellos CIRCUST no arranca, porque el CPCA
   se construye sobre su submatriz.
2. El resto de genes ritmicos, con fase uniforme sobre el ciclo.
3. Genes no ritmicos, ruido puro.

Ruido
-----
De momento gaussiano e independiente por gen, relativo a la amplitud, como en
DCPR. Pendiente de sustituir por un ruido con estructura, con efecto compartido
por muestra y efectos de lote, que es lo que hace realista un conjunto tipo
GTEx y lo que el resumen de la memoria promete.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd

PERIOD_HOURS: float = 24.0

SAMPLING_SCHEMES = ("equiespaciado", "huecos", "concentrado")


# ===========================================================================
# Conocimiento previo: fases de los genes core
# ===========================================================================
#
# Horas transcurridas desde el pico de ARNTL, que se toma como origen. Se
# eligen fases relativas y no absolutas porque lo que la etapa de
# sincronizacion utiliza son las relaciones de orden entre genes, no un
# instante absoluto.
#
# Los siete primeros proceden del material suplementario de Ruben et al.
# (2018), columna ``mean_Phase``, que promedia la fase estimada sobre los
# tejidos en los que cada gen cicla. Es el fichero
# ``data/NIHMS1725947-supplement-File_S2 54 cores.xlsx``.
#
# Los cinco restantes no figuran en ese atlas y se sitan a partir de las
# relaciones cualitativas documentadas en la introduccion de la memoria: CLOCK
# y RORA comparten la fase ROR de ARNTL, NR1D1 pica temprano en el periodo
# activo y en oposicion parcial a ARNTL, CRY1 pica despues de los PER, y STAT3
# tiene ritmicidad debil. Estan marcados como aproximados.
#
# Convenio de orientacion
# -----------------------
# Las horas estan medidas en sentido **retrogrado** respecto al eje del atlas
# de Ruben, es decir, con el signo cambiado. El motivo es que ambos ejes
# describen la misma biologia pero recorren el ciclo en sentidos opuestos.
#
# En el eje de Ruben, DBP alcanza su maximo unas once horas despues de ARNTL.
# La etapa de sincronizacion de CIRCUST rota la escala para situar el pico de
# ARNTL en pi y despues espera encontrar el de DBP en el semicirculo [0, pi),
# esto es, antes de ARNTL. Con el eje de Ruben sin cambiar, DBP caeria en
# 6,12 radianes, dentro de [pi, 2pi), de modo que CIRCUST invertiria el sentido
# del recorrido en todos los escenarios. Peor aun, ese valor dista solo 0,16
# radianes del origen, con lo que se disparara la condicion de aviso
# 1 - cos(t_U) <= 0,1 y la orientacion quedaria decidida por el gen de
# consistencia en lugar de por el de direccion.
#
# Al cambiar el signo, DBP pica unas once horas **antes** que ARNTL y cae en
# [0, pi) tras la rotacion, de modo que la regla de direccion resuelve en el
# sentido correcto. La biologia generada es la misma, solo cambia el sentido en
# que se etiqueta el paso del tiempo.
#
# Observacion sobre la condicion de aviso (no explotada en este trabajo)
# ---------------------------------------------------------------------
# Invertir el eje corrige el sentido pero no evita la condicion de aviso, y no
# puede evitarla ninguna eleccion de convenio. La razon es que DBP y ARNTL
# pican con una separacion de unas 11,4 horas, casi media vuelta, de manera que
# al situar ARNTL en pi el pico de DBP aterriza junto al origen cualquiera que
# sea el sentido de recorrido. Para que la regla del gen de direccion decidiese
# por si sola, esa separacion tendria que ser inferior a unas 10,3 horas.
#
# En consecuencia, con relaciones de fase realistas la condicion de aviso se
# cumple siempre y la orientacion queda decidida por el gen de consistencia,
# CRY1, y no por el de direccion. Eso concuerda con lo que el propio material
# suplementario de CIRCUST deja entrever al declarar ARNTL, DBP y CRY1 como las
# tres claves de sincronizacion y excluirlos de su analisis de robustez frente
# a cambios en la lista de genes semilla, aunque el articulo no describe ni la
# condicion de aviso ni el papel de CRY1, que solo aparecen en el codigo.
#
# Se deja constancia aqui por si en algun momento interesa desarrollarlo. No
# forma parte de los objetivos de este trabajo.
CORE_PHASE_H: dict[str, float] = {
    # --- de Ruben et al. (2018), con el eje invertido ---
    "ARNTL":   0.00,
    "TEF":    -10.68 % 24.0,   # 13.32
    "DBP":    -11.40 % 24.0,   # 12.60
    "PER3":   -11.76 % 24.0,   # 12.24
    "CRY2":   -12.86 % 24.0,   # 11.14
    "PER2":   -13.65 % 24.0,   # 10.35
    "PER1":   -14.19 % 24.0,   #  9.81
    # --- aproximados a partir de relaciones documentadas, mismo eje ---
    "CLOCK":    0.00,
    "RORA":    -1.00 % 24.0,   # 23.00
    "NR1D1":   -6.00 % 24.0,   # 18.00
    "CRY1":   -16.00 % 24.0,   #  8.00
    "STAT3":  -12.00 % 24.0,   # 12.00
}

CORE_PHASE_SOURCE: dict[str, str] = {
    g: ("ruben2018" if i < 7 else "aproximado")
    for i, g in enumerate(CORE_PHASE_H)
}

# Amplitud relativa media de Ruben et al. (2018), columna ``rAmp_mean``. Se usa
# para que los genes core no tengan todos la misma intensidad de oscilacion.
CORE_RAMP: dict[str, float] = {
    "ARNTL": 0.615, "TEF": 0.577, "DBP": 0.809, "PER3": 0.715,
    "CRY2": 0.376, "PER2": 0.556, "PER1": 0.489,
    "CLOCK": 0.50, "RORA": 0.50, "NR1D1": 0.60, "CRY1": 0.50,
    "STAT3": 0.30,     # ritmicidad menos marcada, ver la memoria
}


# ===========================================================================
# Modelo FMM
# ===========================================================================

def fmm_signal(
    t_h:    np.ndarray,
    M:      float,
    A:      float,
    alpha:  float,
    beta:   float,
    omega:  float,
    period: float = PERIOD_HOURS,
) -> np.ndarray:
    r"""Senal esperada del modelo FMM monocomponente.

    .. math::
        \mu(t) = M + A\cos\{\beta + \phi(t;\alpha,\omega)\},
        \qquad
        \phi(t;\alpha,\omega) = 2\arctan\{\omega\tan((t-\alpha)/2)\}

    Se evalua por la forma racional equivalente, sin llamar a ``tan``, de modo
    que desaparece la singularidad en ``t = alpha + pi``. Con
    ``u = (t-alpha)/2`` y ``den = cos^2 u + omega^2 sin^2 u``, que es una suma
    de cuadrados y por tanto positiva:

        cos(phi) = (cos^2 u - omega^2 sin^2 u) / den
        sin(phi) = 2 omega sin u cos u / den

    y se compone con ``cos(beta + phi) = cos b cos phi - sin b sin phi``.
    """
    t_rad = 2.0 * np.pi * np.asarray(t_h, dtype=float) / period
    u = 0.5 * (t_rad - alpha)
    cu, su = np.cos(u), np.sin(u)
    cu2, su2 = cu * cu, su * su
    den = cu2 + (omega * omega) * su2
    den = np.where(den < 1e-300, 1e-300, den)
    cos_phi = (cu2 - (omega * omega) * su2) / den
    sin_phi = (2.0 * omega * su * cu) / den
    return M + A * (np.cos(beta) * cos_phi - np.sin(beta) * sin_phi)


def peak_time_rad(alpha: float, beta: float, omega: float) -> float:
    """Instante de maxima expresion ``t_U`` del modelo FMM, en radianes.

    Es la formula cerrada del modelo, la misma que emplea ``FMMModel.peak_time``
    en el paquete ``circust``.
    """
    return float(
        (alpha + 2.0 * np.arctan2(
            (1.0 / omega) * np.sin(-beta / 2.0),
            np.cos(-beta / 2.0),
        )) % (2.0 * np.pi)
    )


def alpha_for_peak(t_u_rad: float, beta: float, omega: float) -> float:
    """Inversa de :func:`peak_time_rad`: el ``alpha`` que coloca el pico en ``t_u_rad``.

    Es lo que permite generar un gen cuyo maximo de expresion cae exactamente
    donde se quiere, que es como se sitan los genes \\textit{core} en sus fases
    conocidas. Despejar ``alpha`` es inmediato porque solo aparece como sumando.
    """
    return float(
        (t_u_rad - 2.0 * np.arctan2(
            (1.0 / omega) * np.sin(-beta / 2.0),
            np.cos(-beta / 2.0),
        )) % (2.0 * np.pi)
    )


# ===========================================================================
# Esquemas de muestreo
# ===========================================================================

def sampling_times(
    n_samples: int,
    scheme:    str,
    rng:       np.random.Generator,
    period:    float = PERIOD_HOURS,
    gap:       tuple[float, float] = (2.0, 8.0),
    cluster_center: float = 12.0,
    cluster_sd:     float = 1.5,
    cluster_frac:   float = 0.6,
) -> np.ndarray:
    """Instantes de muestreo en horas, ordenados, dentro de ``[0, period)``.

    Los tres esquemas reproducen los tres regimenes que motivan la estimacion
    estabilizada del modelo FMM.

    ``equiespaciado``
        Cobertura uniforme del ciclo. Es el caso favorable y sirve de
        referencia.

    ``huecos``
        Cobertura uniforme salvo una ventana sin observar. Reproduce el diseno
        con lagunas, donde la transicion de la fase puede caer integramente
        dentro del intervalo no observado.

    ``concentrado``
        Una fraccion de las muestras se concentra en torno a un instante y el
        resto se reparte. Reproduce el diseno con instantes casi duplicados,
        donde los incrementos temporales muy pequenos desestabilizan las
        pendientes.
    """
    if scheme == "equiespaciado":
        t = np.linspace(0.0, period, n_samples, endpoint=False)

    elif scheme == "huecos":
        lo, hi = gap
        if not 0.0 <= lo < hi <= period:
            raise ValueError(f"ventana {gap} invalida para un periodo {period}.")
        # Se muestrea directamente sobre el complemento de la ventana, de modo
        # que el procedimiento termina siempre.
        libre = period - (hi - lo)
        u = np.linspace(0.0, libre, n_samples, endpoint=False)
        t = np.where(u < lo, u, u + (hi - lo))

    elif scheme == "concentrado":
        n_clu = int(round(cluster_frac * n_samples))
        clu = rng.normal(cluster_center, cluster_sd, size=n_clu) % period
        uni = np.linspace(0.0, period, n_samples - n_clu, endpoint=False)
        t = np.concatenate([clu, uni])

    else:
        raise ValueError(
            f"scheme debe ser uno de {SAMPLING_SCHEMES}, recibido {scheme!r}."
        )

    return np.sort(np.mod(t, period))


# ===========================================================================
# Generador
# ===========================================================================

@dataclass
class SyntheticDataset:
    """Un escenario sintetico completo."""
    name:         str
    matrix:       pd.DataFrame          # genes x muestras
    true_times_h: np.ndarray            # (n_muestras,) en [0, 24)
    gene_params:  pd.DataFrame          # una fila por gen
    meta:         dict = field(default_factory=dict)


def generate(
    name:            str = "escenario",
    *,
    n_samples:       int   = 60,
    n_genes:         int   = 2000,
    prop_rhythmic:   float = 0.05,
    sampling:        str   = "equiespaciado",
    omega:           float = 1.00,
    omega_jitter:    float = 0.05,
    beta_mode:       str   = "aleatorio",
    snr:             float = 3.0,
    core_phase_sd_h: float = 1.0,
    seed:            int   = 0,
    period:          float = PERIOD_HOURS,
    **sampling_kwargs,
) -> SyntheticDataset:
    """Genera un escenario sintetico.

    Parameters
    ----------
    n_samples, n_genes
        Dimensiones de la matriz.
    prop_rhythmic
        Fraccion de genes con senal circadiana. Por defecto 0,05, que es la
        proporcion que indica el guion del capitulo y el orden de magnitud
        habitual en tejido real.
    sampling
        Uno de :data:`SAMPLING_SCHEMES`.
    omega
        Apuntamiento de la forma de onda. ``1`` es la sinusoide simetrica y
        valores menores producen perfiles progresivamente mas apuntados y
        asimetricos. Es el factor principal del estudio.
    omega_jitter
        Dispersion de ``omega`` entre genes, para que no todos compartan
        exactamente la misma forma.
    beta_mode
        ``"aleatorio"`` sortea ``beta`` uniformemente, de modo que la asimetria
        apunte tanto al maximo como al minimo. ``"cero"`` lo fija en cero, que
        produce un unico tipo de asimetria y sirve como caso controlado.
    snr
        Razon senal-ruido. La desviacion tipica del ruido de un gen es su
        amplitud dividida por este valor, como en DCPR, de modo que todos los
        genes ritmicos tienen la misma dificultad con independencia de su
        amplitud.
    core_phase_sd_h
        Dispersion, en horas, con que se sortea la fase de cada gen
        \\textit{core} alrededor de su valor de referencia. Con cero las fases
        serian exactamente las de :data:`CORE_PHASE_H`, lo que haria trivial la
        etapa de sincronizacion y convertiria su validacion en una tautologia.
    seed
        Semilla del generador de numeros aleatorios.
    """
    if not 0.0 < prop_rhythmic <= 1.0:
        raise ValueError(f"prop_rhythmic debe estar en (0, 1], recibido {prop_rhythmic}.")
    if not 0.0 < omega <= 1.0:
        raise ValueError(f"omega debe estar en (0, 1], recibido {omega}.")

    rng = np.random.default_rng(seed)
    two_pi = 2.0 * np.pi

    n_rhythmic = int(round(prop_rhythmic * n_genes))
    core = list(CORE_PHASE_H)
    if n_rhythmic < len(core):
        raise ValueError(
            f"prop_rhythmic={prop_rhythmic} sobre {n_genes} genes da "
            f"{n_rhythmic} ritmicos, menos que los {len(core)} genes core que "
            f"CIRCUST necesita para arrancar el CPCA."
        )

    t_h = sampling_times(n_samples, sampling, rng, period, **sampling_kwargs)
    n = t_h.size

    # ── Parametros por gen ────────────────────────────────────────────────
    n_extra = n_rhythmic - len(core)
    nombres = core + [f"RIT{i:05d}" for i in range(n_extra)] \
                   + [f"BG{i:05d}" for i in range(n_genes - n_rhythmic)]

    t_u = np.full(n_genes, np.nan)      # pico objetivo, en radianes
    amp = np.full(n_genes, np.nan)
    om  = np.full(n_genes, np.nan)
    bet = np.full(n_genes, np.nan)
    alp = np.full(n_genes, np.nan)
    fuente = ["-"] * n_genes

    # Genes core: fase de referencia mas dispersion, y amplitud del atlas
    for i, g in enumerate(core):
        fase_h = rng.normal(CORE_PHASE_H[g], core_phase_sd_h)
        t_u[i] = (two_pi * fase_h / period) % two_pi
        amp[i] = CORE_RAMP[g]
        fuente[i] = CORE_PHASE_SOURCE[g]

    # Resto de ritmicos: pico uniforme sobre el ciclo
    sl = slice(len(core), n_rhythmic)
    t_u[sl] = rng.uniform(0.0, two_pi, size=n_extra)
    amp[sl] = rng.uniform(0.5, 1.0, size=n_extra)
    for i in range(len(core), n_rhythmic):
        fuente[i] = "uniforme"

    # Forma de onda comun a todos los ritmicos, con dispersion
    r = slice(0, n_rhythmic)
    om[r] = np.clip(rng.normal(omega, omega_jitter, size=n_rhythmic), 0.05, 1.0)
    bet[r] = (rng.uniform(0.0, two_pi, size=n_rhythmic) if beta_mode == "aleatorio"
              else np.zeros(n_rhythmic))
    # alpha se despeja para que el pico caiga en t_u
    for i in range(n_rhythmic):
        alp[i] = alpha_for_peak(t_u[i], bet[i], om[i])

    # Genes de fondo: sin senal, solo amplitud de escala del ruido
    bg = slice(n_rhythmic, n_genes)
    amp[bg] = rng.uniform(0.5, 1.0, size=n_genes - n_rhythmic)

    # ── Matriz ────────────────────────────────────────────────────────────
    valores = np.empty((n_genes, n), dtype=float)
    ruido = rng.standard_normal((n_genes, n))

    for i in range(n_rhythmic):
        valores[i] = fmm_signal(t_h, 0.0, amp[i], alp[i], bet[i], om[i], period)
    valores[bg] = 0.0
    valores += (amp[:, None] / snr) * ruido

    matriz = pd.DataFrame(
        valores,
        index=nombres,
        columns=[f"S{j:03d}" for j in range(n)],
    )
    params = pd.DataFrame({
        "gene":       nombres,
        "ritmico":    [True] * n_rhythmic + [False] * (n_genes - n_rhythmic),
        "fuente_fase": fuente,
        "amplitud":   amp,
        "alpha":      alp,
        "beta":       bet,
        "omega":      om,
        "t_U_rad":    t_u,
        "t_U_h":      t_u * period / two_pi,
    })

    return SyntheticDataset(
        name=name,
        matrix=matriz,
        true_times_h=t_h,
        gene_params=params,
        meta=dict(
            seed=seed, n_samples=n, n_genes=n_genes,
            prop_rhythmic=prop_rhythmic, sampling=sampling,
            omega=omega, omega_jitter=omega_jitter, beta_mode=beta_mode,
            snr=snr, core_phase_sd_h=core_phase_sd_h, period=period,
            ruido="gaussiano independiente, relativo a la amplitud",
        ),
    )
