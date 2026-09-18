"""
Calidad del atlas: cuanto se parece la onda publicada por gen a la verdadera.

Estas son las medidas del bloque de regularizacion. El FMM estabilizado no
puede cambiar el orden de muestras, porque ese orden sale de una SVD sobre
los genes nucleares donde el FMM no interviene en ningun paso. Lo que si
cambia es la estimacion de los parametros gen a gen, y eso es exactamente lo
que el atlas publica. Medir el orden para comparar base contra estabilizado
seria medir donde los dos metodos son identicos por construccion.

Que se guarda y por que
-----------------------
La medida principal es el RMSE de la curva. Es la unica libre de
parametrizacion: ``alpha``, ``beta`` y ``omega`` interactuan, de modo que dos
ternas muy distintas pueden describir casi la misma onda y compararlas una a
una penalizaria diferencias que no existen. La curva no tiene ese problema
porque es el objeto que el atlas afirma haber estimado.

El resto acompana porque cada una retiene entre el 75 y el 84 % de varianza
propia frente al RMSE de la curva (medido sobre 2160 ajustes), asi que
ninguna es redundante:

    pico y valle   en horas, interpretables directamente
    subida         |error del pico - error del valle|, la forma de la onda
    omega          el mecanismo por el que el estabilizado mejora
    R2             lo que el atlas declara de si mismo, sin usar la verdad
    VPP            de los genes del atlas, cuantos son ritmicos de verdad

Todos los errores se guardan CON SIGNO y acompanados de su dispersion, para
que el analisis derive sesgo y varianza sin volver a ejecutar nada.

Por que no se compara la amplitud
---------------------------------
El preprocesado normaliza min-max cada gen al intervalo [-1, 1]. Despues de
eso ``M`` y ``A`` absolutos ya no son recuperables ni en principio: la onda
llega al FMM reescalada por un factor que depende del rango observado, que a
su vez depende del ruido. Lo que ``amplitude`` estima es la fraccion del
rango observado que explica el ritmo, o sea un indicador de senal a ruido, no
el parametro ``A`` del generador. Compararlos daria un error que crece con el
ruido aunque la onda se haya recuperado perfectamente. Por eso las curvas se
comparan estandarizadas, ``(f - M)/A``, que es justo la parte identificable,
y el papel de indicador de calidad lo cumple el ``R2`` que el propio atlas
publica.

El giro global
--------------
Los tiempos de pico del atlas viven en la escala circular estimada, que puede
estar girada respecto de la verdadera. Ese giro es un error de anclaje, ya
medido por ``error_fase_medio_h`` sobre las muestras, y si se dejara dentro dominaria el
error de pico y taparia justo lo que se quiere ver. Con ``quitar_giro`` se
descuenta el giro comun a todos los genes y se registra aparte.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from circust.fitting.fmm import FMMModel
from simulation.generate import fmm_wave
from simulation.truth import TWO_PI

# Malla sobre la que se compara la onda estimada contra la verdadera. 360
# puntos son un grado de resolucion, muy por debajo de cualquier error que
# importe, y el coste es despreciable frente a un ajuste FMM.
_MALLA = np.linspace(0.0, TWO_PI, 360, endpoint=False)

_CLAVES = (
    "atlas_error_onda",
    "atlas_err_pico_h", "atlas_variabilidad_pico_h",
    "atlas_err_valle_h", "atlas_err_valle_sd_h",
    "atlas_err_subida_h",
    "atlas_sesgo_omega", "atlas_omega_sd", "atlas_omega_rmse",
    "atlas_r2_mediana",
    "atlas_giro_h",
)


def _a_horas_con_signo(rad) -> np.ndarray:
    """Diferencia circular con signo, envuelta a ``(-12, 12]`` horas."""
    d = (np.asarray(rad, dtype=float) + np.pi) % TWO_PI - np.pi
    return d * 24.0 / TWO_PI


def _giro_comun(err_rad: np.ndarray) -> float:
    """Media circular de los errores de pico, el giro que comparten todos."""
    if err_rad.size == 0:
        return 0.0
    return float(np.arctan2(np.sin(err_rad).mean(), np.cos(err_rad).mean()))


def _ritmicos(ajuste: pd.DataFrame, verdad: pd.DataFrame) -> np.ndarray:
    """Mascara de los genes del atlas que tienen onda verdadera que comparar."""
    nombres = ajuste["gene_name"].to_numpy()
    faltan = set(nombres) - set(verdad.index)
    if faltan:
        raise KeyError(
            f"{len(faltan)} genes del atlas no estan en la verdad, p.ej. "
            f"{sorted(faltan)[:3]}"
        )
    return verdad.loc[nombres, "kind"].to_numpy() != "fondo"


def per_gene_errors(
    ajuste: pd.DataFrame,
    verdad: pd.DataFrame,
    quitar_giro: bool = False,
) -> pd.DataFrame:
    """
    Una fila por gen evaluable, con su error frente a la onda verdadera.

    Es la tabla de la que ``atlas_metrics`` saca sus resumenes, y la que
    alimenta las figuras por gen del capitulo: sin ella habria que repetir
    este calculo en el codigo de las figuras y las dos copias podrian
    separarse en silencio.

    Los genes de fondo quedan fuera, porque no tienen onda con la que
    comparar. Ver el encabezado del modulo para el resto de decisiones.
    """
    es_ritmico = _ritmicos(ajuste, verdad)
    a = ajuste.loc[es_ritmico].reset_index(drop=True)
    v = verdad.loc[a["gene_name"].to_numpy()]

    # --- pico y valle, con el giro global fuera si se pide ---------------
    pico_hat = np.array([FMMModel.peak_time(al, be, om) for al, be, om
                         in zip(a["alpha"], a["beta"], a["omega"])])
    pico_ref = np.array([FMMModel.peak_time(al, be, om) for al, be, om
                         in zip(v["alpha"], v["beta"], v["omega"])])
    giro = _giro_comun(pico_hat - pico_ref) if quitar_giro else 0.0

    alpha_hat = a["alpha"].to_numpy() - giro
    pico_hat -= giro
    valle_hat = np.array([FMMModel.trough_time(al, be, om) for al, be, om
                          in zip(alpha_hat, a["beta"], a["omega"])])
    valle_ref = np.array([FMMModel.trough_time(al, be, om) for al, be, om
                          in zip(v["alpha"], v["beta"], v["omega"])])

    e_pico = _a_horas_con_signo(pico_hat - pico_ref)
    e_valle = _a_horas_con_signo(valle_hat - valle_ref)

    # --- RMSE de la curva, sobre las ondas estandarizadas -----------------
    # (f - M)/A en las dos, porque el min-max del preprocesado borra M y A.
    # Ver el encabezado del modulo.
    rmse = np.empty(len(v))
    for i in range(len(v)):
        fa = fmm_wave(_MALLA, 0.0, 1.0, alpha_hat[i],
                      a["beta"].iloc[i], a["omega"].iloc[i])
        fv = fmm_wave(_MALLA, 0.0, 1.0, v["alpha"].iloc[i],
                      v["beta"].iloc[i], v["omega"].iloc[i])
        rmse[i] = np.sqrt(np.mean((fa - fv) ** 2))

    e_omega = a["omega"].to_numpy() - v["omega"].to_numpy()
    r2 = (a["r2_fmm"].to_numpy() if "r2_fmm" in a.columns
          else np.full(len(a), np.nan))

    return pd.DataFrame({
        "gene_name":   a["gene_name"].to_numpy(),
        "kind":        v["kind"].to_numpy(),
        "omega_true":  v["omega"].to_numpy(),
        "omega_hat":   a["omega"].to_numpy(),
        "alpha_true":  v["alpha"].to_numpy(),
        "alpha_hat":   alpha_hat % TWO_PI,
        "beta_true":   v["beta"].to_numpy(),
        "beta_hat":    a["beta"].to_numpy(),
        "pico_true_h": pico_ref * 24.0 / TWO_PI,
        "pico_hat_h":  (pico_hat % TWO_PI) * 24.0 / TWO_PI,
        "err_pico_h":  e_pico,
        "err_valle_h": e_valle,
        "err_omega":   e_omega,
        "rmse_curva":  rmse,
        "r2_fmm":      r2,
        "giro_rad":    giro,
    })


def atlas_metrics(
    ajuste: pd.DataFrame,
    verdad: pd.DataFrame,
    quitar_giro: bool = False,
) -> dict[str, float]:
    """
    Compara las ondas del atlas con las verdaderas, gen a gen.

    Parameters
    ----------
    ajuste
        Lo que publica el atlas, una fila por gen seleccionado, con columnas
        ``gene_name, M, amplitude, alpha, beta, omega``.
    verdad
        ``truth.gene_params`` indexado por gen, con ``kind, M, A, alpha,
        beta, omega``. Los genes de fondo no tienen onda verdadera, asi que
        cuentan para el VPP pero no para los errores.
    quitar_giro
        Descuenta el giro comun a todos los genes antes de medir pico, valle
        y curva. Ver el encabezado del modulo.
    """
    n_atlas = int(len(ajuste))
    salida: dict[str, float] = {k: float("nan") for k in _CLAVES}
    salida["atlas_n_atlas"] = float(n_atlas)

    if n_atlas == 0:
        salida["atlas_vpp"] = float("nan")
        salida["atlas_n_evaluados"] = 0.0
        return salida

    es_ritmico = _ritmicos(ajuste, verdad)
    salida["atlas_vpp"] = float(es_ritmico.mean())
    salida["atlas_n_evaluados"] = float(es_ritmico.sum())
    if not es_ritmico.any():
        return salida

    g = per_gene_errors(ajuste, verdad, quitar_giro)
    e_pico, e_valle = g["err_pico_h"].to_numpy(), g["err_valle_h"].to_numpy()
    e_omega, rmse = g["err_omega"].to_numpy(), g["rmse_curva"].to_numpy()
    r2, giro = g["r2_fmm"].to_numpy(), float(g["giro_rad"].iloc[0])

    salida.update({
        # mediana: el RMSE por gen ya es una media cuadratica, y un solo gen
        # mal ajustado no debe decidir la celda
        "atlas_error_onda":   float(np.median(rmse)),
        # media: un sesgo es una media por definicion, y la sd va aparte
        "atlas_err_pico_h":   float(e_pico.mean()),
        "atlas_variabilidad_pico_h": float(e_pico.std(ddof=0)),
        "atlas_err_valle_h":  float(e_valle.mean()),
        "atlas_err_valle_sd_h": float(e_valle.std(ddof=0)),
        "atlas_err_subida_h": float(np.abs(e_pico - e_valle).mean()),
        "atlas_sesgo_omega":  float(e_omega.mean()),
        "atlas_omega_sd":     float(e_omega.std(ddof=0)),
        "atlas_omega_rmse":   float(np.sqrt(np.mean(e_omega ** 2))),
        "atlas_r2_mediana":   float(np.median(r2)),
        "atlas_giro_h":       float(_a_horas_con_signo(giro)),
    })
    return salida


COLUMNAS_ATLAS = tuple(_CLAVES) + ("atlas_vpp", "atlas_n_atlas",
                                   "atlas_n_evaluados")
