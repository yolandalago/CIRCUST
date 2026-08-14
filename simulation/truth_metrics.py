#!/usr/bin/env python3
"""
simulation/truth_metrics.py
===========================
Metricas que solo pueden calcularse sobre datos sinteticos, porque comparan la
salida del algoritmo contra los parametros con los que se genero cada gen.

Que anaden respecto a lo que ya se mide
---------------------------------------
``validation/compare_methods.py`` mide el error de las fases de pico contra una
**referencia estimada**: ajusta el FMM sobre los tiempos reales de recogida y
compara contra eso. Es lo mejor que puede hacerse con datos reales, pero
confunde dos cosas, el error de reconstruir el orden y el error del propio
ajuste FMM sobre tiempos perfectos.

Aqui la comparacion es contra los ``alpha``, ``beta`` y ``omega`` verdaderos,
sin ajuste intermedio, de modo que las dos fuentes de error quedan separadas.
Ninguna validacion previa de CIRCUST hace esto: su estudio de simulacion
compara picos estimados contra picos estimados y no mide ni omega, ni la
calidad del conjunto TOP, ni la sincronizacion.

Las tres familias
-----------------
1. **Fases de pico y sincronizacion.** Error de ``t_U`` frente al verdadero, y
   ademas el giro y la reflexion que hicieron falta para alinear, que son
   precisamente el error de la etapa de sincronizacion.
2. **Omega.** Error y sesgo del apuntamiento estimado. Es la metrica que juzga
   directamente la regularizacion del FMM, y el sesgo con signo importa mas
   que su valor absoluto.
3. **Conjunto TOP.** Que fraccion de los genes seleccionados son de verdad
   ritmicos, y como de bien cubren el ciclo.

Convenio de nombres
-------------------
Todas las claves llevan el prefijo ``verdad_`` para que puedan convivir en la
misma fila con las que ``compare_methods`` calcula contra la referencia
estimada, sin ambiguedad sobre cual es cual.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from validation.metrics import align_phases_by_mad, circular_abs_diff

PERIOD_HOURS: float = 24.0
DOS_PI: float = 2.0 * np.pi


def _centrar(delta_h: float, period: float = PERIOD_HOURS) -> float:
    """Lleva un desplazamiento circular al intervalo ``(-period/2, period/2]``.

    Un giro de 20 horas es el mismo que uno de -4. Reportar 20 sugeriria un
    fallo enorme de sincronizacion donde hay uno pequeno.
    """
    d = float(delta_h) % period
    return d - period if d > period / 2.0 else d


def _comunes(atlas: pd.DataFrame, params: pd.DataFrame) -> pd.DataFrame:
    """Une atlas y verdad por gen, conservando solo los genes de ambos.

    Los genes de fondo no estan en ``params`` a proposito (ver
    ``simulation/cells.py``), de modo que esta union descarta por si sola los
    genes sin senal que el conjunto TOP haya podido colar.
    """
    izq = atlas.rename(columns={"gene_name": "gene"})
    return izq.merge(params, on="gene", how="inner", suffixes=("_est", "_ver"))


# ===========================================================================
# 1. Fases de pico y sincronizacion
# ===========================================================================

def _error_de_sincronizacion(delta: float, sign: int, params: pd.DataFrame,
                             anchor_gene: str, period: float) -> float:
    """Cuanto se desvia el origen del ciclo respecto de donde deberia estar.

    El desplazamiento crudo de la alineacion **no** es el error de
    sincronizacion, porque incluye un desplazamiento debido al convenio de
    anclaje. CIRCUST situa el pico del gen ancla en ``pi``, o sea en la hora
    ``period / 2`` de la escala reconstruida, mientras que en el generador ese
    gen pica donde diga el atlas de referencia. Con ARNTL, que pica en el
    origen, el desfase esperado es de doce horas y un desplazamiento crudo
    cercano a cero significaria justo lo contrario de lo que parece.

    De ahi que el error se mida frente al desplazamiento esperado. Sea ``t_a``
    la fase verdadera del gen ancla. La alineacion aplica
    ``alineado = sign * estimado + delta``, y el estimado del ancla vale
    ``period / 2`` por construccion, luego el delta correcto es
    ``t_a - sign * period / 2``.
    """
    fila = params.loc[params["gene"] == anchor_gene, "t_U_h"]
    if fila.empty:
        return float("nan")
    esperado = float(fila.iloc[0]) - sign * period / 2.0
    return float(circular_abs_diff(np.array([delta]),
                                   np.array([esperado]), period)[0])


def error_picos(atlas: pd.DataFrame, params: pd.DataFrame,
                period: float = PERIOD_HOURS,
                anchor_gene: str = "ARNTL") -> dict:
    """Error de las fases de pico frente a las verdaderas, en horas.

    Las fases estimadas viven sobre la escala circular reconstruida, que solo
    esta definida salvo rotacion y reflexion, asi que hay que alinear antes de
    medir. Ese es el mismo tratamiento que da ``compare_methods.error_de_picos``
    y sin el una estimacion perfecta podria reportar un error arbitrario.

    La novedad es que aqui la transformacion de alineamiento **no se descarta**.
    Sobre datos sinteticos la escala verdadera se conoce, de modo que el giro y
    la reflexion necesarios informan sobre la etapa de sincronizacion, que
    ninguna validacion previa de CIRCUST mide:

    ``verdad_sinc_error_h``
        Desviacion del origen del ciclo respecto de donde deberia estar segun
        el convenio de anclaje. Cero significa que la sincronizacion acerto.
        NaN si el gen ancla no llego al atlas.
    ``verdad_sinc_reflejado``
        1 si el ciclo reconstruido recorre el tiempo al reves. Es un fallo del
        criterio de direccion, no un error graduable.
    ``verdad_sinc_delta_h``
        El desplazamiento crudo de la alineacion, que se conserva por
        trazabilidad. No es interpretable por si solo, ver
        :func:`_error_de_sincronizacion`.
    """
    d = _comunes(atlas, params)
    if len(d) < 3:
        return {"verdad_picos_medae_h": np.nan, "verdad_picos_p75_h": np.nan,
                "verdad_picos_n": len(d), "verdad_sinc_error_h": np.nan,
                "verdad_sinc_reflejado": np.nan, "verdad_sinc_delta_h": np.nan}

    est_h = d["peak_phase_rad"].to_numpy(dtype=float) * period / DOS_PI
    ver_h = d["t_U_h"].to_numpy(dtype=float)

    ali = align_phases_by_mad(est_h, ver_h, period=period)
    err = circular_abs_diff(ali.aligned_pred, ver_h, period)
    return {
        "verdad_picos_medae_h":   float(np.median(err)),
        "verdad_picos_p75_h":     float(np.percentile(err, 75)),
        "verdad_picos_n":         int(len(d)),
        "verdad_sinc_error_h":    _error_de_sincronizacion(
            ali.delta, ali.sign, params, anchor_gene, period),
        "verdad_sinc_reflejado":  int(ali.sign < 0),
        "verdad_sinc_delta_h":    _centrar(ali.delta, period),
    }


# ===========================================================================
# 2. Omega
# ===========================================================================

def error_omega(atlas: pd.DataFrame, params: pd.DataFrame) -> dict:
    """Error y sesgo del apuntamiento estimado.

    ``omega`` es invariante frente a rotacion y a reflexion del eje temporal,
    a diferencia de ``alpha`` y ``beta``, de modo que se compara en bruto sin
    ninguna alineacion previa.

    El **sesgo con signo** es la cifra importante. Un sesgo negativo significa
    subestimar omega, o sea estimar formas mas apuntadas que las reales, que es
    exactamente lo que se observa en babuino al ajustar con doce muestras y el
    FMM base. Un error absoluto pequeno con sesgo nulo y un error absoluto
    pequeno con sesgo sistematico son situaciones muy distintas.

    Devuelve NaN si el atlas no trae la columna, en vez de fallar, para que un
    atlas reducido a las fases de pico siga siendo evaluable en lo demas.
    """
    if "omega" not in atlas.columns:
        return {"verdad_omega_medae": np.nan, "verdad_omega_sesgo": np.nan,
                "verdad_omega_n": 0}

    d = _comunes(atlas, params)
    if d.empty:
        return {"verdad_omega_medae": np.nan, "verdad_omega_sesgo": np.nan,
                "verdad_omega_n": 0}

    err = (d["omega_est"].to_numpy(dtype=float)
           - d["omega_ver"].to_numpy(dtype=float))
    return {
        "verdad_omega_medae": float(np.median(np.abs(err))),
        "verdad_omega_sesgo": float(np.median(err)),
        "verdad_omega_n":     int(len(d)),
    }


def error_omega_por_clase(atlas: pd.DataFrame,
                          params: pd.DataFrame) -> pd.DataFrame:
    """El error de omega desglosado por clase de forma de onda.

    Permite responder la pregunta concreta del bloque R: la regularizacion,
    que ayuda a los genes apuntados, no estropea los sinusoidales.
    """
    d = _comunes(atlas, params)
    filas = []
    for clase, sub in d.groupby("clase", sort=True):
        err = (sub["omega_est"].to_numpy(dtype=float)
               - sub["omega_ver"].to_numpy(dtype=float))
        filas.append({
            "clase":              clase,
            "verdad_omega_medae": float(np.median(np.abs(err))),
            "verdad_omega_sesgo": float(np.median(err)),
            "verdad_omega_n":     int(len(sub)),
        })
    return pd.DataFrame(filas)


# ===========================================================================
# 3. Conjunto TOP
# ===========================================================================

def deteccion_top(genes_top: list[str], params: pd.DataFrame,
                  period: float = PERIOD_HOURS) -> dict:
    """Calidad del conjunto TOP frente a la verdad.

    ``verdad_top_precision``
        Fraccion de los genes seleccionados que son de verdad ritmicos. Los
        que no lo son son ruido que aparento ritmicidad, y arrastran la etapa
        robusta.
    ``verdad_top_hueco_max_h``
        Mayor hueco angular, en horas, entre picos verdaderos consecutivos de
        los genes del TOP. Mide la cobertura del ciclo, que es justo lo que
        persigue el cupo por sector, y es mas informativo que una dispersion
        promedio porque un unico sector vacio ya compromete la reconstruccion.
    """
    vacio = {"verdad_top_n": len(genes_top), "verdad_top_precision": np.nan,
             "verdad_top_n_ritmicos": 0, "verdad_top_hueco_max_h": np.nan}
    if not genes_top:
        return vacio

    ritmicos = params.set_index("gene")
    dentro = [g for g in genes_top if g in ritmicos.index]

    salida = dict(vacio)
    salida["verdad_top_precision"] = len(dentro) / len(genes_top)
    salida["verdad_top_n_ritmicos"] = len(dentro)
    if len(dentro) >= 2:
        fases = np.sort(ritmicos.loc[dentro, "t_U_h"].to_numpy(dtype=float)
                        % period)
        # El circulo se cierra, luego el ultimo hueco va del mayor al primero.
        huecos = np.diff(np.append(fases, fases[0] + period))
        salida["verdad_top_hueco_max_h"] = float(huecos.max())
    return salida


# ===========================================================================
# Agregado
# ===========================================================================

def evaluar_contra_verdad(atlas: pd.DataFrame, params: pd.DataFrame, *,
                          genes_top: list[str] | None = None,
                          period: float = PERIOD_HOURS) -> dict:
    """Las tres familias de metricas en una sola fila.

    ``genes_top`` es opcional porque no todos los usos disponen de la lista.
    Si no se pasa, las metricas del conjunto TOP salen a NaN en vez de fallar.
    """
    fila = {}
    fila.update(error_picos(atlas, params, period))
    fila.update(error_omega(atlas, params))
    fila.update(deteccion_top(genes_top or [], params, period))
    return fila
