"""
Tablas y figuras del capitulo de simulaciones.

Lee el CSV largo que produce ``run.py`` y produce las tablas que van a la
memoria. El grueso trabaja sobre el bloque ``rejilla``, donde cada escenario
y replica tiene las 30 variantes (3 core x 2 FMM x 5 metodos) corridas sobre
la misma matriz.

Las tres preguntas y como se responden
--------------------------------------
1. **Reconstruccion.** Como se comporta CIRCUST tal cual (FMM base, best_k
   con r2_median) a lo largo de los factores. ``tabla_reconstruccion``.
2. **Regularizacion del FMM.** Diferencia PAREADA CircaFMM menos base,
   replica a replica y con el metodo de resumen fijo. ``contraste_fmm``.
3. **Agregacion.** Diferencia PAREADA de cada metodo de resumen contra la
   referencia, con el FMM fijo. ``contraste_metodos``.

Sobre cualquiera de los dos contrastes se aplican las mismas tres lecturas:
``efectos_principales`` (mediana del delta, dispersion y tasa de victorias
por nivel de cada factor), ``interaccion`` (mapa de dos factores) y
``modelo_lineal`` (que factores modulan el delta, con numeros). Es lo que
responde a "recorrer las mejoras a lo largo de los hiperparametros".

Tres decisiones de agregacion, todas tomadas midiendo y no por costumbre
-----------------------------------------------------------------------
1. **Las medidas continuas se agregan con MEDIANA, no con media.** ``rho_circ``
   cambia de signo en torno al 1-2 % de las replicas aunque la reconstruccion
   sea perfecta, porque su direccion media queda mal determinada cuando las
   muestras cubren todo el ciclo. Una sola replica asi arrastra la media de
   una celda. Por eso la correlacion que lidera es ``rho_fl``, la de Fisher y
   Lee, que no necesita direccion media.

2. **La orientacion se agrega como TASA de acierto**, no como mediana, porque
   es una proporcion. Y va en columna propia porque ``error_fase_h`` mezcla dos
   cosas: si el sentido de giro sale invertido da 5 a 7 h aunque el orden sea
   perfecto, con lo que en las celdas donde la orientacion no es decidible esa
   metrica mide suerte y no reconstruccion.

3. **Las filas con fallo se excluyen de las medianas pero su proporcion se
   reporta** como columna. La tasa de fallo por celda es un resultado del
   estudio, es la frontera de aplicabilidad del metodo.

Uso
---
    python -m simulation.analyze --in output/simulacion --out output/simulacion
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

# Factores de generacion y de configuracion que identifican una celda.
FACTORES = ["n_samples", "sampling", "omega", "snr", "n_genes",
            "prop_rhythmic", "n_core"]

# Identidad de una replica: los factores mas el numero de replica. Dos filas
# con la misma clave corrieron sobre la MISMA matriz.
CLAVE_REPLICA = FACTORES + ["rep"]

FMM_REFERENCIA    = "base"
METODO_REFERENCIA = "best_k:r2_median"

# Nombre largo de cada medida, el que ven las figuras y las tablas del
# capitulo. Las columnas del CSV llevan el nombre corto de la izquierda.
ETIQUETAS = {
    "error_orden_h":             "error del orden relativo (h)",
    "error_fase_h":              "error de fase absoluta (h)",
    "error_fase_medio_h":        "error de fase absoluta, media (h)",
    "concordancia":              "concordancia circular",
    "cre":                       "error relativo de concordancia (CRE)",
    "tau_circular":              "tau circular con signo",
    "acierto_sentido":           "acierto del sentido de giro",
    "atlas_error_onda":          "error de la onda del atlas",
    "atlas_variabilidad_pico_h": "variabilidad del error en la hora del pico (h)",
    "atlas_sesgo_omega":         "sesgo de omega",
    "atlas_vpp":                 "fraccion de genes ritmicos en el atlas",
    "rho_fl":                    "correlacion circular de Fisher y Lee",
    "rho_circ":                  "correlacion circular de Jammalamadaka",
    "tau":                       "tau de Kendall con corte fijo",
    "n_top":                     "tamano del conjunto TOP",
    "snr":                       "A/sigma",
    "n_samples":                 "n muestras",
    "sampling":                  "muestreo",
    "omega":                     "omega",
    "n_genes":                   "n genes",
    "prop_rhythmic":             "proporcion de ritmicos",
    "n_core":                    "genes core declarados",
}


def etiqueta(clave: str) -> str:
    """Nombre largo de una medida o factor, o la clave si no tiene."""
    return ETIQUETAS.get(clave, clave)

# Medidas continuas, se agregan con mediana.
#
# Va ``concordancia`` y no ``tau_circular``. ``tau_circular`` cambia de signo cuando el
# ciclo se recorre al reves, de modo que su SIGNO y su MAGNITUD responden a dos
# preguntas distintas: el signo dice la orientacion y la magnitud la calidad
# del orden. ``acierto_sentido`` se define como ``tau_circular > 0``, asi que el
# signo ya esta contado en su propia columna. Resumir con signo mezcla las
# dos: cuatro replicas +0,8 +0,8 -0,8 -0,8 dan mediana 0,0, que se lee como
# "no se recupero nada" cuando la calidad real es 0,8 en las cuatro.
MEDIDAS = ["error_fase_h", "error_fase_medio_h", "error_orden_h", "concordancia",
           "cre", "rho_fl", "rho_circ", "tau",
           "p2", "r2_anchor", "evidencia", "n_top"]

# Proporciones, se agregan con media porque son tasas.
TASAS = ["acierto_sentido"]

# Las cuatro del atlas que van a las tablas del capitulo. Cada una responde a
# una pregunta distinta:
#
#   resumen       atlas_error_onda      libre de parametrizacion
#   interpretable atlas_variabilidad_pico_h   en horas
#   mecanismo     atlas_sesgo_omega     por que CircaFMM mejora
#   practica      atlas_vpp             cuantos genes del atlas son ritmicos
MEDIDAS_ATLAS = ["atlas_error_onda", "atlas_variabilidad_pico_h",
                 "atlas_sesgo_omega", "atlas_vpp"]

# Medidas de orden sobre las que se contrastan los metodos de resumen. El
# error del orden relativo y |tau_circular| miden el orden al margen del anclaje, error_fase_h es el
# resultado conjunto que ve el usuario. El CRE mide lo mismo desde la senal:
# cuanto se aleja el ajuste FMM con las fases estimadas del ajuste con las
# verdaderas.
MEDIDAS_ORDEN = ["error_orden_h", "concordancia", "cre", "error_fase_h"]

# El CRE se calcula con el modelo FMM de cada configuracion, asi que sus
# valores con FMM base y con CircaFMM no son comparables entre si (CircaFMM da
# curvas mas suaves y eso, por si solo, baja el CRE). Vale para contrastar
# metodos de resumen dentro de un modelo, no para el contraste entre modelos.
MEDIDAS_ORDEN_ENTRE_FMM = [m for m in MEDIDAS_ORDEN if m != "cre"]

TIEMPOS = ["t_preprocess_s", "t_cpca_s", "t_sync_s", "t_top_s",
           "t_robust_s", "t_total_s"]


# ═══════════════════════════════════════════════════════════════════════════
# Carga y columnas derivadas
# ═══════════════════════════════════════════════════════════════════════════

def _con_nivel(tabla: pd.DataFrame, columnas: list[str],
               nivel: str) -> pd.DataFrame:
    """Convierte las columnas en un MultiIndex ``(medida, nivel)``."""
    salida = tabla.copy()
    salida.columns = pd.MultiIndex.from_product([columnas, [nivel]])
    return salida


def _etiqueta_metodo(df: pd.DataFrame) -> pd.Series:
    """La misma etiqueta que ``RobustConfig.label``, reconstruida del CSV."""
    agg = df["aggregation"].astype(str)
    esc = df["consensus_scale"].astype(str)
    sel = df["selection_criterion"].astype(str)
    return pd.Series(
        np.where(df["robust_method"] == "best_k", "best_k:" + sel,
                 np.where(df["robust_method"] == "aggregate",
                          "aggregate:" + agg + ":" + esc,
                          df["robust_method"].astype(str))),
        index=df.index,
    )


def _derivadas(df: pd.DataFrame) -> pd.DataFrame:
    """
    Columnas que se calculan de otras y por eso no ocupan sitio en el CSV.

    ``concordancia`` es la calidad del orden al margen del sentido de giro.
    ``metodo`` es la etiqueta del metodo de resumen.
    """
    df = df.copy()
    if "tau_circular" in df.columns and "concordancia" not in df.columns:
        df["concordancia"] = df["tau_circular"].abs()
    if "robust_method" in df.columns and "metodo" not in df.columns:
        df["metodo"] = _etiqueta_metodo(df)
    return df


def load(path: Path | str) -> pd.DataFrame:
    """Carga el CSV largo de resultados."""
    return _derivadas(pd.read_csv(path, low_memory=False))


def _rejilla(df: pd.DataFrame) -> pd.DataFrame:
    r = _derivadas(df)
    r = r[r["bloque"] == "rejilla"]
    if r.empty:
        raise ValueError("no hay filas del bloque rejilla")
    return r


def _referencia(df: pd.DataFrame) -> pd.DataFrame:
    """CIRCUST tal cual: FMM base y best_k con r2_median."""
    r = _rejilla(df)
    return r[(r["fmm_method"] == FMM_REFERENCIA)
             & (r["metodo"] == METODO_REFERENCIA)]


# ═══════════════════════════════════════════════════════════════════════════
# Resumen por niveles
# ═══════════════════════════════════════════════════════════════════════════

def _resumen(df: pd.DataFrame, por: list[str]) -> pd.DataFrame:
    """
    Mediana y dispersion de las medidas por grupo, tasas, y tasa de fallo.

    Las filas con fallo no entran en las medianas pero si en el denominador
    de la tasa de fallo.
    """
    df = _derivadas(df)
    fallo = df.groupby(por)["status"].apply(
        lambda s: float((s == "fallo").mean())
    )
    ok = df[df["status"] == "ok"]
    g  = ok.groupby(por)

    medidas = [m for m in MEDIDAS if m in ok.columns]
    tabla = _con_nivel(g[medidas].median(), medidas, "mediana")
    iqr   = _con_nivel(g[medidas].quantile(0.75) - g[medidas].quantile(0.25),
                       medidas, "iqr")
    tabla = pd.concat([tabla, iqr], axis=1)

    for t in (t for t in TASAS if t in ok.columns):
        tabla[(t, "tasa")] = g[t].mean()
    tabla[("tasa_fallo", "")] = fallo
    tabla[("n", "")] = df.groupby(por).size()
    return tabla.sort_index(axis=1)


def tabla_reconstruccion(df: pd.DataFrame,
                         por: list[str] | None = None) -> pd.DataFrame:
    """
    Pregunta 1: como reconstruye CIRCUST tal cual a lo largo de los factores.

    Por defecto se agrupa por muestreo, forma de onda y senal-ruido, que son
    los factores que mas discriminan. Cualquier otro subconjunto de
    ``FACTORES`` vale.
    """
    ref = _referencia(df)
    if ref.empty:
        raise ValueError("no hay filas de la configuracion de referencia")
    return _resumen(ref, por or ["sampling", "omega", "snr"])


def tabla_por_factor(df: pd.DataFrame) -> pd.DataFrame:
    """
    Efectos principales de la reconstruccion: una fila por factor y nivel,
    promediando sobre todo lo demas. Es la vista de una sola pagina.
    """
    ref = _referencia(df)
    trozos = []
    for f in FACTORES:
        t = _resumen(ref, [f]).reset_index()
        t.insert(0, "factor", f)
        t = t.rename(columns={f: "nivel"})
        trozos.append(t)
    salida = pd.concat(trozos, ignore_index=True)
    salida["nivel"] = salida["nivel"].astype(str)
    return salida.set_index(["factor", "nivel"])


def tabla_sin_ruido(df: pd.DataFrame) -> pd.DataFrame:
    """
    Verificacion de que el pipeline funciona con datos perfectos.

    Cualquier error aqui no puede venir del ruido, asi que aisla los sesgos
    del propio algoritmo. Se agrupa por forma de onda y muestreo, que son los
    factores que discriminan.
    """
    s = _derivadas(df)
    s = s[s["bloque"] == "sin_ruido"]
    if s.empty:
        raise ValueError("no hay filas del bloque sin_ruido")
    por = [c for c in ("omega", "sampling", "n_samples") if c in s.columns]
    return _resumen(s, por[:2] if len(por) > 1 else por)


def tabla_invariancia(df: pd.DataFrame) -> pd.DataFrame:
    """
    Bloque de invariancia: comprobacion de que el pipeline solo ve A/sigma.

    Los dos miembros de cada par tienen el mismo cociente y distinta escala
    absoluta, asi que deben dar metricas identicas.
    """
    inv = _derivadas(df)
    inv = inv[(inv["bloque"] == "invariancia") & (inv["status"] == "ok")]
    if inv.empty:
        raise ValueError("no hay filas del bloque invariancia")
    medidas = [m for m in MEDIDAS if m in inv.columns]
    return inv.groupby(["snr", "sigma"])[medidas].median()


def tabla_atlas(df: pd.DataFrame,
                metodo: str = METODO_REFERENCIA) -> pd.DataFrame:
    """
    Calidad del atlas en NIVELES, base al lado de CircaFMM.

    ``contraste_fmm`` da la diferencia pareada, que es la que tiene potencia,
    pero una diferencia sola no dice si el atlas es utilizable. Esta tabla
    da los valores absolutos para poder decir, por ejemplo, que a A/sigma=4
    el error de pico esta por debajo de la hora con los dos metodos.

    Solo se resumen las filas con la orientacion bien resuelta. Una escala
    reflejada da un atlas especular, con errores de pico repartidos por todo
    el ciclo, y mezclarlos con los aciertos no promedia dos calidades sino
    dos regimenes distintos.
    """
    g = _rejilla(df)
    g = g[(g["metodo"] == metodo) & (g["status"] == "ok")
          & (g["acierto_sentido"] == 1.0)]
    if g.empty:
        raise ValueError("no hay filas con orientacion resuelta en la rejilla")
    medidas = [m for m in MEDIDAS_ATLAS if m in g.columns]
    por = ["snr", "omega", "fmm_method"]
    gr = g.groupby(por)
    tabla = _con_nivel(gr[medidas].median(), medidas, "mediana")
    iqr   = _con_nivel(gr[medidas].quantile(0.75) - gr[medidas].quantile(0.25),
                       medidas, "iqr")
    tabla[("n_celdas", "")] = gr.size()
    return pd.concat([tabla, iqr], axis=1).sort_index(axis=1)


# ═══════════════════════════════════════════════════════════════════════════
# Contrastes pareados
# ═══════════════════════════════════════════════════════════════════════════

def _emparejar(a: pd.DataFrame, b: pd.DataFrame, medidas: list[str],
               sufijos: tuple[str, str]) -> pd.DataFrame:
    """
    Une dos brazos replica a replica y calcula ``delta_<medida>`` como
    ``b - a``. Solo entran las replicas con ``status == ok`` en los dos
    brazos: una replica que fallo en uno no tiene con que emparejarse.
    """
    a = a[a["status"] == "ok"]
    b = b[b["status"] == "ok"]
    if a.empty or b.empty:
        raise ValueError("uno de los dos brazos del contraste esta vacio")
    columnas = CLAVE_REPLICA + [m for m in medidas if m in a.columns]
    par = a[columnas].merge(b[columnas], on=CLAVE_REPLICA, suffixes=sufijos)
    if par.empty:
        raise ValueError("el emparejamiento no cuadra, revisa las semillas")
    for m in medidas:
        ca, cb = f"{m}{sufijos[0]}", f"{m}{sufijos[1]}"
        if ca in par.columns and cb in par.columns:
            par[f"delta_{m}"] = par[cb] - par[ca]
    return par


def contraste_fmm(df: pd.DataFrame, metodo: str = METODO_REFERENCIA,
                  solo_orientadas: bool = True) -> pd.DataFrame:
    """
    Pregunta 2: CircaFMM menos base, pareado, con el metodo de resumen
    fijo.

    El orden de muestras sale de una SVD sobre los genes nucleares donde el
    FMM no interviene, asi que sobre el orden los dos son casi identicos por
    construccion. Donde se diferencian es en el atlas, y por eso las medidas
    del atlas van delante. Las de orden se conservan como comprobacion.

    Con ``solo_orientadas`` se exige que la orientacion este resuelta en los
    DOS brazos: una escala reflejada da un atlas especular, que es otro
    regimen y no una calidad peor.
    """
    r = _rejilla(df)
    r = r[r["metodo"] == metodo]
    medidas = MEDIDAS_ATLAS + MEDIDAS_ORDEN_ENTRE_FMM + ["acierto_sentido", "t_total_s"]
    par = _emparejar(r[r["fmm_method"] == "base"],
                     r[r["fmm_method"] == "circafmm"],
                     medidas, ("_base", "_circafmm"))
    if solo_orientadas:
        par = par[(par["acierto_sentido_base"] == 1.0)
                  & (par["acierto_sentido_circafmm"] == 1.0)]
    return par.reset_index(drop=True)


def contraste_metodos(df: pd.DataFrame, fmm: str = FMM_REFERENCIA,
                      referencia: str = METODO_REFERENCIA) -> pd.DataFrame:
    """
    Pregunta 3: cada metodo de resumen menos la referencia, pareado, con el
    FMM fijo.

    Devuelve un DataFrame largo con una columna ``metodo`` y los
    ``delta_<medida>`` de ese metodo contra la referencia en la misma
    replica. Negativo en el error significa que el metodo gana.
    """
    r = _rejilla(df)
    r = r[r["fmm_method"] == fmm]
    ref = r[r["metodo"] == referencia]
    medidas = MEDIDAS_ORDEN + ["acierto_sentido", "t_robust_s"]
    trozos = []
    for metodo, sub in r.groupby("metodo"):
        if metodo == referencia:
            continue
        par = _emparejar(ref, sub, medidas, ("_ref", "_met"))
        par.insert(0, "metodo", metodo)
        trozos.append(par)
    if not trozos:
        raise ValueError("no hay metodos distintos de la referencia")
    return pd.concat(trozos, ignore_index=True)


# ═══════════════════════════════════════════════════════════════════════════
# Lecturas de un contraste
# ═══════════════════════════════════════════════════════════════════════════

def _resumir_delta(d: pd.Series) -> dict:
    d = d.dropna()
    if d.empty:
        return dict(n=0, mediana=np.nan, iqr=np.nan, tasa_victorias=np.nan)
    return dict(
        n              = int(d.size),
        mediana        = float(d.median()),
        iqr            = float(d.quantile(0.75) - d.quantile(0.25)),
        # la fraccion de replicas donde el brazo nuevo mejora (delta < 0).
        # Un delta mediano pequeno con el 100 % de victorias es una mejora
        # mas solida que uno grande con el 55 %.
        tasa_victorias = float((d < 0).mean()),
    )


def efectos_principales(par: pd.DataFrame, medida: str,
                        factores: list[str] | None = None,
                        por: list[str] | None = None) -> pd.DataFrame:
    """
    Mediana, dispersion y tasa de victorias del delta por nivel de cada
    factor, promediando sobre los demas.

    ``por`` permite separar ademas por una columna extra, por ejemplo
    ``["metodo"]`` en el contraste de metodos.
    """
    col = f"delta_{medida}"
    if col not in par.columns:
        raise KeyError(f"no existe la columna {col}")
    factores = factores or FACTORES
    extra = por or []
    filas = []
    for f in factores:
        for claves, g in par.groupby(extra + [f]):
            claves = claves if isinstance(claves, tuple) else (claves,)
            fila = dict(zip(extra, claves[:-1]))
            fila.update(factor=f, nivel=str(claves[-1]))
            fila.update(_resumir_delta(g[col]))
            filas.append(fila)
    tabla = pd.DataFrame(filas)
    return tabla.set_index(extra + ["factor", "nivel"])


def interaccion(par: pd.DataFrame, medida: str, fila: str, columna: str,
                estadistico: str = "mediana") -> pd.DataFrame:
    """
    Mapa de dos factores: mediana del delta o tasa de victorias.

    Filas ``fila``, columnas ``columna``, promediando sobre el resto.
    """
    col = f"delta_{medida}"
    if estadistico == "mediana":
        return par.pivot_table(index=fila, columns=columna, values=col,
                               aggfunc="median")
    if estadistico == "victorias":
        return par.assign(gana=lambda x: (x[col] < 0).astype(float)) \
                  .pivot_table(index=fila, columns=columna, values="gana")
    raise ValueError(f"estadistico desconocido: {estadistico!r}")


def modelo_lineal(par: pd.DataFrame, medida: str,
                  factores: list[str] | None = None) -> pd.DataFrame:
    """
    Regresion del delta sobre los factores como categoricos, sin
    interacciones. Responde con numeros a que factores modulan la mejora.

    Cada coeficiente es el cambio del delta al pasar del nivel de referencia
    (el menor de cada factor) al nivel indicado, a igualdad del resto. Se
    reporta con su error tipico y su p-valor, y el R2 del modelo dice cuanta
    de la variacion del delta explican los factores en conjunto. Es un
    modelo descriptivo: las replicas de una celda no son independientes del
    todo, asi que los p-valores son orientativos y lo que hay que mirar es
    la magnitud y el signo.
    """
    import statsmodels.formula.api as smf

    col = f"delta_{medida}"
    factores = factores or FACTORES
    datos = par[[col] + factores].dropna().copy()
    if datos.empty:
        raise ValueError(f"no hay datos para el modelo de {col}")
    for f in factores:
        datos[f] = datos[f].astype(str)
    datos = datos.rename(columns={col: "delta"})
    formula = "delta ~ " + " + ".join(f"C({f})" for f in factores)
    ajuste = smf.ols(formula, data=datos).fit()
    tabla = pd.DataFrame({
        "coef": ajuste.params,
        "error_tipico": ajuste.bse,
        "p": ajuste.pvalues,
    })
    tabla.attrs["r2"] = float(ajuste.rsquared)
    tabla.attrs["n"] = int(ajuste.nobs)
    return tabla


# ═══════════════════════════════════════════════════════════════════════════
# Tiempos
# ═══════════════════════════════════════════════════════════════════════════

def tabla_tiempos_python(df: pd.DataFrame) -> pd.DataFrame:
    """
    Coste en Python por variante del FMM, metodo de resumen y tamano.

    Las etapas 0 a 3 son comunes a los metodos de resumen, asi que su tiempo
    se repite en las cinco filas de un grupo. ``t_robust_s`` es lo unico que
    cambia entre metodos.
    """
    r = _rejilla(df)
    r = r[r["status"] == "ok"]
    tiempos = [t for t in TIEMPOS if t in r.columns]
    return r.groupby(["fmm_method", "metodo", "n_genes", "n_samples"])[
        tiempos].median()


def tabla_tiempos(df_tiempos: pd.DataFrame) -> pd.DataFrame:
    """Bloque C: tiempo por implementacion y tamano de matriz, R contra Python."""
    agg   = df_tiempos.groupby(["implementacion", "n_genes", "n_samples"])
    ancho = agg["t_total_s"].median().unstack("implementacion")
    if {"python", "r"} <= set(ancho.columns):
        ancho["factor_r_sobre_python"] = ancho["r"] / ancho["python"]
    return ancho


# ═══════════════════════════════════════════════════════════════════════════
# Calibracion del aviso de orientacion
# ═══════════════════════════════════════════════════════════════════════════

def calibracion_orientacion(df: pd.DataFrame, indicador: str = "evidencia",
                            n_tramos: int = 8) -> pd.DataFrame:
    """
    Relacion entre un indicador **observable sin la verdad** y el acierto de
    la orientacion.

    Es la base del aviso de fase absoluta: en un tejido nuevo no se puede
    medir el error, pero si se puede medir el indicador, y esta tabla dice
    que esperar para cada valor. El umbral no se fija aqui, se elige mirando
    la curva con un criterio explicito de precision y cobertura.
    """
    s = _derivadas(df)
    s = s[(s["status"] == "ok") & s[indicador].notna()].copy()
    if s.empty:
        raise ValueError(f"no hay filas con {indicador} definido")
    s["tramo"] = pd.qcut(s[indicador], n_tramos, duplicates="drop")
    g = s.groupby("tramo", observed=True)
    return pd.DataFrame({
        "n":                  g.size(),
        f"{indicador}_medio": g[indicador].median(),
        "acierto_orientacion": g["acierto_sentido"].mean(),
        "error_fase_h":            g["error_fase_h"].median(),
    })


# ═══════════════════════════════════════════════════════════════════════════
# Figuras
# ═══════════════════════════════════════════════════════════════════════════

def _mpl():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 9, "figure.dpi": 160,
                         "savefig.bbox": "tight"})
    return plt


def _figura_lineas(datos, x, y, hue, titulo, xlabel, ylabel, out: Path) -> None:
    """
    Grafico de lineas con una serie por nivel de ``hue``.

    El eje x lleva marcas exactamente en los niveles del factor y con sus
    valores tal cual (1, 1,5, 2, 4), sin escala logaritmica: con cuatro
    niveles las potencias de 2 solo confunden.
    """
    # Sin esto se escribian PNG vacios cuando el bloque todavia no se habia
    # ejecutado, y un fichero vacio parece una figura valida.
    if len(datos) == 0:
        raise ValueError(f"no hay datos para la figura {out.name}")
    plt = _mpl()
    fig, ax = plt.subplots(figsize=(7.0, 4.2))
    for nivel, sub in datos.groupby(hue):
        sub = sub.sort_values(x)
        ax.plot(sub[x], sub[y], marker="o", label=f"{etiqueta(hue)} = {nivel}")
    niveles = sorted(datos[x].unique())
    ax.set_xticks(niveles, [f"{v:g}" for v in niveles])
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(titulo)
    ax.legend(frameon=False)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out)
    plt.close(fig)


def figura_snr(df: pd.DataFrame, out: Path) -> None:
    """Degradacion del error con la senal-ruido, una serie por muestreo."""
    p = _referencia(df)
    p = p[p["status"] == "ok"]
    datos = p.groupby(["sampling", "snr"])["error_fase_h"].median().reset_index()
    _figura_lineas(datos, "snr", "error_fase_h", "sampling",
                   "Reconstruccion frente a la razon senal-ruido",
                   etiqueta("snr"), "mediana del " + etiqueta("error_fase_h"), out)


def figura_orientacion(df: pd.DataFrame, out: Path) -> None:
    """Tasa de acierto de orientacion, que va aparte del error."""
    p = _referencia(df)
    p = p[p["status"] == "ok"]
    datos = p.groupby(["omega", "snr"])["acierto_sentido"].mean().reset_index()
    _figura_lineas(datos, "snr", "acierto_sentido", "omega",
                   "Acierto del sentido de giro",
                   etiqueta("snr"), "fraccion de replicas con el sentido correcto", out)


def figura_efectos(par: pd.DataFrame, medida: str, titulo: str,
                   out: Path, por: str | None = None) -> None:
    """
    Efectos principales de un contraste: un panel por factor, la mediana del
    delta por nivel con su rango intercuartilico. Con ``por`` se dibuja una
    serie por nivel de esa columna (por ejemplo, por metodo).
    """
    col = f"delta_{medida}"
    if col not in par.columns or par[col].dropna().empty:
        raise ValueError(f"no hay datos para la figura {out.name}")
    plt = _mpl()
    fig, axes = plt.subplots(1, len(FACTORES), figsize=(2.3 * len(FACTORES), 3.2),
                             sharey=True)
    series = [None] if por is None else sorted(par[por].unique())
    for ax, f in zip(axes, FACTORES):
        for i, s in enumerate(series):
            sub = par if s is None else par[par[por] == s]
            g = sub.groupby(f)[col]
            m, q1, q3 = g.median(), g.quantile(0.25), g.quantile(0.75)
            x = np.arange(len(m)) + (0 if s is None else (i - len(series) / 2) * 0.12)
            ax.errorbar(x, m.values, yerr=[m.values - q1.values, q3.values - m.values],
                        fmt="o", ms=4, capsize=2,
                        label=None if s is None else str(s))
        ax.set_xticks(np.arange(len(m)), [f"{v:g}" if isinstance(v, float) else str(v)
                                          for v in m.index], rotation=45)
        ax.axhline(0, color="0.4", lw=0.8, ls="--")
        ax.set_title(etiqueta(f), fontsize=9)
    axes[0].set_ylabel("cambio en el\n" + etiqueta(medida), fontsize=8)
    if por is not None:
        axes[-1].legend(frameon=False, fontsize=7, title=por)
    fig.suptitle(titulo, fontsize=10)
    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out)
    plt.close(fig)


# ═══════════════════════════════════════════════════════════════════════════
# CLI
# ═══════════════════════════════════════════════════════════════════════════

def _escribir(nombre: str, fn, destino: Path) -> None:
    try:
        t = fn()
        t.to_csv(destino / f"{nombre}.csv")
        extra = ""
        if "r2" in t.attrs:
            extra = f" (R2 = {t.attrs['r2']:.3f}, n = {t.attrs['n']})"
        print(f"  escrita {nombre}.csv{extra}")
    except (KeyError, ValueError) as exc:
        print(f"  {nombre} omitida: {exc}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Tablas y figuras del capitulo de simulaciones.",
    )
    parser.add_argument("--in", dest="entrada", type=str,
                        default="output/simulacion")
    parser.add_argument("--out", dest="salida", type=str,
                        default="output/simulacion")
    parser.add_argument("--images", type=str, default="docs/images")
    args = parser.parse_args()

    entrada = Path(args.entrada)
    df = load(entrada / "resultados.csv")

    tablas = Path(args.salida) / "tablas"
    tablas.mkdir(parents=True, exist_ok=True)
    imagenes = Path(args.images)

    # Pregunta 1
    _escribir("tabla_reconstruccion", lambda: tabla_reconstruccion(df), tablas)
    _escribir("tabla_por_factor",     lambda: tabla_por_factor(df), tablas)
    _escribir("tabla_sin_ruido",      lambda: tabla_sin_ruido(df), tablas)
    _escribir("tabla_invariancia",    lambda: tabla_invariancia(df), tablas)
    _escribir("tabla_atlas",          lambda: tabla_atlas(df), tablas)
    _escribir("tabla_tiempos_python", lambda: tabla_tiempos_python(df), tablas)

    # Pregunta 2: CircaFMM
    try:
        fmm = contraste_fmm(df)
        for m in MEDIDAS_ATLAS + MEDIDAS_ORDEN_ENTRE_FMM:
            _escribir(f"fmm_efectos_{m}", lambda m=m: efectos_principales(fmm, m), tablas)
            _escribir(f"fmm_modelo_{m}",  lambda m=m: modelo_lineal(fmm, m), tablas)
        _escribir("fmm_interaccion_n_snr",
                  lambda: interaccion(fmm, "atlas_error_onda", "n_samples", "snr"),
                  tablas)
        _escribir("fmm_interaccion_omega_n",
                  lambda: interaccion(fmm, "atlas_error_onda", "omega", "n_samples"),
                  tablas)
        figura_efectos(fmm, "atlas_error_onda",
                       "CircaFMM menos base (negativo = mejora)",
                       imagenes / "sim-fmm-efectos.png")
        print("  escrita sim-fmm-efectos.png")
    except (KeyError, ValueError) as exc:
        print(f"  contraste FMM omitido: {exc}")

    # Pregunta 3: los metodos de resumen
    try:
        met = contraste_metodos(df)
        for m in MEDIDAS_ORDEN:
            _escribir(f"metodos_efectos_{m}",
                      lambda m=m: efectos_principales(met, m, por=["metodo"]), tablas)
        for metodo, sub in met.groupby("metodo"):
            corto = metodo.replace(":", "_")
            _escribir(f"metodos_modelo_{corto}",
                      lambda sub=sub: modelo_lineal(sub, "error_orden_h"), tablas)
        figura_efectos(met, "error_orden_h",
                       "Metodo de resumen menos best_k:r2_median (negativo = mejora)",
                       imagenes / "sim-metodos-efectos.png", por="metodo")
        print("  escrita sim-metodos-efectos.png")
    except (KeyError, ValueError) as exc:
        print(f"  contraste de metodos omitido: {exc}")

    for indicador in ("evidencia", "r2_anchor"):
        _escribir(f"calibracion_{indicador}",
                  lambda i=indicador: calibracion_orientacion(df, i), tablas)

    ruta_tiempos = entrada / "tiempos.csv"
    if ruta_tiempos.exists():
        _escribir("tabla_tiempos",
                  lambda: tabla_tiempos(pd.read_csv(ruta_tiempos)), tablas)

    for nombre, fn in (("sim-snr", figura_snr),
                       ("sim-orientacion", figura_orientacion)):
        try:
            fn(df, imagenes / f"{nombre}.png")
            print(f"  escrita {nombre}.png")
        except (KeyError, ValueError, IndexError) as exc:
            print(f"  {nombre}.png omitida: {exc}")


if __name__ == "__main__":
    main()
