"""Tests de las extensiones del evaluador.

``comparar_unidad`` es codigo compartido entre los experimentos con datos
reales y el estudio de simulacion, de modo que estos tests vigilan tanto lo
nuevo como que lo anterior siga funcionando igual.

El test de integracion corre el pipeline completo sobre una celda pequena, asi
que tarda unos segundos. Es deliberado: lo que hay que comprobar es que las
piezas encajan de verdad, no que un mock devuelve lo que le han dicho.
"""
import numpy as np
import pytest

from simulation.cells import GENE_COLUMN, read_truth, write_cell
from simulation.generate import generate
from simulation.truth_metrics import evaluar_contra_verdad
from validation.compare_methods import METODOS, comparar_unidad, resolver_metodos


# ===========================================================================
# Seleccion de metodos
# ===========================================================================

def test_sin_seleccion_salen_los_seis_metodos():
    assert resolver_metodos() == list(METODOS)


def test_se_puede_pedir_un_subconjunto_por_etiqueta():
    r = resolver_metodos(["best_k:r2_median", "aggregate:tsp3"])
    assert [e for e, _ in r] == ["best_k:r2_median", "aggregate:tsp3"]
    assert r[0][1]["selection_criterion"] == "r2_median"


def test_una_sola_etiqueta_tambien_vale():
    assert [e for e, _ in resolver_metodos("median")] == ["median"]


def test_se_pueden_pasar_pares_explicitos():
    # Es la via del barrido de tau: el mismo metodo repetido con etiquetas
    # distintas y un parametro cambiado.
    par = ("best_k:n_well_fitted@0.70",
           dict(method="best_k", selection_criterion="n_well_fitted"))
    assert resolver_metodos([par]) == [par]


def test_un_metodo_desconocido_falla_con_mensaje_claro():
    with pytest.raises(ValueError, match="inventado"):
        resolver_metodos(["inventado"])


# ===========================================================================
# Integracion sobre una celda sintetica
# ===========================================================================

@pytest.fixture(scope="module")
def celda(tmp_path_factory):
    destino = tmp_path_factory.mktemp("celda") / "rep00"
    ds = generate("celda", n_samples=24, n_genes=3000, prop_rhythmic=0.05,
                  sampling="aleatorio", waveform_mix="gtex", snr=1.5, seed=7)
    ruta = write_cell(ds, destino)
    return ruta, read_truth(destino)


def test_el_evaluador_lee_parquet_si_se_le_dice_la_columna_de_genes(celda):
    ruta, verdad = celda
    assert ruta.suffix == ".parquet"
    filas = comparar_unidad("c", ruta, verdad.tiempos.to_dict(), [100],
                            ["base"], metodos=["best_k:r2_median"],
                            gene_column=GENE_COLUMN)
    assert len(filas) == 1
    assert filas[0]["metodo"] == "best_k:r2_median"


def test_solo_se_evaluan_los_metodos_pedidos(celda):
    ruta, verdad = celda
    filas = comparar_unidad("c", ruta, verdad.tiempos.to_dict(), [100],
                            ["base"], metodos=["median", "aggregate:tsp3"],
                            gene_column=GENE_COLUMN)
    assert [f["metodo"] for f in filas] == ["median", "aggregate:tsp3"]


def test_el_desacuerdo_entre_repeticiones_se_reporta(celda):
    ruta, verdad = celda
    filas = comparar_unidad("c", ruta, verdad.tiempos.to_dict(), [100],
                            ["base"], metodos=["best_k:r2_median"],
                            gene_column=GENE_COLUMN)
    sd = filas[0]["sd_rho_reps"]
    assert np.isfinite(sd) and sd >= 0.0


def test_el_atlas_de_agregacion_ya_trae_omega(celda):
    # Antes solo devolvia la fase de pico, lo que dejaba a los metodos de
    # agregacion sin poder compararse en apuntamiento.
    ruta, verdad = celda
    capturado = {}

    def guardar(rob, top, atlas, fmm_method):
        capturado["atlas"] = atlas
        return {}

    comparar_unidad("c", ruta, verdad.tiempos.to_dict(), [100], ["base"],
                    metodos=["aggregate:tsp3"], gene_column=GENE_COLUMN,
                    extra=guardar)
    assert "omega" in capturado["atlas"].columns
    assert capturado["atlas"]["omega"].between(0.0, 1.0).all()


def test_el_callback_extra_engancha_las_metricas_contra_la_verdad(celda):
    ruta, verdad = celda

    def contra_verdad(rob, top, atlas, fmm_method):
        return evaluar_contra_verdad(atlas, verdad.params,
                                     genes_top=list(rob.top_gene_names))

    filas = comparar_unidad("c", ruta, verdad.tiempos.to_dict(), [100],
                            ["base"], metodos=["best_k:r2_median"],
                            gene_column=GENE_COLUMN, extra=contra_verdad)
    fila = filas[0]
    assert 0.0 <= fila["verdad_top_precision"] <= 1.0
    assert np.isfinite(fila["verdad_picos_medae_h"])


def test_saltarse_la_referencia_estimada_ahorra_sus_columnas(celda):
    ruta, verdad = celda
    filas = comparar_unidad("c", ruta, verdad.tiempos.to_dict(), [100],
                            ["base"], metodos=["best_k:r2_median"],
                            gene_column=GENE_COLUMN, con_referencia=False)
    assert "picos_medae_h" not in filas[0]
    # Lo que no depende de la referencia sigue estando.
    assert "rho" in filas[0]


def test_robust_kwargs_llega_al_estimador(celda):
    """El umbral tau tiene que poder cambiarse desde fuera.

    Es la via del barrido de calibracion. Para el criterio
    ``n_well_fitted``, ``selection_score`` es el numero de genes del TOP que
    superan el umbral, de modo que subirlo a un valor imposible tiene que
    llevarlo a cero. Si el parametro no llegase, el score no se moveria.
    """
    ruta, verdad = celda
    comun = dict(metodos=["best_k:n_well_fitted"], gene_column=GENE_COLUMN,
                 con_referencia=False)
    capturado = {}

    def guardar(clave):
        def _f(rob, top, atlas, fmm_method):
            capturado[clave] = rob.selection_score
            return {}
        return _f

    comparar_unidad("c", ruta, verdad.tiempos.to_dict(), [100], ["base"],
                    robust_kwargs={"selection_r2_threshold": 0.60},
                    extra=guardar("bajo"), **comun)
    comparar_unidad("c", ruta, verdad.tiempos.to_dict(), [100], ["base"],
                    robust_kwargs={"selection_r2_threshold": 0.999},
                    extra=guardar("imposible"), **comun)
    assert capturado["bajo"] > 0
    assert capturado["imposible"] == 0
