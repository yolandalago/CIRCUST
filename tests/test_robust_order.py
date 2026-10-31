"""
Tests de la correccion de orientaciones de la Etapa 4.

``_fix_orientations`` detecta las repeticiones cuyo sentido de recorrido
salio invertido y las corrige antes de que se usen para nada. Decidia por el
**signo** de la correlacion de Jammalamadaka, y ese signo no es fiable cuando
las muestras cubren todo el ciclo: sobre repeticiones genuinamente invertidas
del pipeline se midieron valores positivos de +0,11 a +0,36, con lo que la
condicion ``corr < 0`` no se cumplia nunca y no se corregia ninguna.

Se usa la de Fisher y Lee, que se construye sobre diferencias por pares y no
necesita ninguna direccion media. Medido sobre 720 escenarios de inversion
conocida: Jammalamadaka acierta el signo el 95,3 % de las veces y Fisher-Lee
el 100 %. Frente a la alternativa por tripletes, que tambien acierta el
100 %, es 49 veces mas barata con n=96 (0,22 ms contra 10,63 ms).
"""
import numpy as np
import pytest

from circust.fitting.fmm import FMMModel
from circust.robust_order import (
    RobustOrderEstimator, _build_theta, _fisher_lee,
)

TWO_PI = 2.0 * np.pi


def _rep(theta_por_muestra):
    """Convierte fases por muestra en (orden, escala) como los produce la etapa."""
    orden = np.argsort(theta_por_muestra % TWO_PI)
    return orden, (theta_por_muestra % TWO_PI)[orden]


def _escenario(n=48, seed=0):
    """Preliminar mas dos reps: la 0 bien orientada y la 1 reflejada."""
    rng = np.random.default_rng(seed)
    theta_pre = np.sort(rng.uniform(0, TWO_PI, n))
    bien = (theta_pre + rng.normal(0, 0.05, n)) % TWO_PI
    mal  = (TWO_PI - bien) % TWO_PI
    o0, e0 = _rep(bien)
    o1, e1 = _rep(mal)
    return theta_pre, np.vstack([o0, o1]), np.vstack([e0, e1])


def test_fisher_lee_identifica_la_rep_invertida():
    theta_pre, ordenes, escalas = _escenario()
    Theta = _build_theta(ordenes, escalas, theta_pre.size)
    assert _fisher_lee(Theta[0], theta_pre) > 0.8
    assert _fisher_lee(Theta[1], theta_pre) < -0.8


def test_fix_orientations_corrige_la_rep_invertida():
    theta_pre, ordenes, escalas = _escenario()
    est = RobustOrderEstimator(n_reps=2, verbose=False)
    o_fix, e_fix, _, _ = est._fix_orientations(
        ordenes, escalas, None, None,
        theta_pre, theta_pre.size)
    Theta = _build_theta(o_fix, e_fix, theta_pre.size)
    assert _fisher_lee(Theta[0], theta_pre) > 0.8
    assert _fisher_lee(Theta[1], theta_pre) > 0.8, (
        "la repeticion reflejada deberia haber quedado corregida"
    )


def test_fix_orientations_no_toca_lo_que_ya_esta_bien():
    rng = np.random.default_rng(1)
    n = 48
    theta_pre = np.sort(rng.uniform(0, TWO_PI, n))
    reps = [_rep((theta_pre + rng.normal(0, 0.05, n)) % TWO_PI) for _ in range(3)]
    ordenes = np.vstack([o for o, _ in reps])
    escalas = np.vstack([e for _, e in reps])

    est = RobustOrderEstimator(n_reps=3, verbose=False)
    o_fix, e_fix, _, _ = est._fix_orientations(
        ordenes, escalas, None, None, theta_pre, n)
    assert np.array_equal(o_fix, ordenes)
    assert np.allclose(e_fix, escalas)


def test_fisher_lee_es_invariante_a_rotacion():
    """El origen de cada repeticion es arbitrario y no debe influir."""
    rng = np.random.default_rng(2)
    theta = np.sort(rng.uniform(0, TWO_PI, 40))
    for delta in (0.0, 1.0, 3.0, 5.5):
        assert _fisher_lee(theta, (theta + delta) % TWO_PI) == pytest.approx(1.0)


def test_fisher_lee_con_una_sola_muestra():
    assert np.isnan(_fisher_lee(np.array([0.0]), np.array([0.0])))


# ---------------------------------------------------------------------------
# Calibracion de lambda de CircaFMM
# ---------------------------------------------------------------------------

def test_lambda_no_se_recalibra_por_repeticion_por_defecto():
    """
    Recalibrar lambda en cada una de las K repeticiones era el valor por
    defecto y costaba el 99 % del sobrecoste de CircaFMM: 27 s por
    calibracion, seis por celda, y el bloque de regularizacion pasaba de 4 a
    20 horas.

    Y no aportaba nada, medido:

    1. La calibracion no determina lambda de forma estable. Sobre la MISMA
       escala, cambiando solo la semilla del Monte Carlo, devuelve 1,23,
       1,23, 4,12, 4,12 y 8,76. Esa dispersion (193 %) es MAYOR que la que
       hay entre repeticiones distintas (150 %), asi que lo que varia entre
       repeticiones es ruido del estimador, no informacion sobre los datos.
    2. El valor tampoco cambia el resultado donde el metodo funciona: con
       omega=0,5 y A/sigma=2, pasar de lambda 1,23 a 8,76 da el mismo orden
       bit a bit, 0,000 h de diferencia de fase.
    """
    est = RobustOrderEstimator(verbose=False)
    assert est.recalibrate_lambda_per_rep is False


def test_run_pipeline_no_recalibra_lambda():
    """El estudio y la herramienta tienen que ejecutar lo mismo."""
    import inspect
    import sys
    from pathlib import Path
    raiz = Path(__file__).resolve().parent.parent
    if str(raiz / "scripts") not in sys.path:
        sys.path.insert(0, str(raiz / "scripts"))
    import run_pipeline
    fuente = inspect.getsource(run_pipeline)
    assert "recalibrate_lambda_per_rep" in fuente, (
        "run_pipeline deberia fijarlo explicitamente, no depender del default"
    )


# ---------------------------------------------------------------------------
# La correccion de orientacion tiene que alcanzar tambien al atlas
# ---------------------------------------------------------------------------

def test_corregir_la_orientacion_refleja_tambien_los_parametros_fmm():
    """
    ``_fix_orientations`` invierte el orden de muestras y refleja la escala de
    las repeticiones que salieron al reves. Los parametros FMM de esas
    repeticiones estan expresados en la escala de ANTES, asi que hay que
    reflejarlos igual o el atlas queda especular respecto del orden que se
    publica junto a el.

    Bajo theta -> 2*pi - theta la onda FMM se transforma como
    ``(alpha, beta, omega) -> (-alpha, -beta, omega)``, de donde pico y valle
    tambien cambian de signo. M, A y R2 no se tocan: la reflexion no cambia
    ni la escala vertical ni la bondad del ajuste.
    """
    import numpy as np
    from circust.robust_order import RobustOrderEstimator

    n_samp, n_top, K = 12, 3, 2
    est = RobustOrderEstimator(n_reps=K, verbose=False)

    # rep 0 concuerda con el preliminar, rep 1 esta invertida
    theta_pre = np.linspace(0.0, 2 * np.pi, n_samp, endpoint=False)
    sample_orders = np.array([np.arange(n_samp), np.arange(n_samp)])
    circ_scales   = np.array([theta_pre.copy(), (2 * np.pi - theta_pre) % (2 * np.pi)])

    stats = np.zeros((K * n_top, 12))
    for k in range(K):
        for i in range(n_top):
            alpha, beta, omega = 1.0 + 0.3 * i, 2.0 + 0.2 * i, 0.6
            stats[k * n_top + i, 0] = 5.0            # M
            stats[k * n_top + i, 1] = 2.0            # A
            stats[k * n_top + i, 2] = alpha
            stats[k * n_top + i, 3] = beta
            stats[k * n_top + i, 4] = omega
            stats[k * n_top + i, 5] = FMMModel.peak_time(alpha, beta, omega)
            stats[k * n_top + i, 6] = FMMModel.trough_time(alpha, beta, omega)
            stats[k * n_top + i, 11] = 0.9           # R2

    so, cs, st, _ = est._fix_orientations(
        sample_orders, circ_scales, stats.copy(), None,
        theta_pre=theta_pre, n_samp=n_samp,
    )

    dos_pi = 2 * np.pi
    for i in range(n_top):
        # la rep 0 no se toca
        np.testing.assert_allclose(st[i, :], stats[i, :])
        # la rep 1 sale reflejada
        f, o = st[n_top + i], stats[n_top + i]
        assert f[2] == pytest.approx((-o[2]) % dos_pi)      # alpha
        assert f[3] == pytest.approx((-o[3]) % dos_pi)      # beta
        assert f[4] == pytest.approx(o[4])                  # omega intacta
        assert f[5] == pytest.approx((-o[5]) % dos_pi)      # pico
        assert f[6] == pytest.approx((-o[6]) % dos_pi)      # valle
        assert f[0] == pytest.approx(o[0])                  # M intacta
        assert f[1] == pytest.approx(o[1])                  # A intacta
        assert f[11] == pytest.approx(o[11])                # R2 intacto
        # pico y valle siguen siendo los de los parametros reflejados
        assert f[5] == pytest.approx(FMMModel.peak_time(f[2], f[3], f[4]))
        assert f[6] == pytest.approx(FMMModel.trough_time(f[2], f[3], f[4]))


# ---------------------------------------------------------------------------
# Atlas del tejido (build_atlas)
# ---------------------------------------------------------------------------

def test_build_atlas_best_k_devuelve_los_parametros_de_la_rep_elegida():
    """Con best_k no se re-ajusta nada: el atlas es ``selected_fmm`` tal cual."""
    import pandas as pd
    from circust.robust_order import build_atlas

    sel = pd.DataFrame({
        "gene_name": ["G1", "G2"], "r2_fmm": [0.6, 0.9], "M": [0.0, 0.1],
        "amplitude": [1.0, 2.0], "alpha": [1.0, 4.0], "beta": [3.0, 2.0],
        "omega": [0.5, 0.8], "peak_phase_rad": [1.2, 4.5],
        "trough_phase_rad": [4.0, 1.0],
    })
    top = pd.DataFrame(np.zeros((2, 5)), index=["G1", "G2"])
    atlas = build_atlas("best_k", sel, top, np.arange(5), np.linspace(0, 6, 5))

    assert atlas["gene_name"].tolist() == ["G2", "G1"]       # por R2 descendente
    assert atlas["classification"].tolist() == ["night", "day"]
    fila = atlas.set_index("gene_name").loc["G1"]
    assert fila["alpha"] == 1.0 and fila["omega"] == 0.5


def test_build_atlas_aggregate_reajusta_sobre_el_orden_de_consenso():
    """Con aggregate cada gen se re-ajusta sobre ``top_matrix`` reordenada por
    ``sample_order``, igual que un ajuste FMM directo sobre ese orden."""
    import pandas as pd
    from circust.robust_order import build_atlas

    rng   = np.random.default_rng(0)
    n     = 40
    fases = np.sort(rng.uniform(0, TWO_PI, n))
    orden = rng.permutation(n)                  # columna orden[j] tiene fase fases[j]
    datos = np.empty((2, n))
    for g, (a, b, w) in enumerate([(1.0, 3.0, 0.4), (4.0, 2.0, 0.8)]):
        onda = np.cos(b + 2 * np.arctan(w * np.tan((fases - a) / 2)))
        datos[g, orden] = onda + rng.normal(0, 0.05, n)
    top = pd.DataFrame(datos, index=["G1", "G2"])

    atlas = build_atlas("aggregate", pd.DataFrame(), top, orden, fases, n_jobs=1)
    for g in ["G1", "G2"]:
        directo = FMMModel().fit(top.loc[g].values[orden], fases)
        fila    = atlas.set_index("gene_name").loc[g]
        assert fila["r2_fmm"] == pytest.approx(directo.r2)
        assert fila["alpha"] == pytest.approx(directo.params["alpha"])
        assert fila["peak_phase_rad"] == pytest.approx(directo.peak_time)
