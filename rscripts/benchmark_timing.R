# =============================================================================
# rscripts/benchmark_timing.R
# =============================================================================
# Benchmark de tiempos de la implementacion CIRCUST **antigua en R** sobre la
# matriz canonica matrixIn (la misma que el pipeline Python consume desde
# matrixIn.parquet). Cronometra por etapas (proc.time) y repite N veces.
# Es el lado R del Experimento 1 (R vs Python) de la validacion del TFG.
#
# Etapas cronometradas (mismas de runCIRCUST.R):
#   preprocess         giveMatIniNP_v3_cores   (preproceso + CPCA + fits core)
#   preliminary_order  basicPreOder_cores + basicOder_cores
#   top_genes          computeNP + minSetCand3 + computeMinSet + compRandomSel + obtainTop
#   robust_order       robustEst_v3_cores + robustSincroDBP_v3_cores
#
# Dependencias C: la ruta NP usa busquedaMejorC (.Call). Este script compila
# busquedaMejor.c -> busquedaMejor.so si no existe, lo precarga, y neutraliza
# los dyn.load("*.dll") (Windows) del codigo fuente con un stub. El paquete
# 'circular' (solo plotting) se sustituye por stubs, como en validate_export.R.
#
# Uso (desde la RAIZ del repo):
#   Rscript rscripts/benchmark_timing.R                 # full matrixIn, 3 reps
#   Rscript rscripts/benchmark_timing.R --reps 5
#   Rscript rscripts/benchmark_timing.R --reps 1 --max-genes 800   # smoke rapido
#
# Salida:  validation/reference/timing_R.csv  (columnas: rep, stage, seconds)
# =============================================================================

# --------------------------------------------------------------------------
# Argumentos
# --------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
get_opt <- function(flag, default) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) return(default)
  args[i + 1]
}
N_REPS    <- as.integer(get_opt("--reps", "3"))
MAX_GENES <- suppressWarnings(as.integer(get_opt("--max-genes", "NA")))  # NA = todas

REPO_ROOT  <- normalizePath(".")
RSCRIPTS   <- file.path(REPO_ROOT, "rscripts")
RDATA_PATH <- file.path(REPO_ROOT, "data", "matrixIn.RData")
OUT_DIR    <- file.path(REPO_ROOT, "validation", "reference")

# Matriz alternativa en CSV, con los genes como primera columna. Si se pasa,
# sustituye a matrixIn.RData. La usa el bloque computacional del estudio de
# simulacion (simulation/rbench.py), que exporta matriz_r.csv con los genes
# core ya renombrados a los nombres que este script espera en coreG.
MATRIX_CSV <- get_opt("--matrix", NA)
# Directorio donde volcar el resultado estimado, no solo los tiempos. Sirve
# para comparar la precision de R contra la de Python sobre la misma matriz.
DUMP_DIR   <- get_opt("--dump", NA)
OUT_CSV    <- get_opt("--out", file.path(OUT_DIR, "timing_R.csv"))
LABEL      <- get_opt("--label", "matrixIn")
if (!is.na(MATRIX_CSV)) MATRIX_CSV <- normalizePath(MATRIX_CSV)
if (!grepl("^/", OUT_CSV)) OUT_CSV <- file.path(REPO_ROOT, OUT_CSV)
dir.create(dirname(OUT_CSV), showWarnings = FALSE, recursive = TRUE)

coreG <- c("PER1", "PER2", "PER3", "CRY1", "CRY2", "ARNTL", "CLOCK",
           "NR1D1", "RORA", "DBP", "TEF", "STAT3")
K <- 5
nameTissue <- "Example"

# --------------------------------------------------------------------------
# 1) Compilar y precargar la libreria C (busquedaMejor); stub de dyn.load .dll
# --------------------------------------------------------------------------
old_wd <- getwd()
setwd(RSCRIPTS)

if (!file.exists("busquedaMejor.so")) {
  cat("Compilando busquedaMejor.c -> busquedaMejor.so ...\n")
  st <- system("R CMD SHLIB busquedaMejor.c")
  if (st != 0) stop("Fallo al compilar busquedaMejor.c")
}
base::dyn.load("busquedaMejor.so")
cat("  busquedaMejor.so cargada.\n")

# Stub: las fuentes hacen dyn.load("pava.dll")/dyn.load("busquedaMejor.dll").
# En Linux no existen; como ya cargamos el .so, las hacemos no-op.
.real_dynload <- base::dyn.load
dyn.load <- function(x, ...) {
  if (grepl("\\.dll$", x)) return(invisible(NULL))
  .real_dynload(x, ...)
}

# Dispositivo grafico nulo: el pipeline R dibuja plantillas circulares
# (plantillaCircular_cores, etc.). En batch las mandamos a un device nulo para
# no generar Rplots.pdf; el paquete 'circular' (ya instalado) hace el resto.
# Es el modo headless estandar de R, no altera ningun calculo.
pdf(NULL)

# --------------------------------------------------------------------------
# 2) Cargar el codigo fuente R (manteniendo el NP)
# --------------------------------------------------------------------------
cat("Cargando fuentes R ...\n")
source("FMM.R")

# OJO: no usar 'lines' como variable, colisiona con el stub grafico lines().
src_lines <- readLines("functionGTEX_cores.R")
skip_patterns <- c(
  "install.packages"
)
keep <- rep(TRUE, length(src_lines))
for (pat in skip_patterns) keep <- keep & !grepl(pat, src_lines, fixed = TRUE)
eval(parse(text = paste(src_lines[keep], collapse = "\n")), envir = .GlobalEnv)
setwd(old_wd)
cat("  Fuentes cargadas.\n")

# --------------------------------------------------------------------------
# 3) Cargar matrixIn (opcionalmente submuestreado para smoke)
# --------------------------------------------------------------------------
if (!is.na(MATRIX_CSV)) {
  cat(sprintf("Cargando %s ...\n", MATRIX_CSV))
  df       <- read.csv(MATRIX_CSV, row.names = 1, check.names = FALSE)
  matrixIn <- as.matrix(df)
  faltan   <- setdiff(coreG, rownames(matrixIn))
  if (length(faltan) > 0) {
    stop(sprintf("Faltan genes core en la matriz: %s",
                 paste(faltan, collapse = ", ")))
  }
} else {
  cat(sprintf("Cargando %s ...\n", RDATA_PATH))
  load(RDATA_PATH)                     # define matrixIn
  if (!is.na(MAX_GENES) && MAX_GENES < nrow(matrixIn)) {
    set.seed(1)
    core_idx  <- which(rownames(matrixIn) %in% coreG)
    other_idx <- setdiff(seq_len(nrow(matrixIn)), core_idx)
    keep_idx  <- sort(c(core_idx, sample(other_idx, MAX_GENES - length(core_idx))))
    matrixIn  <- matrixIn[keep_idx, ]
    cat(sprintf("  Submuestreo smoke: %d genes (incluye core).\n", nrow(matrixIn)))
  }
}
cat(sprintf("  matrixIn: %d genes x %d muestras\n", nrow(matrixIn), ncol(matrixIn)))

# --------------------------------------------------------------------------
# 4) Una pasada cronometrada del pipeline R completo
# --------------------------------------------------------------------------
run_once <- function(inputExample) {
  tt <- list()
  el <- function(t0) as.numeric((proc.time() - t0)["elapsed"])

  t0 <- proc.time()
  giveMat <- giveMatIniNP_v3_cores(inputExample, nameTissue, coreG)
  tt$preprocess <- el(t0)

  t0 <- proc.time()
  preOrd  <- basicPreOder_cores(giveMat, nameTissue, coreG)
  basicO  <- basicOder_cores(nameTissue, preOrd, coreG)
  tt$preliminary_order <- el(t0)

  t0 <- proc.time()
  obtainNP   <- computeNP(preOrd[[3]])
  outSetCand <- minSetCand3_v2_cores(obtainNP[[5]], nameTissue, preOrd[[3]],
                                     preOrd[[2]], tam = 50, basicO[[5]])
  refSet <- computeMinSet_v3_cores(outSetCand[[1]], outSetCand[[2]], outSetCand[[3]],
                                   outSetCand[[4]], nameTissue, outSetCand[[7]],
                                   preOrd[[1]], preOrd[[2]], preOrd[[3]], preOrd[[5]],
                                   preOrd[[4]], obtainNP[[5]], coreG)
  selRefG <- compRandomSel_v3(K, ceiling(2/3 * length(refSet[[3]])), refSet[[1]],
                              refSet[[2]], refSet[[3]], refSet[[4]], nameTissue, TRUE,
                              refSet[[18]], refSet[[13]], refSet[[12]], preOrd[[3]],
                              refSet[[19]])
  topRefG <- obtainTop_v2_cores(coreG, refSet[[1]], preOrd[[3]], refSet[[11]],
                                refSet[[12]], refSet[[19]])
  tt$top_genes <- el(t0)

  t0 <- proc.time()
  robEst   <- robustEst_v3_cores(K, selRefG[[2]], refSet[[11]], topRefG, nameTissue,
                                 refSet[[19]], refSet[[13]], refSet[[12]], coreG)
  robSincro <- robustSincroDBP_v3_cores(K, topRefG, robEst, nameTissue, coreG)
  tt$robust_order <- el(t0)

  tt$elapsed_total <- tt$preprocess + tt$preliminary_order + tt$top_genes + tt$robust_order

  # Fase estimada por muestra. Es la salida de muestras que reporta CIRCUST:
  # el orden preliminar, porque la etapa robusta resume las K repeticiones con
  # una mediana de parametros FMM y no las colapsa en un orden de muestras.
  #
  # OJO: cuando basicOder_cores invierte la orientacion reconstruye la matriz
  # desde cero y solo le asigna rownames, de modo que los nombres de columna
  # se pierden. En ese caso se recuperan invirtiendo los de la etapa anterior,
  # que es exactamente lo que hace con los datos.
  nombres <- colnames(basicO[[3]])
  if (is.null(nombres)) nombres <- rev(colnames(preOrd[[3]]))
  fases <- data.frame(sample = nombres, phase = as.numeric(basicO[[2]]))

  # Atlas: parametros FMM por gen y repeticion, tal como los publica
  # runCIRCUST.R. La mediana sobre las K se hace al analizar.
  atlas <- data.frame(
    TOP   = rep(rownames(topRefG), each = K),
    K     = rep(1:K, times = nrow(topRefG)),
    robSincro[[1]][, c(1:7, 12)]
  )
  colnames(atlas) <- c("TOP", "K", "M", "A", "alpha", "beta", "omega",
                       "t_U", "t_L", "R2")

  list(tiempos = tt, fases = fases, atlas = atlas)
}

# --------------------------------------------------------------------------
# 5) N repeticiones
# --------------------------------------------------------------------------
# El pipeline R antiguo escribe decenas de PNG/TXT (plantillas, TOP cands) en
# el cwd via png()/write.table. Para no ensuciar el repo, ejecutamos en un
# directorio scratch temporal (OUT_DIR es absoluto, asi que timing_R.csv se
# guarda bien igualmente). No es un truco de computo: solo aisla los efectos
# secundarios de escritura.
scratch <- file.path(tempdir(), "circust_r_bench")
dir.create(scratch, showWarnings = FALSE, recursive = TRUE)
setwd(scratch)

rows <- list()
for (r in seq_len(N_REPS)) {
  cat(sprintf("\n=== Repeticion %d / %d ===\n", r, N_REPS))
  salida <- run_once(matrixIn)
  tt <- salida$tiempos
  for (stage in names(tt)) {
    rows[[length(rows) + 1]] <- data.frame(rep = r, stage = stage,
                                           seconds = tt[[stage]])
  }
  cat(sprintf("  total: %.1f s\n", tt$elapsed_total))

  if (!is.na(DUMP_DIR)) {
    dir.create(DUMP_DIR, showWarnings = FALSE, recursive = TRUE)
    write.csv(salida$fases, file.path(DUMP_DIR, "r_fases.csv"), row.names = FALSE)
    write.csv(salida$atlas, file.path(DUMP_DIR, "r_atlas.csv"), row.names = FALSE)
    cat(sprintf("  volcado en %s\n", DUMP_DIR))
  }
}

setwd(REPO_ROOT)
unlink(scratch, recursive = TRUE)   # limpia los PNG/TXT del scratch

timing_R <- do.call(rbind, rows)
timing_R$implementacion <- "r"
timing_R$celda <- LABEL
out_csv <- OUT_CSV
write.csv(timing_R, out_csv, row.names = FALSE)
cat(sprintf("\nGuardado: %s  (%d filas)\n", out_csv, nrow(timing_R)))
