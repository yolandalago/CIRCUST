# =============================================================================
# rscripts/volcado_circust.R
# =============================================================================
# Ejecuta el pipeline R de CIRCUST (runCIRCUST.R) SIN modificar su logica y
# vuelca a CSV todo lo que calcula internamente, para poder compararlo con la
# version Python "CIRCUST original".
#
# Lo que el codigo R publica de forma natural (runCIRCUST.R) es solo el data
# frame `outs` con K filas por gen TOP. Aqui se exporta ademas lo que el
# codigo R ya tiene calculado pero no escribe: el orden preliminar, los K
# ordenes, las K matrices ordenadas, los K ajustes y la mediana de los
# parametros.
#
# Uso (desde la RAIZ del repo):
#   Rscript rscripts/volcado_circust.R --dump salida_R
#   Rscript rscripts/volcado_circust.R --dump salida_R --seed 1 --full-matrix
#   Rscript rscripts/volcado_circust.R --dump salida_R --matrix datos.csv
#
# Argumentos:
#   --dump DIR        carpeta de salida (obligatorio)
#   --seed N          set.seed(N) antes de ejecutar. El R original NO fija
#                     semilla, asi que sin esto dos ejecuciones difieren.
#   --matrix CSV      matriz alternativa (genes en filas, primera columna =
#                     nombres). Por defecto data/matrixIn.RData.
#   --max-genes N     submuestreo para pruebas rapidas
#   --full-matrix     escribe tambien la matriz COMPLETA ordenada por cada k.
#                     Pesa (todos los genes x muestras x K), por eso es opcional.
#   --K N             numero de ordenes (por defecto 5, como runCIRCUST.R)
#   --parametros-paper  carga el codigo R con los MISMOS parametros que usa
#                     circust_original en Python, sin tocar los ficheros R:
#                       outliers: ceil(0.05 n) candidatas en vez de 8
#                       etapa 3:  R2_ORI >= 0.5, R2_FMM > 0.5 (omega > 0.1),
#                                 se para al reunir `tam` genes por region
#                     Quedan sin igualar la segunda sincronizacion de la Etapa 3
#                     y las condiciones del sorteo de la Etapa 4, que son
#                     estructura del codigo R y no parametros.
#
# Salidas en DIR:
#   r_resumen.txt              genes/muestras, outliers, orientacion, tiempos
#   r_preliminar.csv           muestra, fase (orden preliminar, Etapa 2)
#   r_outliers.csv             muestras eliminadas por la deteccion de outliers
#   r_top_genes.csv            genes TOP, su pico y R2 en el orden preliminar,
#                              y el sector al que pertenecen (Etapa 3)
#   r_orden_k<k>.csv           muestra, posicion, fase del orden k (Etapa 4)
#   r_top_matriz_k<k>.csv      matriz TOP ordenada segun el orden k
#   r_top_ajustes_k<k>.csv     valores ajustados FMM de los genes TOP, orden k
#   r_matriz_completa_k<k>.csv matriz completa ordenada segun k (--full-matrix)
#   r_parametros_K.csv         `outs` de runCIRCUST.R: K filas por gen TOP
#   r_atlas_mediana.csv        mediana de los K parametros por gen TOP
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
get_opt <- function(flag, default) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) return(default)
  args[i + 1]
}
DUMP_DIR  <- get_opt("--dump", NA)
if (is.na(DUMP_DIR)) stop("Falta --dump DIR")
SEED      <- suppressWarnings(as.integer(get_opt("--seed", "NA")))
MATRIX_CSV <- get_opt("--matrix", NA)
MAX_GENES <- suppressWarnings(as.integer(get_opt("--max-genes", "NA")))
FULL      <- "--full-matrix" %in% args
K         <- as.integer(get_opt("--K", "5"))
PAPER     <- "--parametros-paper" %in% args

REPO_ROOT  <- normalizePath(".")
RSCRIPTS   <- file.path(REPO_ROOT, "rscripts")
RDATA_PATH <- file.path(REPO_ROOT, "data", "matrixIn.RData")
if (!grepl("^/", DUMP_DIR)) DUMP_DIR <- file.path(REPO_ROOT, DUMP_DIR)
dir.create(DUMP_DIR, showWarnings = FALSE, recursive = TRUE)
if (!is.na(MATRIX_CSV)) MATRIX_CSV <- normalizePath(MATRIX_CSV)

coreG <- c("PER1", "PER2", "PER3", "CRY1", "CRY2", "ARNTL", "CLOCK",
           "NR1D1", "RORA", "DBP", "TEF", "STAT3")
nameTissue <- "Example"

# --- preparacion identica a benchmark_timing.R (solo plataforma, no logica) --
old_wd <- getwd()
setwd(RSCRIPTS)
if (!file.exists("busquedaMejor.so")) {
  st <- system("R CMD SHLIB busquedaMejor.c")
  if (st != 0) stop("Fallo al compilar busquedaMejor.c")
}
base::dyn.load("busquedaMejor.so")
.real_dynload <- base::dyn.load
dyn.load <- function(x, ...) {
  if (grepl("\\.dll$", x)) return(invisible(NULL))
  .real_dynload(x, ...)
}
pdf(NULL)
source("FMM.R")
src_lines <- readLines("functionGTEX_cores.R")
keep <- !grepl("install.packages", src_lines, fixed = TRUE)
src_txt <- paste(src_lines[keep], collapse = "\n")

# Sustituye `pat` por `rep` en el texto fuente y falla si el numero de
# coincidencias no es el esperado, para que un parche no aplicado no pase
# desapercibido.
parchea <- function(txt, pat, rep, esperadas, fixed = TRUE, nombre = pat) {
  n <- if (fixed) lengths(regmatches(txt, gregexpr(pat, txt, fixed = TRUE)))
       else lengths(regmatches(txt, gregexpr(pat, txt, perl = TRUE)))
  if (n != esperadas) stop(sprintf("Parche '%s': %d coincidencias, se esperaban %d", nombre, n, esperadas))
  if (fixed) gsub(pat, rep, txt, fixed = TRUE) else gsub(pat, rep, txt, perl = TRUE)
}
if (PAPER) {
  cat("Aplicando parametros del paper al codigo R (en memoria) ...\n")
  # outliers: ceil(0.05 n) muestras examinadas en lugar de 8
  src_txt <- parchea(src_txt, "initialTissue<-obtainCPCA13(mTissueCoreGNorm,nameT,coreG,8,FALSE)",
                     "initialTissue<-obtainCPCA13(mTissueCoreGNorm,nameT,coreG,ceiling(0.05*ncol(mTissueCoreGNorm)),FALSE)", 1,
                     nombre = "outliers (1a llamada)")
  src_txt <- parchea(src_txt, "obtainCPCA13(mTissueCoreGNorm,nameT,coreG,8,TRUE)",
                     "obtainCPCA13(mTissueCoreGNorm,nameT,coreG,ceiling(0.05*ncol(mTissueCoreGNorm)),TRUE)", 1,
                     nombre = "outliers (2a llamada)")
  # etapa 3.1: R2_ORI >= 0.5 en lugar de mayor que la mediana
  src_txt <- parchea(src_txt, "which(R2Tissue>R2NP50)", "which(R2Tissue>=0.5)", 1, nombre = "filtro ORI")
  src_txt <- parchea(src_txt, "which(R2NPTissue>median(R2NPTissue))", "which(R2NPTissue>=0.5)", 1,
                     nombre = "filtro ORI (2a pasada)")
  # etapa 3.3: R2_FMM > 0.5 en lugar de > min(0.4, Q1 de los core)
  src_txt <- parchea(src_txt, "ifelse(quantile(R2ParcoreG)[2]>0.4,0.4,quantile(R2ParcoreG)[2])", "0.5", 8,
                     nombre = "R2_FMM > 0.5")
  # etapa 3.3: parar solo al reunir la cuota y mientras R2_ORI >= 0.5
  src_txt <- parchea(src_txt, "\\(in([1-8])<tam \\| \\(in\\1>=tam & r2Reg\\1Ord\\[rec\\1\\]<quantile\\(distR2Tissue\\)\\[4\\]\\) \\)",
                     "(in\\1<tam)", 8, fixed = FALSE, nombre = "parada por cuota")
  src_txt <- parchea(src_txt, "r2Reg([1-8])Ord\\[rec\\1\\]>max\\(quantile\\(r2Reg\\1Ord,probs=0.25\\),0.4\\)",
                     "r2Reg\\1Ord[rec\\1]>=0.5", 8, fixed = FALSE, nombre = "R2_ORI >= 0.5 en el recorrido")
}
eval(parse(text = src_txt), envir = .GlobalEnv)
setwd(old_wd)

# --- datos ------------------------------------------------------------------
if (!is.na(MATRIX_CSV)) {
  # Sin row.names=1: matrixIn tiene nombres de gen duplicados (R se queda con
  # el de mayor MAD en giveMatIniNP_v3_cores) y read.csv no los admite.
  df <- read.csv(MATRIX_CSV, check.names = FALSE)
  matrixIn <- as.matrix(df[, -1])
  rownames(matrixIn) <- df[[1]]
} else {
  load(RDATA_PATH)
}
if (!is.na(MAX_GENES) && MAX_GENES < nrow(matrixIn)) {
  set.seed(1)
  core_idx  <- which(rownames(matrixIn) %in% coreG)
  other_idx <- setdiff(seq_len(nrow(matrixIn)), core_idx)
  matrixIn  <- matrixIn[sort(c(core_idx, sample(other_idx, MAX_GENES - length(core_idx)))), ]
}
sample_names <- colnames(matrixIn)
cat(sprintf("matrixIn: %d genes x %d muestras\n", nrow(matrixIn), ncol(matrixIn)))
if (!is.na(SEED)) set.seed(SEED)

# --- pipeline R, llamadas identicas a runCIRCUST.R ---------------------------
scratch <- file.path(tempdir(), "circust_r_volcado")
dir.create(scratch, showWarnings = FALSE, recursive = TRUE)
setwd(scratch)          # el codigo R escribe PNG/TXT en el cwd
tt <- list(); el <- function(t0) as.numeric((proc.time() - t0)["elapsed"])

t0 <- proc.time()
giveMat <- giveMatIniNP_v3_cores(matrixIn, nameTissue, coreG)
tt$preprocess <- el(t0)

t0 <- proc.time()
preOrd <- basicPreOder_cores(giveMat, nameTissue, coreG)
basicO <- basicOder_cores(nameTissue, preOrd, coreG)
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
robEst    <- robustEst_v3_cores(K, selRefG[[2]], refSet[[11]], topRefG, nameTissue,
                                refSet[[19]], refSet[[13]], refSet[[12]], coreG)
robSincro <- robustSincroDBP_v3_cores(K, topRefG, robEst, nameTissue, coreG)
tt$robust_order <- el(t0)
setwd(old_wd)

# --- volcado -----------------------------------------------------------------
w <- function(x, f, ...) write.csv(x, file.path(DUMP_DIR, f), row.names = FALSE, ...)

# Orden preliminar. Cuando basicOder_cores invierte la orientacion pierde los
# nombres de columna, asi que se reconstruyen desde el indice de columnas.
# preOrd[[1]] ordena las columnas de giveMat[[3]] (matriz completa normalizada
# tras quitar outliers) y basicO[[1]] las de preOrd[[3]].
drop_idx   <- giveMat[[8]]
nm_clean   <- if (length(drop_idx) > 0) sample_names[-drop_idx] else sample_names
nm_prelim1 <- nm_clean[giveMat[[1]]]            # tras CPCA preliminar
nm_prelim  <- nm_prelim1[preOrd[[1]]]           # tras rotacion/direccion 2.1
nm_prelim  <- nm_prelim[basicO[[8]]]            # tras verificacion 2.2
w(data.frame(sample = nm_prelim, phase = as.numeric(basicO[[2]])), "r_preliminar.csv")
w(data.frame(sample = if (length(drop_idx) > 0) sample_names[drop_idx] else character(0)),
  "r_outliers.csv")

# Genes TOP con su pico en el orden preliminar y su sector (Etapa 3)
cand_names <- refSet[[1]]
w(data.frame(gene = cand_names, peak_prelim = refSet[[2]], r2_fmm = refSet[[3]],
             sector = refSet[[4]],
             en_top = cand_names %in% rownames(topRefG)), "r_candidatos_etapa3.csv")
w(data.frame(gene = rownames(topRefG),
             es_core = rownames(topRefG) %in% coreG), "r_top_genes.csv")

# Columnas de `top` (orden sincronizado de la Etapa 3) respecto a nm_prelim
col_top <- nm_prelim[refSet[[12]]]
stopifnot(ncol(topRefG) == length(col_top))

for (k in seq_len(K)) {
  perm_k   <- robEst[[3]][[k]][[1]][ robEst[[2]][[k]][[1]] ]   # indices en columnas de top
  nm_k     <- col_top[perm_k]
  escala_k <- robSincro[[6]][k, ]
  # comprobacion interna: reordenar top por perm_k reproduce la matriz de R
  chk <- max(abs(topRefG[, perm_k] - robSincro[[2]][[k]]))
  if (chk > 1e-10) warning(sprintf("k=%d: la permutacion reconstruida no reproduce mF (dif=%g)", k, chk))
  w(data.frame(sample = nm_k, position = seq_along(nm_k), circular_phase = as.numeric(escala_k)),
    sprintf("r_orden_k%d.csv", k))
  mk <- robSincro[[2]][[k]]; rownames(mk) <- rownames(topRefG); colnames(mk) <- nm_k
  write.csv(mk, file.path(DUMP_DIR, sprintf("r_top_matriz_k%d.csv", k)))
  ak <- robSincro[[3]][[k]]; rownames(ak) <- rownames(topRefG); colnames(ak) <- nm_k
  write.csv(ak, file.path(DUMP_DIR, sprintf("r_top_ajustes_k%d.csv", k)))
  if (FULL) {
    # matriz completa en el orden preliminar sincronizado, reordenada por k
    full <- preOrd[[3]][, refSet[[12]]][, perm_k]
    colnames(full) <- nm_k
    write.table(full, file.path(DUMP_DIR, sprintf("r_matriz_completa_k%d.csv", k)),
                sep = ",", quote = FALSE, col.names = NA)
  }
}

# `outs` de runCIRCUST.R: K filas por gen TOP
outs <- data.frame(rep(rownames(topRefG), each = K), rep(1:K, times = nrow(topRefG)),
                   robSincro[[1]][, c(1:7, 12)])
colnames(outs) <- c("TOP", "K", "M", "A", "alpha", "beta", "omega", "t_U", "t_L", "R2")
w(outs, "r_parametros_K.csv")

# Mediana de los K parametros por gen. Los angulos (alpha, beta, t_U, t_L) se
# resumen con la mediana CIRCULAR, porque la ordinaria se rompe en el corte
# 0/2pi. R no hace este paso: runCIRCUST.R lo deja al usuario.
med_circ <- function(x) as.numeric(median.circular(circular(x)) %% (2 * pi))
genes <- unique(outs$TOP)
atlas <- do.call(rbind, lapply(genes, function(g) {
  d <- outs[outs$TOP == g, ]
  data.frame(TOP = g, M = median(d$M), A = median(d$A),
             alpha = med_circ(d$alpha), beta = med_circ(d$beta),
             omega = median(d$omega), t_U = med_circ(d$t_U), t_L = med_circ(d$t_L),
             R2 = median(d$R2))
}))
w(atlas, "r_atlas_mediana.csv")

writeLines(c(
  sprintf("genes_entrada=%d", nrow(matrixIn)), sprintf("muestras_entrada=%d", ncol(matrixIn)),
  sprintf("outliers=%d", length(drop_idx)), sprintf("genes_TOP=%d", nrow(topRefG)),
  sprintf("orientacion_cambiada_prelim=%s", as.character(basicO[[7]])),
  sprintf("K=%d", K), sprintf("seed=%s", ifelse(is.na(SEED), "ninguna", SEED)),
  sprintf("t_%s=%.1f", names(tt), unlist(tt))), file.path(DUMP_DIR, "r_resumen.txt"))
cat(sprintf("Volcado en %s\n", DUMP_DIR))
