# CIRCUST
CIRCular robUST statistical methodology to analyse human molecular rhythms from post-mortem samples.

# How to use 
CIRCUST is achieved in R and is easy to use. 
The code provided in this GitHub replicated the four steps described in CIRCUST paper (https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1011510).

Run the R script named runCIRCURST.R to ,conduct the methodology.
The file matrixIn.RData is loaded on the R script and serves as example of unorderd post-mortem gene expression matrix.
matrixIn: gene expression matrix of size 56200X479 (genes X unordered samples/individuals) at a given tissue. Note that it is a real example, so it requires a reasonable amount of computational time.

INPUTS, see CIRCUST paper for details: 
  - Expression matrix (matrixIn).
  - Number of random selections of the genes at the TOP (K).
  - Name tissue (nameTissue).
  - Core clock gene names (coreG).

0. Preparatory work.
When you use our tool, you should source R source functions and install some R packages detailed. Run in the Rscrip:
  source("functionGTEX_cores.R")

1. Preprocessing.
Run the code line under this name in runCIRCURST.R to clean and normalize the data.

2. Preliminary order.
Run the code lines under this name in runCIRCURST.R to obtain a preliminary order based on the core clok gene. Note that the second element of the list obtained as output from basicOder_cores function provides CIRCUST estimates for the circadian phase time for each sample. In the code we usually refer to these values using the term "esc".

3. TOP Rhythmic orderings.
Run the code lines under this name in runCIRCURST.R to derive the tissue-specific TOP gene list and K circular orders based on K random selections. 

4. Robust Estimation.
Run the code lines under this name in runCIRCURST.R to compute FMM predictions as functions of K circular ordering for the TOP genes. 

OUTPUT:
  - Data frame with the FMM parameter estimations for TOP genes, k=1,...,K (outs).
