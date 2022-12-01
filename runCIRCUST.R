
##############################################################
#		0.-	Preliminary: Load source code			 #
##############################################################

source("functionGTEX_cores.R")

##############################################################
#		INPUTS			     				 #
##############################################################

# Load data example: Tissue-specific gene expression data from GTEx (dowloaded from:)
load(matrixIn.RData) # Provided in this GitHub

# Expression matrix
inputExample<-matrixIn 

# Tissue name used for files
nameTissue<-"Example"

# Name of the 12 core clock genes used in CIRCUST paper
coreG<-c("PER1","PER2","PER3","CRY1","CRY2","ARNTL","CLOCK","NR1D1",
		"RORA","DBP","TEF","STAT3") 

# Number of random selection of the genes at TOP, see CIRCUST paper
K<-5 

##############################################################
#		1.-	Preprocessing        				 #
##############################################################

giveMatIniRefG<-giveMatIniNP_v3_cores(inputExample,nameTissue,coreG)
	
##############################################################
#		2.-	Preliminary Order      				 #
##############################################################

# (2.1)
preOrdRefG2<-basicPreOder_cores(giveMatIniRefG,nameTissue,coreG)
# (2.2) & (2.3)
basicOrdRefG2<-basicOder_cores(nameTissue,preOrdRefG2,coreG)
	
##############################################################
#		3.-	TOP Rhythmic Orderings   	      	 #
##############################################################

# (3.1)
obtainNPRefG<-computeNP(preOrdRefG2[[3]])

# (3.2)
outSetCandRefG<-minSetCand3_v2_cores(obtainNPRefG[[5]],nameTissue,preOrdRefG2[[3]],preOrdRefG2[[2]],tam=50,basicOrdRefG2[[5]])

# (3.3) & (3.4)
refSetRefG<-computeMinSet_v3_cores(outSetCandRefG[[1]],outSetCandRefG[[2]],outSetCandRefG[[3]],outSetCandRefG[[4]],nameTissue,outSetCandRefG[[7]],preOrdRefG2[[1]],preOrdRefG2[[2]],preOrdRefG2[[3]],
	preOrdRefG2[[5]],preOrdRefG2[[4]],obtainNPRefG[[5]],coreG)

# (3.5) & (3.6)
selRefG<-compRandomSel_v3(K,ceiling(2/3*length(refSetRefG[[3]])),refSetRefG[[1]],refSetRefG[[2]],
	refSetRefG[[3]],refSetRefG[[4]],nameTissue,TRUE,refSetRefG[[18]],refSetRefG[[13]],
	refSetRefG[[12]],preOrdRefG2[[3]],refSetRefG[[19]])

# (3.7)
topRefG<-obtainTop_v2_cores(coreG,refSetRefG[[1]],preOrdRefG2[[3]],refSetRefG[[11]],refSetRefG[[12]],refSetRefG[[19]])

##############################################################
#		4.-	Robust Estimation        	           	 #
##############################################################

robEstRefG<-robustEst_v3_cores(K,selRefG[[2]],refSetRefG[[11]],topRefG,nameTissue,
	refSetRefG[[19]],refSetRefG[[13]],refSetRefG[[12]],coreG)
robustSincroRefG<-robustSincroDBP_v3_cores(K,topRefG,robEstRefG,nameTissue,coreG)

##############################################################
#		OUTPUTS		        	           	       #
##############################################################

#data.frame with FMM parameters to compute medians
outs<-data.frame(rep(rownames(topRefG),each=K),rep(1:K,times=nrow(topRefG)),robustSincroRefG[[1]][,c(1:7,12)])
colnames(outs)<-c("TOP","K","M","A","alpha","beta","omega","t_U","t_L","R2")
