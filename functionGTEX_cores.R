require("Iso")
require("CircStats")
require("circular")
require("FMM.R")
source("upDownUp_NP_Code_modif.R")
source("upDownUp_NP_Code_NoParalelizado.R")
source("upDownUp_NP_Code_Paralelizado.R")
source("NucleoComun.R)

#datos<-mFullNormRefG[[numTissue]]

computeNP<-function(datos){
	require('Iso')
	fitNP<-matrix(0,nrow(datos),ncol(datos))
	flat<-matrix(0,nrow(datos),ncol(datos))
	msseNP<-c()
	msseFlat<-c()
	R2<-c()
	for(i in 1:nrow(datos)){
		print(paste0(i,"/",nrow(datos)))
  		fitNP[i,]<-function1Local_modif(datos[i,], parallel = "auto")[[1]]#function1Local(datos[i,])[[1]]
  		flat[i,]<-rep(mean(datos[i,]),ncol(datos))
  		msseNP[i]<-sum((datos[i,]-fitNP[i,])^2)/ncol(datos)
  		msseFlat[i]<-sum((datos[i,]-flat[i,])^2)/ncol(datos)
  		R2[i]<-1-msseNP[i]/msseFlat[i]
  
	}
	return(list(fitNP,flat,msseNP,msseFlat,R2))
}




normalice2<-function(z){
  a <- min(z)
  b <- max(z)
  if(a==b){
    rr<-rep(0,length(z))
  }else{
    rr<-2*(((z-a)/(b-a))-0.5)
  }
  
  return (rr)
}

normalice<-function(z){
        a <- min(z)
        b <- max(z)
        return (2*(((z-a)/(b-a))-0.5))
}
secTimes=function(peri){
  newTime<-seq(0,2*pi,by=2*pi/peri)
  newTime<-newTime[-length(newTime)]
  return(newTime)
}
#dat<-mNorm8
centrado<-function(dat){
dat2<-matrix(0,nrow(dat),ncol(dat))
for(i in 1:nrow(dat2)){
dat2[i,]<-(dat[i,]-mean(dat[i,],na.rm=TRUE))
}
return(dat2)
}


outF2<-function(dat,param,time){
  return(param[1] + param[2] * cos(param[4]+2*atan2(param[5]*sin((time-param[3])/2),cos((time-param[3])/2))))
}

outC<-function(param,time){
return(param[1]+param[2]*cos(time+param[3]))
}

compLL<-function(al,be,om){
return((((al+2*atan2(1/om*sin((pi-be)/2),cos((pi-be)/2))))%%(2*pi)))
}

compUU<-function(al,be,om){
return((((al+2*atan2(1/om*sin(-be/2),cos(-be/2))))%%(2*pi)))
}

funcionCosinor<-function(datos,time,periodo){#time entre 0y2pi divididos por el periodo

  xx<-cos(time)
  zz<-sin(time)

  fit<-lm((datos)~xx+zz)
  Mest<-fit$coefficients[1]
  bb<-fit$coefficients[2]
  gg<-fit$coefficients[3]
  phiEst<-atan2(-gg,bb)%%(2*pi)
  Aest<-sqrt(bb^2+gg^2)

  return(list(Mest+Aest*cos(time+phiEst),Mest,Aest,(time+phiEst)%%(2*pi),phiEst))
}

#coreG=coreG
#mName=namesRefRefG[[numTissue]]
#mFullNorm=mFullNormRefG[[numTissue]]
#mMinCandNormOrd=mFullSincroRefG[[numTissue]]
#indOrd=indSincroRefG[[numTissue]]
#genAdd=addArtGene[[numTissue]]

#coreG=coreG2
#mName="Epidermis"
#mFullNorm=mFullNormRefGSkin
#mMinCandNormOrd=mFullSincroRefGSkin
#indOrd=indSincroRefGSkin
#genAdd=addArtGeneSkin


obtainTop_v2_cores<-function(coreG,mName,mFullNorm,mMinCandNormOrd,indOrd,genAdd){
	top<-c()
	if(length(genAdd)>0){
		if(!is.na(match(genAdd[1],mName)))mName=mName[-match(genAdd[1],mName)]
		if(!is.na(match(genAdd[2],mName)))mName=mName[-match(genAdd[2],mName)]
	}
		if(sum(is.na(match(coreG,mName)))==length(coreG)){
			coreGenes<-mFullNorm[match(coreG,rownames(mFullNorm)),indOrd]
			top<-rbind(coreGenes,mMinCandNormOrd[mName,])
			rownames(top)<-c(coreG,mName)
			
			
		}else{
			mm<-match(coreG,mName)[!is.na(match(coreG,mName))]
			coreGenes<-mFullNorm[match(coreG,rownames(mFullNorm)),indOrd]
			top<-rbind(coreGenes,mMinCandNormOrd[match((mName)[-mm],rownames(mMinCandNormOrd)),])
			rownames(top)<-c(coreG,(mName)[-mm])
		}
		png(filename="coreGenesSynchroInitialOrder.png",width=800,height=400)
		par(mfrow=c(compPerfSq(1,length(coreG)),compPerfSq(1,length(coreG))))
		for(i in 1:length(coreG)){
			plot(indOrd,coreGenes[i,],type="b",main=paste(coreG[i]," synchro",sep=""))
		}
		dev.off()
	return(top)
}

reduceSectors_v2<-function(peakS,R2S,namesS,nSamples){
R2PhaseS75<-quantile(R2S,probs=0.75)#[4]

condR2PhaseS75<-length(peakS)>25 | length(R2S[which(R2S>R2PhaseS75)])>25
n5<-ceiling(nSamples)
nS<-length(peakS)
dropS<-c()
while(condR2PhaseS75){
	dropi<-which.min(R2S)
	peakS<-peakS[-dropi]
	R2S<-R2S[-dropi]
	namesS<-namesS[-dropi]
	dropS<-c(dropS,dropi)
	if(min(R2S)<R2PhaseS75){
		nS<-nS-1
		condR2PhaseS75<-TRUE
	}else{
		condR2PhaseS75<-FALSE
	}
}
return(list(peakS,R2S,namesS))
}

#nRep=K
#tam=ceiling(2/3*length(r2RefRefGSkin))
#namesRedMuscle=namesRefRefGSkin
#peakParRedMuscle=peaksRefRefGSkin
#R2RedMuscle=r2RefRefGSkin
#sectorBelongingMuscle=belongRefRefGSkin
#nameTissue="Epidermis"
#plotting=TRUE
#peakCosSin=peaksCosRefRefGSkin
#escSin=escSincroRefGSkin
#indSin=indSincroRefGSkin
#mFull=mFullNormRefGSkin
#genesAdd=addArtGeneSkin

compRandomSel_v2<-function(nRep,tam,namesRedMuscle,peakParRedMuscle,R2RedMuscle,sectorBelongingMuscle,nameTissue,plotting,peakCosSin,escSin,indSin,mFull,genesAdd){
	k<-0
	dosSect=FALSE
	mSel<-matrix(0,nRep,tam)
	mName<-matrix(0,nRep,tam)
	mPeak<-matrix(0,nRep,tam)
	mCosPeak<-matrix(0,nRep,tam)
	mR2<-matrix(0,nRep,tam)
	sectorBelongingMuscle2=c()
	for(i in 1:length(sectorBelongingMuscle)){
		if(sectorBelongingMuscle[i]==1 | sectorBelongingMuscle[i]==2)sectorBelongingMuscle2[i]=1
		if(sectorBelongingMuscle[i]==3 | sectorBelongingMuscle[i]==4)sectorBelongingMuscle2[i]=2
		if(sectorBelongingMuscle[i]==5 | sectorBelongingMuscle[i]==6)sectorBelongingMuscle2[i]=3
		if(sectorBelongingMuscle[i]==7 | sectorBelongingMuscle[i]==8)sectorBelongingMuscle2[i]=4
	}
	veces=1
contando=0
storeSel=c()
storeMaxDistEquis8=c()
	while(k<nRep){
		print(paste("m.a.s. ",veces,sep=""))
		if(length(genesAdd)>0){
			if(!is.na(match(genesAdd[1],namesRedMuscle))){
				R2RedMuscle=R2RedMuscle[-match(genesAdd[1],namesRedMuscle)]
				peakParRedMuscle=peakParRedMuscle[-match(genesAdd[1],namesRedMuscle)]
				sectorBelongingMuscle2=sectorBelongingMuscle2[-match(genesAdd[1],namesRedMuscle)]
				peakCosSin=peakCosSin[-match(genesAdd[1],namesRedMuscle)]
				namesRedMuscle=namesRedMuscle[-match(genesAdd[1],namesRedMuscle)]
			}
			if(!is.na(match(genesAdd[2],namesRedMuscle))){
				R2RedMuscle=R2RedMuscle[-match(genesAdd[2],namesRedMuscle)]
				peakParRedMuscle=peakParRedMuscle[-match(genesAdd[2],namesRedMuscle)]
				sectorBelongingMuscle2=sectorBelongingMuscle2[-match(genesAdd[2],namesRedMuscle)]
				peakCosSin=peakCosSin[-match(genesAdd[2],namesRedMuscle)]
				namesRedMuscle=namesRedMuscle[-match(genesAdd[2],namesRedMuscle)]
			}
		}
		if(length(unique(sectorBelongingMuscle2))<=2)dosSect=TRUE
		sel<-sample(1:length(namesRedMuscle),tam)
		highR2sel<-min(R2RedMuscle[sel])>0.5#tampoco va a hacer nada, ya lo hemos asegurado esto
		distCond=TRUE
		if(highR2sel==FALSE & (namesRedMuscle[sel][which.min(R2RedMuscle[sel])]=="ARNTL" |
			namesRedMuscle[sel][which.min(R2RedMuscle[sel])]=="DBP"))highR2sel=TRUE
		atLeast2Sectors<-length(unique(sectorBelongingMuscle2[sel]))>2
		atLeast1Sectors<-length(unique(sectorBelongingMuscle2[sel]))>1
		#CONDICION DISTANCIA ENTRE LOSPEAKS DE COSINOR
		distPeakCos=c()
		peakCosSinSelSort=sort(peakParRedMuscle[sel])
		#(sectorBelongingMuscle2[sel])[order(peakParRedMuscle[sel])]
		for(i in 2:length(sel)){
			distPeakCos[i-1]=1-cos(peakCosSinSelSort[i-1]-peakCosSinSelSort[i])
		}
		
		distPeakCos[length(sel)]=1-cos(peakCosSinSelSort[length(sel)]-peakCosSinSelSort[1])
		distS1=distPeakCos[which((sectorBelongingMuscle2[sel])[order(peakParRedMuscle[sel])]==1)]
		distS2=distPeakCos[which((sectorBelongingMuscle2[sel])[order(peakParRedMuscle[sel])]==2)]
		distS3=distPeakCos[which((sectorBelongingMuscle2[sel])[order(peakParRedMuscle[sel])]==3)]
		distS4=distPeakCos[which((sectorBelongingMuscle2[sel])[order(peakParRedMuscle[sel])]==4)]
		condS1=TRUE;condS2=TRUE;condS3=TRUE;condS4=TRUE
		if(length(distS1)>0)condS1=min(distS1)>0.00001 | (min(distS1)<0.00001 & sum(distS1<0.00001)<=ceiling(0.1*length(distS1)))
		if(length(distS2)>0)condS2=min(distS2)>0.00001 | (min(distS2)<0.00001 & sum(distS2<0.00001)<=ceiling(0.1*length(distS2)))
		if(length(distS3)>0)condS3=min(distS3)>0.00001 | (min(distS3)<0.00001 & sum(distS3<0.00001)<=ceiling(0.1*length(distS3)))
		if(length(distS4)>0)condS4=min(distS4)>0.00001 | (min(distS4)<0.00001 & sum(distS4<0.00001)<=ceiling(0.1*length(distS4)))
		#condicion para que la distancia max entre muestras de la rep, 
		#no sea mayor la dist max que hay en el orden inicail entre las muestras cuando esta sea alta (mas de 2 mean+sd)
		timEuc=escSin/(2*pi)*length(escSin)
		distEquis=c()
		for(i in 2:length(escSin)){
			distEquis[i-1]=abs(timEuc[i]-timEuc[i-1])
		}
		distEquis[length(escSin)]=length(escSin)-timEuc[length(escSin)]+timEuc[1]
		if(max(distEquis)>1.5*IQR(distEquis)){#(mean(distEquis)+2*sd(distEquis))){
				contando=contando+1
				if(length(genesAdd)>0){
				mAux2=mFull[match(namesRedMuscle[sel],rownames(mFull)),indSin]
				if(!is.na(match(genesAdd[1],rownames(mAux2)))){
					mAux3=mAux2[-match(genesAdd[1],rownames(mAux2)),]
				}else{
					mAux3=mAux2
				}
				if(!is.na(match(genesAdd[2],rownames(mAux3)))){
					mAux4=mAux3[-match(genesAdd[2],rownames(mAux3)),]
				}else{
					mAux4=mAux3
				}
				mAux=mAux4
			}else{
				mAux=mFull[match(namesRedMuscle[sel],rownames(mFull)),indSin]
			}
			cp8<-prcomp(centrado(mAux), scale. = TRUE, center = FALSE)
			varPer8<-c(cp8[[1]][1]^2/sum(cp8[[1]]^2),
			cp8[[1]][2]^2/sum(cp8[[1]]^2),
			cp8[[1]][3]^2/sum(cp8[[1]]^2))#varianza explicada
			eigen18<-cp8$rotation[,1]
			eigen28<-cp8$rotation[,2]
			eigen38<-cp8$rotation[,3]
			xi8<-eigen18/sqrt((eigen18^2+ eigen28^2))
			yi8<-eigen28/sqrt((eigen18^2+ eigen28^2))
			phi8<-(atan2(yi8,xi8))%%(2*pi)
			orderCPCA8<-order(phi8)
			escalaPhi8<-phi8[orderCPCA8]
			timEuc8=escalaPhi8/(2*pi)*length(escalaPhi8)
			distEquis8=c()
			for(i in 2:length(escalaPhi8)){
				distEquis8[i-1]=abs(timEuc8[i]-timEuc8[i-1])
			}
			distEquis8[length(escalaPhi8)]=length(escalaPhi8)-timEuc8[length(escalaPhi8)]+timEuc8[1]
			if(max(distEquis8)>max(distEquis))distCond=FALSE
			if(distCond==FALSE)storeMaxDistEquis8=c(storeMaxDistEquis8,max(distEquis8))
			if(distCond==FALSE)storeSel=rbind(storeSel,sel)
		}
#x11()
#plot(escalaPhi8,mAux[match(namesRedMuscle[3],rownames(mAux)),orderCPCA8],type="b",main="despues STAT2")
		if(((atLeast2Sectors & !dosSect) | (atLeast1Sectors & dosSect)) & highR2sel & 
			condS1 & condS2 & condS3 & condS4 & distCond){
			k<-k+1
			print(paste("Seleccion ",k," de ",veces," m.a.s.",sep=""))
			mSel[k,]<-sel
			mName[k,]<-namesRedMuscle[sel]
			mPeak[k,]<-peakParRedMuscle[sel]
			mCosPeak[k,]<-peakCosSin[sel]
			mR2[k,]<-R2RedMuscle[sel]

			if(plotting){
				
				#x11()
				png(file=paste("randomSel_Rep ",k," from ",nameTissue,"CosRefG.RData",sep=""))
				par(mfrow=c(1,1))
				par(mar=c(1,1,1,1))
				plot(as.circular(mPeak[k,]),main=paste("Rep ",k," from Cos ",nameTissue,sep=""))
				axis.circular(at=circular(mPeak[k,]), 
					labels=mName[k,],tcl.text=-0.07,tick=FALSE)
				dev.off()
			}

		}
		if(veces>100000){
			if(k==0){
				mSel<-storeSel[(order(storeMaxDistEquis8))[1:nRep],]
				for(kkk in 1:nRep){
					mName[kkk,]<-namesRedMuscle[mSel[kkk,]]
					mPeak[kkk,]<-peakParRedMuscle[mSel[kkk,]]
					mCosPeak[kkk,]<-peakCosSin[mSel[kkk,]]
					mR2[kkk,]<-R2RedMuscle[mSel[kkk,]]
				}
				if(plotting){
				
					for(kkk in 1:nRep){
						png(file=paste("randomSel_Rep ",kkk," from ",nameTissue,"CosRefG.RData",sep=""))
						par(mfrow=c(1,1))
						par(mar=c(1,1,1,1))
						plot(as.circular(mPeak[kkk,]),main=paste("Rep ",kkk," from Cos ",nameTissue,sep=""))
						axis.circular(at=circular(mPeak[kkk,]), 
							labels=mName[kkk,],tcl.text=-0.07,tick=FALSE)
						dev.off()
					}
				}
			}else{
				mSel<-storeSel[(order(storeMaxDistEquis8))[1:(nRep-k)],]
				for(kkk in (k+1):nRep){
					mName[kkk,]<-namesRedMuscle[mSel[kkk-k,]]
					mPeak[kkk,]<-peakParRedMuscle[mSel[kkk-k,]]
					mCosPeak[kkk,]<-peakCosSin[mSel[kkk-k,]]
					mR2[kkk,]<-R2RedMuscle[mSel[kkk-k,]]
				}
				if(plotting){
				
					for(kkk in (k+1):nRep){
						png(file=paste("randomSel_Rep ",kkk," from ",nameTissue,"CosRefG.RData",sep=""))
						par(mfrow=c(1,1))
						par(mar=c(1,1,1,1))
						plot(as.circular(mPeak[kkk-k,]),main=paste("Rep ",kkk," from Cos ",nameTissue,sep=""))
						axis.circular(at=circular(mPeak[kkk-k,]), 
							labels=mName[kkk-k,],tcl.text=-0.07,tick=FALSE)
						dev.off()
					}
				}
			}
			k=99999
			
		}
	veces=veces+1
	}
	#selRefG<-list(mSel,mName,mPeak,mR2,mCosPeak,dosSect,distPeakCos,distEquis,distEquis8,distCond)
	return(list(mSel,mName,mPeak,mR2,mCosPeak,dosSect,distPeakCos,distEquis,distEquis8,distCond))
}

compRandomSel_v3<-function(nRep,tam,namesRedMuscle,peakParRedMuscle,R2RedMuscle,sectorBelongingMuscle,nameTissue,plotting,peakCosSin,escSin,indSin,mFull,genesAdd){
  k<-0
  dosSect=FALSE
  mSel<-matrix(0,nRep,tam)
  mName<-matrix(0,nRep,tam)
  mPeak<-matrix(0,nRep,tam)
  mCosPeak<-matrix(0,nRep,tam)
  mR2<-matrix(0,nRep,tam)
  sectorBelongingMuscle2=c()
  for(i in 1:length(sectorBelongingMuscle)){
    if(sectorBelongingMuscle[i]==1 | sectorBelongingMuscle[i]==2)sectorBelongingMuscle2[i]=1
    if(sectorBelongingMuscle[i]==3 | sectorBelongingMuscle[i]==4)sectorBelongingMuscle2[i]=2
    if(sectorBelongingMuscle[i]==5 | sectorBelongingMuscle[i]==6)sectorBelongingMuscle2[i]=3
    if(sectorBelongingMuscle[i]==7 | sectorBelongingMuscle[i]==8)sectorBelongingMuscle2[i]=4
  }
  veces=1
  contando=0
  storeSel=c()
  storeMaxDistEquis8=c()
  while(k<nRep){
    print(paste("m.a.s. ",veces,sep=""))
    if(length(genesAdd)>0){
      if(!is.na(match(genesAdd[1],namesRedMuscle))){
        R2RedMuscle=R2RedMuscle[-match(genesAdd[1],namesRedMuscle)]
        peakParRedMuscle=peakParRedMuscle[-match(genesAdd[1],namesRedMuscle)]
        sectorBelongingMuscle2=sectorBelongingMuscle2[-match(genesAdd[1],namesRedMuscle)]
        peakCosSin=peakCosSin[-match(genesAdd[1],namesRedMuscle)]
        namesRedMuscle=namesRedMuscle[-match(genesAdd[1],namesRedMuscle)]
      }
      if(!is.na(match(genesAdd[2],namesRedMuscle))){
        R2RedMuscle=R2RedMuscle[-match(genesAdd[2],namesRedMuscle)]
        peakParRedMuscle=peakParRedMuscle[-match(genesAdd[2],namesRedMuscle)]
        sectorBelongingMuscle2=sectorBelongingMuscle2[-match(genesAdd[2],namesRedMuscle)]
        peakCosSin=peakCosSin[-match(genesAdd[2],namesRedMuscle)]
        namesRedMuscle=namesRedMuscle[-match(genesAdd[2],namesRedMuscle)]
      }
    }
    if(length(unique(sectorBelongingMuscle2))<=2)dosSect=TRUE
    sel<-sample(1:length(namesRedMuscle),tam)
    highR2sel<-min(R2RedMuscle[sel])>0.5#tampoco va a hacer nada, ya lo hemos asegurado esto
    distCond=TRUE
    if(highR2sel==FALSE & (namesRedMuscle[sel][which.min(R2RedMuscle[sel])]=="ARNTL" |
                           namesRedMuscle[sel][which.min(R2RedMuscle[sel])]=="DBP"))highR2sel=TRUE
    atLeast2Sectors<-length(unique(sectorBelongingMuscle2[sel]))>2
    atLeast1Sectors<-length(unique(sectorBelongingMuscle2[sel]))>1
    #CONDICION DISTANCIA ENTRE LOSPEAKS DE COSINOR
    distPeakCos=c()
    peakCosSinSelSort=sort(peakParRedMuscle[sel])
    #(sectorBelongingMuscle2[sel])[order(peakParRedMuscle[sel])]
    for(i in 2:length(sel)){
      distPeakCos[i-1]=1-cos(peakCosSinSelSort[i-1]-peakCosSinSelSort[i])
    }
    
    distPeakCos[length(sel)]=1-cos(peakCosSinSelSort[length(sel)]-peakCosSinSelSort[1])
    distS1=distPeakCos[which((sectorBelongingMuscle2[sel])[order(peakParRedMuscle[sel])]==1)]
    distS2=distPeakCos[which((sectorBelongingMuscle2[sel])[order(peakParRedMuscle[sel])]==2)]
    distS3=distPeakCos[which((sectorBelongingMuscle2[sel])[order(peakParRedMuscle[sel])]==3)]
    distS4=distPeakCos[which((sectorBelongingMuscle2[sel])[order(peakParRedMuscle[sel])]==4)]
    condS1=TRUE;condS2=TRUE;condS3=TRUE;condS4=TRUE
    if(length(distS1)>0)condS1=min(distS1)>0.00001 | (min(distS1)<0.00001 & sum(distS1<0.00001)<=ceiling(0.1*length(distS1)))
    if(length(distS2)>0)condS2=min(distS2)>0.00001 | (min(distS2)<0.00001 & sum(distS2<0.00001)<=ceiling(0.1*length(distS2)))
    if(length(distS3)>0)condS3=min(distS3)>0.00001 | (min(distS3)<0.00001 & sum(distS3<0.00001)<=ceiling(0.1*length(distS3)))
    if(length(distS4)>0)condS4=min(distS4)>0.00001 | (min(distS4)<0.00001 & sum(distS4<0.00001)<=ceiling(0.1*length(distS4)))
    #condicion para que la distancia max entre muestras de la rep, 
    #no sea mayor la dist max que hay en el orden inicail entre las muestras cuando esta sea alta (mas de 2 mean+sd)
    timEuc=escSin/(2*pi)*length(escSin)
    distEquis=c()
    for(i in 2:length(escSin)){
      distEquis[i-1]=abs(timEuc[i]-timEuc[i-1])
    }
    distEquis[length(escSin)]=length(escSin)-timEuc[length(escSin)]+timEuc[1]
    if(max(distEquis)>1.5*IQR(distEquis)){#(mean(distEquis)+2*sd(distEquis))){
      contando=contando+1
      if(length(genesAdd)>0){
        mAux2=mFull[match(namesRedMuscle[sel],rownames(mFull)),indSin]
        if(!is.na(match(genesAdd[1],rownames(mAux2)))){
          mAux3=mAux2[-match(genesAdd[1],rownames(mAux2)),]
        }else{
          mAux3=mAux2
        }
        if(!is.na(match(genesAdd[2],rownames(mAux3)))){
          mAux4=mAux3[-match(genesAdd[2],rownames(mAux3)),]
        }else{
          mAux4=mAux3
        }
        mAux=mAux4
      }else{
        mAux=mFull[match(namesRedMuscle[sel],rownames(mFull)),indSin]
      }
      cp8<-prcomp(centrado(mAux), scale. = TRUE, center = FALSE)
      varPer8<-c(cp8[[1]][1]^2/sum(cp8[[1]]^2),
                 cp8[[1]][2]^2/sum(cp8[[1]]^2),
                 cp8[[1]][3]^2/sum(cp8[[1]]^2))#varianza explicada
      eigen18<-cp8$rotation[,1]
      eigen28<-cp8$rotation[,2]
      eigen38<-cp8$rotation[,3]
      xi8<-eigen18/sqrt((eigen18^2+ eigen28^2))
      yi8<-eigen28/sqrt((eigen18^2+ eigen28^2))
      phi8<-(atan2(yi8,xi8))%%(2*pi)
      orderCPCA8<-order(phi8)
      escalaPhi8<-phi8[orderCPCA8]
      timEuc8=escalaPhi8/(2*pi)*length(escalaPhi8)
      distEquis8=c()
      for(i in 2:length(escalaPhi8)){
        distEquis8[i-1]=abs(timEuc8[i]-timEuc8[i-1])
      }
      distEquis8[length(escalaPhi8)]=length(escalaPhi8)-timEuc8[length(escalaPhi8)]+timEuc8[1]
      if(max(distEquis8)>max(distEquis))distCond=FALSE
      if(distCond==FALSE)storeMaxDistEquis8=c(storeMaxDistEquis8,max(distEquis8))
      if(distCond==FALSE)storeSel=rbind(storeSel,sel)
    }
    #x11()
    #plot(escalaPhi8,mAux[match(namesRedMuscle[3],rownames(mAux)),orderCPCA8],type="b",main="despues STAT2")
    if(((atLeast2Sectors & !dosSect) | (atLeast1Sectors & dosSect)) & highR2sel & 
       condS1 & condS2 & condS3 & condS4 & distCond){
      k<-k+1
      print(paste("Seleccion ",k," de ",veces," m.a.s.",sep=""))
      mSel[k,]<-sel
      mName[k,]<-namesRedMuscle[sel]
      mPeak[k,]<-peakParRedMuscle[sel]
      mCosPeak[k,]<-peakCosSin[sel]
      mR2[k,]<-R2RedMuscle[sel]
      
      if(plotting){
        
        #x11()
        png(file=paste("randomSel_Rep ",k," from ",nameTissue,"CosRefG.RData",sep=""))
        par(mfrow=c(1,1))
        par(mar=c(1,1,1,1))
        plot(as.circular(mPeak[k,]),main=paste("Rep ",k," from Cos ",nameTissue,sep=""))
        axis.circular(at=circular(mPeak[k,]), 
                      labels=mName[k,],tcl.text=-0.07,tick=FALSE)
        dev.off()
      }
      
    }
    if(veces>5000){
      if(k==0){
        mSel<-storeSel[(order(storeMaxDistEquis8))[1:nRep],]
        for(kkk in 1:nRep){
          mName[kkk,]<-namesRedMuscle[mSel[kkk,]]
          mPeak[kkk,]<-peakParRedMuscle[mSel[kkk,]]
          mCosPeak[kkk,]<-peakCosSin[mSel[kkk,]]
          mR2[kkk,]<-R2RedMuscle[mSel[kkk,]]
        }
        if(plotting){
          
          for(kkk in 1:nRep){
            png(file=paste("randomSel_Rep ",kkk," from ",nameTissue,"CosRefG.RData",sep=""))
            par(mfrow=c(1,1))
            par(mar=c(1,1,1,1))
            plot(as.circular(mPeak[kkk,]),main=paste("Rep ",kkk," from Cos ",nameTissue,sep=""))
            axis.circular(at=circular(mPeak[kkk,]), 
                          labels=mName[kkk,],tcl.text=-0.07,tick=FALSE)
            dev.off()
          }
        }
      }else{
        mSel<-storeSel[(order(storeMaxDistEquis8))[1:(nRep-k)],]
        for(kkk in (k+1):nRep){
          mName[kkk,]<-namesRedMuscle[mSel[kkk-k,]]
          mPeak[kkk,]<-peakParRedMuscle[mSel[kkk-k,]]
          mCosPeak[kkk,]<-peakCosSin[mSel[kkk-k,]]
          mR2[kkk,]<-R2RedMuscle[mSel[kkk-k,]]
        }
        if(plotting){
          
          for(kkk in (k+1):nRep){
            png(file=paste("randomSel_Rep ",kkk," from ",nameTissue,"CosRefG.RData",sep=""))
            par(mfrow=c(1,1))
            par(mar=c(1,1,1,1))
            plot(as.circular(mPeak[kkk-k,]),main=paste("Rep ",kkk," from Cos ",nameTissue,sep=""))
            axis.circular(at=circular(mPeak[kkk-k,]), 
                          labels=mName[kkk-k,],tcl.text=-0.07,tick=FALSE)
            dev.off()
          }
        }
      }
      k=4999
      
    }
    veces=veces+1
  }
  #selRefG<-list(mSel,mName,mPeak,mR2,mCosPeak,dosSect,distPeakCos,distEquis,distEquis8,distCond)
  return(list(mSel,mName,mPeak,mR2,mCosPeak,dosSect,distPeakCos,distEquis,distEquis8,distCond))
}
compRandomSel<-function(nRep,tam,namesRedMuscle,peakParRedMuscle,R2RedMuscle,sectorBelongingMuscle,nameTissue,plotting){
	k<-0
	dosSect=FALSE
	mSel<-matrix(0,nRep,tam)
	mName<-matrix(0,nRep,tam)
	mPeak<-matrix(0,nRep,tam)
	mR2<-matrix(0,nRep,tam)
	if(length(unique(sectorBelongingMuscle))<=2)dosSect=TRUE
	while(k<nRep){
		sel<-sample(1:length(namesRedMuscle),tam)
		highR2sel<-min(R2RedMuscle[sel])>0.5#tampoco va a hacer nada, ya lo hemos asegurado esto
		if(highR2sel==FALSE & namesRedMuscle[sel[which.min(R2RedMuscle[sel])]]=="ARNTL")highR2sel=TRUE
		atLeast2Sectors<-length(unique(sectorBelongingMuscle[sel]))>2
		atLeast1Sectors<-length(unique(sectorBelongingMuscle[sel]))>1
		if(((atLeast2Sectors & !dosSect) | (atLeast1Sectors & dosSect))& highR2sel){
			k<-k+1
			mSel[k,]<-sel
			mName[k,]<-namesRedMuscle[sel]
			mPeak[k,]<-peakParRedMuscle[sel]
			mR2[k,]<-R2RedMuscle[sel]
			
			


			if(plotting){
				#x11()
				png(file=paste("randomSel_Rep ",k," from ",nameTissue,"RefG.RData",sep=""))
				par(mfrow=c(1,1))
				par(mar=c(1,1,1,1))
				plot(as.circular(mPeak[k,]),main=paste("Rep ",k," from ",nameTissue,sep=""))
				axis.circular(at=circular(mPeak[k,]), 
					labels=mName[k,],tcl.text=-0.07,tick=FALSE)
				dev.off()
			}

		}
	}

	return(list(mSel,mName,mPeak,mR2))
}

giveR2<-function(dat,fit){
		return(1-mean((dat-fit)^2)/mean((dat-rep(mean(dat),length(dat)))^2))
	}

#nRep=K
#nameSel=mNameSelRefG[[numTissue]]
#mFullNorm=mFullSincroRefG[[numTissue]]
#load(file=paste("top",nameTissue,"RefG.RData",sep=""))

#top=topRefG
#nameTissue=tissues40[1]
#genAd=addArtGene[[numTissue]]
#escS=escSincroRefG[[numTissue]]
#indS=indSincroRefG[[numTissue]]

#load(file=paste("robEst",nameTissue,"RefG.RData",sep=""))

#nRep=K
#nameSel=mNameSelRefG[[numTissue]]
#mFullNorm=mFullSincroRefG[[numTissue]]
#top=topRefG
#nameTissue=nameTissue
#genAd=addArtGene[[numTissue]]
#escS=escSincroRefG[[numTissue]]
#indS=indSincroRefG[[numTissue]]
#coreG=coreG12

robustEst_v3_cores<-function(nRep,nameSel,mFullNorm,top,nameTissue,genAd,escS,indS,coreG){
giveCPCA<-list()
reFitSamples<-list()
tabSamples<-matrix(0,nRep*nrow(top),25)
mTopKNoSincro<-list()
mTopCoreKNoSincro<-list()
storeParRaw<-array(0,dim=c(length(coreG),7,nRep))
allFmmRaw<-list()
escKReorderNoSincro<-list()
listaArg<-list()
objOrdPre<-list()
objOrd<-list()
fmmSincroTop<-list()
cosSincroTop<-list()
npSincroTop<-list()
for(k in 1:nRep){
	giveCPCA[[k]]<-obtainCPCAMany(mFullNorm[match(nameSel[k,],rownames(mFullNorm)),],
	paste(nameTissue," Rep ",k,sep=""),nameSel[k,],8,FALSE)

	##############################################################################
	#a aprtir de este nuevo orden, oredenar la submatriz de cores y ajustar FMM
	#############################################################################

	#esc de la  ueva ordenacion 
	escKReorderNoSincroAux<-giveCPCA[[k]][[2]]
	escKReorderNoSincro[[k]]<-escKReorderNoSincroAux

	#matrix top reordenada sin sincro para cada k
	mTopKNoSincroAux<-top[,giveCPCA[[k]][[1]]]
	rownames(mTopKNoSincroAux)<-rownames(top)
	mTopKNoSincro[[k]]<-mTopKNoSincroAux
	
	#matrix top cores reordenada sin sincro para cada k
	mTopCoreKNoSincroAux<-top[match(coreG,rownames(top)),giveCPCA[[k]][[1]]]
	rownames(mTopCoreKNoSincroAux)<-coreG
	mTopCoreKNoSincro[[k]]<-mTopCoreKNoSincroAux
	
	#plot(1:length(mTopCoreKNoSincroAux[match("DBP",coreG),]),mTopCoreKNoSincroAux[match("DBP",coreG),],type="b")

	allFmmAux<-list()
	for(i in 1:length(coreG)){
		allFmmAux[[i]]<-fitFMM_Par(mTopCoreKNoSincro[[k]][i,],giveCPCA[[k]][[2]])
		storeParRaw[i,,k]<-c(allFmmAux[[i]][[5]],allFmmAux[[i]][[6]],allFmmAux[[i]][[2]],
						allFmmAux[[i]][[3]],allFmmAux[[i]][[4]],
						compUU(allFmmAux[[i]][[2]],allFmmAux[[i]][[3]],allFmmAux[[i]][[4]]),
						giveR2(mTopCoreKNoSincro[[k]][i,],allFmmAux[[i]][[1]]))
	}	
	allFmmRaw[[k]]<-allFmmAux
	
	#necesito una lista con los argumenttos que se usan en basicPreOder:peaks,esc,mFull,R2,alpha
	listaArgAux<-list(storeParRaw[,6,k],escKReorderNoSincro[[k]],mTopKNoSincro[[k]],
				storeParRaw[,7,k],storeParRaw[,1,k],storeParRaw[,2,k],storeParRaw[,3,k],
				storeParRaw[,4,k],storeParRaw[,5,k])
	listaArg[[k]]<-listaArgAux

	#preorder
	objOrdPreAux<-basicPreOder_v3_cores(listaArg[[k]],k,nameTissue,coreG)
	#objOrdPreAux<-list(oNew,escNew,matNew,peakEpidermis2,listaStep1[[4]],parCore,coreDay,r2Day,namesDay,coreNight,r2Night,namesNight)
	objOrdPre[[k]]<-objOrdPreAux
	#plantillaCircular_cores("ADISPOSE PRUEBA k=1 Pre",objOrdPreAux[[4]],coreG,objOrdPreAux[[5]])

	#order con orientacion
	objOrdAux<-basicOder_v4_cores(k,nameTissue,objOrdPre[[k]],coreG)
	#objOrdAux<-list(oNewNew,escNewNew,matNewNew,peaksPre,objPrep[[5]],pars,cambiaOri,indNewNew)
	objOrd[[k]]<-objOrdAux
	#plantillaCircular_cores("Epidermis k=2 Pre",objOrdAux[[4]],coreG,objOrdAux[[5]])

	#necesito sacar si hay cambio e orientacion 

	#preparamos salidas iguales a lo antiguo

	#En la matriz reordenada ajustes de FMM, COSINOR Y NP
	mSincroK<-objOrdAux[[3]]
	escSincroK<-objOrdAux[[2]]
	fmmSincroKTop<-list()
	cosSincroKTop<-list()
	npSincroKTop<-list()
	mFitP3N<-matrix(0,nrow(mSincroK),ncol(mSincroK))
	mFitC3N<-matrix(0,nrow(mSincroK),ncol(mSincroK))
	mFitNP3N<-matrix(0,nrow(mSincroK),ncol(mSincroK))
	mStats<-c()
	for(i in 1:nrow(mSincroK)){
	
		#esc
		escE=(giveCPCA[[k]][[2]]-listaArg[[k]][[1]][[6]]+pi)%%(2*pi)
		
		#datos
		vvv<-mTopKNoSincro[[k]][i,]#mSincroK[i,]
		#plot(1:length(vvv),vvv,type="b")
		#FMM
		fmmSincroKTopAux<-fitFMM_Par(vvv,giveCPCA[[k]][[2]])#escSincroK)
		fmmSincroKTop[[i]]<-fmmSincroKTopAux#no sirve para mucho

		#1.-par y luego fitOut
		al<-(fmmSincroKTopAux[[2]]-listaArg[[k]][[1]][match("ARNTL",coreG)]+pi)%%(2*pi)
		be<-fmmSincroKTopAux[[3]]
		om<-fmmSincroKTopAux[[4]]
		em<-fmmSincroKTopAux[[5]]
		aa<-fmmSincroKTopAux[[6]]

		mFitP3N[i,]<-fmmSincroKTopAux[[1]]
		
		if(objOrdPre[[k]][[13]]==objOrd[[k]][[7]]){#si la orientacion no cambia
			pars<-c(em,aa,al,be,om)
			peaks<-c(compUU(pars[3],pars[4],pars[5]),compLL(pars[3],pars[4],pars[5]),
				compUU(pars[3],pars[4],pars[5])/(2*pi)*100,compLL(pars[3],pars[4],pars[5])/(2*pi)*100)
			sFMM<-sum((vvv-mFitP3N[i,])^2)/(length(mFitP3N[i,])-5)
			mseFMM<-sum((vvv-mFitP3N[i,])^2)/(length(mFitP3N[i,]))
			r2FMM<-giveR2(vvv,mFitP3N[i,])
			statFMM<-c(pars,peaks,sFMM,mseFMM,r2FMM)
			mFitP3N[i,]<-mFitP3N[i,order(escE)]

		}else{#si cambia la orintacion
			pars<-c(em,aa,(2*pi-al)%%(2*pi),(2*pi-be)%%(2*pi),om)
			peaks<-c(compUU(pars[3],pars[4],pars[5]),compLL(pars[3],pars[4],pars[5]),
				compUU(pars[3],pars[4],pars[5])/(2*pi)*100,compLL(pars[3],pars[4],pars[5])/(2*pi)*100)
			sFMM<-sum((vvv-mFitP3N[i,])^2)/(length(mFitP3N[i,])-5)
			mseFMM<-sum((vvv-mFitP3N[i,])^2)/(length(mFitP3N[i,]))
			r2FMM<-giveR2(vvv,mFitP3N[i,])
			statFMM<-c(pars,peaks,sFMM,mseFMM,r2FMM)
			mFitP3N[i,]<-rev(mFitP3N[i,order(escE)])
		}
		#Cosinor
		cosSincroKTopAux<-funcionCosinor(vvv,giveCPCA[[k]][[2]],length(giveCPCA[[k]][[2]]))#escSincroK,length(escSincroK))
		cosSincroKTop[[i]]<-cosSincroKTopAux

		emC<-cosSincroKTopAux[[2]]
		aC<-cosSincroKTopAux[[3]]	
		phiC<-(cosSincroKTopAux[[5]]-listaArg[[k]][[1]][[6]]+pi)%%(2*pi)


		mFitC3N[i,]<-cosSincroKTopAux[[1]]
		if(objOrdPre[[k]][[13]]==objOrd[[k]][[7]]){#si la orientacion no cambia
			pars<-c(emC,aC,phiC)
			peaks<-c((-phiC)%%(2*pi),(pi-phiC)%%(2*pi),
				((-phiC)%%(2*pi))/(2*pi)*100,((pi-phiC)%%(2*pi))/(2*pi)*100)
			sCos<-sum((vvv-mFitC3N[i,])^2)/(length(mFitP3N[i,])-3)
			mseCos<-sum((vvv-mFitC3N[i,])^2)/(length(mFitC3N[i,]))
			r2Cos<-giveR2(vvv,mFitC3N[i,])
			statCos<-c(pars,peaks,sCos,mseCos,r2Cos)
			mFitC3N[i,]<-mFitC3N[i,order(escE)]
		}else{
			pars<-c(emC,aC,(2*pi-phiC)%%(2*pi))
			peaks<-c((-cosSincroKTopAux[[5]])%%(2*pi),(pi-cosSincroKTopAux[[5]])%%(2*pi),
				((-cosSincroKTopAux[[5]])%%(2*pi))/(2*pi)*100,((pi-cosSincroKTopAux[[5]])%%(2*pi))/(2*pi)*100)
			sCos<-sum((vvv-mFitC3N[i,])^2)/(length(mFitP3N[i,])-3)
			mseCos<-sum((vvv-mFitC3N[i,])^2)/(length(mFitC3N[i,]))
			r2Cos<-giveR2(vvv,mFitC3N[i,])
			statCos<-c(pars,peaks,sCos,mseCos,r2Cos)
			mFitC3N[i,]<-rev(mFitC3N[i,order(escE)])
		}
		#plot(objOrd[[k]][[2]],mFitP3N[i,],type="b")
		#lines(objOrd[[k]][[2]],objOrd[[k]][[3]][i,],col=2)
		#lines(objOrd[[k]][[2]],mFitC3N[i,],col=3)
		#lines(objOrd[[k]][[2]],mFitNP3N[i,],col=4)
		#NP
		if(objOrdPre[[k]][[13]]==objOrd[[k]][[7]]){
			npSincroKTopAux<-function1Local(vvv[order(escE)])#function1Local(mTopKNoSincro[[k]][i,order(escE)])
			mFitNP3N[i,]<-npSincroKTopAux[[1]]
			sNp<-sum((vvv[order(escE)]-mFitNP3N[i,])^2)/(length(mFitNP3N[i,])-3)
			mseNp<-sum((vvv[order(escE)]-mFitNP3N[i,])^2)/(length(mFitNP3N[i,]))
			r2Np<-giveR2(vvv[order(escE)],mFitNP3N[i,])
			statNp<-c(sNp,mseNp,r2Np)
			npSincroKTop[[i]]<-npSincroKTopAux
		}else{
			npSincroKTopAux<-function1Local(rev(vvv[order(escE)]))
			mFitNP3N[i,]<-npSincroKTopAux[[1]]
			sNp<-sum((rev(vvv[order(escE)])-mFitNP3N[i,])^2)/(length(mFitNP3N[i,])-3)
			mseNp<-sum((rev(vvv[order(escE)])-mFitNP3N[i,])^2)/(length(mFitNP3N[i,]))
			r2Np<-giveR2(rev(vvv[order(escE)]),mFitNP3N[i,])
			statNp<-c(sNp,mseNp,r2Np)
			npSincroKTop[[i]]<-npSincroKTopAux
		}
		if(objOrdPre[[k]][[13]]!=objOrd[[k]][[7]])mFitNP3N[i,]<-rev(npSincroKTopAux[[1]][order(escE)])

		

		statVAux<-c(statFMM,statCos,statNp)
		mStats<-rbind(mStats,statVAux)
	}
	reFitAux<-list(objOrdAux[[1]],objOrdAux[[2]],objOrdAux[[8]],mFitP3N,mFitC3N,mFitNP3N,mStats)
	fmmSincroTop[[k]]<-fmmSincroKTop
	cosSincroTop[[k]]<-cosSincroKTop
	npSincroTop[[k]]<-npSincroKTop
	reFitSamples[[k]]<-reFitAux
	tabSamples[seq(k,nrow(top)*nRep,nRep),]<-reFitSamples[[k]][[7]]
	
}
return(list(tabSamples,reFitSamples,giveCPCA,objOrdAux[[7]],fmmSincroTop,cosSincroTop,npSincroTop,
		mTopKNoSincro,storeParRaw,allFmmRaw,escKReorderNoSincro,listaArg,objOrdPre,objOrd))
}



#listaStep1<-listaArg[[k]]
#repe<-k
#tejidoN<-"Epidermis"

basicPreOder_v3_cores<-function(listaStep1,repe,tejidoN,coreG){
peaksEpidermis2<-listaStep1[[1]]


peakRefEpidermis2<-peaksEpidermis2[match("ARNTL",coreG)]
peakDBPEpidermis2<-peaksEpidermis2[match("DBP",coreG)]
escTEpidermis<-order((listaStep1[[2]]-peakRefEpidermis2+pi)%%(2*pi))
esc2<-((listaStep1[[2]]-peakRefEpidermis2+pi)%%(2*pi))[escTEpidermis]
#x11()
png(filename=paste0("12CorePreRep",repe,"_",tejidoN,".png"))
par(mfrow=c(compPerfSq(1,length(coreG)),compPerfSq(1,length(coreG))))
par(mar=c(2,2,2,2))
peakEpidermis2<-c()
a=TRUE
parCore<-matrix(0,length(coreG),5)
for(i in 1:length(coreG)){
	dat_i<-listaStep1[[3]][match(coreG[i],rownames(listaStep1[[3]])),]
	plot(0,0,type="n",main=paste0(tejidoN," ",coreG[i]," R2=",round(listaStep1[[4]][i],2)),xaxt="n",yaxt="n",ylab="Exp",xlab="Time",ylim=c(-1,1),xlim=c(0,2*pi))
		newAlphaEpidermis<-(listaStep1[[7]][i]-peakRefEpidermis2+pi)%%(2*pi)
		if((peakDBPEpidermis2-peakRefEpidermis2+pi)%%(2*pi)<pi){
		#if(!a){
			lines(esc2,dat_i[escTEpidermis],col=rainbowColor(coreG)[i])
			points(esc2,dat_i[escTEpidermis],col=rainbowColor(coreG)[i])
			lines(seq(0,2*pi,length.out=length(dat_i)),fitMob(listaStep1[[5]][i],listaStep1[[6]][i],newAlphaEpidermis,
				listaStep1[[8]][i],listaStep1[[9]][i],seq(0,2*pi,length.out=length(dat_i))),col=1)
			peakEpidermis2[i]<-compUU(newAlphaEpidermis,listaStep1[[8]][i],listaStep1[[9]][i])
			parCore[i,]<-c(listaStep1[[5]][i],listaStep1[[6]][i],
				newAlphaEpidermis,listaStep1[[8]][i],listaStep1[[9]][i])
		}else{
			lines(2*pi-rev(esc2),rev(dat_i[escTEpidermis]),col=rainbowColor(coreG)[i])
			points(2*pi-rev(esc2),rev(dat_i[escTEpidermis]),col=rainbowColor(coreG)[i])
			dd<-outF2(rev(dat_i),c(listaStep1[[5]][i],listaStep1[[6]][i],
				2*pi-newAlphaEpidermis,2*pi-listaStep1[[8]][i],listaStep1[[9]][i]),2*pi-rev(esc2))
			lines(2*pi-rev(esc2),dd,col=1)
			peakEpidermis2[i]<-compUU(2*pi-newAlphaEpidermis,2*pi-listaStep1[[8]][i],listaStep1[[9]][i])
			parCore[i,]<-c(listaStep1[[5]][i],listaStep1[[6]][i],
				2*pi-newAlphaEpidermis,2*pi-listaStep1[[8]][i],listaStep1[[9]][i])
		}
}
matNew<-matrix(0,nrow(listaStep1[[3]]),ncol(listaStep1[[3]]))
cambioOriPre<-FALSE
if((peakDBPEpidermis2-peakRefEpidermis2+pi)%%(2*pi)<pi){
	oNew<-escTEpidermis
	escNew<-esc2
	matNew<-listaStep1[[3]][,escTEpidermis]
	
}else{
	oNew<-rev(escTEpidermis)
	escNew<-2*pi-rev(esc2)
	for(i in 1:nrow(listaStep1[[3]])){
		matNew[i,]<-rev(listaStep1[[3]][i,escTEpidermis])
	}
	cambioOriPre<-TRUE
}
rownames(matNew)<-rownames(listaStep1[[3]])
colnames(matNew)<-colnames(listaStep1[[3]])
dev.off()
coreDay<-c();namesDay<-c();r2Day<-c()
coreNight<-c();namesNight<-c();r2Night<-c()
for(i in 1:length(peakEpidermis2)){
	if(coreG[i]!="ARNTL" & coreG[i]!="DBP"){
		if(0<=peakEpidermis2[i] & peakEpidermis2[i]<pi ){
			coreDay<-c(coreDay,peakEpidermis2[i])
			r2Day<-c(r2Day,listaStep1[[4]][i])
			namesDay<-c(namesDay,coreG[i])
		}else{
			coreNight<-c(coreNight,peakEpidermis2[i])
			r2Night<-c(r2Night,listaStep1[[4]][i])
			namesNight<-c(namesNight,coreG[i])
		}
	}
}
return(list(oNew,escNew,matNew,peakEpidermis2,listaStep1[[4]],parCore,coreDay,r2Day,namesDay,coreNight,r2Night,namesNight,cambioOriPre))
}

#repi=k
#tejidoN=nameTissue
#objPrep=objOrdPre[[k]]

basicOder_v3_cores<-function(repi,tejidoN,objPrep,coreG){
	peaksPre<-objPrep[[4]]
	aviso<-FALSE
	peakD<-0
	p6am<-1-cos(peaksPre[match("DBP",coreG)])
	p6pm<-1-cos(peaksPre[match("DBP",coreG)]-pi)
	for(i in 1:length(coreG)){
		if(i!=match("ARNTL",coreG) & i!=match("DBP",coreG) & peaksPre[i]>0 & peaksPre[i]<pi)peakD<-peakD+1
	}
	mitad<-floor(length(coreG)-2)/2
	if(p6am<=0.1 | p6pm<=0.1 | peakD<mitad)aviso<-TRUE
	cambiaOri<-FALSE
	peaksPreNew<-c()
	if(aviso & !(peaksPre[match("DBP",coreG)]<peaksPre[match("CRY1",coreG)] & peaksPre[match("CRY1",coreG)]<peaksPre[match("ARNTL",coreG)])){
		peaksPreNew<-(2*pi-peaksPre)%%(2*pi)
		print(tejidoN)
		oNewNew<-rev(objPrep[[1]])
		escNewNew<-2*pi-rev(objPrep[[2]])
		matNewNew<-matrix(0,nrow(objPrep[[3]]),ncol(objPrep[[3]]))
		pars<-matrix(0,length(coreG),5)
		for( i in 1:length(coreG)){
			pars[i,]<-c(objPrep[[6]][i,1],objPrep[[6]][i,2],2*pi-objPrep[[6]][i,3],2*pi-objPrep[[6]][i,4],objPrep[[6]][i,5])
		}
		for(i in 1:nrow(objPrep[[3]])){
			matNewNew[i,]<-rev(objPrep[[3]][i,])
		}
		rownames(matNewNew)<-rownames(objPrep[[3]])
		indNewNew<-rev(1:length(oNewNew))
		cambiaOri<-TRUE
	}else{
		peaksPreNew<-peaksPre
		oNewNew<-objPrep[[1]]
		escNewNew<-objPrep[[2]]
		matNewNew<-objPrep[[3]]
		rownames(matNewNew)<-rownames(objPrep[[3]])
		indNewNew<-1:length(oNewNew)
		pars<-matrix(0,length(coreG),5)
		for( i in 1:length(coreG)){
			pars[i,]<-c(objPrep[[6]][i,1],objPrep[[6]][i,2],objPrep[[6]][i,3],objPrep[[6]][i,4],objPrep[[6]][i,5])
		}
	}

	#x11()
if(!cambiaOri)png(filename=paste0("12CoreF_Rep",repi,"_",tejidoN,".png"))
if(cambiaOri)png(filename=paste0("!!12CoreF_Rep",repi,"_",tejidoN,".png"))
par(mfrow=c(compPerfSq(1,length(coreG)),compPerfSq(1,length(coreG))))
par(mar=c(2,2,2,2))
peakEpidermis2<-c()
for(i in 1:length(coreG)){
	dat_i<-matNewNew[match(coreG[i],rownames(matNewNew)),]
	plot(0,0,type="n",main=paste0("Rep",repi,". ",tejidoN," ",coreG[i]," R2=",round(objPrep[[5]][i],2)),xaxt="n",yaxt="n",ylab="Exp",xlab="Time",ylim=c(-1,1),xlim=c(0,2*pi))
	lines(escNewNew,dat_i,col=rainbowColor(coreG)[i])
	points(escNewNew,dat_i,col=rainbowColor(coreG)[i])
	ss<-seq(0,2*pi,length.out=length(dat_i))
	lines(ss,fitMob(pars[i,1],pars[i,2],pars[i,3],pars[i,4],pars[i,5],seq(0,2*pi,length.out=length(dat_i))),col=1)
	axis(1,at=c(ss[1],ss[ceiling(length(ss)/2)]),labels=c("6AM","6PM"))
}
dev.off()

if(!cambiaOri)png(filename=paste0("12PeaksRep",repi,"_",tejidoN,".png"))
if(cambiaOri)png(filename=paste0("!!12PeaksRep",repi,"_",tejidoN,".png"))

plantillaCircular_cores(tejidoN,peaksPreNew,coreG,objPrep[[5]])
dev.off()


	return(list(oNewNew,escNewNew,matNewNew,peaksPreNew,objPrep[[5]],pars,cambiaOri,indNewNew))
}
#repi=k
#tejidoN=nameTissue
#objPrep=objOrdPre[[k]]
#coreG=coreG12
basicOder_v4_cores<-function(repi,tejidoN,objPrep,coreG){
	peaksPre<-objPrep[[4]]
	aviso1<-FALSE;aviso2<-FALSE;aviso3<-FALSE;cond31<-FALSE;cond32<-FALSE
	peakD<-0
	p6am<-1-cos(peaksPre[match("DBP",coreG)])
	p6pm<-1-cos(peaksPre[match("DBP",coreG)]-pi)
	for(i in 1:length(coreG)){
		if(i!=match("ARNTL",coreG) & i!=match("DBP",coreG) & peaksPre[i]>0 & peaksPre[i]<pi)peakD<-peakD+1
	}
	mitad<-floor(length(coreG)-2)/2
	if(p6am<=0.1 )aviso1<-TRUE# | peakD<mitad)aviso1<-TRUE
	if(peaksPre[match("DBP",coreG)]>pi/2 & p6pm>0.1)aviso2<-TRUE
	if(p6pm<=0.1 )aviso3<-TRUE
	cond1<-(aviso1 & peakD<=mitad &
	!(peaksPre[match("DBP",coreG)]<peaksPre[match("CRY1",coreG)] & peaksPre[match("CRY1",coreG)]<peaksPre[match("ARNTL",coreG)]) )|
	(aviso1 & peakD<mitad &
	(peaksPre[match("DBP",coreG)]<peaksPre[match("CRY1",coreG)] & peaksPre[match("CRY1",coreG)]<peaksPre[match("ARNTL",coreG)]) )
	cond2<-aviso2 & peakD<mitad &
		!(peaksPre[match("DBP",coreG)]<peaksPre[match("ARNTL",coreG)] & peaksPre[match("ARNTL",coreG)]<peaksPre[match("CRY1",coreG)])
	if(peaksPre[match("TEF",coreG)]>pi & (peaksPre[match("CRY1",coreG)]>peaksPre[match("TEF",coreG)] | peaksPre[match("CRY1",coreG)]<pi))cond31<-TRUE
	if(peaksPre[match("TEF",coreG)]<pi & peaksPre[match("CRY1",coreG)]>peaksPre[match("TEF",coreG)] & peaksPre[match("CRY1",coreG)]<pi)cond32<-TRUE
	cond3<- aviso3 & !cond31 & !cond32 
	cond4<-!aviso1 & !aviso2 &  (1-cos(peaksPre[match("CRY1",coreG)]-pi))<0.1 & peakD<mitad
	cambiaOri<-FALSE
	peaksPreNew<-c()
	if(cond1 | cond2 |cond3 | cond4){
		peaksPreNew<-(2*pi-peaksPre)%%(2*pi)
		print(tejidoN)
		oNewNew<-rev(objPrep[[1]])
		escNewNew<-2*pi-rev(objPrep[[2]])
		matNewNew<-matrix(0,nrow(objPrep[[3]]),ncol(objPrep[[3]]))
		pars<-matrix(0,length(coreG),5)
		for( i in 1:length(coreG)){
			pars[i,]<-c(objPrep[[6]][i,1],objPrep[[6]][i,2],2*pi-objPrep[[6]][i,3],2*pi-objPrep[[6]][i,4],objPrep[[6]][i,5])
		}
		for(i in 1:nrow(objPrep[[3]])){
			matNewNew[i,]<-rev(objPrep[[3]][i,])
		}
		rownames(matNewNew)<-rownames(objPrep[[3]])
		indNewNew<-rev(1:length(oNewNew))
		cambiaOri<-TRUE
	}else{
		peaksPreNew<-peaksPre
		oNewNew<-objPrep[[1]]
		escNewNew<-objPrep[[2]]
		matNewNew<-objPrep[[3]]
		rownames(matNewNew)<-rownames(objPrep[[3]])
		indNewNew<-1:length(oNewNew)
		pars<-matrix(0,length(coreG),5)
		for( i in 1:length(coreG)){
			pars[i,]<-c(objPrep[[6]][i,1],objPrep[[6]][i,2],objPrep[[6]][i,3],objPrep[[6]][i,4],objPrep[[6]][i,5])
		}
	}

	#x11()
if(!cambiaOri)png(filename=paste0("12CoreF_Rep",repi,"_",tejidoN,".png"))
if(cambiaOri)png(filename=paste0("!!12CoreF_Rep",repi,"_",tejidoN,".png"))
par(mfrow=c(compPerfSq(1,length(coreG)),compPerfSq(1,length(coreG))))
par(mar=c(2,2,2,2))
peakEpidermis2<-c()
for(i in 1:length(coreG)){
	dat_i<-matNewNew[match(coreG[i],rownames(matNewNew)),]
	plot(0,0,type="n",main=paste0("Rep",repi,". ",tejidoN," ",coreG[i]," R2=",round(objPrep[[5]][i],2)),xaxt="n",yaxt="n",ylab="Exp",xlab="Time",ylim=c(-1,1),xlim=c(0,2*pi))
	lines(escNewNew,dat_i,col=rainbowColor(coreG)[i])
	points(escNewNew,dat_i,col=rainbowColor(coreG)[i])
	ss<-seq(0,2*pi,length.out=length(dat_i))
	lines(ss,fitMob(pars[i,1],pars[i,2],pars[i,3],pars[i,4],pars[i,5],seq(0,2*pi,length.out=length(dat_i))),col=1)
	axis(1,at=c(ss[1],ss[ceiling(length(ss)/2)]),labels=c("6AM","6PM"))
}
dev.off()

if(!cambiaOri)png(filename=paste0("12PeaksRep",repi,"_",tejidoN,".png"))
if(cambiaOri)png(filename=paste0("!!12PeaksRep",repi,"_",tejidoN,".png"))

plantillaCircular_cores(tejidoN,peaksPreNew,coreG,objPrep[[5]])
dev.off()


	return(list(oNewNew,escNewNew,matNewNew,peaksPreNew,objPrep[[5]],pars,cambiaOri,indNewNew))
}


robustEst<-function(nRep,nameSel,mFullNorm,top,nameTissue,genAd,escS,indS){
giveCPCA<-list()
reFitSamples<-list()
tabSamples<-matrix(0,nRep*nrow(top),25)
for(k in 1:nRep){
	giveCPCA[[k]]<-obtainCPCAMany(mFullNorm[match(nameSel[k,],rownames(mFullNorm)),],
	paste(nameTissue," Rep ",k,sep=""),nameSel[k,],8,FALSE)
	reFitSamples[[k]]<-reajustarPlotsNewReferences_vInicial("ARNTL","DBP",top[,giveCPCA[[k]][[1]]],#OJO
		giveCPCA[[k]][[1]],giveCPCA[[k]][[2]],rownames(top))#OJO
	#reFitSamples[[1]]<-list(oN3,pN3,indN3,mFitP3N,mFitC3N,mFitNP3N,mStatistics)
	tabSamples[seq(k,nrow(top)*nRep,nRep),]<-reFitSamples[[k]][[7]]
	
}
return(list(tabSamples,reFitSamples,giveCPCA))
}




#nameTissue=allTissues[1]
#load(file=paste("D://gTEX//dataBase//infoInicial//rdata//R2_NPRefG_",ageMaleTissues[[ageGroup]][1],".RData",sep=""))
#R2Tissue=R2_NPRefGaux
#load(file=paste("D://gTEX//dataBase//tissueOutputsNew2//giveMatIniRefG",allTissues[1],".RData",sep=""))
#mFullTissueNorm=giveMatIniRefG[[3]]
#load(file=paste("D://gTEX//dataBase//tissueOutputsNew2//outSetCandRefG",allTissues[1],".RData",sep=""))
#namesAnal=outSetCandRefG[[1]]
#datAnal=outSetCandRefG[[2]]
#parAnal=outSetCandRefG[[4]]
#peakAnal=outSetCandRefG[[7]]
#r2Anal=outSetCandRefG[[9]]
#escParTissue=giveMatIniRefG[[2]]

ritGen=function(nameTissue,R2Tissue,mFullTissueNorm,namesAnal,datAnal,parAnal,peakAnal,r2Anal,escParTissue,descartes){
#primero la distribucion de los R2
distR2Tissue<-R2Tissue[order(R2Tissue,decreasing=TRUE)]
mDistR2Tissue<-matrix(distR2Tissue,length(distR2Tissue),1)
rownames(mDistR2Tissue)<-rownames(mFullTissueNorm)[order(R2Tissue,decreasing=TRUE)]
colnames(mDistR2Tissue)<-"RankR2_NP"

genesR2Tissue<-rownames(mDistR2Tissue)[1:(which(mDistR2Tissue<0.5)[1]-1)]
entra=1
datosParTissue<-c()
paramParamTissue<-c()
peaksIn<-c()
r2ParIn<-c()
namesIn<-c()
ss=0
ss1=0
storeF=list()
sigue=TRUE
while(sigue==TRUE & entra<=length(genesR2Tissue)){
	if(is.na(match(genesR2Tissue[entra],namesAnal))){
		if(is.na(match(genesR2Tissue[entra],descartes))){
			datParTissue<-mFullTissueNorm[match(genesR2Tissue[entra],rownames(mFullTissueNorm)),]
			fParTissue<-fitFMM_Par(datParTissue,escParTissue)
			wTissue<-fParTissue[[4]]
			peakTissue<-compUU(fParTissue[[2]],fParTissue[[3]],fParTissue[[4]])
      		r2ParAux<-1-mean((datParTissue-fParTissue[[1]])^2)/mean((datParTissue-mean(datParTissue))^2)
			print(paste(" Analysing gene ",rownames(mDistR2Tissue)[entra]," R2NP=",round(mDistR2Tissue[entra],3)," R2Par=",round(r2ParAux,3),sep=""))
			if(r2ParAux<0.6 & wTissue>0.1 ){
				sigue=FALSE
			}else{
				if(r2ParAux>0.6 & wTissue>0.1){
					ss=ss+1
					print(paste("En ",nameTissue," hay ",ss," genes ritmicos",sep=""))
					datosParTissue<-rbind(datosParTissue,datParTissue)
					paramParamTissue<-rbind(paramParamTissue,c(fParTissue[[5]],fParTissue[[6]],fParTissue[[2]],fParTissue[[3]],fParTissue[[4]]))
					peaksIn<-c(peaksIn,peakTissue)
					r2ParIn<-c(r2ParIn,r2ParAux)
					namesIn<-c(namesIn,rownames(mDistR2Tissue)[entra])
					storeF[[ss]]=fParTissue
				}
			}
		}
	}else{
		if(parAnal[match(genesR2Tissue[entra],rownames(datAnal)),5]>0.6){
			ss1=ss1+1
			print(paste("En ",nameTissue," hay ",ss+ss1," genes ritmicos",sep=""))
			#nameTissue,R2Tissue,mFullTissueNorm,namesAnal,datAnal,parAnal,peakAnal,r2Anal
			datosParTissue<-rbind(datosParTissue,datAnal[match(genesR2Tissue[entra],rownames(datAnal)),])
			paramParamTissue<-rbind(paramParamTissue,parAnal[match(genesR2Tissue[entra],rownames(datAnal)),])
			peaksIn<-c(peaksIn,peakAnal[match(genesR2Tissue[entra],rownames(datAnal))])
			r2ParIn<-c(r2ParIn,r2Anal[match(genesR2Tissue[entra],rownames(datAnal))])
			namesIn<-c(namesIn,genesR2Tissue[entra])
		}
	}
	entra=entra+1
}
	total=ss+ss1
	return(list(total,datosParTissue,paramParamTissue,peaksIn,r2ParIn,namesIn,storeF))

}












#R2Tissue=R2_NPRefG[[numTissue]]
#nameTissue=allTissues[numTissue]
#mFullTissueNorm=mFullNormRefG[[numTissue]]
#escParTissue=escRefG[[numTissue]]
ritTiss=function(R2Tissue,nameTissue,mFullTissueNorm,escParTissue){
#primero la distribucion de los R2
distR2Tissue<-R2Tissue[order(R2Tissue,decreasing=TRUE)]
mDistR2Tissue<-matrix(distR2Tissue,length(distR2Tissue),1)
rownames(mDistR2Tissue)<-rownames(mFullTissueNorm)[order(R2Tissue,decreasing=TRUE)]
colnames(mDistR2Tissue)<-"RankR2_NP"

genesR2Tissue<-rownames(mFullTissueNorm)[which(R2Tissue>0.5)]
cosPeakTissueAll<-c()
for(i in 1:length(genesR2Tissue)){
	ff<-mFullTissueNorm[match(genesR2Tissue[i],rownames(mFullTissueNorm)),]
	cosPeakTissueAll[i]<-(-funcionCosinor(ff,escParTissue,length(escParTissue))[[5]])%%(2*pi)
}

#seran candidatos al conjunto de genes minimos no mas de los primeros 500 genes 
#con R2_NP mayor del percentil 50, con w>0.05 en el paramertirvo y r2 par >0.5
in1=0
entra<-1
hasta=which(mDistR2Tissue[,1]<0.5)[1]
sigue=TRUE
adjCosTissue<-c()
datosParTissue<-c()
paramCosTissue<-c()
paramParamTissue<-c()
adjParamTissue<-c()
r2ParIn<-c()
peaksIn=c()
namesIn=c()
while(entra<=hasta & sigue==TRUE){# & in1<25){
	datParTissue<-mFullTissueNorm[match(rownames(mDistR2Tissue)[entra],rownames(mFullTissueNorm)),]
	#original#fParTissue<-fitFMM_Par2(datParTissue,escParTissue)
	fParTissue<-fitFMM_Par(datParTissue,escParTissue)
	wTissue<-fParTissue[[4]]
	peakTissue<-compUU(fParTissue[[2]],fParTissue[[3]],fParTissue[[4]])
      r2ParAux<-1-mean((datParTissue-fParTissue[[1]])^2)/mean((datParTissue-mean(datParTissue))^2)
	print(paste(" Analysing gene ",rownames(mDistR2Tissue)[entra]," R2NP=",round(mDistR2Tissue[entra],3)," R2Par=",round(r2ParAux,3),sep=""))
	if(r2ParAux<0.7 & wTissue>0.1 ){
		sigue=FALSE
	}else{
		if(r2ParAux>0.7 & wTissue>0.1 & in1<25){
			in1<-in1+1
			print(paste("There are  ",entra," pot rhyth genes",sep=""))
			print(paste("There are  ",in1," rhythmic genes",sep=""))
			datosParTissue<-rbind(datosParTissue,datParTissue)
			fCosTissue<-funcionCosinor(datParTissue,escParTissue,length(escParTissue))
			paramCosTissue<-rbind(paramCosTissue,c(fCosTissue[[2]],fCosTissue[[3]],fCosTissue[[5]]))
			adjCosTissue<-rbind(adjCosTissue,fCosTissue[[1]])
			adjParamTissue<-rbind(adjParamTissue,fParTissue[[1]])
			paramParamTissue<-rbind(paramParamTissue,c(fParTissue[[5]],fParTissue[[6]],fParTissue[[2]],fParTissue[[3]],fParTissue[[4]]))
			peaksIn<-c(peaksIn,peakTissue)
			r2ParIn<-c(r2ParIn,r2ParAux)
			namesIn<-c(namesIn,rownames(mDistR2Tissue)[entra])
			save(datosParTissue,file=paste("datosParTissue_",nameTissue,".RData",sep=""))
			save(adjParamTissue,file=paste("adjParamTissue_",nameTissue,".RData",sep=""))
			save(paramParamTissue,file=paste("paramParamTissue_",nameTissue,".RData",sep=""))
			save(adjCosTissue,file=paste("adjCosTissue_",nameTissue,".RData",sep=""))
			save(paramCosTissue,file=paste("paramCosTissue_",nameTissue,".RData",sep=""))
			save(peaksIn,file=paste("peaksIn_",nameTissue,".RData",sep=""))
			save(r2ParIn,file=paste("r2ParIn_",nameTissue,".RData",sep=""))
			save(namesIn,file=paste("namesIn_",nameTissue,".RData",sep=""))
		}
	}
	entra=entra+1
}
rownames(datosParTissue)<-namesIn
m<-cbind(r2ParIn,peaksIn)
rownames(m)<-namesIn
colnames(m)<-c("R2_Par","Peak_Par")
write.table(m,file=paste("mCand",nameTissue,".txt",sep=""),sep="\t",dec=",",col.names=TRUE,row.names=TRUE)
return(list(namesIn,datosParTissue,adjParamTissue,paramParamTissue,adjCosTissue,paramCosTissue,peaksIn,r2ParIn,m,mDistR2Tissue,entra))
}

#R2Tissue=R2_NPRefG[[numTissue]]
#nameTissue=nameTissue
#mFullTissueNorm=mFullNormRefG[[numTissue]]
#escParTissue=escRefG[[numTissue]]
#tam=2

#R2Tissue=R2_NPRefGSkin
#nameTissue="Epidermis"
#mFullTissueNorm=mFullNormRefGSkin
#escParTissue=escRefGSkin
#tam=2
#R2ParcoreG=basicOrdRefSkin[[5]]

minSetCand3_v2_cores<-function(R2Tissue,nameTissue,mFullTissueNorm,escParTissue,tam,R2ParcoreG){
#primero la distribucion de los R2
distR2Tissue<-R2Tissue[order(R2Tissue,decreasing=TRUE)]
mDistR2Tissue<-matrix(distR2Tissue,length(distR2Tissue),1)
rownames(mDistR2Tissue)<-rownames(mFullTissueNorm)[order(R2Tissue,decreasing=TRUE)]
colnames(mDistR2Tissue)<-"RankR2_NP"
#write.table(mDistR2Tissue,file=paste("mDistR2Tissue",nameTissue,".txt",sep=""),sep="\t",dec=",",row.names=TRUE)

if(quantile(distR2Tissue)[3]>0.5){
	R2NP50<-quantile(distR2Tissue)[3]
}else{
	R2NP50<-quantile(distR2Tissue)[3]
	print("Se ha puesto 0.5 en r2 de NP")
}
genesR2Tissue<-rownames(mFullTissueNorm)[which(R2Tissue>R2NP50)]
cosPeakTissueAll<-c()
for(i in 1:length(genesR2Tissue)){
	ff<-mFullTissueNorm[match(genesR2Tissue[i],rownames(mFullTissueNorm)),]
	cosPeakTissueAll[i]<-(-funcionCosinor(ff,escParTissue,length(escParTissue))[[5]])%%(2*pi)
}

#seran candidatos al conjunto de genes minimos no mas de los primeros 500 genes 
#con R2_NP mayor del percentil 50, con w>0.05 en el paramertirvo y r2 par >0.5
entra<-0

adjCosTissue<-c()
namGenesInTissue<-c()
datosParTissue<-c()
paramCosTissue<-c()
paramParamTissue<-c()
adjParamTissue<-c()
regIn<-c()
r2ParIn<-c()
peaksIn<-c()
namesIn<-c()
seguir<-TRUE
entra<-0
peakRef<-median.circular(as.circular(cosPeakTissueAll))%%(2*pi)
m1<-peakRef
m2<-(peakRef-pi/4)%%(2*pi)
m3<-(peakRef-pi/2)%%(2*pi)
m4<-(peakRef-3*pi/4)%%(2*pi)
m5<-(peakRef-pi)%%(2*pi)
m6<-(peakRef-5*pi/4)%%(2*pi)
m7<-(peakRef-3*pi/2)%%(2*pi)
m8<-(peakRef-7*pi/4)%%(2*pi)
r10<-0
r20<-(m2-peakRef)%%(2*pi)
r30<-(m3-peakRef)%%(2*pi)
r40<-(m4-peakRef)%%(2*pi)
r50<-(m5-peakRef)%%(2*pi)
r60<-(m6-peakRef)%%(2*pi)
r70<-(m7-peakRef)%%(2*pi)
r80<-(m8-peakRef)%%(2*pi)
peakRot<-c()
peakReg<-c()
genesReg2<-c();genesReg3<-c();genesReg4<-c();genesReg1<-c()
genesReg5<-c();genesReg6<-c();genesReg7<-c();genesReg8<-c()
for(i in 1:length(cosPeakTissueAll)){
	peakRot[i]<-(cosPeakTissueAll[i]-peakRef)%%(2*pi)
	if(peakRot[i]>r20){
		peakReg[i]<-1
		genesReg1<-c(genesReg1,genesR2Tissue[i])
	}
	if(peakRot[i]>r30 & peakRot[i]<r20 ){
		peakReg[i]<-2
		genesReg2<-c(genesReg2,genesR2Tissue[i])
	}
	if(peakRot[i]>r40 & peakRot[i]<r30 ){
		peakReg[i]<-3
		genesReg3<-c(genesReg3,genesR2Tissue[i])
	}
	if(peakRot[i]>r50 & peakRot[i]<r40 ){
		peakReg[i]<-4
		genesReg4<-c(genesReg4,genesR2Tissue[i])
	}
	if(peakRot[i]>r60 & peakRot[i]<r50 ){
		peakReg[i]<-5
		genesReg5<-c(genesReg5,genesR2Tissue[i])
	}
	if(peakRot[i]>r70 & peakRot[i]<r60 ){
		peakReg[i]<-6
		genesReg6<-c(genesReg6,genesR2Tissue[i])
	}
	if(peakRot[i]>r80 & peakRot[i]<r70 ){
		peakReg[i]<-7
		genesReg7<-c(genesReg7,genesR2Tissue[i])
	}
	if(peakRot[i]<r80 ){
		peakReg[i]<-8
		genesReg8<-c(genesReg8,genesR2Tissue[i])
	}
}
l1<-length(genesReg1)>0;l2<-length(genesReg2)>0;l3<-length(genesReg3)>0;l4<-length(genesReg4)>0
l5<-length(genesReg5)>0;l6<-length(genesReg6)>0;l7<-length(genesReg7)>0;l8<-length(genesReg8)>0

nameDiscard=c()

r2Reg1<-R2Tissue[match(genesReg1,rownames(mFullTissueNorm))]
r2Reg1Ord<-r2Reg1[order(r2Reg1,decreasing=TRUE)]
genesReg1Ord<-genesReg1[order(r2Reg1,decreasing=TRUE)]
in1<-0
rec1<-1
while(l1 & (in1<tam | (in1>=tam & r2Reg1Ord[rec1]<quantile(distR2Tissue)[4]) ) & rec1<=length(genesReg1Ord) & r2Reg1Ord[rec1]>max(quantile(r2Reg1Ord,probs=0.25),0.4) ){
	datParTissue<-mFullTissueNorm[match(genesReg1Ord[rec1],rownames(mFullTissueNorm)),]
	#original#fParTissue<-fitFMM_Par2(datParTissue,escParTissue)
	fParTissue<-fitFMM_Par(datParTissue,escParTissue)
	wTissue<-fParTissue[[4]]
	peakTissue<-compUU(fParTissue[[2]],fParTissue[[3]],fParTissue[[4]])
        r2ParAux<-1-mean((datParTissue-fParTissue[[1]])^2)/mean((datParTissue-mean(datParTissue))^2)
	print(paste(" Analysing gene ",genesReg1Ord[rec1]," R2NP=",round(r2Reg1Ord[rec1],3)," R2Par=",round(r2ParAux,3)," omega=",round(wTissue,3),sep=""))
	if(wTissue>0.1 & r2ParAux > ifelse(quantile(R2ParcoreG)[2]>0.4,0.4,quantile(R2ParcoreG)[2]) ){
		in1<-in1+1
		entra<-entra+1
		print(paste("There are  ",entra," genes",sep=""))
		print(paste("There are  ",in1," genes in reg1",sep=""))
		datosParTissue<-rbind(datosParTissue,datParTissue)
		fCosTissue<-funcionCosinor(datParTissue,escParTissue,length(escParTissue))
		paramCosTissue<-rbind(paramCosTissue,c(fCosTissue[[2]],fCosTissue[[3]],fCosTissue[[5]]))
		adjCosTissue<-rbind(adjCosTissue,fCosTissue[[1]])
		adjParamTissue<-rbind(adjParamTissue,fParTissue[[1]])
		paramParamTissue<-rbind(paramParamTissue,c(fParTissue[[5]],fParTissue[[6]],fParTissue[[2]],fParTissue[[3]],fParTissue[[4]]))
		peaksIn<-c(peaksIn,peakTissue)
		regIn<-c(regIn,1)
		r2ParIn<-c(r2ParIn,r2ParAux)
		namesIn<-c(namesIn,genesReg1Ord[rec1])
		mAux<-cbind(regIn,r2ParIn,peaksIn)
		rownames(mAux)<-namesIn
		write.table(mAux,file=paste("mAux",nameTissue,".txt",sep=""),sep="\t",dec=",",col.names=TRUE,row.names=TRUE)
		save(datosParTissue,file=paste("datosParTissue_",nameTissue,".RData",sep=""))
		save(adjParamTissue,file=paste("adjParamTissue_",nameTissue,".RData",sep=""))
		save(paramParamTissue,file=paste("paramParamTissue_",nameTissue,".RData",sep=""))
		save(adjCosTissue,file=paste("adjCosTissue_",nameTissue,".RData",sep=""))
		save(paramCosTissue,file=paste("paramCosTissue_",nameTissue,".RData",sep=""))
	}else{
		nameDiscard=c(nameDiscard,genesReg1Ord[rec1])
	}
	rec1<-rec1+1
}
#region 2
r2Reg2<-R2Tissue[match(genesReg2,rownames(mFullTissueNorm))]
r2Reg2Ord<-r2Reg2[order(r2Reg2,decreasing=TRUE)]
genesReg2Ord<-genesReg2[order(r2Reg2,decreasing=TRUE)]
in2<-0
rec2<-1
while(l2 & (in2<tam | (in2>=tam & r2Reg2Ord[rec2]<quantile(distR2Tissue)[4]) ) & rec2<=length(genesReg2Ord)  & r2Reg2Ord[rec2]>max(quantile(r2Reg2Ord,probs=0.25),0.4) ){
	datParTissue<-mFullTissueNorm[match(genesReg2Ord[rec2],rownames(mFullTissueNorm)),]
	#original#fParTissue<-fitFMM_Par2(datParTissue,escParTissue)
	fParTissue<-fitFMM_Par(datParTissue,escParTissue)
	wTissue<-fParTissue[[4]]
	peakTissue<-compUU(fParTissue[[2]],fParTissue[[3]],fParTissue[[4]])
        r2ParAux<-1-mean((datParTissue-fParTissue[[1]])^2)/mean((datParTissue-mean(datParTissue))^2)
	print(paste(" Analysing gene ",genesReg2Ord[rec2]," R2NP=",round(r2Reg2Ord[rec2],3)," R2Par=",round(r2ParAux,3),sep=""))
	if(wTissue>0.1 & r2ParAux >ifelse(quantile(R2ParcoreG)[2]>0.4,0.4,quantile(R2ParcoreG)[2])  ){
		in2<-in2+1
		entra<-entra+1
		print(paste("There are  ",entra," genes",sep=""))
		print(paste("There are  ",in2," genes in reg2",sep=""))
		datosParTissue<-rbind(datosParTissue,datParTissue)
		fCosTissue<-funcionCosinor(datParTissue,escParTissue,length(escParTissue))
		paramCosTissue<-rbind(paramCosTissue,c(fCosTissue[[2]],fCosTissue[[3]],fCosTissue[[5]]))
		adjCosTissue<-rbind(adjCosTissue,fCosTissue[[1]])
		adjParamTissue<-rbind(adjParamTissue,fParTissue[[1]])
		paramParamTissue<-rbind(paramParamTissue,c(fParTissue[[5]],fParTissue[[6]],fParTissue[[2]],fParTissue[[3]],fParTissue[[4]]))
		peaksIn<-c(peaksIn,peakTissue)
		regIn<-c(regIn,2)
		r2ParIn<-c(r2ParIn,r2ParAux)
		namesIn<-c(namesIn,genesReg2Ord[rec2])
		mAux<-cbind(regIn,r2ParIn,peaksIn)
		rownames(mAux)<-namesIn
		write.table(mAux,file=paste("mAux",nameTissue,".txt",sep=""),sep="\t",dec=",",col.names=TRUE,row.names=TRUE)
		save(datosParTissue,file=paste("datosParTissue_",nameTissue,".RData",sep=""))
		save(adjParamTissue,file=paste("adjParamTissue_",nameTissue,".RData",sep=""))
		save(paramParamTissue,file=paste("paramParamTissue_",nameTissue,".RData",sep=""))
		save(adjCosTissue,file=paste("adjCosTissue_",nameTissue,".RData",sep=""))
		save(paramCosTissue,file=paste("paramCosTissue_",nameTissue,".RData",sep=""))
	}else{
		nameDiscard=c(nameDiscard,genesReg2Ord[rec2])
	}
	rec2<-rec2+1
	
}

#region 3
r2Reg3<-R2Tissue[match(genesReg3,rownames(mFullTissueNorm))]
r2Reg3Ord<-r2Reg3[order(r2Reg3,decreasing=TRUE)]
genesReg3Ord<-genesReg3[order(r2Reg3,decreasing=TRUE)]
in3<-0
rec3<-1
while(l3 & (in3<tam | (in3>=tam & r2Reg3Ord[rec3]<quantile(distR2Tissue)[4]) ) & rec3<=length(genesReg3Ord)  & r2Reg3Ord[rec3]>max(quantile(r2Reg3Ord,probs=0.25),0.4) ){
	datParTissue<-mFullTissueNorm[match(genesReg3Ord[rec3],rownames(mFullTissueNorm)),]
	#original#fParTissue<-fitFMM_Par2(datParTissue,escParTissue)
	fParTissue<-fitFMM_Par(datParTissue,escParTissue)
	wTissue<-fParTissue[[4]]
	peakTissue<-compUU(fParTissue[[2]],fParTissue[[3]],fParTissue[[4]])
        r2ParAux<-1-mean((datParTissue-fParTissue[[1]])^2)/mean((datParTissue-mean(datParTissue))^2)
	print(paste(" Analysing gene ",genesReg3Ord[rec3]," R2NP=",round(r2Reg3Ord[rec3],3)," R2Par=",round(r2ParAux,3),sep=""))
	if(wTissue>0.1 & r2ParAux >ifelse(quantile(R2ParcoreG)[2]>0.4,0.4,quantile(R2ParcoreG)[2])  ){
		in3<-in3+1
		entra<-entra+1
		print(paste("There are  ",entra," genes",sep=""))
		print(paste("There are  ",in3," genes in reg3",sep=""))
		datosParTissue<-rbind(datosParTissue,datParTissue)
		fCosTissue<-funcionCosinor(datParTissue,escParTissue,length(escParTissue))
		paramCosTissue<-rbind(paramCosTissue,c(fCosTissue[[2]],fCosTissue[[3]],fCosTissue[[5]]))
		adjCosTissue<-rbind(adjCosTissue,fCosTissue[[1]])
		adjParamTissue<-rbind(adjParamTissue,fParTissue[[1]])
		paramParamTissue<-rbind(paramParamTissue,c(fParTissue[[5]],fParTissue[[6]],fParTissue[[2]],fParTissue[[3]],fParTissue[[4]]))
		peaksIn<-c(peaksIn,peakTissue)
		regIn<-c(regIn,3)
		r2ParIn<-c(r2ParIn,r2ParAux)
		namesIn<-c(namesIn,genesReg3Ord[rec3])
		mAux<-cbind(regIn,r2ParIn,peaksIn)
		rownames(mAux)<-namesIn
		write.table(mAux,file=paste("mAux",nameTissue,".txt",sep=""),sep="\t",dec=",",col.names=TRUE,row.names=TRUE)
		save(datosParTissue,file=paste("datosParTissue_",nameTissue,".RData",sep=""))
		save(adjParamTissue,file=paste("adjParamTissue_",nameTissue,".RData",sep=""))
		save(paramParamTissue,file=paste("paramParamTissue_",nameTissue,".RData",sep=""))
		save(adjCosTissue,file=paste("adjCosTissue_",nameTissue,".RData",sep=""))
		save(paramCosTissue,file=paste("paramCosTissue_",nameTissue,".RData",sep=""))
	}else{
		nameDiscard=c(nameDiscard,genesReg3Ord[rec3])
	}
	rec3<-rec3+1
	
}

#region 4
r2Reg4<-R2Tissue[match(genesReg4,rownames(mFullTissueNorm))]
r2Reg4Ord<-r2Reg4[order(r2Reg4,decreasing=TRUE)]
genesReg4Ord<-genesReg4[order(r2Reg4,decreasing=TRUE)]
in4<-0
rec4<-1
while(l4 & (in4<tam | (in4>=tam & r2Reg4Ord[rec4]<quantile(distR2Tissue)[4]) ) & rec4<=length(genesReg4Ord)  & r2Reg4Ord[rec4]>max(quantile(r2Reg4Ord,probs=0.25),0.4) ){
	datParTissue<-mFullTissueNorm[match(genesReg4Ord[rec4],rownames(mFullTissueNorm)),]
	#original#fParTissue<-fitFMM_Par2(datParTissue,escParTissue)
	fParTissue<-fitFMM_Par(datParTissue,escParTissue)
	wTissue<-fParTissue[[4]]
	peakTissue<-compUU(fParTissue[[2]],fParTissue[[3]],fParTissue[[4]])
        r2ParAux<-1-mean((datParTissue-fParTissue[[1]])^2)/mean((datParTissue-mean(datParTissue))^2)
	print(paste(" Analysing gene ",genesReg4Ord[rec4]," R2NP=",round(r2Reg4Ord[rec4],3)," R2Par=",round(r2ParAux,3),sep=""))
	if(wTissue>0.1 & r2ParAux >ifelse(quantile(R2ParcoreG)[2]>0.4,0.4,quantile(R2ParcoreG)[2])  ){
		in4<-in4+1
		entra<-entra+1
		print(paste("There are  ",entra," genes",sep=""))
		print(paste("There are  ",in4," genes in reg4",sep=""))
		datosParTissue<-rbind(datosParTissue,datParTissue)
		fCosTissue<-funcionCosinor(datParTissue,escParTissue,length(escParTissue))
		paramCosTissue<-rbind(paramCosTissue,c(fCosTissue[[2]],fCosTissue[[3]],fCosTissue[[5]]))
		adjCosTissue<-rbind(adjCosTissue,fCosTissue[[1]])
		adjParamTissue<-rbind(adjParamTissue,fParTissue[[1]])
		paramParamTissue<-rbind(paramParamTissue,c(fParTissue[[5]],fParTissue[[6]],fParTissue[[2]],fParTissue[[3]],fParTissue[[4]]))
		peaksIn<-c(peaksIn,peakTissue)
		regIn<-c(regIn,4)
		r2ParIn<-c(r2ParIn,r2ParAux)
		namesIn<-c(namesIn,genesReg4Ord[rec4])
		mAux<-cbind(regIn,r2ParIn,peaksIn)
		rownames(mAux)<-namesIn
		write.table(mAux,file=paste("mAux",nameTissue,".txt",sep=""),sep="\t",dec=",",col.names=TRUE,row.names=TRUE)
		save(datosParTissue,file=paste("datosParTissue_",nameTissue,".RData",sep=""))
		save(adjParamTissue,file=paste("adjParamTissue_",nameTissue,".RData",sep=""))
		save(paramParamTissue,file=paste("paramParamTissue_",nameTissue,".RData",sep=""))
		save(adjCosTissue,file=paste("adjCosTissue_",nameTissue,".RData",sep=""))
		save(paramCosTissue,file=paste("paramCosTissue_",nameTissue,".RData",sep=""))
	}else{
		nameDiscard=c(nameDiscard,genesReg4Ord[rec4])
	}
	rec4<-rec4+1
}


#region 5
r2Reg5<-R2Tissue[match(genesReg5,rownames(mFullTissueNorm))]
r2Reg5Ord<-r2Reg5[order(r2Reg5,decreasing=TRUE)]
genesReg5Ord<-genesReg5[order(r2Reg5,decreasing=TRUE)]
in5<-0
rec5<-1
while(l5 & (in5<tam | (in5>=tam & r2Reg5Ord[rec5]<quantile(distR2Tissue)[4]) ) & rec5<=length(genesReg5Ord)  & r2Reg5Ord[rec5]>max(quantile(r2Reg5Ord,probs=0.25),0.4) ){
	datParTissue<-mFullTissueNorm[match(genesReg5Ord[rec5],rownames(mFullTissueNorm)),]
	#original#fParTissue<-fitFMM_Par2(datParTissue,escParTissue)
	fParTissue<-fitFMM_Par(datParTissue,escParTissue)
	wTissue<-fParTissue[[4]]
	peakTissue<-compUU(fParTissue[[2]],fParTissue[[3]],fParTissue[[4]])
        r2ParAux<-1-mean((datParTissue-fParTissue[[1]])^2)/mean((datParTissue-mean(datParTissue))^2)
	print(paste(" Analysing gene ",genesReg5Ord[rec5]," R2NP=",round(r2Reg5Ord[rec5],3)," R2Par=",round(r2ParAux,3),sep=""))
	if(wTissue>0.1 & r2ParAux >ifelse(quantile(R2ParcoreG)[2]>0.4,0.4,quantile(R2ParcoreG)[2])  ){
		in5<-in5+1
		entra<-entra+1
		print(paste("There are  ",entra," genes",sep=""))
		print(paste("There are  ",in5," genes in reg5",sep=""))
		datosParTissue<-rbind(datosParTissue,datParTissue)
		fCosTissue<-funcionCosinor(datParTissue,escParTissue,length(escParTissue))
		paramCosTissue<-rbind(paramCosTissue,c(fCosTissue[[2]],fCosTissue[[3]],fCosTissue[[5]]))
		adjCosTissue<-rbind(adjCosTissue,fCosTissue[[1]])
		adjParamTissue<-rbind(adjParamTissue,fParTissue[[1]])
		paramParamTissue<-rbind(paramParamTissue,c(fParTissue[[5]],fParTissue[[6]],fParTissue[[2]],fParTissue[[3]],fParTissue[[4]]))
		peaksIn<-c(peaksIn,peakTissue)
		regIn<-c(regIn,5)
		r2ParIn<-c(r2ParIn,r2ParAux)
		namesIn<-c(namesIn,genesReg5Ord[rec5])
		mAux<-cbind(regIn,r2ParIn,peaksIn)
		rownames(mAux)<-namesIn
		write.table(mAux,file=paste("mAux",nameTissue,".txt",sep=""),sep="\t",dec=",",col.names=TRUE,row.names=TRUE)
		save(datosParTissue,file=paste("datosParTissue_",nameTissue,".RData",sep=""))
		save(adjParamTissue,file=paste("adjParamTissue_",nameTissue,".RData",sep=""))
		save(paramParamTissue,file=paste("paramParamTissue_",nameTissue,".RData",sep=""))
		save(adjCosTissue,file=paste("adjCosTissue_",nameTissue,".RData",sep=""))
		save(paramCosTissue,file=paste("paramCosTissue_",nameTissue,".RData",sep=""))
	}else{
		nameDiscard=c(nameDiscard,genesReg5Ord[rec5])
	}
	rec5<-rec5+1
}

#region 6
r2Reg6<-R2Tissue[match(genesReg6,rownames(mFullTissueNorm))]
r2Reg6Ord<-r2Reg6[order(r2Reg6,decreasing=TRUE)]
genesReg6Ord<-genesReg6[order(r2Reg6,decreasing=TRUE)]
in6<-0
rec6<-1
while(l6 & (in6<tam | (in6>=tam & r2Reg6Ord[rec6]<quantile(distR2Tissue)[4]) ) & rec6<=length(genesReg6Ord)  & r2Reg6Ord[rec6]>max(quantile(r2Reg6Ord,probs=0.25),0.4) ){
	datParTissue<-mFullTissueNorm[match(genesReg6Ord[rec6],rownames(mFullTissueNorm)),]
	#original#fParTissue<-fitFMM_Par2(datParTissue,escParTissue)
	fParTissue<-fitFMM_Par(datParTissue,escParTissue)
	wTissue<-fParTissue[[4]]
	peakTissue<-compUU(fParTissue[[2]],fParTissue[[3]],fParTissue[[4]])
        r2ParAux<-1-mean((datParTissue-fParTissue[[1]])^2)/mean((datParTissue-mean(datParTissue))^2)
	print(paste(" Analysing gene ",genesReg6Ord[rec6]," R2NP=",round(r2Reg6Ord[rec6],3)," R2Par=",round(r2ParAux,3),sep=""))
	if(wTissue>0.1 & r2ParAux >ifelse(quantile(R2ParcoreG)[2]>0.4,0.4,quantile(R2ParcoreG)[2])  ){
		in6<-in6+1
		entra<-entra+1
		print(paste("There are  ",entra," genes",sep=""))
		print(paste("There are  ",in6," genes in reg6",sep=""))
		datosParTissue<-rbind(datosParTissue,datParTissue)
		fCosTissue<-funcionCosinor(datParTissue,escParTissue,length(escParTissue))
		paramCosTissue<-rbind(paramCosTissue,c(fCosTissue[[2]],fCosTissue[[3]],fCosTissue[[5]]))
		adjCosTissue<-rbind(adjCosTissue,fCosTissue[[1]])
		adjParamTissue<-rbind(adjParamTissue,fParTissue[[1]])
		paramParamTissue<-rbind(paramParamTissue,c(fParTissue[[5]],fParTissue[[6]],fParTissue[[2]],fParTissue[[3]],fParTissue[[4]]))
		peaksIn<-c(peaksIn,peakTissue)
		regIn<-c(regIn,6)
		r2ParIn<-c(r2ParIn,r2ParAux)
		namesIn<-c(namesIn,genesReg6Ord[rec6])
		mAux<-cbind(regIn,r2ParIn,peaksIn)
		rownames(mAux)<-namesIn
		write.table(mAux,file=paste("mAux",nameTissue,".txt",sep=""),sep="\t",dec=",",col.names=TRUE,row.names=TRUE)
		save(datosParTissue,file=paste("datosParTissue_",nameTissue,".RData",sep=""))
		save(adjParamTissue,file=paste("adjParamTissue_",nameTissue,".RData",sep=""))
		save(paramParamTissue,file=paste("paramParamTissue_",nameTissue,".RData",sep=""))
		save(adjCosTissue,file=paste("adjCosTissue_",nameTissue,".RData",sep=""))
		save(paramCosTissue,file=paste("paramCosTissue_",nameTissue,".RData",sep=""))
	}else{
		nameDiscard=c(nameDiscard,genesReg6Ord[rec6])
	}
	rec6<-rec6+1
}

#region 7
r2Reg7<-R2Tissue[match(genesReg7,rownames(mFullTissueNorm))]
r2Reg7Ord<-r2Reg7[order(r2Reg7,decreasing=TRUE)]
genesReg7Ord<-genesReg7[order(r2Reg7,decreasing=TRUE)]
in7<-0
rec7<-1
while(l7 & (in7<tam | (in7>=tam & r2Reg7Ord[rec7]<quantile(distR2Tissue)[4]) ) & rec7<=length(genesReg7Ord)  & r2Reg7Ord[rec7]>max(quantile(r2Reg7Ord,probs=0.25),0.4) ){
	datParTissue<-mFullTissueNorm[match(genesReg7Ord[rec7],rownames(mFullTissueNorm)),]
	#original#fParTissue<-fitFMM_Par2(datParTissue,escParTissue)
	fParTissue<-fitFMM_Par(datParTissue,escParTissue)
	wTissue<-fParTissue[[4]]
	peakTissue<-compUU(fParTissue[[2]],fParTissue[[3]],fParTissue[[4]])
        r2ParAux<-1-mean((datParTissue-fParTissue[[1]])^2)/mean((datParTissue-mean(datParTissue))^2)
	print(paste(" Analysing gene ",genesReg7Ord[rec7]," R2NP=",round(r2Reg7Ord[rec7],3)," R2Par=",round(r2ParAux,3),sep=""))
	if(wTissue>0.1 & r2ParAux >ifelse(quantile(R2ParcoreG)[2]>0.4,0.4,quantile(R2ParcoreG)[2])  ){
		in7<-in7+1
		entra<-entra+1
		print(paste("There are  ",entra," genes",sep=""))
		print(paste("There are  ",in7," genes in reg7",sep=""))
		datosParTissue<-rbind(datosParTissue,datParTissue)
		fCosTissue<-funcionCosinor(datParTissue,escParTissue,length(escParTissue))
		paramCosTissue<-rbind(paramCosTissue,c(fCosTissue[[2]],fCosTissue[[3]],fCosTissue[[5]]))
		adjCosTissue<-rbind(adjCosTissue,fCosTissue[[1]])
		adjParamTissue<-rbind(adjParamTissue,fParTissue[[1]])
		paramParamTissue<-rbind(paramParamTissue,c(fParTissue[[5]],fParTissue[[6]],fParTissue[[2]],fParTissue[[3]],fParTissue[[4]]))
		peaksIn<-c(peaksIn,peakTissue)
		regIn<-c(regIn,7)
		r2ParIn<-c(r2ParIn,r2ParAux)
		namesIn<-c(namesIn,genesReg7Ord[rec7])
		mAux<-cbind(regIn,r2ParIn,peaksIn)
		rownames(mAux)<-namesIn
		write.table(mAux,file=paste("mAux",nameTissue,".txt",sep=""),sep="\t",dec=",",col.names=TRUE,row.names=TRUE)
		save(datosParTissue,file=paste("datosParTissue_",nameTissue,".RData",sep=""))
		save(adjParamTissue,file=paste("adjParamTissue_",nameTissue,".RData",sep=""))
		save(paramParamTissue,file=paste("paramParamTissue_",nameTissue,".RData",sep=""))
		save(adjCosTissue,file=paste("adjCosTissue_",nameTissue,".RData",sep=""))
		save(paramCosTissue,file=paste("paramCosTissue_",nameTissue,".RData",sep=""))
	}else{
		nameDiscard=c(nameDiscard,genesReg7Ord[rec7])
	}
	rec7<-rec7+1
}

#region 8
r2Reg8<-R2Tissue[match(genesReg8,rownames(mFullTissueNorm))]
r2Reg8Ord<-r2Reg8[order(r2Reg8,decreasing=TRUE)]
genesReg8Ord<-genesReg8[order(r2Reg8,decreasing=TRUE)]
in8<-0
rec8<-1
while(l8 & (in8<tam | (in8>=tam & r2Reg8Ord[rec8]<quantile(distR2Tissue)[4]) ) & rec8<=length(genesReg8Ord)  & r2Reg8Ord[rec8]>max(quantile(r2Reg8Ord,probs=0.25),0.4) ){
	datParTissue<-mFullTissueNorm[match(genesReg8Ord[rec8],rownames(mFullTissueNorm)),]
	#original#fParTissue<-fitFMM_Par2(datParTissue,escParTissue)
	fParTissue<-fitFMM_Par(datParTissue,escParTissue)
	wTissue<-fParTissue[[4]]
	peakTissue<-compUU(fParTissue[[2]],fParTissue[[3]],fParTissue[[4]])
        r2ParAux<-1-mean((datParTissue-fParTissue[[1]])^2)/mean((datParTissue-mean(datParTissue))^2)
	print(paste(" Analysing gene ",genesReg8Ord[rec8]," R2NP=",round(r2Reg8Ord[rec8],3)," R2Par=",round(r2ParAux,3),sep=""))
	if(wTissue>0.1 & r2ParAux >ifelse(quantile(R2ParcoreG)[2]>0.4,0.4,quantile(R2ParcoreG)[2])  ){
		in8<-in8+1
		entra<-entra+1
		print(paste("There are  ",entra," genes",sep=""))
		print(paste("There are  ",in8," genes in reg8",sep=""))
		datosParTissue<-rbind(datosParTissue,datParTissue)
		fCosTissue<-funcionCosinor(datParTissue,escParTissue,length(escParTissue))
		paramCosTissue<-rbind(paramCosTissue,c(fCosTissue[[2]],fCosTissue[[3]],fCosTissue[[5]]))
		adjCosTissue<-rbind(adjCosTissue,fCosTissue[[1]])
		adjParamTissue<-rbind(adjParamTissue,fParTissue[[1]])
		paramParamTissue<-rbind(paramParamTissue,c(fParTissue[[5]],fParTissue[[6]],fParTissue[[2]],fParTissue[[3]],fParTissue[[4]]))
		peaksIn<-c(peaksIn,peakTissue)
		regIn<-c(regIn,8)
		r2ParIn<-c(r2ParIn,r2ParAux)
		namesIn<-c(namesIn,genesReg8Ord[rec8])
		mAux<-cbind(regIn,r2ParIn,peaksIn)
		rownames(mAux)<-namesIn
		write.table(mAux,file=paste("mAux",nameTissue,".txt",sep=""),sep="\t",dec=",",col.names=TRUE,row.names=TRUE)
		save(datosParTissue,file=paste("datosParTissue_",nameTissue,".RData",sep=""))
		save(adjParamTissue,file=paste("adjParamTissue_",nameTissue,".RData",sep=""))
		save(paramParamTissue,file=paste("paramParamTissue_",nameTissue,".RData",sep=""))
		save(adjCosTissue,file=paste("adjCosTissue_",nameTissue,".RData",sep=""))
		save(paramCosTissue,file=paste("paramCosTissue_",nameTissue,".RData",sep=""))
	}else{
		nameDiscard=c(nameDiscard,genesReg8Ord[rec8])
	}
	rec8<-rec8+1
}

rownames(datosParTissue)<-namesIn
m<-cbind(regIn,r2ParIn,peaksIn)
rownames(m)<-namesIn
colnames(m)<-c("Region","R2_Par","Peak_Par")
write.table(m,file=paste("mCand",nameTissue,".txt",sep=""),sep="\t",dec=",",col.names=TRUE,row.names=TRUE)
return(list(namesIn,datosParTissue,adjParamTissue,paramParamTissue,adjCosTissue,paramCosTissue,peaksIn,regIn,r2ParIn,m,mDistR2Tissue,
r10,r20,r30,r40,r50,r60,r70,r80,m1,m2,m3,m4,m5,m6,m7,m8,nameDiscard))
}
#ptosRotIni=c(r10,r20,r30,r40,r50,r60,r70,r80)
#ptosMIni=c(m1,m2,m3,m4,m5,m6,m7,m8)

#gen1="ARNTL"
#gen2="DBP"
#M=top[,giveCPCA[[k]][[1]]]
#o=giveCPCA[[k]][[1]]
#esc=giveCPCA[[k]][[2]]
#genes=rownames(top)

reajustarPlotsNewReferences_vInicial<-function(gen1,gen2,M,o,esc,genes){
mStatistics<-matrix(0,nrow(M),25)
mFitP3N<-matrix(0,nrow(M),ncol(M))
mFitC3N<-matrix(0,nrow(M),ncol(M))
mFitNP3N<-matrix(0,nrow(M),ncol(M))
vvv<-M[match(gen1,rownames(M)),]
#rrr<-M[match(gen2,rownames(M))]
muevoR<-FALSE;muevoL<-FALSE
#gen 1 de referencia ppicar medio
#Fijandonos en el ajuste para metrico 
#por analogia trabajamso con L y K sin restar esc[1] al haver fitFMM_Par2
fittingRef<-fitFMM_Par(vvv,esc)#M_Par2(M[match(P,cyc),o],p)
#original#fittingRef<-fitFMM_Par2(vvv,esc)#M_Par2(M[match(P,cyc),o],p)
#fittingRefL<-fittingRef
#fittingRefL0<-fittingRef
#fittingRefK0<-fittingRef
#fittingRefK<-fittingRef
#fittingRef2<-fittingRefK
eme<-fittingRef[[5]]
aa<-fittingRef[[6]]
al<-fittingRef[[2]]
be<-fittingRef[[3]]
om<-fittingRef[[4]]
UNC<-(((al+2*atan2(1/om*sin(-be/2),cos(-be/2))))%%(2*pi))
alT=(al-UNC+pi)%%(2*pi)
UNC2<-(((alT+2*atan2(1/om*sin(-be/2),cos(-be/2))))%%(2*pi))

#mseOld=sum((vvv-outF2(vvv,c(eme,aa,al,be,om),esc))^2)/length(vvv)
#nuevo punto medio
medio<-pi#(escT[length(escT)]+escT[1])/2#12.5*(escT[length(escT)]-escT[1])/25#

#cuanto me muevo y hacia donde
if(UNC<medio){

escE=(esc-UNC+pi)%%(2*pi)
vvv2<-c(vvv,vvv)
ooo2<-c(o,o)
oPart1<-o[order(escE)]
vpart1<-vvv[order(escE)]
escBis<-escE[order(escE)]
oN3<-oPart1
vN3<-vpart1
pN3<-escBis
indN3<-order(escE)



}else{#el peak dek gen1 esta a la derecha del pi#UNC<medio


escE=(esc-UNC+pi)%%(2*pi)
vvv2<-c(vvv,vvv)
ooo2<-c(o,o)
oPart1<-o[order(escE)]
vpart1<-vvv[order(escE)]
escBis<-escE[order(escE)]
oN3<-oPart1
vN3<-vpart1
pN3<-escBis
indN3<-order(escE)



}

#original#fit3N<-fitFMM_Par2(rev(vN3),rev(pN3))
fit3N<-fitFMM_Par(vvv,esc)#fitFMM_Par(vN3,pN3)
eme3<-fit3N[[5]]
aa3<-fit3N[[6]]
al3<-(fit3N[[2]]-UNC+pi)%%(2*pi)
be3<-fit3N[[3]]
om3<-fit3N[[4]]
UNBis<-(((al3+2*atan2(1/om3*sin(-be3/2),cos(-be3/2))))%%(2*pi))





for(i in 1:nrow(M)){
vi<-M[i,indN3]
vvv<-M[i,]
#x11()
#plot(pN3,vi,type="b")
#lines(pN3,fit3N[[1]],col=4)
#abline(v=pi,col=3)
mseFlat<-sum((vvv-rep(mean(vvv),length(vvv)))^2)/length(vvv)#sum((vi-rep(mean(vi),length(vi)))^2)/length(vi)

#if( i!=match(gen1,genes) | (i==match(gen1,genes) & UNBis2>UNBis))fitAuxN<-fitFMM_Par2(vi,pN3)
#if(i==match(gen1,genes) & UNBis2<UNBis)fitAuxN<-fit3N
#original#if( i!=match(gen1,genes))fitAuxN<-fitFMM_Par2(vi,pN3)
if( i!=match(gen1,genes))fitAuxN<-fitFMM_Par(vvv,esc)#fitFMM_Par(vi,pN3)

if(i==match(gen1,genes))fitAuxN<-fit3N
fitAllN<-fitAuxN[[1]]
al<-(fitAuxN[[2]]-UNC+pi)%%(2*pi);be<-fitAuxN[[3]];om<-fitAuxN[[4]];MM<-fitAuxN[[5]];A<-fitAuxN[[6]]
peakPar<-compUU(al,be,om);troughPar<-compLL(al,be,om);peakRelPar<-peakPar/(2*pi)*100;troughRelPar<-troughPar/(2*pi)*100
sigmaPar<-sum((vvv-fitAllN)^2)/(length(vvv)-5)#sum((vi-fitAllN)^2)/(length(vi)-5)
msePar<-sum((vvv-fitAllN)^2)/length(vvv)#sum((vi-fitAllN)^2)/length(vi)
R2Par<-1-msePar/mseFlat
mEstPar<-c(MM,A,al,be,om,peakPar,troughPar,peakRelPar,troughRelPar,sigmaPar,msePar,R2Par)
mFitP3N[i,]<-fitAllN[order(escE)]

#x11()
#plot(escE[order(escE)],vvv[order(escE)],type="b")
#lines(seq(0,2*pi,length.out=100), outF2(vvv,c(MM,A,al,be,om),seq(0,2*pi,length.out=100)),col=4)
#abline(v=pi,col=3)

#x11()
#plot(escE[order(escE)],vvv[order(escE)],type="b")
#lines(escE[order(escE)], fitAllN[order(escE)],col=4)
#abline(v=pi,col=3)


fCos<-funcionCosinor(vvv,esc,length(esc))#funcionCosinor(vi,pN3,length(pN3))
mFitC3N[i,]<-fCos[[1]]
eme<-fCos[[2]]
a<-fCos[[3]]
phi<-(fCos[[5]]-UNC+pi)%%(2*pi)
mseCos<-sum((vvv-mFitC3N[i,])^2)/length(vvv)#sum((vi-mFitC3N[i,])^2)/length(vi)
R2Cos<-1-mseCos/mseFlat
peakCos<-(-phi)%%(2*pi)#(-fCos[[5]])%%(2*pi)
troughCos<-(pi-phi)%%(2*pi)#(pi-fCos[[5]])%%(2*pi)
peakRelCos<-peakCos/(2*pi)*100
troughRelCos<-troughCos/(2*pi)*100
sigmaCos<-sum((vvv-mFitC3N[i,])^2)/(length(vvv)-3)#sum((vi-mFitC3N[i,])^2)/(length(vi)-3)
statCos<-c(eme,a,phi,peakCos,troughCos,peakRelCos,troughRelCos,sigmaCos,mseCos,R2Cos)
mFitC3N[i,]<-mFitC3N[i,order(escE)]

mFitNP3N[i,]<-function1Local(vi)[[1]]
mseNP<-sum((vi-mFitNP3N[i,])^2)/length(vi)
sigmaNP<-sum((vi-mFitNP3N[i,])^2)/(length(vi)-length(unique(mFitNP3N[i,])))
R2NP<-1-mseNP/mseFlat
statNP<-c(sigmaNP,mseNP,R2NP)

mStatistics[i,]<-c(mEstPar,statCos,statNP)
}
colnames(mStatistics)<-c("M","A","al","be","om","peakPar","troughPar","peakRelPar","troughRelPar","sigmaPar","msePar","R2Par",
"eme","a","phi","peakCos","troughCos","peakRelCos","troughRelCos","sigmaCos","mseCos","R2Cos",
"sigmaNP","mseNP","R2NP")
return(list(oN3,pN3,indN3,mFitP3N,mFitC3N,mFitNP3N,mStatistics))
}



#reFitSamples[[k]]=list(oN3,pN3,indN3,mFitP3N,mFitC3N,mFitNP3N,mStatistics)

plotTop<-function(plotIni,plotFin,repi,esci,mFini,paramFi,nameTissue,r2){
	#x11()
	side<-compPerfSq(plotIni,plotFin)
	par(mfrow=c(side,side))
	#par(mfrow=c(7,8))
	par(mar=c(1,1,1,1))
	esci<-esci[repi,]
	mFini<-mFini[[repi]]
	paramFi<-paramFi[[repi]]
	r2Ord<-order(r2[-c(1:12)],decreasing=TRUE)+12
	ii<-1
	for(j in plotIni:plotFin){
		if(j<=12){
			plot(esci,mFini[j,],type="b",main=paste(nameTissue," k=",repi," ",rownames(mFini)[j]," R2=",round(r2[j],3),sep=""),
				xaxt="n",yaxt="n")
			tes<-seq(esci[1],esci[length(esci)],length.out=100)
			lines(tes,outF2(mFini[j,],paramFi[j,1:5],tes),col=4)
			#if(rownames(mFini)[j]=="PER1")abline(v=pi,col=3)
		}else{
			plot(esci,mFini[r2Ord[ii],],type="b",main=paste(nameTissue," k=",repi," ",rownames(mFini)[r2Ord[ii]]," R2=",round(r2[r2Ord[ii]],3),sep=""),
				xaxt="n",yaxt="n")
			tes<-seq(esci[1],esci[length(esci)],length.out=100)
			lines(tes,outF2(mFini[r2Ord[ii],],paramFi[r2Ord[ii],1:5],tes),col=4)
			#if(rownames(mFini)[r2Ord[ii]]=="PER1")abline(v=pi,col=3)
			ii<-ii+1
		}
	}
}
compPerfSq<-function(i,f){
	t<-length(i:f)
	j<-0
	sq<-0
	while(sq<t){
		j<-j+1
		sq<-j*j
	}
	return(j)
}

#nReps<-1
#top<-topMuscleRefG
#robEstTissue<-robEstMuscleRefG
#nameTissue<-"Muscle"


robustSincro<-function(nReps=3,top,robEstTissue,nameTissue){
robEstMuscle<-robEstTissue
topMuscle<-top
paramK<-array(0,dim=c(nrow(topMuscle),25,nReps))
for(k in 1:nReps){
	paramK[,,k]<-robEstMuscle[[2]][[k]][[7]]
}
#dim(paramK1)
#p1<-paramK[,,1];p2<-paramK[,,2];p3<-paramK[,,3];
#rownames(p1)<-rownames(topMuscle)
#rownames(p2)<-rownames(topMuscle)
#rownames(p3)<-rownames(topMuscle)
#vamos a trabajar aqui con la sincronizacion
sincroGenes<-c("ARNTL","CLOCK","DBP","PER1","PER2","PER3","CRY1","RORA","STAT3")
R2sincroGenes<-matrix(0,nReps,length(sincroGenes))
PEAKsincroGenes<-matrix(0,nReps,length(sincroGenes))
for(k in 1:nReps){
R2sincroGenes[k,]<-paramK[match(sincroGenes,rownames(topMuscle)),12,k]
PEAKsincroGenes[k,]<-paramK[match(sincroGenes,rownames(topMuscle)),6,k]
}
#para lung
#quito t1<-c("ARNTL","DBP", "PER2");t2<-c("ARNTL","DBP", "CRY1")
#pongo t1<-c("STAT3","DBP", "PER2");t2<-c("STAT3","DBP", "CRY1")
t1<-c("ARNTL","DBP", "PER2");t2<-c("ARNTL","DBP", "CRY1");t3<-c("CLOCK","DBP", "PER2");t4<-c("CLOCK","DBP","CRY1")
t5<-c("STAT3","DBP", "PER2");t6<-c("STAT3","DBP", "CRY1");

t7<-c("ARNTL","DBP", "CRY1");t8<-c("ARNTL","PER1", "CRY1");t9<-c("ARNTL","PER2", "CRY1");t10<-c("ARNTL","PER3", "CRY1")
t11<-c("CLOCK","DBP", "CRY1");t12<-c("CLOCK","PER1", "CRY1");t13<-c("CLOCK","PER2", "CRY1");t14<-c("CLOCK","PER3", "CRY1")
t15<-c("STAT3","DBP", "CRY1");t16<-c("STAT3","PER1", "CRY1");t17<-c("STAT3","PER2", "CRY1");t18<-c("STAT3","PER3", "CRY1")


t19<-c("ARNTL","DBP", "RORA");t20<-c("ARNTL","PER2", "RORA");t21<-c("ARNTL","CRY1", "RORA")
t22<-c("CLOCK","DBP", "RORA");t23<-c("CLOCK","PER2", "RORA");t24<-c("CLOCK","CRY1", "RORA")
t25<-c("STAT3","DBP", "RORA");t26<-c("STAT3","PER2", "RORA");t27<-c("STAT3","CRY1", "RORA")
mAllTrip<-rbind(t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,
t21,t22,t23,t24,t25,t26,t27)
mAllTrip<-cbind(mAllTrip,1:nrow(mAllTrip))

dropTrip<-c()
mInvolvedTrip<-mAllTrip
for(i in 1:length(sincroGenes)){
	
	if(median(R2sincroGenes[,i])<0.5){
		
		for(j in 1:nrow(mAllTrip)){
			if(!is.na(match(sincroGenes[i],mAllTrip[j,]))){ 
				if(nrow(R2sincroGenes)==1){
					if(min(R2sincroGenes[,match(mAllTrip[j,-c(match(sincroGenes[i],mAllTrip[j,]),4)],sincroGenes)])<0.4){
						dropTrip<-c(dropTrip,j)
					}
				}else{
					if(min(apply(R2sincroGenes[,match(mAllTrip[j,-c(match(sincroGenes[i],mAllTrip[j,]),4)],sincroGenes)],2,median))<0.4){
						dropTrip<-c(dropTrip,j)
					}
				}
			}
		}
		
	}
}
if(length(dropTrip)>0)mInvolvedTrip<-mAllTrip[-dropTrip,]
#obtenemos los peaks
mPeaksTrip<-matrix(0,nrow(mInvolvedTrip),ncol(mInvolvedTrip)-1)
for(i in 1:nrow(mInvolvedTrip)){
	for(j in 1:(ncol(mInvolvedTrip)-1)){
		mPeaksTrip[i,j]<-(median.circular(PEAKsincroGenes[,match(mInvolvedTrip[i,j],sincroGenes)]))%%(2*pi)
	}
}
if(length(dropTrip)>0){
	mPeaksTrip<-cbind(mPeaksTrip,(1:nrow(mAllTrip))[-dropTrip])
}else{
	mPeaksTrip<-cbind(mPeaksTrip,1:nrow(mAllTrip))
}
mPeaksTripRot<-matrix(0,nrow(mInvolvedTrip),ncol(mInvolvedTrip))
mPeaksTripRot1<-matrix(0,nrow(mInvolvedTrip),ncol(mInvolvedTrip))
mPeaksTripRot2<-matrix(0,nrow(mInvolvedTrip),ncol(mInvolvedTrip))
mPeaksTripEnd<-c();dista<-FALSE
mDistsTripEnd<-c()
for(i in 1:nrow(mPeaksTripRot)){
	mPeaksTripRot1[i,1:3]<-(mPeaksTrip[i,1:3]-mPeaksTrip[i,1])%%(2*pi)
	mPeaksTripRot2[i,1:3]<-(mPeaksTrip[i,1:3]-mPeaksTrip[i,3])%%(2*pi)
	a<-mPeaksTripRot1[i,1]<=mPeaksTripRot1[i,2] & mPeaksTripRot1[i,2]<=mPeaksTripRot1[i,3]
	b<-mPeaksTripRot2[i,1]>=mPeaksTripRot2[i,2] & mPeaksTripRot2[i,2]>=mPeaksTripRot2[i,3]
	if(a | b){
		mPeaksTripEnd<-rbind(mPeaksTripEnd,mPeaksTrip[i,])
		if(a)mPeaksTripRot[i,]<-mPeaksTripRot1[i,]
		if(b)mPeaksTripRot[i,]<-mPeaksTripRot2[i,]
		dista<-TRUE
	}else{
		dista<-FALSE
	}
	if(dista & a)mDistsTripEnd<-rbind(mDistsTripEnd,
		c(mPeaksTripRot[i,2]-mPeaksTripRot[i,1],mPeaksTripRot[i,3]-mPeaksTripRot[i,2],2*pi-mPeaksTripRot[i,3],
			mPeaksTrip[i,ncol(mPeaksTrip)]))
	if(dista & b)mDistsTripEnd<-rbind(mDistsTripEnd,
		c(mPeaksTripRot[i,1]-mPeaksTripRot[i,2],mPeaksTripRot[i,2]-mPeaksTripRot[i,3],2*pi-mPeaksTripRot[i,1],
			mPeaksTrip[i,ncol(mPeaksTrip)]))

}

firstTrip<-mPeaksTripEnd[1,ncol(mPeaksTripEnd)]
setFilas1<-match(1:6,mDistsTripEnd[,ncol(mDistsTripEnd)])[!is.na(match(1:6,mDistsTripEnd[,ncol(mDistsTripEnd)]))]
setFilas2<-match(7:18,mDistsTripEnd[,ncol(mDistsTripEnd)])[!is.na(match(7:18,mDistsTripEnd[,ncol(mDistsTripEnd)]))]
setFilas3<-match(19:27,mDistsTripEnd[,ncol(mDistsTripEnd)])[!is.na(match(19:27,mDistsTripEnd[,ncol(mDistsTripEnd)]))]
saveMean<-0;saveDisp<-0
if(firstTrip<=6  ){
	setFilas<-setFilas1
	if(min(mDistsTripEnd[setFilas,1:3])<0.1)print("Peak distance lower than 0.1 rad")
	subPeaksTrip<-mPeaksTripEnd[setFilas,]
	subDistsTrip<-mDistsTripEnd[setFilas,]
	if(length(setFilas)==1){
		subPeaksTrip<-matrix(subPeaksTrip,1,4)
		subDistsTrip<-matrix(subDistsTrip,1,4)
	}
	for(i in 1:nrow(subDistsTrip)){
		if((median.circular(subDistsTrip[i,1:3]))%%(2*pi)>saveMean & (circ.disp(subDistsTrip[i,1:3])[4])%%(2*pi)>saveDisp){
			saveMean<-(median.circular(subDistsTrip[i,1:3]))%%(2*pi)
			saveDisp<-(circ.disp(subDistsTrip[i,1:3])[4])%%(2*pi)
			tripleta<-subDistsTrip[i,4]
		}	
	}
}else{
	if(firstTrip>6 & firstTrip<=18){
		setFilas<-setFilas2
		if(min(mDistsTripEnd[setFilas,1:3])<0.1)print("Peak distance lower than 0.1 rad")
		subPeaksTrip<-mPeaksTripEnd[setFilas,]
		subDistsTrip<-mDistsTripEnd[setFilas,]
		if(length(setFilas)==1){
			subPeaksTrip<-matrix(subPeaksTrip,1,4)
			subDistsTrip<-matrix(subDistsTrip,1,4)
		}
		for(i in 1:nrow(subDistsTrip)){
			if((median.circular(subDistsTrip[i,1:3]))%%(2*pi)>saveMean & (circ.disp(subDistsTrip[i,1:3])[4])%%(2*pi)>saveDisp){
				saveMean<-(median.circular(subDistsTrip[i,1:3]))%%(2*pi)
				saveDisp<-(circ.disp(subDistsTrip[i,1:3])[4])%%(2*pi)
				tripleta<-subDistsTrip[i,4]
			}	
		}
	}else{
		setFilas<-setFilas3
		if(min(mDistsTripEnd[setFilas,1:3])<0.1)print("Peak distance lower than 0.1 rad")
		subPeaksTrip<-mPeaksTripEnd[setFilas,]
		subDistsTrip<-mDistsTripEnd[setFilas,]
		if(length(setFilas)==1){
			subPeaksTrip<-matrix(subPeaksTrip,1,4)
			subDistsTrip<-matrix(subDistsTrip,1,4)
		}
		for(i in 1:nrow(subDistsTrip)){
			if((median.circular(subDistsTrip[i,1:3]))%%(2*pi)>saveMean & (circ.disp(subDistsTrip[i,1:3])[4])%%(2*pi)>saveDisp){
				saveMean<-(median.circular(subDistsTrip[i,1:3]))%%(2*pi)
				saveDisp<-(circ.disp(subDistsTrip[i,1:3])[4])%%(2*pi)
				tripleta<-subDistsTrip[i,4]
			}	
		}
	}
}
selTrip<-mAllTrip[tripleta,1:3]
tripRef<-matrix(0,nReps,3)
oF<-c();pF<-c();statF<-c();indF<-c()
tabSamples<-matrix(0,nReps*nrow(topMuscle),25)
mF<-list();adjF<-list()
paramF<-list()
for(k in 1:nReps){
	peaksK<-c(paramK[match(selTrip[1],rownames(topMuscle)),6,k],
		paramK[match(selTrip[2],rownames(topMuscle)),6,k],paramK[match(selTrip[3],rownames(topMuscle)),6,k])
	peaksKRot<-(peaksK-peaksK[1])%%(2*pi)
	if(peaksKRot[1]<=peaksKRot[2] & peaksKRot[2]<=peaksKRot[3] ){
		tripRef[k,]<-peaksK
		oF<-rbind(oF,robEstMuscle[[2]][[k]][[1]])
		pF<-rbind(pF,robEstMuscle[[2]][[k]][[2]])
		indF<-rbind(indF,robEstMuscle[[2]][[k]][[3]])
		statF<-rbind(statF,robEstMuscle[[2]][[k]][[7]])
		tabSamples[seq(k,nrow(topMuscle)*nReps,nReps),]<-robEstMuscle[[2]][[k]][[7]]
		adjFAux<-matrix(0,nrow(topMuscle),ncol(topMuscle))
		paramAuxF<-matrix(0,nrow(topMuscle),5)
		for(i in 1:nrow(topMuscle)){
			adjFAux[i,]<-outF2(topMuscle[i,oF[k,]],robEstMuscle[[2]][[k]][[7]][i,1:5],pF[k,])
			paramAuxF[i,]<-robEstMuscle[[2]][[k]][[7]][i,1:5]
		}
		paramF[[k]]<-paramAuxF
		mF[[k]]<-topMuscle[,oF[k,]]
		adjF[[k]]<-adjFAux
		print("Direct order")
	}else{
		tripRef[k,]<-rev(peaksK)
		oF<-rbind(oF,rev(robEstMuscle[[2]][[k]][[1]]))
		pF<-rbind(pF,2*pi-rev(robEstMuscle[[2]][[k]][[2]]))
		indF<-rbind(indF,rev(robEstMuscle[[2]][[k]][[3]]))
		statsAux<-robEstMuscle[[2]][[k]][[7]]
		statsAux[,c(3,4,6,7)]<-(2*pi-statsAux[,c(3,4,6,7)])%%(2*pi)
		statsAux[,8]<-statsAux[,6]/(2*pi)
		statsAux[,9]<-statsAux[,7]/(2*pi)
		statF<-statsAux
		tabSamples[seq(k,nrow(topMuscle)*nReps,nReps),]<-statF
		paramAuxF<-matrix(0,nrow(topMuscle),5)
		adjFAux<-matrix(0,nrow(topMuscle),ncol(topMuscle))
		for(i in 1:nrow(topMuscle)){
			adjFAux[i,]<-outF2(topMuscle[i,oF[k,]],statsAux[i,1:5],pF[k,])
			paramAuxF[i,]<-statsAux[i,1:5]
		}
		paramF[[k]]<-paramAuxF
		mF[[k]]<-topMuscle[,oF[k,]]
		adjF[[k]]<-adjFAux
		print("Inverse order")
	}
	print(selTrip)
}


for(k in 1:nReps){
	ll=1
	newPlot=FALSE
	for(j in 1:nrow(mF[[k]])){
		if(ll==1 | newPlot){
			png(filename=paste(nameTissue,". Rep ",k," Top_",ll,".png",sep=""),width=1200,height=800)
			#x11()
			par(mfrow=c(3,4))
			par(mar=c(1,1,1,1))
			newPlot=FALSE
		}
		plot(pF[k,],mF[[k]][j,],type="b",main=paste(nameTissue,". Rep ",k," ",rownames(topMuscle)[j],sep=""),
			xaxt="n",yaxt="n")
		tes<-seq(pF[k,1],pF[k,length(pF[k,])],length.out=100)
		lines(tes,outF2(topMuscle[j,oF[k,]],paramF[[k]][j,1:5],tes),col=4)
		ll=ll+1
		if(j%%12==0 | j==nrow(mF[[k]])){
			dev.off()
			newPlot=TRUE
		}
	}
}
#poner esto bien y lo de la tabla, dejar una funcion para pintar, comprobar en muscle y kidney, poner nuevos fix en blood y muscle
return(list(tabSamples,mF,adjF,paramF,oF,pF,indF,selTrip,tripRef))
}
#list(tabSamples,reFitSamples,giveCPCA)
robustSincroDBP<-function(nReps=3,top,robEstTissue,nameTissue){
robEstMuscle<-robEstTissue
topMuscle<-top
paramK<-array(0,dim=c(nrow(topMuscle),25,nReps))
for(k in 1:nReps){
	paramK[,,k]<-robEstMuscle[[2]][[k]][[7]]
}
#dim(paramK1)
#p1<-paramK[,,1];p2<-paramK[,,2];p3<-paramK[,,3];
#rownames(p1)<-rownames(topMuscle)
#rownames(p2)<-rownames(topMuscle)
#rownames(p3)<-rownames(topMuscle)
#vamos a trabajar aqui con la sincronizacion
sincroGenes<-c("ARNTL","CLOCK","DBP","PER1","PER2","PER3","CRY1","RORA","STAT3")
R2sincroGenes<-matrix(0,nReps,length(sincroGenes))
PEAKsincroGenes<-matrix(0,nReps,length(sincroGenes))
for(k in 1:nReps){
R2sincroGenes[k,]<-paramK[match(sincroGenes,rownames(topMuscle)),12,k]
PEAKsincroGenes[k,]<-paramK[match(sincroGenes,rownames(topMuscle)),6,k]
}
#para lung
#quito t1<-c("ARNTL","DBP", "PER2");t2<-c("ARNTL","DBP", "CRY1")
#pongo t1<-c("STAT3","DBP", "PER2");t2<-c("STAT3","DBP", "CRY1")
t1<-c("ARNTL","DBP", "PER2");t2<-c("ARNTL","DBP", "CRY1");t3<-c("CLOCK","DBP", "PER2");t4<-c("CLOCK","DBP","CRY1")
t5<-c("STAT3","DBP", "PER2");t6<-c("STAT3","DBP", "CRY1");

t7<-c("ARNTL","DBP", "CRY1");t8<-c("ARNTL","PER1", "CRY1");t9<-c("ARNTL","PER2", "CRY1");t10<-c("ARNTL","PER3", "CRY1")
t11<-c("CLOCK","DBP", "CRY1");t12<-c("CLOCK","PER1", "CRY1");t13<-c("CLOCK","PER2", "CRY1");t14<-c("CLOCK","PER3", "CRY1")
t15<-c("STAT3","DBP", "CRY1");t16<-c("STAT3","PER1", "CRY1");t17<-c("STAT3","PER2", "CRY1");t18<-c("STAT3","PER3", "CRY1")


t19<-c("ARNTL","DBP", "RORA");t20<-c("ARNTL","PER2", "RORA");t21<-c("ARNTL","CRY1", "RORA")
t22<-c("CLOCK","DBP", "RORA");t23<-c("CLOCK","PER2", "RORA");t24<-c("CLOCK","CRY1", "RORA")
t25<-c("STAT3","DBP", "RORA");t26<-c("STAT3","PER2", "RORA");t27<-c("STAT3","CRY1", "RORA")
mAllTrip<-rbind(t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,
t21,t22,t23,t24,t25,t26,t27)
mAllTrip<-cbind(mAllTrip,1:nrow(mAllTrip))

dropTrip<-c()
mInvolvedTrip<-mAllTrip
for(i in 1:length(sincroGenes)){
	
	if(median(R2sincroGenes[,i])<0.5){
		
		for(j in 1:nrow(mAllTrip)){
			if(!is.na(match(sincroGenes[i],mAllTrip[j,]))){ 
				if(nrow(R2sincroGenes)==1){
					if(min(R2sincroGenes[,match(mAllTrip[j,-c(match(sincroGenes[i],mAllTrip[j,]),4)],sincroGenes)])<0.4){
						dropTrip<-c(dropTrip,j)
					}
				}else{
					if(min(apply(R2sincroGenes[,match(mAllTrip[j,-c(match(sincroGenes[i],mAllTrip[j,]),4)],sincroGenes)],2,median))<0.4){
						dropTrip<-c(dropTrip,j)
					}
				}
			}
		}
		
	}
}
if(length(dropTrip)>0)mInvolvedTrip<-mAllTrip[-dropTrip,]
#obtenemos los peaks

oF<-c();pF<-c();statF<-c();indF<-c()
tabSamples<-matrix(0,nReps*nrow(topMuscle),25)
mF<-list();adjF<-list()
paramF<-list()
for(k in 1:nReps){
	peaksK<-c(paramK[match("DBP",rownames(topMuscle)),6,k])
	peaksKRot<-peaksK
	if(peaksK<pi ){
		#tripRef[k,]<-peaksK
		oF<-rbind(oF,robEstMuscle[[2]][[k]][[1]])
		pF<-rbind(pF,robEstMuscle[[2]][[k]][[2]])
		indF<-rbind(indF,robEstMuscle[[2]][[k]][[3]])
		statF<-rbind(statF,robEstMuscle[[2]][[k]][[7]])
		tabSamples[seq(k,nrow(topMuscle)*nReps,nReps),]<-robEstMuscle[[2]][[k]][[7]]
		adjFAux<-matrix(0,nrow(topMuscle),ncol(topMuscle))
		paramAuxF<-matrix(0,nrow(topMuscle),5)
		for(i in 1:nrow(topMuscle)){
			adjFAux[i,]<-outF2(topMuscle[i,oF[k,]],robEstMuscle[[2]][[k]][[7]][i,1:5],pF[k,])
			paramAuxF[i,]<-robEstMuscle[[2]][[k]][[7]][i,1:5]
		}
		paramF[[k]]<-paramAuxF
		mF[[k]]<-topMuscle[,oF[k,]]
		adjF[[k]]<-adjFAux
		print("Direct order")
	}else{
		#tripRef[k,]<-rev(peaksK)
		oF<-rbind(oF,rev(robEstMuscle[[2]][[k]][[1]]))
		pF<-rbind(pF,2*pi-rev(robEstMuscle[[2]][[k]][[2]]))
		indF<-rbind(indF,rev(robEstMuscle[[2]][[k]][[3]]))
		statsAux<-robEstMuscle[[2]][[k]][[7]]
		statsAux[,c(3,4,6,7)]<-(2*pi-statsAux[,c(3,4,6,7)])%%(2*pi)
		statsAux[,8]<-statsAux[,6]/(2*pi)
		statsAux[,9]<-statsAux[,7]/(2*pi)
		statF<-statsAux
		tabSamples[seq(k,nrow(topMuscle)*nReps,nReps),]<-statF
		paramAuxF<-matrix(0,nrow(topMuscle),5)
		adjFAux<-matrix(0,nrow(topMuscle),ncol(topMuscle))
		for(i in 1:nrow(topMuscle)){
			adjFAux[i,]<-outF2(topMuscle[i,oF[k,]],statsAux[i,1:5],pF[k,])
			paramAuxF[i,]<-statsAux[i,1:5]
		}
		paramF[[k]]<-paramAuxF
		mF[[k]]<-topMuscle[,oF[k,]]
		adjF[[k]]<-adjFAux
		print("Inverse order")
	}
}


for(k in 1:nReps){
	ll=1
	newPlot=FALSE
	for(j in 1:nrow(mF[[k]])){
		if(ll==1 | newPlot){
			png(filename=paste(nameTissue,". Rep ",k," Top_",ll,".png",sep=""),width=1200,height=800)
			#x11()
			par(mfrow=c(3,4))
			par(mar=c(1,1,1,1))
			newPlot=FALSE
		}
		plot(pF[k,],mF[[k]][j,],type="b",main=paste(nameTissue,". Rep ",k," ",rownames(topMuscle)[j],sep=""),
			xaxt="n",yaxt="n")
		tes<-seq(pF[k,1],pF[k,length(pF[k,])],length.out=100)
		lines(tes,outF2(topMuscle[j,oF[k,]],paramF[[k]][j,1:5],tes),col=4)
		ll=ll+1
		if(j%%12==0 | j==nrow(mF[[k]])){
			dev.off()
			newPlot=TRUE
		}
	}
}
#poner esto bien y lo de la tabla, dejar una funcion para pintar, comprobar en muscle y kidney, poner nuevos fix en blood y muscle
return(list(tabSamples,mF,adjF,paramF,oF,pF,indF))
}

#nReps=K
#top=topRefG
#robEstTissue=robEstRefG
#nameTissue=nameTissue

#nReps=K
#top=topRefGSkin
#robEstTissue=robEstRefG
#nameTissue="Epidermis"
#coreG=coreG2

robustSincroDBP_v3_cores<-function(nReps=3,top,robEstTissue,nameTissue,coreG){


robEstMuscle<-robEstTissue
topMuscle<-top
paramK<-array(0,dim=c(nrow(topMuscle),25,nReps))
for(k in 1:nReps){
	paramK[,,k]<-robEstMuscle[[2]][[k]][[7]]
}



oF<-c();pF<-c();statF<-c();indF<-c()
tabSamples<-matrix(0,nReps*nrow(topMuscle),25)
mF<-list();adjF<-list()
paramF<-list()
for(k in 1:nReps){

		oF<-rbind(oF,robEstMuscle[[2]][[k]][[1]])
		pF<-rbind(pF,robEstMuscle[[14]][[k]][[2]])
		indF<-rbind(indF,robEstMuscle[[2]][[k]][[3]])
		statF<-rbind(statF,robEstMuscle[[2]][[k]][[7]])
		tabSamples[seq(k,nrow(topMuscle)*nReps,nReps),]<-robEstMuscle[[2]][[k]][[7]]
		adjFAux<-matrix(0,nrow(topMuscle),ncol(topMuscle))
		paramAuxF<-matrix(0,nrow(topMuscle),5)
		for(i in 1:nrow(topMuscle)){
			adjFAux[i,]<-outF2(robEstMuscle[[14]][[k]][[3]][i,],robEstMuscle[[2]][[k]][[7]][i,1:5],pF[k,])
			paramAuxF[i,]<-robEstMuscle[[2]][[k]][[7]][i,1:5]
		}
		paramF[[k]]<-paramAuxF
		mF[[k]]<-robEstMuscle[[14]][[k]][[3]]
		adjF[[k]]<-adjFAux
}

for(k in 1:nReps){
	ll=1
	newPlot=FALSE
	for(j in 1:nrow(mF[[k]])){
		if(ll==1 | newPlot){
			png(filename=paste(nameTissue,". Rep ",k," Top_",ll,".png",sep=""),width=1200,height=800)
			#x11()
			par(mfrow=c(compPerfSq(1,length(coreG)),compPerfSq(1,length(coreG))))
			par(mar=c(1,1,1,1))
			newPlot=FALSE
		}
		plot(pF[k,],robEstMuscle[[14]][[k]][[3]][j,],type="b",main=paste(nameTissue,". Rep ",k," ",rownames(topMuscle)[j],sep=""),
			xaxt="n",yaxt="n")
		#tes<-seq(pF[k,1],pF[k,length(pF[k,])],length.out=length(pF[k,]))
		tes<-seq(0,2*pi,length.out=length(pF[k,]))
		abline(v=pi,col=3)
		lines(tes,fitMob(paramF[[k]][j,1],paramF[[k]][j,2],paramF[[k]][j,3],paramF[[k]][j,4],
				paramF[[k]][j,5],tes),col=4)
		ll=ll+1
		if(j%%12==0 | j==nrow(mF[[k]])){
			dev.off()
			newPlot=TRUE
		}
	}
}
#poner esto bien y lo de la tabla, dejar una funcion para pintar, comprobar en muscle y kidney, poner nuevos fix en blood y muscle
return(list(tabSamples,mF,adjF,paramF,oF,pF,indF))
}



#namesIn=mNamRefGSkin
#datosIn=mCandRefGSkin
#adjIn=mAdjRefGSkin
#parIn=mParRefGSkin
#nameTissue="Epidermis"
#peaksIn=mPeaksRefGSkin
#orIn=iniOrdRefGSkin
#escIn=escRefGSkin
#mFullN=mFullNormRefGSkin
#R2noSincroCores=R2FMMNoSincroRefGSkin
#PEAKnoSincroCores=peakFMMNoSincroRefGSkin
#R2NPTissue=R2_NPRefGSkin
#coreG=coreG2

computeMinSet_v3_cores<-function(namesIn,datosIn,adjIn,parIn,nameTissue,peaksIn,orIn,escIn,mFullN,
			R2noSincroCores,PEAKnoSincroCores,R2NPTissue,coreG){

#obtain R2 for parametric
mseParIn<-c()
mseFlatIn<-c()
R2ParIn<-c()
for(i in 1:nrow(datosIn)){
	mseParIn[i]<-mean((datosIn[i,]-adjIn[i,])^2)
	mseFlatIn[i]<-mean((datosIn[i,]-mean(datosIn[i,]))^2)
	R2ParIn[i]<-1-mseParIn[i]/mseFlatIn[i]
}

###VAMOS A HACER EL PEAK SINCRONIZADO
step1SincroCand<-reajustarPlotsNewReferences_vInicial_v3("ARNTL","DBP",datosIn,orIn,escIn,namesIn,mFullN)
step2SincroCand<-robustSincro2DBP_v3_cores(nReps=1,top=step1SincroCand[[8]],
	param25=step1SincroCand,nameTissue=nameTissue,R2noSincroCores,PEAKnoSincroCores,
	mFullN,step1SincroCand[[2]],step1SincroCand[[3]],coreG)
	

#####
#tabSamples,mF,adjF,paramF,oF,pF,indF,selTrip,tripRef
###
namGenesInMuscle<-step1SincroCand[[9]]
paramParMuscleIn<-step2SincroCand[[4]][[1]][,1:5]
R2ParMuscle<-step2SincroCand[[1]][,12]
peakParMucle<-step2SincroCand[[1]][,6]#peaksIn
mMinCandNormOrd<-step2SincroCand[[2]][[1]]#matriz con ARNTL Y SINCRO
rownames(mMinCandNormOrd)<-namGenesInMuscle
indMinCandOrd<-step2SincroCand[[7]][1,]#indice para sincro la mFllNorm
escMinCandOrd<-step2SincroCand[[6]][1,]#esc para sincro la mFllNorm

#R2ParMuscle,peakParMucle,namGenesInMuscle
peakCosMucle=c()
for(i in 1:nrow(mMinCandNormOrd)){
	peakCosMucle[i]<-(-funcionCosinor(mMinCandNormOrd[i,],escMinCandOrd,length(escMinCandOrd))[[5]])%%(2*pi)
}


png(paste("Whole Cand Par Sincro Peaks in ",nameTissue,".png",sep=""),width=480,height=480)
	par(mfrow=c(1,1))
	par(mar=c(1,1,1,1))
	plot(as.circular(peakParMucle),main=paste("Whole Cand Par Sincro Peaks in ",nameTissue))
	axis.circular(at=circular(peakParMucle), 
              labels=namGenesInMuscle,tcl.text=-0.07,tick=FALSE)
dev.off()

png(paste("Whole Cand Par Sincro Peaks Cos in ",nameTissue,".png",sep=""),width=480,height=480)
	par(mfrow=c(1,1))
	par(mar=c(1,1,1,1))
	plot(as.circular(peakCosMucle),main=paste("Whole Cand Par Sincro Peaks  Cos in ",nameTissue))
	axis.circular(at=circular(peakCosMucle), 
              labels=rownames(mMinCandNormOrd),tcl.text=-0.07,tick=FALSE)
dev.off()

#R2ParMuscle,peakParMucle,namGenesInMuscle
#calculo la mediana


genesR2Tissue<-rownames(mFullN)[which(R2NPTissue>median(R2NPTissue))]
cosPeakTissueAll<-c()
for(i in 1:length(genesR2Tissue)){
	ff<-mFullN[match(genesR2Tissue[i],rownames(mFullN)),step2SincroCand[[7]][1,]]
	cosPeakTissueAll[i]<-(-funcionCosinor(ff,escMinCandOrd,length(escMinCandOrd))[[5]])%%(2*pi)
}
med<-(median.circular(as.circular(cosPeakTissueAll)))%%(2*pi)
#med<-(median.circular(as.circular(peakParMucle)))%%(2*pi)
m1<-med
m2<-(med-pi/4)%%(2*pi)
m3<-(med-pi/2)%%(2*pi)
m4<-(med-3*pi/4)%%(2*pi)
m5<-(med-pi)%%(2*pi)
m6<-(med-5*pi/4)%%(2*pi)
m7<-(med-3*pi/2)%%(2*pi)
m8<-(med-7*pi/4)%%(2*pi)


if(m1>m2){
	peakS1<-peakCosMucle[peakCosMucle<m1 & peakCosMucle>m2]
	R2S1<-R2ParMuscle[peakCosMucle<m1 & peakCosMucle>m2]
	namesS1<-namGenesInMuscle[peakCosMucle<m1 & peakCosMucle>m2]
}
if(m1<m2){
	peakS1<-peakCosMucle[peakCosMucle<m1 | peakCosMucle>m2]
	R2S1<-R2ParMuscle[peakCosMucle<m1 | peakCosMucle>m2]
	namesS1<-namGenesInMuscle[peakCosMucle<m1 | peakCosMucle>m2]
}
if(m2>m3){
	peakS2<-peakCosMucle[peakCosMucle<m2 & peakCosMucle>m3]
	R2S2<-R2ParMuscle[peakCosMucle<m2 & peakCosMucle>m3]
	namesS2<-namGenesInMuscle[peakCosMucle<m2 & peakCosMucle>m3]
}
if(m2<m3){
	peakS2<-peakCosMucle[peakCosMucle<m2 | peakCosMucle>m3]
	R2S2<-R2ParMuscle[peakCosMucle<m2 | peakCosMucle>m3]
	namesS2<-namGenesInMuscle[peakCosMucle<m2 | peakCosMucle>m3]
}
if(m3>m4){
	peakS3<-peakCosMucle[peakCosMucle<m3 & peakCosMucle>m4]
	R2S3<-R2ParMuscle[peakCosMucle<m3 & peakCosMucle>m4]
	namesS3<-namGenesInMuscle[peakCosMucle<m3 & peakCosMucle>m4]
}
if(m3<m4){
	peakS3<-peakCosMucle[peakCosMucle<m3 | peakCosMucle>m4]
	R2S3<-R2ParMuscle[peakCosMucle<m3 | peakCosMucle>m4]
	namesS3<-namGenesInMuscle[peakCosMucle<m3 | peakCosMucle>m4]
}
if(m4>m5){
	peakS4<-peakCosMucle[peakCosMucle<m4 & peakCosMucle>m5]
	R2S4<-R2ParMuscle[peakCosMucle<m4 & peakCosMucle>m5]
	namesS4<-namGenesInMuscle[peakCosMucle<m4 & peakCosMucle>m5]
}
if(m4<m5){
	peakS4<-peakCosMucle[peakCosMucle<m4 | peakCosMucle>m5]
	R2S4<-R2ParMuscle[peakCosMucle<m4 | peakCosMucle>m5]
	namesS4<-namGenesInMuscle[peakCosMucle<m4 | peakCosMucle>m5]
}
if(m5>m6){
	peakS5<-peakCosMucle[peakCosMucle<m5 & peakCosMucle>m6]
	R2S5<-R2ParMuscle[peakCosMucle<m5 & peakCosMucle>m6]
	namesS5<-namGenesInMuscle[peakCosMucle<m5 & peakCosMucle>m6]
}
if(m5<m6){
	peakS5<-peakCosMucle[peakCosMucle<m5 | peakCosMucle>m6]
	R2S5<-R2ParMuscle[peakCosMucle<m5 | peakCosMucle>m6]
	namesS5<-namGenesInMuscle[peakCosMucle<m5 | peakCosMucle>m6]
}
if(m6>m7){
	peakS6<-peakCosMucle[peakCosMucle<m6 & peakCosMucle>m7]
	R2S6<-R2ParMuscle[peakCosMucle<m6 & peakCosMucle>m7]
	namesS6<-namGenesInMuscle[peakCosMucle<m6 & peakCosMucle>m7]
}
if(m6<m7){
	peakS6<-peakCosMucle[peakCosMucle<m6 | peakCosMucle>m7]
	R2S6<-R2ParMuscle[peakCosMucle<m6 | peakCosMucle>m7]
	namesS6<-namGenesInMuscle[peakCosMucle<m6 | peakCosMucle>m7]
}
if(m7>m8){
	peakS7<-peakCosMucle[peakCosMucle<m7 & peakCosMucle>m8]
	R2S7<-R2ParMuscle[peakCosMucle<m7 & peakCosMucle>m8]
	namesS7<-namGenesInMuscle[peakCosMucle<m7 & peakCosMucle>m8]
}
if(m7<m8){
	peakS7<-peakCosMucle[peakCosMucle<m7 | peakCosMucle>m8]
	R2S7<-R2ParMuscle[peakCosMucle<m7 | peakCosMucle>m8]
	namesS7<-namGenesInMuscle[peakCosMucle<m7 | peakCosMucle>m8]
}
if(m8>m1){
	peakS8<-peakCosMucle[peakCosMucle<m8 & peakCosMucle>m1]
	R2S8<-R2ParMuscle[peakCosMucle<m8 & peakCosMucle>m1]
	namesS8<-namGenesInMuscle[peakCosMucle<m8 & peakCosMucle>m1]
}
if(m8<m1){
	peakS8<-peakCosMucle[peakCosMucle<m8 | peakCosMucle>m1]
	R2S8<-R2ParMuscle[peakCosMucle<m8 | peakCosMucle>m1]
	namesS8<-namGenesInMuscle[peakCosMucle<m8 | peakCosMucle>m1]
}
outS1<-reduceSectors_v2(peakS1,R2S1,namesS1,ncol(paramParMuscleIn))
outS2<-reduceSectors_v2(peakS2,R2S2,namesS2,ncol(paramParMuscleIn))
outS3<-reduceSectors_v2(peakS3,R2S3,namesS3,ncol(paramParMuscleIn))
outS4<-reduceSectors_v2(peakS4,R2S4,namesS4,ncol(paramParMuscleIn))
outS5<-reduceSectors_v2(peakS5,R2S5,namesS5,ncol(paramParMuscleIn))
outS6<-reduceSectors_v2(peakS6,R2S6,namesS6,ncol(paramParMuscleIn))
outS7<-reduceSectors_v2(peakS7,R2S7,namesS7,ncol(paramParMuscleIn))
outS8<-reduceSectors_v2(peakS8,R2S8,namesS8,ncol(paramParMuscleIn))
#Si en alguno de los sectores hay ams diferencia que la longitu del arco
#en los sectores involucrados incluir todos lso genes
#esto lo activo solo si fuera necesario para cuando hay mucha separacioin por ejemplo en cell cultured
#if(length(outS8)>0){
#	pmax<-max(outS8[[1]])
#	if(length(outS7)>0){
#		pmin<-min(outS7[[1]])
#		dista<-1-cos(pmax-pmin)
#		if(dista>(1-cos(m8-m7))){
#			outS8<-list(peakS8,R2S8,namesS8)
#			outS7<-list(peakS7,R2S7,namesS7)
#		}
#	}
#}
#if(length(outS7)>0){
#	pmax<-max(outS6[[1]])
#	if(length(outS6)>0){
#		pmin<-min(outS6[[1]])
#		dista<-1-cos(pmax-pmin)
#		if(dista>(1-cos(m7-m6))){
#			outS7<-list(peakS7,R2S7,namesS7)
#			outS6<-list(peakS6,R2S6,namesS6)
#		}
#	}
#}
#if(length(outS6)>0){
#	pmax<-max(outS5[[1]])
#	if(length(outS5)>0){
#		pmin<-min(outS5[[1]])
#		dista<-1-cos(pmax-pmin)
#		if(dista>(1-cos(m6-m5))){
#			outS6<-list(peakS6,R2S6,namesS6)
#			outS5<-list(peakS5,R2S5,namesS5)
#		}
#	}
#}
#if(length(outS5)>0){
#	pmax<-max(outS4[[1]])
#	if(length(outS4)>0){
#		pmin<-min(outS4[[1]])
#		dista<-1-cos(pmax-pmin)
#		if(dista>(1-cos(m5-m4))){
#			outS5<-list(peakS5,R2S5,namesS5)
#			outS4<-list(peakS4,R2S4,namesS4)
#		}
#	}
#}
#if(length(outS4)>0){
#	pmax<-max(outS3[[1]])
#	if(length(outS3)>0){
#		pmin<-min(outS3[[1]])
#		dista<-1-cos(pmax-pmin)
#		if(dista>(1-cos(m4-m3))){
#			outS4<-list(peakS4,R2S4,namesS4)
#			outS3<-list(peakS3,R2S3,namesS3)
#		}
#	}
#}
#if(length(outS3)>0){
#	pmax<-max(outS2[[1]])
#	if(length(outS2)>0){
#		pmin<-min(outS2[[1]])
#		dista<-1-cos(pmax-pmin)
#		if(dista>(1-cos(m3-m2))){
#			outS3<-list(peakS3,R2S3,namesS3)
#			outS2<-list(peakS2,R2S2,namesS2)
#		}
#	}
#}
#if(length(outS2)>0){
#	pmax<-max(outS1[[1]])
#	if(length(outS1)>0){
#		pmin<-min(outS1[[1]])
#		dista<-1-cos(pmax-pmin)
#		if(dista>(1-cos(m2-m1))){
#			outS2<-list(peakS2,R2S2,namesS2)
#			outS1<-list(peakS1,R2S1,namesS1)
#		}
#	}
#}
#if(length(outS1)>0){
#	pmax<-max(outS8[[1]])
#	if(length(outS8)>0){
#		pmin<-min(outS8[[1]])
#		dista<-1-cos(pmax-pmin)
#		if(dista>(1-cos(m1-m8))){
#			outS1<-list(peakS1,R2S1,namesS1)
#			outS8<-list(peakS8,R2S8,namesS8)
#		}
#	}
#}

sectorBelongingMuscle<-c(rep(1,length(outS1[[1]])),rep(2,length(outS2[[1]])),rep(3,length(outS3[[1]])),rep(4,length(outS4[[1]])),
rep(5,length(outS5[[1]])),rep(6,length(outS6[[1]])),rep(7,length(outS7[[1]])),rep(8,length(outS8[[1]])))
peakParRedMuscle<-c(outS1[[1]],outS2[[1]],outS3[[1]],outS4[[1]],outS5[[1]],outS6[[1]],outS7[[1]],outS8[[1]])
R2RedMuscle<-c(outS1[[2]],outS2[[2]],outS3[[2]],outS4[[2]],outS5[[2]],outS6[[2]],outS7[[2]],outS8[[2]])
namesRedMuscle<-c(outS1[[3]],outS2[[3]],outS3[[3]],outS4[[3]],outS5[[3]],outS6[[3]],outS7[[3]],outS8[[3]])
png(paste("Min Set of Par Sincro Peaks in ",nameTissue,".png",sep=""),width=480,height=480)
par(mfrow=c(1,1))
par(mar=c(1,1,1,1))
plot(as.circular(peakParRedMuscle),main=paste("Min Set of Par Sincro Peaks in ",nameTissue,sep=""))
axis.circular(at=circular(peakParRedMuscle), 
	labels=namesRedMuscle,tcl.text=-0.07,tick=FALSE)
axis.circular(at=circular(c(m1,m2,m3,m4,m5,m6,m7,m8)), 
	labels=c("m1","m2","m3","m4","m5","m6","m7","m8"),tcl.text=-0.07,tick=TRUE,col=2)
dev.off()
genesArt=step1SincroCand[[10]]
return(list(namesRedMuscle,peakParRedMuscle,R2RedMuscle,sectorBelongingMuscle,m1,m2,m3,m4,
	step1SincroCand,step2SincroCand,mMinCandNormOrd,indMinCandOrd,escMinCandOrd,m5,m6,m7,m8,peakCosMucle,genesArt))

}

computeMinSet_v2<-function(namesIn,datosIn,adjIn,parIn,nameTissue,peaksIn,orIn,escIn,mFullN,R2noSincroCores,PEAKnoSincroCores,R2NPTissue){

#obtain R2 for parametric
mseParIn<-c()
mseFlatIn<-c()
R2ParIn<-c()
for(i in 1:nrow(datosIn)){
	mseParIn[i]<-mean((datosIn[i,]-adjIn[i,])^2)
	mseFlatIn[i]<-mean((datosIn[i,]-mean(datosIn[i,]))^2)
	R2ParIn[i]<-1-mseParIn[i]/mseFlatIn[i]
}

###VAMOS A HACER EL PEAK SINCRONIZADO
step1SincroCand<-reajustarPlotsNewReferences_vInicial_v2("ARNTL","DBP",datosIn,orIn,escIn,namesIn,mFullN)
#step1SincroCand<-list(oN3,pN3,indN3,mFitP3N,mFitC3N,mFitNP3N,mStatistics,M,genes)
save(step1SincroCand,file=paste("step1SincroCand_",nameTissue,".RData",sep=""))
step2SincroCand<-robustSincro2DBP(nReps=1,top=step1SincroCand[[8]],param25=step1SincroCand,nameTissue=nameTissue,R2noSincroCores,PEAKnoSincroCores,
mFullN,step1SincroCand[[2]],step1SincroCand[[3]])
#step2SincroCand<-list(tabSamples,mF,adjF,paramF,oF,pF,indF,selTrip,tripRef,fitCoresOutList)
save(step1SincroCand,file=paste("step2SincroCand_",nameTissue,".RData",sep=""))
	

#####
#tabSamples,mF,adjF,paramF,oF,pF,indF,selTrip,tripRef
###
namGenesInMuscle<-step1SincroCand[[9]]
paramParMuscleIn<-step2SincroCand[[4]][[1]][,1:5]
R2ParMuscle<-step2SincroCand[[1]][,12]
peakParMucle<-step2SincroCand[[1]][,6]#peaksIn
mMinCandNormOrd<-step2SincroCand[[2]][[1]]#matriz con ARNTL Y SINCRO
rownames(mMinCandNormOrd)<-namGenesInMuscle
indMinCandOrd<-step2SincroCand[[7]][1,]#indice para sincro la mFllNorm
escMinCandOrd<-step2SincroCand[[6]][1,]#esc para sincro la mFllNorm

#R2ParMuscle,peakParMucle,namGenesInMuscle
peakCosMucle=c()
for(i in 1:nrow(mMinCandNormOrd)){
	peakCosMucle[i]<-(-funcionCosinor(mMinCandNormOrd[i,],escMinCandOrd,length(escMinCandOrd))[[5]])%%(2*pi)
}


png(paste("Whole Cand Par Sincro Peaks in ",nameTissue,".png",sep=""),width=480,height=480)
	par(mfrow=c(1,1))
	par(mar=c(1,1,1,1))
	plot(as.circular(peakParMucle),main=paste("Whole Cand Par Sincro Peaks in ",nameTissue))
	axis.circular(at=circular(peakParMucle), 
              labels=namGenesInMuscle,tcl.text=-0.07,tick=FALSE)
dev.off()

png(paste("Whole Cand Par Sincro Peaks Cos in ",nameTissue,".png",sep=""),width=480,height=480)
	par(mfrow=c(1,1))
	par(mar=c(1,1,1,1))
	plot(as.circular(peakCosMucle),main=paste("Whole Cand Par Sincro Peaks  Cos in ",nameTissue))
	axis.circular(at=circular(peakCosMucle), 
              labels=rownames(mMinCandNormOrd),tcl.text=-0.07,tick=FALSE)
dev.off()

#R2ParMuscle,peakParMucle,namGenesInMuscle
#calculo la mediana


genesR2Tissue<-rownames(mFullN)[which(R2NPTissue>0.5)]
cosPeakTissueAll<-c()
for(i in 1:length(genesR2Tissue)){
	ff<-mFullN[match(genesR2Tissue[i],rownames(mFullN)),step2SincroCand[[7]][1,]]
	cosPeakTissueAll[i]<-(-funcionCosinor(ff,escMinCandOrd,length(escMinCandOrd))[[5]])%%(2*pi)
}
med<-(median.circular(as.circular(cosPeakTissueAll)))%%(2*pi)
#med<-(median.circular(as.circular(peakParMucle)))%%(2*pi)
m1<-med
m2<-(med-pi/4)%%(2*pi)
m3<-(med-pi/2)%%(2*pi)
m4<-(med-3*pi/4)%%(2*pi)
m5<-(med-pi)%%(2*pi)
m6<-(med-5*pi/4)%%(2*pi)
m7<-(med-3*pi/2)%%(2*pi)
m8<-(med-7*pi/4)%%(2*pi)


if(m1>m2){
	peakS1<-peakCosMucle[peakCosMucle<m1 & peakCosMucle>m2]
	R2S1<-R2ParMuscle[peakCosMucle<m1 & peakCosMucle>m2]
	namesS1<-namGenesInMuscle[peakCosMucle<m1 & peakCosMucle>m2]
}
if(m1<m2){
	peakS1<-peakCosMucle[peakCosMucle<m1 | peakCosMucle>m2]
	R2S1<-R2ParMuscle[peakCosMucle<m1 | peakCosMucle>m2]
	namesS1<-namGenesInMuscle[peakCosMucle<m1 | peakCosMucle>m2]
}
if(m2>m3){
	peakS2<-peakCosMucle[peakCosMucle<m2 & peakCosMucle>m3]
	R2S2<-R2ParMuscle[peakCosMucle<m2 & peakCosMucle>m3]
	namesS2<-namGenesInMuscle[peakCosMucle<m2 & peakCosMucle>m3]
}
if(m2<m3){
	peakS2<-peakCosMucle[peakCosMucle<m2 | peakCosMucle>m3]
	R2S2<-R2ParMuscle[peakCosMucle<m2 | peakCosMucle>m3]
	namesS2<-namGenesInMuscle[peakCosMucle<m2 | peakCosMucle>m3]
}
if(m3>m4){
	peakS3<-peakCosMucle[peakCosMucle<m3 & peakCosMucle>m4]
	R2S3<-R2ParMuscle[peakCosMucle<m3 & peakCosMucle>m4]
	namesS3<-namGenesInMuscle[peakCosMucle<m3 & peakCosMucle>m4]
}
if(m3<m4){
	peakS3<-peakCosMucle[peakCosMucle<m3 | peakCosMucle>m4]
	R2S3<-R2ParMuscle[peakCosMucle<m3 | peakCosMucle>m4]
	namesS3<-namGenesInMuscle[peakCosMucle<m3 | peakCosMucle>m4]
}
if(m4>m5){
	peakS4<-peakCosMucle[peakCosMucle<m4 & peakCosMucle>m5]
	R2S4<-R2ParMuscle[peakCosMucle<m4 & peakCosMucle>m5]
	namesS4<-namGenesInMuscle[peakCosMucle<m4 & peakCosMucle>m5]
}
if(m4<m5){
	peakS4<-peakCosMucle[peakCosMucle<m4 | peakCosMucle>m5]
	R2S4<-R2ParMuscle[peakCosMucle<m4 | peakCosMucle>m5]
	namesS4<-namGenesInMuscle[peakCosMucle<m4 | peakCosMucle>m5]
}
if(m5>m6){
	peakS5<-peakCosMucle[peakCosMucle<m5 & peakCosMucle>m6]
	R2S5<-R2ParMuscle[peakCosMucle<m5 & peakCosMucle>m6]
	namesS5<-namGenesInMuscle[peakCosMucle<m5 & peakCosMucle>m6]
}
if(m5<m6){
	peakS5<-peakCosMucle[peakCosMucle<m5 | peakCosMucle>m6]
	R2S5<-R2ParMuscle[peakCosMucle<m5 | peakCosMucle>m6]
	namesS5<-namGenesInMuscle[peakCosMucle<m5 | peakCosMucle>m6]
}
if(m6>m7){
	peakS6<-peakCosMucle[peakCosMucle<m6 & peakCosMucle>m7]
	R2S6<-R2ParMuscle[peakCosMucle<m6 & peakCosMucle>m7]
	namesS6<-namGenesInMuscle[peakCosMucle<m6 & peakCosMucle>m7]
}
if(m6<m7){
	peakS6<-peakCosMucle[peakCosMucle<m6 | peakCosMucle>m7]
	R2S6<-R2ParMuscle[peakCosMucle<m6 | peakCosMucle>m7]
	namesS6<-namGenesInMuscle[peakCosMucle<m6 | peakCosMucle>m7]
}
if(m7>m8){
	peakS7<-peakCosMucle[peakCosMucle<m7 & peakCosMucle>m8]
	R2S7<-R2ParMuscle[peakCosMucle<m7 & peakCosMucle>m8]
	namesS7<-namGenesInMuscle[peakCosMucle<m7 & peakCosMucle>m8]
}
if(m7<m8){
	peakS7<-peakCosMucle[peakCosMucle<m7 | peakCosMucle>m8]
	R2S7<-R2ParMuscle[peakCosMucle<m7 | peakCosMucle>m8]
	namesS7<-namGenesInMuscle[peakCosMucle<m7 | peakCosMucle>m8]
}
if(m8>m1){
	peakS8<-peakCosMucle[peakCosMucle<m8 & peakCosMucle>m1]
	R2S8<-R2ParMuscle[peakCosMucle<m8 & peakCosMucle>m1]
	namesS8<-namGenesInMuscle[peakCosMucle<m8 & peakCosMucle>m1]
}
if(m8<m1){
	peakS8<-peakCosMucle[peakCosMucle<m8 | peakCosMucle>m1]
	R2S8<-R2ParMuscle[peakCosMucle<m8 | peakCosMucle>m1]
	namesS8<-namGenesInMuscle[peakCosMucle<m8 | peakCosMucle>m1]
}
outS1<-reduceSectors_v2(peakS1,R2S1,namesS1,ncol(paramParMuscleIn))
outS2<-reduceSectors_v2(peakS2,R2S2,namesS2,ncol(paramParMuscleIn))
outS3<-reduceSectors_v2(peakS3,R2S3,namesS3,ncol(paramParMuscleIn))
outS4<-reduceSectors_v2(peakS4,R2S4,namesS4,ncol(paramParMuscleIn))
outS5<-reduceSectors_v2(peakS5,R2S5,namesS5,ncol(paramParMuscleIn))
outS6<-reduceSectors_v2(peakS6,R2S6,namesS6,ncol(paramParMuscleIn))
outS7<-reduceSectors_v2(peakS7,R2S7,namesS7,ncol(paramParMuscleIn))
outS8<-reduceSectors_v2(peakS8,R2S8,namesS8,ncol(paramParMuscleIn))


sectorBelongingMuscle<-c(rep(1,length(outS1[[1]])),rep(2,length(outS2[[1]])),rep(3,length(outS3[[1]])),rep(4,length(outS4[[1]])),
rep(5,length(outS5[[1]])),rep(6,length(outS6[[1]])),rep(7,length(outS7[[1]])),rep(8,length(outS8[[1]])))
peakParRedMuscle<-c(outS1[[1]],outS2[[1]],outS3[[1]],outS4[[1]],outS5[[1]],outS6[[1]],outS7[[1]],outS8[[1]])
R2RedMuscle<-c(outS1[[2]],outS2[[2]],outS3[[2]],outS4[[2]],outS5[[2]],outS6[[2]],outS7[[2]],outS8[[2]])
namesRedMuscle<-c(outS1[[3]],outS2[[3]],outS3[[3]],outS4[[3]],outS5[[3]],outS6[[3]],outS7[[3]],outS8[[3]])
png(paste("Min Set of Par Sincro Peaks in ",nameTissue,".png",sep=""),width=480,height=480)
par(mfrow=c(1,1))
par(mar=c(1,1,1,1))
plot(as.circular(peakParRedMuscle),main=paste("Min Set of Par Sincro Peaks in ",nameTissue,sep=""))
axis.circular(at=circular(peakParRedMuscle), 
	labels=namesRedMuscle,tcl.text=-0.07,tick=FALSE)
axis.circular(at=circular(c(m1,m2,m3,m4,m5,m6,m7,m8)), 
	labels=c("m1","m2","m3","m4","m5","m6","m7","m8"),tcl.text=-0.07,tick=TRUE,col=2)
dev.off()
genesArt=step1SincroCand[[10]]
return(list(namesRedMuscle,peakParRedMuscle,R2RedMuscle,sectorBelongingMuscle,m1,m2,m3,m4,
	step1SincroCand,step2SincroCand,mMinCandNormOrd,indMinCandOrd,escMinCandOrd,m5,m6,m7,m8,peakCosMucle,genesArt))

}


reajustarPlotsNewReferences_vInicial_v2<-function(gen1,gen2,M,o,esc,genes,mNormInput){
#EL gen1 tiene que estar en la matriz, de no ser asi, lo debemos agnadir
addGenes=c()
if(!is.na(match(gen1,rownames(M))) & !is.na(match(gen2,rownames(M)))){
	vvv<-M[match(gen1,rownames(M)),]
}else{
	if(is.na(match(gen1,rownames(M))) & is.na(match(gen2,rownames(M)))){
		mG1=match(c(gen1),rownames(mNormInput))
		mG2=match(c(gen2),rownames(mNormInput))
		vG1<-mNormInput[mG1,]
		vG2<-mNormInput[mG2,]
		M=rbind(vG1,vG2,M)
		genes=c(gen1,gen2,genes)
		rownames(M)=genes
		vvv<-M[match(gen1,rownames(M)),]
		addGenes=c(addGenes,gen1,gen2)
	}else{
		if(!is.na(match(gen1,rownames(M))) & is.na(match(gen2,rownames(M)))){
			mG2=match(c(gen2),rownames(mNormInput))
			vG2<-mNormInput[mG2,]
			M=rbind(vG2,M)
			genes=c(gen2,genes)
			rownames(M)=genes
			vvv<-M[match(gen1,rownames(M)),]
			addGenes=c(addGenes,gen2)
		}else{
			if(is.na(match(gen1,rownames(M))) & !is.na(match(gen2,rownames(M)))){
				mG1=match(c(gen1),rownames(mNormInput))
				vG1<-mNormInput[mG1,]
				M=rbind(vG1,M)
				genes=c(gen1,genes)
				rownames(M)=genes
				vvv<-M[match(gen1,rownames(M)),]
				addGenes=c(addGenes,gen1)
			}
		}
	}
	
}
mStatistics<-matrix(0,nrow(M),25)
mFitP3N<-matrix(0,nrow(M),ncol(M))
mFitC3N<-matrix(0,nrow(M),ncol(M))
mFitNP3N<-matrix(0,nrow(M),ncol(M))


#rrr<-M[match(gen2,rownames(M))]
muevoR<-FALSE;muevoL<-FALSE
#gen 1 de referencia ppicar medio
#Fijandonos en el ajuste para metrico 
#por analogia trabajamso con L y K sin restar esc[1] al haver fitFMM_Par2
fittingRef<-fitFMM_Par(vvv,esc)#M_Par2(M[match(P,cyc),o],p)
#original#fittingRef<-fitFMM_Par2(vvv,esc)#M_Par2(M[match(P,cyc),o],p)
#fittingRefL<-fittingRef
#fittingRefL0<-fittingRef
#fittingRefK0<-fittingRef
#fittingRefK<-fittingRef
#fittingRef2<-fittingRefK
eme<-fittingRef[[5]]
aa<-fittingRef[[6]]
al<-fittingRef[[2]]
be<-fittingRef[[3]]
om<-fittingRef[[4]]
UNC<-(((al+2*atan2(1/om*sin(-be/2),cos(-be/2))))%%(2*pi))
alT=(al-UNC+pi)%%(2*pi)
UNC2<-(((alT+2*atan2(1/om*sin(-be/2),cos(-be/2))))%%(2*pi))

#mseOld=sum((vvv-outF2(vvv,c(eme,aa,al,be,om),esc))^2)/length(vvv)
#nuevo punto medio
medio<-pi#(escT[length(escT)]+escT[1])/2#12.5*(escT[length(escT)]-escT[1])/25#

#cuanto me muevo y hacia donde
if(UNC<medio){

escE=(esc-UNC+pi)%%(2*pi)
vvv2<-c(vvv,vvv)
ooo2<-c(o,o)
oPart1<-o[order(escE)]
vpart1<-vvv[order(escE)]
escBis<-escE[order(escE)]
oN3<-oPart1
vN3<-vpart1
pN3<-escBis
indN3<-order(escE)



}else{#el peak dek gen1 esta a la derecha del pi#UNC<medio


escE=(esc-UNC+pi)%%(2*pi)
vvv2<-c(vvv,vvv)
ooo2<-c(o,o)
oPart1<-o[order(escE)]
vpart1<-vvv[order(escE)]
escBis<-escE[order(escE)]
oN3<-oPart1
vN3<-vpart1
pN3<-escBis
indN3<-order(escE)



}

#original#fit3N<-fitFMM_Par2(rev(vN3),rev(pN3))
fit3N<-fitFMM_Par(vvv,esc)#fitFMM_Par(vN3,pN3)
eme3<-fit3N[[5]]
aa3<-fit3N[[6]]
al3<-(fit3N[[2]]-UNC+pi)%%(2*pi)
be3<-fit3N[[3]]
om3<-fit3N[[4]]
UNBis<-(((al3+2*atan2(1/om3*sin(-be3/2),cos(-be3/2))))%%(2*pi))





for(i in 1:nrow(M)){
vi<-M[i,indN3]
vvv<-M[i,]
#x11()
#plot(pN3,vi,type="b")
#lines(pN3,fit3N[[1]],col=4)
#abline(v=pi,col=3)
mseFlat<-sum((vvv-rep(mean(vvv),length(vvv)))^2)/length(vvv)#sum((vi-rep(mean(vi),length(vi)))^2)/length(vi)

#if( i!=match(gen1,genes) | (i==match(gen1,genes) & UNBis2>UNBis))fitAuxN<-fitFMM_Par2(vi,pN3)
#if(i==match(gen1,genes) & UNBis2<UNBis)fitAuxN<-fit3N
#original#if( i!=match(gen1,genes))fitAuxN<-fitFMM_Par2(vi,pN3)
if( i!=match(gen1,genes))fitAuxN<-fitFMM_Par(vvv,esc)#fitFMM_Par(vi,pN3)

if(i==match(gen1,genes))fitAuxN<-fit3N
fitAllN<-fitAuxN[[1]]
al<-(fitAuxN[[2]]-UNC+pi)%%(2*pi);be<-fitAuxN[[3]];om<-fitAuxN[[4]];MM<-fitAuxN[[5]];A<-fitAuxN[[6]]
peakPar<-compUU(al,be,om);troughPar<-compLL(al,be,om);peakRelPar<-peakPar/(2*pi)*100;troughRelPar<-troughPar/(2*pi)*100
sigmaPar<-sum((vvv-fitAllN)^2)/(length(vvv)-5)#sum((vi-fitAllN)^2)/(length(vi)-5)
msePar<-sum((vvv-fitAllN)^2)/length(vvv)#sum((vi-fitAllN)^2)/length(vi)
R2Par<-1-msePar/mseFlat
mEstPar<-c(MM,A,al,be,om,peakPar,troughPar,peakRelPar,troughRelPar,sigmaPar,msePar,R2Par)
mFitP3N[i,]<-fitAllN[order(escE)]

#x11()
#plot(escE[order(escE)],vvv[order(escE)],type="b")
#lines(seq(0,2*pi,length.out=100), outF2(vvv,c(MM,A,al,be,om),seq(0,2*pi,length.out=100)),col=4)
#abline(v=pi,col=3)

#x11()
#plot(escE[order(escE)],vvv[order(escE)],type="b")
#lines(escE[order(escE)], fitAllN[order(escE)],col=4)
#abline(v=pi,col=3)


fCos<-funcionCosinor(vvv,esc,length(esc))#funcionCosinor(vi,pN3,length(pN3))
mFitC3N[i,]<-fCos[[1]]
eme<-fCos[[2]]
a<-fCos[[3]]
phi<-(fCos[[5]]-UNC+pi)%%(2*pi)
mseCos<-sum((vvv-mFitC3N[i,])^2)/length(vvv)#sum((vi-mFitC3N[i,])^2)/length(vi)
R2Cos<-1-mseCos/mseFlat
peakCos<-(-phi)%%(2*pi)#(-fCos[[5]])%%(2*pi)
troughCos<-(pi-phi)%%(2*pi)#(pi-fCos[[5]])%%(2*pi)
peakRelCos<-peakCos/(2*pi)*100
troughRelCos<-troughCos/(2*pi)*100
sigmaCos<-sum((vvv-mFitC3N[i,])^2)/(length(vvv)-3)#sum((vi-mFitC3N[i,])^2)/(length(vi)-3)
statCos<-c(eme,a,phi,peakCos,troughCos,peakRelCos,troughRelCos,sigmaCos,mseCos,R2Cos)
mFitC3N[i,]<-mFitC3N[i,order(escE)]

mFitNP3N[i,]<-function1Local(vi)[[1]]
mseNP<-sum((vi-mFitNP3N[i,])^2)/length(vi)
sigmaNP<-sum((vi-mFitNP3N[i,])^2)/(length(vi)-length(unique(mFitNP3N[i,])))
R2NP<-1-mseNP/mseFlat
statNP<-c(sigmaNP,mseNP,R2NP)

mStatistics[i,]<-c(mEstPar,statCos,statNP)
}

colnames(mStatistics)<-c("M","A","al","be","om","peakPar","troughPar","peakRelPar","troughRelPar","sigmaPar","msePar","R2Par",
"eme","a","phi","peakCos","troughCos","peakRelCos","troughRelCos","sigmaCos","mseCos","R2Cos",
"sigmaNP","mseNP","R2NP")
return(list(oN3,pN3,indN3,mFitP3N,mFitC3N,mFitNP3N,mStatistics,M,genes,addGenes))
}


reajustarPlotsNewReferences_vInicial_v3<-function(gen1,gen2,M,o,esc,genes,mNormInput){
#EL gen1 tiene que estar en la matriz, de no ser asi, lo debemos agnadir
addGenes=c()
if(!is.na(match(gen1,rownames(M))) & !is.na(match(gen2,rownames(M)))){
	vvv<-M[match(gen1,rownames(M)),]
}else{
	if(is.na(match(gen1,rownames(M))) & is.na(match(gen2,rownames(M)))){
		mG1=match(c(gen1),rownames(mNormInput))
		mG2=match(c(gen2),rownames(mNormInput))
		vG1<-mNormInput[mG1,]
		vG2<-mNormInput[mG2,]
		M=rbind(vG1,vG2,M)
		genes=c(gen1,gen2,genes)
		rownames(M)=genes
		vvv<-M[match(gen1,rownames(M)),]
		addGenes=c(addGenes,gen1,gen2)
	}else{
		if(!is.na(match(gen1,rownames(M))) & is.na(match(gen2,rownames(M)))){
			mG2=match(c(gen2),rownames(mNormInput))
			vG2<-mNormInput[mG2,]
			M=rbind(vG2,M)
			genes=c(gen2,genes)
			rownames(M)=genes
			vvv<-M[match(gen1,rownames(M)),]
			addGenes=c(addGenes,gen2)
		}else{
			if(is.na(match(gen1,rownames(M))) & !is.na(match(gen2,rownames(M)))){
				mG1=match(c(gen1),rownames(mNormInput))
				vG1<-mNormInput[mG1,]
				M=rbind(vG1,M)
				genes=c(gen1,genes)
				rownames(M)=genes
				vvv<-M[match(gen1,rownames(M)),]
				addGenes=c(addGenes,gen1)
			}
		}
	}
	
}
mStatistics<-matrix(0,nrow(M),25)
mFitP3N<-matrix(0,nrow(M),ncol(M))
mFitC3N<-matrix(0,nrow(M),ncol(M))
mFitNP3N<-matrix(0,nrow(M),ncol(M))

oN3<-1:length(vvv)
pN3<-esc
indN3<-1:length(vvv)

for(i in 1:nrow(M)){
vi<-M[i,]
vvv<-M[i,]
mseFlat<-sum((vvv-rep(mean(vvv),length(vvv)))^2)/length(vvv)
fitAuxN<-fitFMM_Par(vvv,esc)

fitAllN<-fitAuxN[[1]]
al<-(fitAuxN[[2]])%%(2*pi);be<-fitAuxN[[3]];om<-fitAuxN[[4]];MM<-fitAuxN[[5]];A<-fitAuxN[[6]]
peakPar<-compUU(al,be,om);troughPar<-compLL(al,be,om);peakRelPar<-peakPar/(2*pi)*100;troughRelPar<-troughPar/(2*pi)*100
sigmaPar<-sum((vvv-fitAllN)^2)/(length(vvv)-5)#sum((vi-fitAllN)^2)/(length(vi)-5)
msePar<-sum((vvv-fitAllN)^2)/length(vvv)#sum((vi-fitAllN)^2)/length(vi)
R2Par<-1-msePar/mseFlat
mEstPar<-c(MM,A,al,be,om,peakPar,troughPar,peakRelPar,troughRelPar,sigmaPar,msePar,R2Par)
mFitP3N[i,]<-fitAllN



fCos<-funcionCosinor(vvv,esc,length(vvv))#funcionCosinor(vi,pN3,length(pN3))
mFitC3N[i,]<-fCos[[1]]
eme<-fCos[[2]]
a<-fCos[[3]]
phi<-(fCos[[5]])%%(2*pi)
mseCos<-sum((vvv-mFitC3N[i,])^2)/length(vvv)#sum((vi-mFitC3N[i,])^2)/length(vi)
R2Cos<-1-mseCos/mseFlat
peakCos<-(-phi)%%(2*pi)#(-fCos[[5]])%%(2*pi)
troughCos<-(pi-phi)%%(2*pi)#(pi-fCos[[5]])%%(2*pi)
peakRelCos<-peakCos/(2*pi)*100
troughRelCos<-troughCos/(2*pi)*100
sigmaCos<-sum((vvv-mFitC3N[i,])^2)/(length(vvv)-3)#sum((vi-mFitC3N[i,])^2)/(length(vi)-3)
statCos<-c(eme,a,phi,peakCos,troughCos,peakRelCos,troughRelCos,sigmaCos,mseCos,R2Cos)
mFitC3N[i,]<-mFitC3N[i,]

mFitNP3N[i,]<-function1Local(vi)[[1]]
mseNP<-sum((vi-mFitNP3N[i,])^2)/length(vi)
sigmaNP<-sum((vi-mFitNP3N[i,])^2)/(length(vi)-length(unique(mFitNP3N[i,])))
R2NP<-1-mseNP/mseFlat
statNP<-c(sigmaNP,mseNP,R2NP)

mStatistics[i,]<-c(mEstPar,statCos,statNP)
}

colnames(mStatistics)<-c("M","A","al","be","om","peakPar","troughPar","peakRelPar","troughRelPar","sigmaPar","msePar","R2Par",
"eme","a","phi","peakCos","troughCos","peakRelCos","troughRelCos","sigmaCos","mseCos","R2Cos",
"sigmaNP","mseNP","R2NP")
return(list(oN3,pN3,indN3,mFitP3N,mFitC3N,mFitNP3N,mStatistics,M,genes,addGenes))
}




robustSincro2<-function(nReps=1,top,param25,nameTissue,R2CoresNoSincro,peakCoresNoSincro,
mNormInput,ppp3,iii3){
topMuscle<-top
paramK<-matrix(0,nrow(topMuscle),25)

paramK<-param25[[7]]

#vamos a trabajar aqui con la sincronizacion
sincroGenes<-c("ARNTL","CLOCK","DBP","PER1","PER2","PER3","CRY1","RORA","STAT3")
R2sincroGenes<-matrix(0,nReps,length(sincroGenes))
PEAKsincroGenes<-matrix(0,nReps,length(sincroGenes))
k=nReps
R2sincroGenes[k,]<-paramK[match(sincroGenes,rownames(topMuscle)),12]
PEAKsincroGenes[k,]<-paramK[match(sincroGenes,rownames(topMuscle)),6]

#para lung
#quito t1<-c("ARNTL","DBP", "PER2");t2<-c("ARNTL","DBP", "CRY1")
#pongo t1<-c("STAT3","DBP", "PER2");t2<-c("STAT3","DBP", "CRY1")
t1<-c("ARNTL","DBP", "PER2");t2<-c("ARNTL","DBP", "CRY1");t3<-c("CLOCK","DBP", "PER2");t4<-c("CLOCK","DBP","CRY1")
t5<-c("STAT3","DBP", "PER2");t6<-c("STAT3","DBP", "CRY1");

t7<-c("ARNTL","DBP", "CRY1");t8<-c("ARNTL","PER1", "CRY1");t9<-c("ARNTL","PER2", "CRY1");t10<-c("ARNTL","PER3", "CRY1")
t11<-c("CLOCK","DBP", "CRY1");t12<-c("CLOCK","PER1", "CRY1");t13<-c("CLOCK","PER2", "CRY1");t14<-c("CLOCK","PER3", "CRY1")
t15<-c("STAT3","DBP", "CRY1");t16<-c("STAT3","PER1", "CRY1");t17<-c("STAT3","PER2", "CRY1");t18<-c("STAT3","PER3", "CRY1")


t19<-c("ARNTL","DBP", "RORA");t20<-c("ARNTL","PER2", "RORA");t21<-c("ARNTL","CRY1", "RORA")
t22<-c("CLOCK","DBP", "RORA");t23<-c("CLOCK","PER2", "RORA");t24<-c("CLOCK","CRY1", "RORA")
t25<-c("STAT3","DBP", "RORA");t26<-c("STAT3","PER2", "RORA");t27<-c("STAT3","CRY1", "RORA")
mAllTrip<-rbind(t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,
t21,t22,t23,t24,t25,t26,t27)
mAllTrip<-cbind(mAllTrip,1:nrow(mAllTrip))

contCoresOut=0
fitCoresOutList=list()
if(sum(is.na(R2sincroGenes))>0){
	for(i in 1:length(sincroGenes)){
		if(is.na(R2sincroGenes[i])){
			vvi=mNormInput[match(sincroGenes[i],rownames(mNormInput)),iii3]
			fitCoresOut=fitFMM_Par(vvi,ppp3)
			contCoresOut=contCoresOut+1
			fitCoresOutList[[contCoresOut]]=fitCoresOut
			PEAKsincroGenes[i]=compUU(fitCoresOut[[2]],fitCoresOut[[3]],fitCoresOut[[4]])
			r2fmmAuxi=sum((vvi-fitCoresOut[[1]])^2)/length(vvi)
			r2flatAuxi=sum((vvi-mean(vvi))^2)/length(vvi)
			R2sincroGenes[i]=1-r2fmmAuxi/r2flatAuxi
			
		}
	}
}

dropTrip<-c()
mInvolvedTrip<-mAllTrip
for(i in 1:length(sincroGenes)){
	if(!is.na(R2sincroGenes[,i])){
		if(median(R2sincroGenes[,i])<0.5){
		
			for(j in 1:nrow(mAllTrip)){
				if(!is.na(match(sincroGenes[i],mAllTrip[j,]))){ 
					if(sum(is.na(R2sincroGenes[,match(mAllTrip[j,-c(match(sincroGenes[i],mAllTrip[j,]),4)],sincroGenes)]))==0){
						if(min(R2sincroGenes[,match(mAllTrip[j,-c(match(sincroGenes[i],mAllTrip[j,]),4)],sincroGenes)])<0.4){
							dropTrip<-c(dropTrip,j)
						}
					}
				}
			}
		
		}
	}
}
if(length(dropTrip)>0)mInvolvedTrip<-mAllTrip[-dropTrip,]
#obtenemos los peaks

mPeaksTrip<-matrix(0,nrow(mInvolvedTrip),ncol(mInvolvedTrip)-1)
for(i in 1:nrow(mInvolvedTrip)){
	for(j in 1:(ncol(mInvolvedTrip)-1)){
		mPeaksTrip[i,j]<-(median.circular(PEAKsincroGenes[,match(mInvolvedTrip[i,j],sincroGenes)]))%%(2*pi)
	}
}
if(length(dropTrip)>0){
	mPeaksTrip<-cbind(mPeaksTrip,(1:nrow(mAllTrip))[-dropTrip])
}else{
	mPeaksTrip<-cbind(mPeaksTrip,1:nrow(mAllTrip))
}
mPeaksTripRot<-matrix(0,nrow(mInvolvedTrip),ncol(mInvolvedTrip))
mPeaksTripRot1<-matrix(0,nrow(mInvolvedTrip),ncol(mInvolvedTrip))
mPeaksTripRot2<-matrix(0,nrow(mInvolvedTrip),ncol(mInvolvedTrip))
mPeaksTripEnd<-c();dista<-FALSE
mDistsTripEnd<-c()
for(i in 1:nrow(mPeaksTripRot)){
	mPeaksTripRot1[i,1:3]<-(mPeaksTrip[i,1:3]-mPeaksTrip[i,1])%%(2*pi)
	mPeaksTripRot2[i,1:3]<-(mPeaksTrip[i,1:3]-mPeaksTrip[i,3])%%(2*pi)
	a<-mPeaksTripRot1[i,1]<=mPeaksTripRot1[i,2] & mPeaksTripRot1[i,2]<=mPeaksTripRot1[i,3]
	b<-mPeaksTripRot2[i,1]>=mPeaksTripRot2[i,2] & mPeaksTripRot2[i,2]>=mPeaksTripRot2[i,3]
	if(a | b){
		mPeaksTripEnd<-rbind(mPeaksTripEnd,mPeaksTrip[i,])
		if(a)mPeaksTripRot[i,]<-mPeaksTripRot1[i,]
		if(b)mPeaksTripRot[i,]<-mPeaksTripRot2[i,]
		dista<-TRUE
	}else{
		dista<-FALSE
	}
	if(dista & a)mDistsTripEnd<-rbind(mDistsTripEnd,
		c(mPeaksTripRot[i,2]-mPeaksTripRot[i,1],mPeaksTripRot[i,3]-mPeaksTripRot[i,2],2*pi-mPeaksTripRot[i,3],
			mPeaksTrip[i,ncol(mPeaksTrip)]))
	if(dista & b)mDistsTripEnd<-rbind(mDistsTripEnd,
		c(mPeaksTripRot[i,1]-mPeaksTripRot[i,2],mPeaksTripRot[i,2]-mPeaksTripRot[i,3],2*pi-mPeaksTripRot[i,1],
			mPeaksTrip[i,ncol(mPeaksTrip)]))

}
require(CircStats)
firstTrip<-mPeaksTripEnd[1,ncol(mPeaksTripEnd)]
setFilas1<-match(1:6,mDistsTripEnd[,ncol(mDistsTripEnd)])[!is.na(match(1:6,mDistsTripEnd[,ncol(mDistsTripEnd)]))]
setFilas2<-match(7:18,mDistsTripEnd[,ncol(mDistsTripEnd)])[!is.na(match(7:18,mDistsTripEnd[,ncol(mDistsTripEnd)]))]
setFilas3<-match(19:27,mDistsTripEnd[,ncol(mDistsTripEnd)])[!is.na(match(19:27,mDistsTripEnd[,ncol(mDistsTripEnd)]))]
saveMean<-0;saveDisp<-0
if(firstTrip<=6  ){
	setFilas<-setFilas1
	if(min(mDistsTripEnd[setFilas,1:3])<0.1)print("Peak distance lower than 0.1 rad")
	subPeaksTrip<-mPeaksTripEnd[setFilas,]
	subDistsTrip<-mDistsTripEnd[setFilas,]
	if(length(setFilas)==1){
		subPeaksTrip<-matrix(subPeaksTrip,1,4)
		subDistsTrip<-matrix(subDistsTrip,1,4)
	}
	for(i in 1:nrow(subDistsTrip)){
		if((median.circular(subDistsTrip[i,1:3]))%%(2*pi)>saveMean & (circ.disp(subDistsTrip[i,1:3])[4])%%(2*pi)>saveDisp){
			saveMean<-(median.circular(subDistsTrip[i,1:3]))%%(2*pi)
			saveDisp<-(circ.disp(subDistsTrip[i,1:3])[4])%%(2*pi)
			tripleta<-subDistsTrip[i,4]
		}	
	}
}else{
	if(firstTrip>6 & firstTrip<=18){
		setFilas<-setFilas2
		if(min(mDistsTripEnd[setFilas,1:3])<0.1)print("Peak distance lower than 0.1 rad")
		subPeaksTrip<-mPeaksTripEnd[setFilas,]
		subDistsTrip<-mDistsTripEnd[setFilas,]
		if(length(setFilas)==1){
			subPeaksTrip<-matrix(subPeaksTrip,1,4)
			subDistsTrip<-matrix(subDistsTrip,1,4)
		}
		for(i in 1:nrow(subDistsTrip)){
			if((median.circular(subDistsTrip[i,1:3]))%%(2*pi)>saveMean & (circ.disp(subDistsTrip[i,1:3])[4])%%(2*pi)>saveDisp){
				saveMean<-(median.circular(subDistsTrip[i,1:3]))%%(2*pi)
				saveDisp<-(circ.disp(subDistsTrip[i,1:3])[4])%%(2*pi)
				tripleta<-subDistsTrip[i,4]
			}	
		}
	}else{
		setFilas<-setFilas3
		if(min(mDistsTripEnd[setFilas,1:3])<0.1)print("Peak distance lower than 0.1 rad")
		subPeaksTrip<-mPeaksTripEnd[setFilas,]
		subDistsTrip<-mDistsTripEnd[setFilas,]
		if(length(setFilas)==1){
			subPeaksTrip<-matrix(subPeaksTrip,1,4)
			subDistsTrip<-matrix(subDistsTrip,1,4)
		}
		for(i in 1:nrow(subDistsTrip)){
			if((median.circular(subDistsTrip[i,1:3]))%%(2*pi)>saveMean & (circ.disp(subDistsTrip[i,1:3])[4])%%(2*pi)>saveDisp){
				saveMean<-(median.circular(subDistsTrip[i,1:3]))%%(2*pi)
				saveDisp<-(circ.disp(subDistsTrip[i,1:3])[4])%%(2*pi)
				tripleta<-subDistsTrip[i,4]
			}	
		}
	}
}
selTrip<-mAllTrip[tripleta,1:3]
tripRef<-matrix(0,nReps,3)
oF<-c();pF<-c();statF<-c();indF<-c()
tabSamples<-matrix(0,nReps*nrow(topMuscle),25)
mF<-list();adjF<-list()
paramF<-list()
for(k in 1:nReps){
	peaksK<-c(PEAKsincroGenes[match(selTrip[1],sincroGenes)],
		PEAKsincroGenes[match(selTrip[2],sincroGenes)],PEAKsincroGenes[match(selTrip[3],sincroGenes)])
	peaksKRot<-(peaksK-peaksK[1])%%(2*pi)
	if(peaksKRot[1]<=peaksKRot[2] & peaksKRot[2]<=peaksKRot[3] ){
		tripRef[k,]<-peaksK
		oF<-rbind(oF,param25[[1]])
		pF<-rbind(pF,param25[[2]])
		indF<-rbind(indF,param25[[3]])
		statF<-rbind(statF,param25[[7]])
		tabSamples[seq(k,nrow(topMuscle)*nReps,nReps),]<-param25[[7]]
		adjFAux<-matrix(0,nrow(topMuscle),ncol(topMuscle))
		paramAuxF<-matrix(0,nrow(topMuscle),5)
		for(i in 1:nrow(topMuscle)){
			adjFAux[i,]<-outF2(topMuscle[i,indF[k,]],param25[[7]][i,1:5],pF[k,])
			paramAuxF[i,]<-param25[[7]][i,1:5]
		}
		paramF[[k]]<-paramAuxF
		#mF[[k]]<-topMuscle[,oF[k,]]
		mF[[k]]<-topMuscle[,indF[k,]]
		adjF[[k]]<-adjFAux
		print("Direct order")
	}else{
		tripRef[k,]<-rev(peaksK)
		oF<-rbind(oF,rev(param25[[1]]))
		pF<-rbind(pF,2*pi-rev(param25[[2]]))
		indF<-rbind(indF,rev(param25[[3]]))
		statsAux<-param25[[7]]
		statsAux[,c(3,4,6,7)]<-(2*pi-statsAux[,c(3,4,6,7)])%%(2*pi)
		statsAux[,8]<-statsAux[,6]/(2*pi)
		statsAux[,9]<-statsAux[,7]/(2*pi)
		statF<-statsAux
		tabSamples[seq(k,nrow(topMuscle)*nReps,nReps),]<-statF
		paramAuxF<-matrix(0,nrow(topMuscle),5)
		adjFAux<-matrix(0,nrow(topMuscle),ncol(topMuscle))
		for(i in 1:nrow(topMuscle)){
			adjFAux[i,]<-outF2(topMuscle[i,indF[k,]],statsAux[i,1:5],pF[k,])
			paramAuxF[i,]<-statsAux[i,1:5]
		}
		paramF[[k]]<-paramAuxF
		#mF[[k]]<-topMuscle[,oF[k,]]
		mF[[k]]<-topMuscle[,indF[k,]]
		adjF[[k]]<-adjFAux
		print("Inverse order")
	}
	print(selTrip)
}

for(k in 1:nReps){
	#x11()
	#nn=compPerfSq(1,nrow(topMuscle))
	#par(mfrow=c(3,4))
	#par(mar=c(1,1,1,1))
	ll=1
	newPlot=FALSE
	for(j in 1:nrow(topMuscle)){
		if(ll==1 | newPlot){
			png(filename=paste(nameTissue," TopCand ",rownames(topMuscle)[j],"_",ll,".png",sep=""),width=1200,height=800)
			par(mfrow=c(3,4))
			par(mar=c(1,1,1,1))
			newPlot=FALSE
		}
		plot(pF[k,],mF[[k]][j,],type="b",main=paste(nameTissue," TopCand ",rownames(topMuscle)[j],sep=""),
			xaxt="n",yaxt="n")
		tes<-seq(pF[k,1],pF[k,length(pF[k,])],length.out=100)
		lines(tes,outF2(topMuscle[j,indF[k,]],paramF[[k]][j,1:5],tes),col=4)
		ll=ll+1
		if(j%%12==0 | j==nrow(topMuscle)){
			dev.off()
			newPlot=TRUE
		}
		#if(rownames(topMuscle)[j]=="PER1")abline(v=pi,col=3)
	}
}
#poner esto bien y lo de la tabla, dejar una funcion para pintar, comprobar en muscle y kidney, poner nuevos fix en blood y muscle
return(list(tabSamples,mF,adjF,paramF,oF,pF,indF,selTrip,tripRef,fitCoresOutList))
}

#nReps=1
##top=step1SincroCand[[8]]
#param25=step1SincroCand
#nameTissue=nameTissue
#R2CoresNoSincro=R2noSincroCores
#peakCoresNoSincro=PEAKnoSincroCores
#mNormInput=mFullN
#ppp3=step1SincroCand[[2]]
#iii3=step1SincroCand[[3]]

#nReps=1
#top=step1SincroCand[[8]]
#param25=step1SincroCand
#nameTissue=nameTissue
#R2noSincroCores
#PEAKnoSincroCores
#mFullN
#step1SincroCand[[2]]
#step1SincroCand[[3]]


robustSincro2DBP<-function(nReps=1,top,param25,nameTissue,R2CoresNoSincro,peakCoresNoSincro,
mNormInput,ppp3,iii3){
topMuscle<-top
paramK<-matrix(0,nrow(topMuscle),25)

paramK<-param25[[7]]

#vamos a trabajar aqui con la sincronizacion
sincroGenes<-c("ARNTL","CLOCK","DBP","PER1","PER2","PER3","CRY1","RORA","STAT3")
R2sincroGenes<-matrix(0,nReps,length(sincroGenes))
PEAKsincroGenes<-matrix(0,nReps,length(sincroGenes))
k=nReps
R2sincroGenes[k,]<-paramK[match(sincroGenes,rownames(topMuscle)),12]
PEAKsincroGenes[k,]<-paramK[match(sincroGenes,rownames(topMuscle)),6]

#para lung
#quito t1<-c("ARNTL","DBP", "PER2");t2<-c("ARNTL","DBP", "CRY1")
#pongo t1<-c("STAT3","DBP", "PER2");t2<-c("STAT3","DBP", "CRY1")
t1<-c("ARNTL","DBP", "PER2");t2<-c("ARNTL","DBP", "CRY1");t3<-c("CLOCK","DBP", "PER2");t4<-c("CLOCK","DBP","CRY1")
t5<-c("STAT3","DBP", "PER2");t6<-c("STAT3","DBP", "CRY1");

t7<-c("ARNTL","DBP", "CRY1");t8<-c("ARNTL","PER1", "CRY1");t9<-c("ARNTL","PER2", "CRY1");t10<-c("ARNTL","PER3", "CRY1")
t11<-c("CLOCK","DBP", "CRY1");t12<-c("CLOCK","PER1", "CRY1");t13<-c("CLOCK","PER2", "CRY1");t14<-c("CLOCK","PER3", "CRY1")
t15<-c("STAT3","DBP", "CRY1");t16<-c("STAT3","PER1", "CRY1");t17<-c("STAT3","PER2", "CRY1");t18<-c("STAT3","PER3", "CRY1")


t19<-c("ARNTL","DBP", "RORA");t20<-c("ARNTL","PER2", "RORA");t21<-c("ARNTL","CRY1", "RORA")
t22<-c("CLOCK","DBP", "RORA");t23<-c("CLOCK","PER2", "RORA");t24<-c("CLOCK","CRY1", "RORA")
t25<-c("STAT3","DBP", "RORA");t26<-c("STAT3","PER2", "RORA");t27<-c("STAT3","CRY1", "RORA")
mAllTrip<-rbind(t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,
t21,t22,t23,t24,t25,t26,t27)
mAllTrip<-cbind(mAllTrip,1:nrow(mAllTrip))

contCoresOut=0
fitCoresOutList=list()
if(sum(is.na(R2sincroGenes))>0){
	for(i in 1:length(sincroGenes)){
		if(is.na(R2sincroGenes[i])){
			vvi=mNormInput[match(sincroGenes[i],rownames(mNormInput)),iii3]
			fitCoresOut=fitFMM_Par(vvi,ppp3)
			contCoresOut=contCoresOut+1
			fitCoresOutList[[contCoresOut]]=fitCoresOut
			PEAKsincroGenes[i]=compUU(fitCoresOut[[2]],fitCoresOut[[3]],fitCoresOut[[4]])
			r2fmmAuxi=sum((vvi-fitCoresOut[[1]])^2)/length(vvi)
			r2flatAuxi=sum((vvi-mean(vvi))^2)/length(vvi)
			R2sincroGenes[i]=1-r2fmmAuxi/r2flatAuxi
			
		}
	}
}

dropTrip<-c()
mInvolvedTrip<-mAllTrip
for(i in 1:length(sincroGenes)){
	if(!is.na(R2sincroGenes[,i])){
		if(median(R2sincroGenes[,i])<0.5){
		
			for(j in 1:nrow(mAllTrip)){
				if(!is.na(match(sincroGenes[i],mAllTrip[j,]))){ 
					if(sum(is.na(R2sincroGenes[,match(mAllTrip[j,-c(match(sincroGenes[i],mAllTrip[j,]),4)],sincroGenes)]))==0){
						if(min(R2sincroGenes[,match(mAllTrip[j,-c(match(sincroGenes[i],mAllTrip[j,]),4)],sincroGenes)])<0.4){
							dropTrip<-c(dropTrip,j)
						}
					}
				}
			}
		
		}
	}
}
if(length(dropTrip)>0)mInvolvedTrip<-mAllTrip[-dropTrip,]
#obtenemos los peaks
oF<-c();pF<-c();statF<-c();indF<-c()
tabSamples<-matrix(0,nReps*nrow(topMuscle),25)
mF<-list();adjF<-list()
paramF<-list()
for(k in 1:nReps){
	peaksK<-c(PEAKsincroGenes[match("DBP",sincroGenes)])
	peaksKRot<-peaksK
	if(peaksK<pi){
		#tripRef[k,]<-peaksK
		oF<-rbind(oF,param25[[1]])
		pF<-rbind(pF,param25[[2]])
		indF<-rbind(indF,param25[[3]])
		statF<-rbind(statF,param25[[7]])
		tabSamples[seq(k,nrow(topMuscle)*nReps,nReps),]<-param25[[7]]
		adjFAux<-matrix(0,nrow(topMuscle),ncol(topMuscle))
		paramAuxF<-matrix(0,nrow(topMuscle),5)
		for(i in 1:nrow(topMuscle)){
			adjFAux[i,]<-outF2(topMuscle[i,indF[k,]],param25[[7]][i,1:5],pF[k,])
			paramAuxF[i,]<-param25[[7]][i,1:5]
		}
		paramF[[k]]<-paramAuxF
		#mF[[k]]<-topMuscle[,oF[k,]]
		mF[[k]]<-topMuscle[,indF[k,]]
		adjF[[k]]<-adjFAux
		print("Direct order")
	}else{
		#tripRef[k,]<-rev(peaksK)
		oF<-rbind(oF,rev(param25[[1]]))
		pF<-rbind(pF,2*pi-rev(param25[[2]]))
		indF<-rbind(indF,rev(param25[[3]]))
		statsAux<-param25[[7]]
		statsAux[,c(3,4,6,7)]<-(2*pi-statsAux[,c(3,4,6,7)])%%(2*pi)
		statsAux[,8]<-statsAux[,6]/(2*pi)
		statsAux[,9]<-statsAux[,7]/(2*pi)
		statF<-statsAux
		tabSamples[seq(k,nrow(topMuscle)*nReps,nReps),]<-statF
		paramAuxF<-matrix(0,nrow(topMuscle),5)
		adjFAux<-matrix(0,nrow(topMuscle),ncol(topMuscle))
		for(i in 1:nrow(topMuscle)){
			adjFAux[i,]<-outF2(topMuscle[i,indF[k,]],statsAux[i,1:5],pF[k,])
			paramAuxF[i,]<-statsAux[i,1:5]
		}
		paramF[[k]]<-paramAuxF
		#mF[[k]]<-topMuscle[,oF[k,]]
		mF[[k]]<-topMuscle[,indF[k,]]
		adjF[[k]]<-adjFAux
		print("Inverse order")
	}
}

for(k in 1:nReps){
	#x11()
	#nn=compPerfSq(1,nrow(topMuscle))
	#par(mfrow=c(3,4))
	#par(mar=c(1,1,1,1))
	ll=1
	newPlot=FALSE
	for(j in 1:nrow(topMuscle)){
		if(ll==1 | newPlot){
			png(filename=paste(nameTissue," TopCand ",rownames(topMuscle)[j],"_",ll,".png",sep=""),width=1200,height=800)
			par(mfrow=c(3,4))
			par(mar=c(1,1,1,1))
			newPlot=FALSE
		}
		plot(pF[k,],mF[[k]][j,],type="b",main=paste(nameTissue," TopCand ",rownames(topMuscle)[j],sep=""),
			xaxt="n",yaxt="n")
		tes<-seq(pF[k,1],pF[k,length(pF[k,])],length.out=100)
		lines(tes,outF2(topMuscle[j,indF[k,]],paramF[[k]][j,1:5],tes),col=4)
		ll=ll+1
		if(j%%12==0 | j==nrow(topMuscle)){
			dev.off()
			newPlot=TRUE
		}
	}
}
#poner esto bien y lo de la tabla, dejar una funcion para pintar, comprobar en muscle y kidney, poner nuevos fix en blood y muscle
return(list(tabSamples,mF,adjF,paramF,oF,pF,indF,fitCoresOutList))
}





robustSincro2DBP_v3_cores<-function(nReps=1,top,param25,nameTissue,R2CoresNoSincro,peakCoresNoSincro,
mNormInput,ppp3,iii3,coreG){
topMuscle<-top
paramK<-matrix(0,nrow(topMuscle),25)

paramK<-param25[[7]]



#obtenemos los peaks
oF<-c();pF<-c();statF<-c();indF<-c()
tabSamples<-matrix(0,nReps*nrow(topMuscle),25)
mF<-list();adjF<-list()
paramF<-list()
for(k in 1:nReps){

		#tripRef[k,]<-peaksK
		oF<-rbind(oF,param25[[1]])
		pF<-rbind(pF,param25[[2]])
		indF<-rbind(indF,param25[[3]])
		statF<-rbind(statF,param25[[7]])
		tabSamples[seq(k,nrow(topMuscle)*nReps,nReps),]<-param25[[7]]
		adjFAux<-matrix(0,nrow(topMuscle),ncol(topMuscle))
		paramAuxF<-matrix(0,nrow(topMuscle),5)
		for(i in 1:nrow(topMuscle)){
			adjFAux[i,]<-outF2(topMuscle[i,indF[k,]],param25[[7]][i,1:5],pF[k,])
			paramAuxF[i,]<-param25[[7]][i,1:5]
		}
		paramF[[k]]<-paramAuxF
		#mF[[k]]<-topMuscle[,oF[k,]]
		mF[[k]]<-topMuscle[,indF[k,]]
		adjF[[k]]<-adjFAux
		print("Direct order")
	
}

for(k in 1:nReps){

	ll=1
	newPlot=FALSE
	for(j in 1:nrow(topMuscle)){
		if(ll==1 | newPlot){
			png(filename=paste(nameTissue," TopCand ",rownames(topMuscle)[j],"_",ll,".png",sep=""),width=1200,height=800)
			par(mfrow=c(compPerfSq(1,length(coreG)),compPerfSq(1,length(coreG))))
			par(mar=c(1,1,1,1))
			newPlot=FALSE
		}
		plot(pF[k,],mF[[k]][j,],type="b",main=paste(nameTissue," TopCand ",rownames(topMuscle)[j],sep=""),
			xaxt="n",yaxt="n")
		tes<-seq(pF[k,1],pF[k,length(pF[k,])],length.out=100)
		lines(tes,outF2(topMuscle[j,indF[k,]],paramF[[k]][j,1:5],tes),col=4)
		ll=ll+1
		if(j%%12==0 | j==nrow(topMuscle)){
			dev.off()
			newPlot=TRUE
		}
	}
}
#poner esto bien y lo de la tabla, dejar una funcion para pintar, comprobar en muscle y kidney, poner nuevos fix en blood y muscle
return(list(tabSamples,mF,adjF,paramF,oF,pF,indF))
}





#funcion para leer muestras por tejido, sexo, raza y edad
#dbCyc: dbCyc<-db[match(ensgSel,rowNames1),] 
#unique(covDB[,"SMTSD"])
#mIn<-covDB
#tisNam<-ageMaleTissues[[5]][i]
#sexG=1
#ageG=5
#raceG=3
#mIn=covDB
#tisNam=ageMaleTissues[[ageGroup]][numTissue]
#sexG=1
#ageG=ageGroup
#raceG=3
	
giveSampleId<-function(mIn,tisNam,sexG,ageG,raceG){
	selTis<-which(mIn[,match("SMTSD",colnames(mIn))]==tisNam,arr.ind=TRUE)
	selSex<-which(mIn[,match("SEX",colnames(mIn))]==sexG,arr.ind=TRUE)
	selAge<-which(mIn[,match("AGE",colnames(mIn))]==ageG,arr.ind=TRUE)
	selRace<-which(mIn[,match("RACE",colnames(mIn))]==raceG,arr.ind=TRUE)
	selRow<-intersect(selTis,intersect(selSex,intersect(selAge,selRace)))
	selSample<-mIn[selRow,1]
	return(selSample)
}

giveSampleId_v2<-function(mIn,tisNam,raceG){
	selTis<-which(mIn[,match("SMTSD",colnames(mIn))]==tisNam,arr.ind=TRUE)
	selRace<-which(mIn[,match("RACE",colnames(mIn))]==raceG,arr.ind=TRUE)
	selRow<-intersect(selTis,selRace)
	selSample<-mIn[selRow,1]
	return(selSample)
}





obtainCPCAMany<-function(mNorm8,tissue,coreG,nOuts,printing){
#CPCA a partir de los 8 Original, non-equispaced
cp8<-prcomp(centrado(mNorm8), scale. = TRUE, center = FALSE)
varPer8<-c(cp8[[1]][1]^2/sum(cp8[[1]]^2),
cp8[[1]][2]^2/sum(cp8[[1]]^2),
cp8[[1]][3]^2/sum(cp8[[1]]^2))#varianza explicada
eigen18<-cp8$rotation[,1]
eigen28<-cp8$rotation[,2]
eigen38<-cp8$rotation[,3]
xi8<-eigen18/sqrt((eigen18^2+ eigen28^2))
yi8<-eigen28/sqrt((eigen18^2+ eigen28^2))
phi8<-(atan2(yi8,xi8))%%(2*pi)
orderCPCA8<-order(phi8)
escalaPhi8<-phi8[orderCPCA8]
rojo8<-FALSE
maxi<-round(max(abs(eigen18),abs(eigen28)),1)
d8<-c()
for(i in 1:length(eigen18)){
d8[i]<-sqrt(eigen18[i]^2+eigen28[i]^2)
}
obs8<-order(d8)[1:nOuts]
if(printing==TRUE)x11(width=15,height=15)
if(printing==TRUE)par(mfrow=c(1,1))
if(printing==TRUE)plot(eigen18,eigen28,xlim=c(-maxi,maxi),ylim=c(-maxi,maxi),
main=paste(tissue,'. Eigen1 (',round(varPer8[1]*100,2),'%) vs Eigen2 (',round(varPer8[2]*100,2),'%). CoreG',sep=''))

s8<-0
for(i in 1:length(obs8)){
if(d8[obs8[i]]<=0.1){
if(printing==TRUE)points(eigen18[obs8[i]],eigen28[obs8[i]],pch=0,col=i)
s8<-s8+1
}
}
if(s8==0){
for(i in 1:length(obs8)){
if(d8[obs8[i]]<=0.15){
if(printing==TRUE)points(eigen18[obs8[i]],eigen28[obs8[i]],pch=0,col=i)
s8<-s8+1
rojo8<-TRUE
}
}
}
ss8<-s8
if(printing==TRUE)lines(cos(secTimes(length(eigen28)))*0.1,sin(secTimes(length(eigen28)))*0.1,lty=2,col=2)
if(printing==TRUE)lines(cos(secTimes(length(eigen28)))*0.15,sin(secTimes(length(eigen28)))*0.15,lty=2)
if(rojo8 & printing==TRUE)lines(printing==TRUE & cos(secTimes(length(eigen28)))*0.15,sin(secTimes(length(eigen28)))*0.15,lty=2,col=2)
if(printing==TRUE)lines(cos(secTimes(length(eigen28)))*0.2,sin(secTimes(length(eigen28)))*0.2,lty=2)
if(printing==TRUE & maxi>=0.3)lines(cos(secTimes(length(eigen28)))*0.3,sin(secTimes(length(eigen28)))*0.3,lty=2)
if(printing==TRUE & maxi>=0.4)lines(cos(secTimes(length(eigen28)))*0.4,sin(secTimes(length(eigen28)))*0.4,lty=2)


mOutliers<-matrix(c(obs8,ss8),1,nOuts+1)
rownames(mOutliers)<-c("orderCPCACoreG")
colnames(mOutliers)<-c(paste0("obsOut",1:nOuts),"nOuts")

tableOuts<-mOutliers
tableOrders<-orderCPCA8

outs<-match((tableOuts)[1,1:(tableOuts)[1,nOuts+1]],tableOrders)


if(printing==TRUE)x11(width=16,height=12)
if(printing==TRUE)par(mfrow=c(1,3))
if(printing==TRUE)par(mar=c(1,1,1,1))
#coreG genes con el orden CPCA54 y outliers para el orden CPCA54
if(printing==TRUE)plot(escalaPhi8,eigen18[tableOrders],type="b",xaxt="n",yaxt="n",
main=paste("1st Eigenegene",". %Var= ",round(varPer8[1],3),sep=""))
if(printing==TRUE)points(escalaPhi8[outs],(eigen18[tableOrders])[outs],col=c(1:(tableOuts)[1,nOuts+1]),pch=0,cex=1.4)
if(printing==TRUE)plot(escalaPhi8,eigen28[tableOrders],type="b",xaxt="n",yaxt="n",
main=paste("2nd Eigenegene",". %Var= ",round(varPer8[2],3),sep=""))
if(printing==TRUE)points(escalaPhi8[outs],(eigen28[tableOrders])[outs],col=c(1:(tableOuts)[1,nOuts+1]),pch=0,cex=1.4)
if(printing==TRUE)plot(escalaPhi8,eigen38[tableOrders],type="b",xaxt="n",yaxt="n",
main=paste("3rd Eigenegene",". %Var= ",round(varPer8[3],3),sep=""))
if(printing==TRUE)points(escalaPhi8[outs],(eigen38[tableOrders])[outs],col=c(1:(tableOuts)[1,nOuts+1]),pch=0,cex=1.4)

return(list(orderCPCA8,escalaPhi8,mOutliers,outs,varPer8,eigen18,eigen28,eigen38))
}




#mNorm8=mTissueCoreGNorm
#tissue=nameT
#coreG=coreG
#nOuts=8
#printing=FALSE

obtainCPCA12<-function(mNorm8,tissue,coreG,nOuts,printing){
#CPCA a partir de los 8 Original, non-equispaced
cp8<-prcomp(centrado(mNorm8), scale. = TRUE, center = FALSE)
varPer8<-c(cp8[[1]][1]^2/sum(cp8[[1]]^2),
cp8[[1]][2]^2/sum(cp8[[1]]^2),
cp8[[1]][3]^2/sum(cp8[[1]]^2))#varianza explicada
eigen18<-cp8$rotation[,1]
eigen28<-cp8$rotation[,2]
eigen38<-cp8$rotation[,3]
xi8<-eigen18/sqrt((eigen18^2+ eigen28^2))
yi8<-eigen28/sqrt((eigen18^2+ eigen28^2))
phi8<-(atan2(yi8,xi8))%%(2*pi)
orderCPCA8<-order(phi8)
escalaPhi8<-phi8[orderCPCA8]
rojo8<-FALSE
maxi<-round(max(abs(eigen18),abs(eigen28)),1)
d8<-c()
for(i in 1:length(eigen18)){
d8[i]<-sqrt(eigen18[i]^2+eigen28[i]^2)
}
obs8<-order(d8)[1:nOuts]

if(printing==TRUE)x11(width=15,height=15)
if(printing==TRUE)par(mfrow=c(1,1))
if(printing==TRUE)plot(eigen18,eigen28,xlim=c(-maxi,maxi),ylim=c(-maxi,maxi),
main=paste(tissue,'. Eigen1 (',round(varPer8[1]*100,2),'%) vs Eigen2 (',round(varPer8[2]*100,2),'%). CoreG',sep=''))

s8<-0
for(i in 1:length(obs8)){
if(d8[obs8[i]]<=0.1){
if(printing==TRUE)points(eigen18[obs8[i]],eigen28[obs8[i]],pch=0,col=i)
s8<-s8+1
}
}
if(s8==0){
for(i in 1:length(obs8)){
if(d8[obs8[i]]<=0.15){
if(printing==TRUE)points(eigen18[obs8[i]],eigen28[obs8[i]],pch=0,col=i)
s8<-s8+1
rojo8<-TRUE
}
}
}
ss8<-s8
if(printing==TRUE)lines(cos(secTimes(length(eigen28)))*0.1,sin(secTimes(length(eigen28)))*0.1,lty=2,col=2)
if(printing==TRUE)lines(cos(secTimes(length(eigen28)))*0.15,sin(secTimes(length(eigen28)))*0.15,lty=2)
if(rojo8 & printing==TRUE)lines(cos(secTimes(length(eigen28)))*0.15,sin(secTimes(length(eigen28)))*0.15,lty=2,col=2)
if(printing==TRUE)lines(cos(secTimes(length(eigen28)))*0.2,sin(secTimes(length(eigen28)))*0.2,lty=2)
if(maxi>=0.3 & printing==TRUE)lines(cos(secTimes(length(eigen28)))*0.3,sin(secTimes(length(eigen28)))*0.3,lty=2)
if(maxi>=0.4 & printing==TRUE)lines(cos(secTimes(length(eigen28)))*0.4,sin(secTimes(length(eigen28)))*0.4,lty=2)


mOutliers<-matrix(c(obs8,ss8),1,nOuts+1)
rownames(mOutliers)<-c("orderCPCACoreG")
colnames(mOutliers)<-c(paste0("obsOut",1:nOuts),"nOuts")

tableOuts<-mOutliers
tableOrders<-orderCPCA8

outs<-match((tableOuts)[1,1:(tableOuts)[1,nOuts+1]],tableOrders)

if(printing==TRUE)x11(width=16,height=12)
if(printing==TRUE)par(mfrow=c(compPerfSq(1,length(coreG)),compPerfSq(1,length(coreG))))
if(printing==TRUE)par(mar=c(1,1,1,1))
#coreG genes con el orden CPCA54 y outliers para el orden CPCA54
for(i in 1:length(coreG)){
if(printing==TRUE)plot(escalaPhi8,mNorm8[i,tableOrders],type="b",xlab="",ylab="",xaxt="n",yaxt="n",main=coreG[i])

if(printing==TRUE)points(escalaPhi8[outs],mNorm8[i,tableOrders[outs]],col=c(1:(tableOuts)[1,nOuts+1]),pch=0,cex=1.4)
}
if(printing==TRUE)plot(escalaPhi8,eigen18[tableOrders],type="b",xaxt="n",yaxt="n",
main=paste("1st Eigenegene",". %Var= ",round(varPer8[1],3),sep=""))
if(printing==TRUE)points(escalaPhi8[outs],(eigen18[tableOrders])[outs],col=c(1:(tableOuts)[1,nOuts+1]),pch=0,cex=1.4)
if(printing==TRUE)plot(escalaPhi8,eigen28[tableOrders],type="b",xaxt="n",yaxt="n",
main=paste("2nd Eigenegene",". %Var= ",round(varPer8[2],3),sep=""))
if(printing==TRUE)points(escalaPhi8[outs],(eigen28[tableOrders])[outs],col=c(1:(tableOuts)[1,nOuts+1]),pch=0,cex=1.4)
if(printing==TRUE)plot(escalaPhi8,eigen38[tableOrders],type="b",xaxt="n",yaxt="n",
main=paste("3rd Eigenegene",". %Var= ",round(varPer8[3],3),sep=""))
if(printing==TRUE)points(escalaPhi8[outs],(eigen38[tableOrders])[outs],col=c(1:(tableOuts)[1,nOuts+1]),pch=0,cex=1.4)

return(list(orderCPCA8,escalaPhi8,mOutliers,outs,varPer8,eigen18,eigen28,eigen38))
}

obtainCPCA13<-function(mNorm8,tissue,coreG,nOuts,printing){
#CPCA a partir de los 8 Original, non-equispaced
cp8<-prcomp(centrado(mNorm8), scale. = TRUE, center = FALSE)
varPer8<-c(cp8[[1]][1]^2/sum(cp8[[1]]^2),
cp8[[1]][2]^2/sum(cp8[[1]]^2),
cp8[[1]][3]^2/sum(cp8[[1]]^2))#varianza explicada
eigen18<-cp8$rotation[,1]
eigen28<-cp8$rotation[,2]
eigen38<-cp8$rotation[,3]
xi8<-eigen18/sqrt((eigen18^2+ eigen28^2))
yi8<-eigen28/sqrt((eigen18^2+ eigen28^2))
phi8<-(atan2(yi8,xi8))%%(2*pi)
orderCPCA8<-order(phi8)
escalaPhi8<-phi8[orderCPCA8]
rojo8<-FALSE
maxi<-round(max(abs(eigen18),abs(eigen28)),1)
d8<-c()
for(i in 1:length(eigen18)){
d8[i]<-sqrt(eigen18[i]^2+eigen28[i]^2)
}
obs8<-order(d8)[1:nOuts]

if(printing==TRUE)x11(width=15,height=15)
if(printing==TRUE)par(mfrow=c(1,1))
if(printing==TRUE)plot(eigen18,eigen28,xlim=c(-maxi,maxi),ylim=c(-maxi,maxi),
main=paste(tissue,'. Eigen1 (',round(varPer8[1]*100,2),'%) vs Eigen2 (',round(varPer8[2]*100,2),'%). CoreG',sep=''))

s8<-0
for(i in 1:length(obs8)){
if(d8[obs8[i]]<=0.1){
if(printing==TRUE)points(eigen18[obs8[i]],eigen28[obs8[i]],pch=0,col=i)
s8<-s8+1
}
}
if(s8==0){
for(i in 1:length(obs8)){
if(d8[obs8[i]]<=0.15){
if(printing==TRUE)points(eigen18[obs8[i]],eigen28[obs8[i]],pch=0,col=i)
s8<-s8+1
rojo8<-TRUE
}
}
}
ss8<-s8
fCos1<-funcionCosinor(eigen18[orderCPCA8],escalaPhi8,length(escalaPhi8))[[1]]
fCos2<-funcionCosinor(eigen28[orderCPCA8],escalaPhi8,length(escalaPhi8))[[1]]
fNP1<-function1Local(eigen18[orderCPCA8])[[1]]
fNP2<-function1Local(eigen28[orderCPCA8])[[1]]
if(printing==TRUE)lines(fCos1,fCos2,col=3)
if(printing==TRUE)lines(c(fCos1[1],fCos1[length(fCos1)]),c(fCos2[1],fCos2[length(fCos2)]),col=3)
if(printing==TRUE)lines(fNP1,fNP2,col=4)
if(printing==TRUE)lines(c(fNP1[1],fNP1[length(fNP1)]),c(fNP2[1],fNP2[length(fNP2)]),col=2)

#lines(cos(secTimes(length(eigen28)))*0.1,sin(secTimes(length(eigen28)))*0.1,lty=2,col=2)
#lines(cos(secTimes(length(eigen28)))*0.15,sin(secTimes(length(eigen28)))*0.15,lty=2)
#if(rojo8)lines(cos(secTimes(length(eigen28)))*0.15,sin(secTimes(length(eigen28)))*0.15,lty=2,col=2)
#lines(cos(secTimes(length(eigen28)))*0.2,sin(secTimes(length(eigen28)))*0.2,lty=2)
#if(maxi>=0.3)lines(cos(secTimes(length(eigen28)))*0.3,sin(secTimes(length(eigen28)))*0.3,lty=2)
#if(maxi>=0.4)lines(cos(secTimes(length(eigen28)))*0.4,sin(secTimes(length(eigen28)))*0.4,lty=2)


mOutliers<-matrix(c(obs8,ss8),1,nOuts+1)
rownames(mOutliers)<-c("orderCPCACoreG")
colnames(mOutliers)<-c(paste0("obsOut",1:nOuts),"nOuts")

tableOuts<-mOutliers
tableOrders<-orderCPCA8

outs<-match((tableOuts)[1,1:(tableOuts)[1,nOuts+1]],tableOrders)

if(printing==TRUE)x11(width=16,height=12)
if(printing==TRUE)par(mfrow=c(compPerfSq(1,length(coreG)),compPerfSq(1,length(coreG))))
if(printing==TRUE)par(mar=c(1,1,1,1))
#coreG genes con el orden CPCA54 y outliers para el orden CPCA54
for(i in 1:length(coreG)){
if(printing==TRUE)plot(escalaPhi8,mNorm8[i,tableOrders],type="b",xlab="",ylab="",xaxt="n",yaxt="n",main=coreG[i])

if(printing==TRUE)points(escalaPhi8[outs],mNorm8[i,tableOrders[outs]],col=c(1:(tableOuts)[1,nOuts+1]),pch=0,cex=1.4)
}
if(printing==TRUE)plot(escalaPhi8,eigen18[tableOrders],type="b",xaxt="n",yaxt="n",
main=paste("1st Eigenegene",". %Var= ",round(varPer8[1],3),sep=""))
if(printing==TRUE)points(escalaPhi8[outs],(eigen18[tableOrders])[outs],col=c(1:(tableOuts)[1,nOuts+1]),pch=0,cex=1.4)
if(printing==TRUE)plot(escalaPhi8,eigen28[tableOrders],type="b",xaxt="n",yaxt="n",
main=paste("2nd Eigenegene",". %Var= ",round(varPer8[2],3),sep=""))
if(printing==TRUE)points(escalaPhi8[outs],(eigen28[tableOrders])[outs],col=c(1:(tableOuts)[1,nOuts+1]),pch=0,cex=1.4)
if(printing==TRUE)plot(escalaPhi8,eigen38[tableOrders],type="b",xaxt="n",yaxt="n",
main=paste("3rd Eigenegene",". %Var= ",round(varPer8[3],3),sep=""))
if(printing==TRUE)points(escalaPhi8[outs],(eigen38[tableOrders])[outs],col=c(1:(tableOuts)[1,nOuts+1]),pch=0,cex=1.4)

return(list(orderCPCA8,escalaPhi8,mOutliers,outs,varPer8,eigen18,eigen28,eigen38))
}





#FUNCION PARA PINTAR OUTLIERS Y GENES CON ORDEN, PERO SOLO CON JUSTE np

#requiere db, fullGenes, dbCoreG, normalice, coreG, obtainCPCA12
#sampIdTissue=giveSampleRefG
#nameT=nameTissue


#matrixAllSamples<-epidermis
#nameT<-"Epidermis"
#coreG<-coreG2

giveMatIniNP_v3_cores<-function(matrixAllSamples,nameT,coreG){

mFull0Tissue<-matrixAllSamples
indNoNames<-which(is.na(rownames(mFull0Tissue)),arr.ind=TRUE)
if(length(indNoNames)>0)mFull0Tissue<-mFull0Tissue[-indNoNames,]

#quirtar de mFull0Tissue, qaqueeloos en los que mas del 30 ciento sean ceros
drop<-c()
for(i in 1:nrow(mFull0Tissue)){
  if(length(which(mFull0Tissue[i,]==0,arr.ind=TRUE))>0.30*length(mFull0Tissue[i,])
		| sum(is.na(mFull0Tissue[i,]))>0.30*length(mFull0Tissue[i,])){
    drop<-c(drop,i)
  }
}
print(paste("Deleted ",length(drop)," genes with more than 30% of zeros"))
#length(drop)
if(length(drop)>0){
	mFull1Tissue<-mFull0Tissue[-drop,]
	rownames(mFull1Tissue)<-rownames(mFull0Tissue)[-drop]
}else{
	mFull1Tissue<-mFull0Tissue
	rownames(mFull1Tissue)<-rownames(mFull0Tissue)
}


namesRep<-names(which(sort(table(rownames(mFull1Tissue)))>1))
#entre los genes con el mismo nombre seleccionar aquel que tenga mayor MAD
if(length(namesRep)>0){
filasOut<-c()
for(i in 1:length(namesRep)){
	filas<-which(namesRep[i]==rownames(mFull1Tissue),arr.ind=TRUE)
	mads<-apply(mFull1Tissue[filas,],1,mad,constant=1)
    	filasOut<-c(filasOut,filas[-which.max(mads)])
}


mFull1Tissue<-mFull1Tissue[-filasOut,]
}
mFullTissueNorm1<-t(apply(mFull1Tissue,1,normalice))
rownames(mFullTissueNorm1)<-rownames(mFull1Tissue)

#coreG genes
dropTissueOut<-c()
mTissueCoreG<-mFull1Tissue[match(coreG,rownames(mFull1Tissue)),]
mTissueCoreGNorm<-mFullTissueNorm1[match(coreG,rownames(mFullTissueNorm1)),]
#decidir si se elimina algun gene de los core porque no verifican el requisito de las vareianz
if(obtainCPCA12(mTissueCoreGNorm,nameT,coreG,8,FALSE)[[5]][1]*0.4<obtainCPCA12(mTissueCoreGNorm,nameT,coreG,8,FALSE)[[5]][1]){
	print("ok")


	#decidir si eliminamos outliers. Usar modelo parametrico
	initialTissue<-obtainCPCA13(mTissueCoreGNorm,nameT,coreG,8,FALSE)

	resTissue<-matrix(0,length(coreG)+3,ncol(mTissueCoreGNorm))
	resStTissue<-matrix(0,length(coreG)+3,ncol(mTissueCoreGNorm))
	resParTissue<-matrix(0,length(coreG)+3,ncol(mTissueCoreGNorm))
	resParStTissue<-matrix(0,length(coreG)+3,ncol(mTissueCoreGNorm))
	fitCosCore<-matrix(0,length(coreG)+3,ncol(mTissueCoreGNorm))
	fitParCore<-matrix(0,length(coreG)+3,ncol(mTissueCoreGNorm))
	datCore<-matrix(0,length(coreG)+3,ncol(mTissueCoreGNorm))
	FMMParCoreG<-matrix(0,length(coreG)+3,5)
	CosParCoreG<-matrix(0,length(coreG)+3,3)
	phisCoreG<-c()
	peaksCoreG<-c()
	for(i in 1:(length(coreG)+3)){
		#data
		if(i<=length(coreG))datCore[i,]<-mTissueCoreGNorm[i,initialTissue[[1]]]
		if(i==(length(coreG)+1))datCore[i,]<-initialTissue[[6]][initialTissue[[1]]]
		if(i==(length(coreG)+2))datCore[i,]<-initialTissue[[7]][initialTissue[[1]]]
		if(i==(length(coreG)+3))datCore[i,]<-initialTissue[[8]][initialTissue[[1]]]
		#fits
		if(i<=length(coreG))funCos<-funcionCosinor(mTissueCoreGNorm[i,initialTissue[[1]]],initialTissue[[2]],length(initialTissue[[1]]))
		if(i<=length(coreG))funPar<-fitFMM_Par(mTissueCoreGNorm[i,initialTissue[[1]]],initialTissue[[2]])
		if(i<=length(coreG))resTissue[i,]<-datCore[i,]-funCos[[1]]
		if(i<=length(coreG))resParTissue[i,]<-datCore[i,]-funPar[[1]]
		if(i==(length(coreG)+1))funCos<-funcionCosinor(datCore[i,],initialTissue[[2]],length(initialTissue[[1]]))
		if(i==(length(coreG)+1))funPar<-fitFMM_Par(datCore[i,],initialTissue[[2]])
		if(i==(length(coreG)+1))resTissue[i,]<-datCore[i,]-funCos[[1]]
		if(i==(length(coreG)+1))resParTissue[i,]<-datCore[i,]-funPar[[1]]
		if(i==(length(coreG)+2))funCos<-funcionCosinor(datCore[i,],initialTissue[[2]],length(initialTissue[[1]]))
		if(i==(length(coreG)+2))funPar<-fitFMM_Par(datCore[i,],initialTissue[[2]])
		if(i==(length(coreG)+2))resTissue[i,]<-datCore[i,]-funCos[[1]]
		if(i==(length(coreG)+2))resParTissue[i,]<-datCore[i,]-funPar[[1]]
		if(i==(length(coreG)+3))funCos<-funcionCosinor(datCore[i,],initialTissue[[2]],length(initialTissue[[1]]))
		if(i==(length(coreG)+3))funPar<-fitFMM_Par(datCore[i,],initialTissue[[2]])
		if(i==(length(coreG)+3))resTissue[i,]<-datCore[i,]-funCos[[1]]
		if(i==(length(coreG)+3))resParTissue[i,]<-datCore[i,]-funPar[[1]]
		fitCosCore[i,]<-funCos[[1]]
		fitParCore[i,]<-funPar[[1]]
	
		#res
		resStTissue[i,]<-(resTissue[i,]-mean(resTissue[i,]))/sd(resTissue[i,])
		resParStTissue[i,]<-(resParTissue[i,]-mean(resParTissue[i,]))/sd(resParTissue[i,])

		#params
		phisCoreG[i]<-(-funCos[[5]])%%(2*pi)
		FMMParCoreG[i,]<-c(funPar[[2]],funPar[[3]],funPar[[4]],funPar[[5]],funPar[[6]])
		CosParCoreG[i,]<-c(funCos[[2]],funCos[[3]],funCos[[5]])
		peaksCoreG[i]<-compUU(FMMParCoreG[i,1],FMMParCoreG[i,2],FMMParCoreG[i,3])
		

	}

	
	
	#PLOT DE LOS CORE CLOCK GENES CON PUNTOS OUT Y AJUSTE PAR
	#x11(width=16,height=12)
	png(paste("outPar_",nameT,".png",sep=""),width=960,height=480)
	par(mfrow=c(compPerfSq(1,length(coreG)+3),compPerfSq(1,length(coreG)+3)))
	par(mar=c(1,1,1,1))
	for(i in 1:(length(coreG)+3)){
		if(i<=(length(coreG)))plot(initialTissue[[2]],datCore[i,],type="b",xaxt="n",yaxt="n",
			main=paste(rownames(mTissueCoreGNorm)[i],sep=""),ylim=c(-1.15,1.15))
		if(i>=(length(coreG)+1))plot(initialTissue[[2]],datCore[i,],type="b",xaxt="n",yaxt="n",
			main=paste(nameT," Eigen ",i-12,". Var=",round(initialTissue[[5]][i-12],3),"%",sep=""))
		points(initialTissue[[2]][initialTissue[[4]]],datCore[i,][initialTissue[[4]]],col=1:8,pch=rep(0,8))
		lines(seq(initialTissue[[2]][1],initialTissue[[2]][length(initialTissue[[2]])],length.out=100),
				outF2(datCore[i,],FMMParCoreG[i,c(4,5,1,2,3)],
				seq(initialTissue[[2]][1],initialTissue[[2]][length(initialTissue[[2]])],length.out=100)),col=4)
	}
	dev.off()
	#PLOT DE LOS RESIDUOS CORE CLOCK GENES CON PUNTOS OUT Y AJUSTE PAR
	#x11(width=16,height=12)
	png(paste("resOutPar_",nameT,".png",sep=""),width=960,height=480)
	par(mfrow=c(compPerfSq(1,length(coreG)+3),compPerfSq(1,length(coreG)+3)))
	par(mar=c(1,1,1,1))
	for(i in 1:(length(coreG)+3)){
		if(i<=(length(coreG)))plot(initialTissue[[2]],resParStTissue[i,],type="b",xaxt="n",yaxt="n",ylim=c(-10,10),
			main=paste(rownames(mTissueCoreGNorm)[i],sep=""))
		if(i>=(length(coreG)+1))plot(initialTissue[[2]],resParStTissue[i,],type="b",xaxt="n",yaxt="n",ylim=c(-10,10),
			main=paste(nameT," Eigen ",i-12,". Var=",round(initialTissue[[5]][i-12],3),"%",sep=""))
		points(initialTissue[[2]][initialTissue[[4]]],resParStTissue[i,][initialTissue[[4]]],col=1:8,pch=rep(0,8))
		abline(h=c(-3,3),col=2,lty=2)
		abline(h=c(-4,4),col=2,lty=2)
		abline(h=c(-5,5),col=2,lty=2)
	}
	dev.off()
	#identificar todos los outliers: fuera de -3,3 y marcadaos on eigen
	outsMult<-c()
	outsMultG<-c()
	outsUni<-c()
	outsUniG<-c()
	mOutsUni<-c()
	mOutsMult<-c()
	nOutsTot<-0
	for(i in 1:length(coreG)){
		if(sum(abs(resParStTissue[i,])>3)>0 &  nOutsTot<=ceiling(0.05*ncol(resParStTissue))){
			cuantos<-which(abs(resParStTissue[i,])>3,arr.ind=TRUE)
			for(j in 1:length(cuantos)){
				if(!is.na(match(initialTissue[[1]][cuantos[j]],initialTissue[[3]])) &
					match(initialTissue[[1]][cuantos[j]],initialTissue[[3]])<=initialTissue[[3]][1,length(initialTissue[[3]])]){
					outsMult<-c(outsMult,initialTissue[[1]][cuantos[j]])
					outsMultG<-c(outsMultG,rep(coreG[i],length(initialTissue[[1]][cuantos[j]])))
					nOutsTot<-nOutsTot+1
				}
			}
		}
		
		if(sum(abs(resParStTissue[i,])>4)>0 &  nOutsTot<=ceiling(0.05*ncol(resParStTissue))){
			#print(i)
			outsUni<-c(outsUni,initialTissue[[1]][which(abs(resParStTissue[i,])>4,arr.ind=TRUE)])
			outsUniG<-c(outsUniG,rep(coreG[i],length(initialTissue[[1]][which(abs(resParStTissue[i,])>4,arr.ind=TRUE)])))
			nOutsTot<-nOutsTot+1
		}
	}
	if(length(outsMult)>0)mOutsMult<-matrix(outsMult,1,length(outsMult))
	if(length(outsMult)>0)colnames(mOutsMult)<-outsMultG
	if(length(outsUni)>0)mOutsUni<-matrix(outsUni,1,length(outsUni))
	if(length(outsUni)>0)colnames(mOutsUni)<-outsUniG

	print(paste(length(outsUni)," univariate outliers were deleted",sep=""))
	print(paste(length(outsMult)," multivariate outliers were deleted",sep=""))
	dropTissueOut<-union(mOutsUni,mOutsMult)#initialTissue[[3]][1,c(5)]
	if(length(dropTissueOut)>0)mTissueCoreG<-mTissueCoreG[,-dropTissueOut]
	if(length(dropTissueOut)>0)mTissueCoreGNorm<-t(apply(mTissueCoreG,1,normalice))
	if(length(dropTissueOut)>0){
		initial2Tissue<-obtainCPCA13(mTissueCoreGNorm,nameT,coreG,8,TRUE)
	}else{
		initial2Tissue<-initialTissue
	}
	#con este orden y para todos genes ver con NP cuales son ritmicos, criterio R2
	initialOrdTissue<-initial2Tissue[[1]]
	initialEscTissue<-initial2Tissue[[2]]
	if(length(dropTissueOut)>0){
		mFullTissue<-(mFull1Tissue[,-dropTissueOut])[,initialOrdTissue]
	}else{
		mFullTissue<-mFull1Tissue[,initialOrdTissue]
	}
	mFullTissueNorm<-t(apply(mFullTissue,1,normalice))

	rownames(mFullTissueNorm)<-rownames(mFullTissue)


	#plot de los doce cores con ajustes FMM, COS y NP
	#x11()
	png(paste("coreOr12NP",nameT,".png",sep=""),height=480,width=960)
	par(mfrow=c(compPerfSq(1,length(coreG)),compPerfSq(1,length(coreG))))
	par(mar=c(1,1,1,1))
	allParAfter<-list();allCosAfter<-list();allNPAfter<-list()
	fitParAfter<-c();fitCosAfter<-c();fitNPAfter<-c()
	r2ParAfter<-c();r2CosAfter<-c();r2NPAfter<-c();phisCosAfter<-c();phisFMMAfter<-c()
	for(i in 1:length(coreG)){
		datAfter<-mFullTissueNorm[match(coreG[i],rownames(mFullTissueNorm)),]
		#par
		#original#allParAfter[[i]]<-fitFMM_Par2(datAfter,initialEscTissue)
		allParAfter[[i]]<-fitFMM_Par(datAfter,initialEscTissue)
		fitParAfter<-c(fitParAfter,allParAfter[[i]][[1]])
		phisFMMAfter=c(phisFMMAfter,compUU(allParAfter[[i]][[2]],allParAfter[[i]][[3]],allParAfter[[i]][[4]]))
		r2ParAfter[i]<-1-mean((datAfter-allParAfter[[i]][[1]])^2)/mean((datAfter-mean(datAfter))^2)
		#cos
		allCosAfter[[i]]<-funcionCosinor(datAfter,initialEscTissue,length(initialEscTissue))
		fitCosAfter<-c(fitCosAfter,allCosAfter[[i]][[1]])
		phisCosAfter<-c(phisCosAfter,(-allCosAfter[[i]][[5]])%%(2*pi))
		r2CosAfter[i]<-1-mean((datAfter-allCosAfter[[i]][[1]])^2)/mean((datAfter-mean(datAfter))^2)
		#NP
		allNPAfter[[i]]<-function1Local(datAfter)
		fitNPAfter<-c(fitNPAfter,allNPAfter[[i]][[1]])
		r2NPAfter[i]<-1-mean((datAfter-allNPAfter[[i]][[1]])^2)/mean((datAfter-mean(datAfter))^2)
		plot(initialEscTissue,datAfter,type="b",xaxt="n",yaxt="n",main=paste(coreG[i]," R2NP_Aft=",round(r2NPAfter[i],3),sep=""))
		equis<-seq(initialEscTissue[1],initialEscTissue[length(initialEscTissue)],length.out=100)
		lines(equis,outC(c(allCosAfter[[i]][[2]],allCosAfter[[i]][[3]],allCosAfter[[i]][[5]]),equis),col=3)
		lines(initialEscTissue,allNPAfter[[i]][[1]],col=2)
		params<-c(allParAfter[[i]][[5]],allParAfter[[i]][[6]],allParAfter[[i]][[2]],allParAfter[[i]][[3]],allParAfter[[i]][[4]])
		lines(equis,outF2(allParAfter[[i]][[1]],params,equis),col=4)
	}
	dev.off()
	#tablaR2 PARA LOS TREWS MODEDLOS
	m3r2<-cbind(r2ParAfter,r2CosAfter,r2NPAfter)
	rownames(m3r2)<-coreG
	colnames(m3r2)<-paste0(c("r2Par_","r2Cos_","r2NP_"),nameT)
	print(m3r2)
	
	#PLOT DE LOS PHIS DE FMM
	require('circular')
	#x11()
	png(paste("peaksFMMCoresAfter_",nameT,".png",sep=""),width=480,height=480)
	par(mfrow=c(1,1))
	par(mar=c(1,1,1,1))
	plot(as.circular(phisFMMAfter),main=paste(nameT,". FMM Peaks in Clock Core Genes",sep=""),col=3)
	axis.circular(at=as.circular(phisFMMAfter),labels=coreG)
	dev.off()

	#PLOT DE LOS PHIS DE COSINOR
	require('circular')
	#x11()
	png(paste("peaksCosCoresAfter_",nameT,".png",sep=""),width=480,height=480)
	par(mfrow=c(1,1))
	par(mar=c(1,1,1,1))
	plot(as.circular(phisCosAfter),main=paste(nameT,". Cosinor Phi's in Clock Core Genes",sep=""),col=3)
	axis.circular(at=as.circular(phisCosAfter),labels=coreG)
	dev.off()
}else{
	print("programar funcion")
}



return(list(initialOrdTissue,initialEscTissue,mFullTissueNorm,
		mFullTissue,initial2Tissue,mTissueCoreGNorm,mTissueCoreG,
		dropTissueOut,drop,fitCosCore,fitParCore,CosParCoreG,FMMParCoreG,
		peaksCoreG,phisCoreG,resStTissue,resParStTissue,initialTissue,
		allParAfter,allCosAfter,allNPAfter,m3r2,mFullTissueNorm1,phisCosAfter,
		phisFMMAfter,outsUni,outsMult,mOutsUni,mOutsMult))

}


fitMob<-function(M,A,alpha,beta,omega,t){
  return(M+A*cos(beta+2*atan2(omega*sin((t-alpha)/2),cos((t-alpha)/2))))
}


#inputs: giveMatIniRefG2[[kk]], tejidoN
#kk<-1
#tejidoN<-tissues40[kk]
#listaStep1<-giveMatIniRefG

#tejidoN<-tejChar2[i]
#coreG<-coreG
#listaStep1<-giveMatIniRefG

basicPreOder_cores<-function(listaStep1,tejidoN,coreG){
peaksEpidermis2<-listaStep1[[25]]


peakRefEpidermis2<-peaksEpidermis2[match("ARNTL",coreG)]
peakDBPEpidermis2<-peaksEpidermis2[match("DBP",coreG)]
escTEpidermis<-order((listaStep1[[2]]-peakRefEpidermis2+pi)%%(2*pi))
esc2<-((listaStep1[[2]]-peakRefEpidermis2+pi)%%(2*pi))[escTEpidermis]
#x11()
png(filename=paste0("12CorePre_",tejidoN,".png"))
par(mfrow=c(compPerfSq(1,length(coreG)),compPerfSq(1,length(coreG))))
par(mar=c(2,2,2,2))
peakEpidermis2<-c()
a=TRUE
parCore<-matrix(0,length(coreG),5)
for(i in 1:length(coreG)){
	dat_i<-listaStep1[[3]][match(coreG[i],rownames(listaStep1[[3]])),]
	plot(0,0,type="n",main=paste0(tejidoN," ",coreG[i]," R2=",round(listaStep1[[22]][i,1],2)),xaxt="n",yaxt="n",ylab="Exp",xlab="Time",ylim=c(-1.2,1.2),xlim=c(0,2*pi))
		newAlphaEpidermis<-(listaStep1[[19]][[i]][[2]]-peakRefEpidermis2+pi)%%(2*pi)
		if((peakDBPEpidermis2-peakRefEpidermis2+pi)%%(2*pi)<pi){
		#if(!a){
			lines(esc2,dat_i[escTEpidermis],col=rainbowColor(coreG)[i])
			points(esc2,dat_i[escTEpidermis],col=rainbowColor(coreG)[i])
			lines(seq(0,2*pi,length.out=length(dat_i)),fitMob(listaStep1[[19]][[i]][[5]],listaStep1[[19]][[i]][[6]],newAlphaEpidermis,
				listaStep1[[19]][[i]][[3]],listaStep1[[19]][[i]][[4]],seq(0,2*pi,length.out=length(dat_i))),col=1)
			peakEpidermis2[i]<-compUU(newAlphaEpidermis,listaStep1[[19]][[i]][[3]],listaStep1[[19]][[i]][[4]])
			parCore[i,]<-c(listaStep1[[19]][[i]][[5]],listaStep1[[19]][[i]][[6]],
				newAlphaEpidermis,listaStep1[[19]][[i]][[3]],listaStep1[[19]][[i]][[4]])
		}else{
			lines(2*pi-rev(esc2),rev(dat_i[escTEpidermis]),col=rainbowColor(coreG)[i])
			points(2*pi-rev(esc2),rev(dat_i[escTEpidermis]),col=rainbowColor(coreG)[i])
			dd<-outF2(rev(dat_i),c(listaStep1[[19]][[i]][[5]],listaStep1[[19]][[i]][[6]],
				2*pi-newAlphaEpidermis,2*pi-listaStep1[[19]][[i]][[3]],listaStep1[[19]][[i]][[4]]),seq(0,2*pi,length.out=100))
			lines(seq(0,2*pi,length.out=100),dd,col=1)
			peakEpidermis2[i]<-compUU(2*pi-newAlphaEpidermis,2*pi-listaStep1[[19]][[i]][[3]],listaStep1[[19]][[i]][[4]])
			parCore[i,]<-c(listaStep1[[19]][[i]][[5]],listaStep1[[19]][[i]][[6]],
				2*pi-newAlphaEpidermis,2*pi-listaStep1[[19]][[i]][[3]],listaStep1[[19]][[i]][[4]])
		}
}
matNew<-matrix(0,nrow(listaStep1[[3]]),ncol(listaStep1[[3]]))
if((peakDBPEpidermis2-peakRefEpidermis2+pi)%%(2*pi)<pi){
	oNew<-escTEpidermis
	escNew<-esc2
	matNew<-listaStep1[[3]][,escTEpidermis]
	reverse<-FALSE
}else{
	oNew<-rev(escTEpidermis)
	escNew<-2*pi-rev(esc2)
	for(i in 1:nrow(listaStep1[[3]])){
		matNew[i,]<-rev(listaStep1[[3]][i,escTEpidermis])
	}
	reverse<-TRUE
}
rownames(matNew)<-rownames(listaStep1[[3]])
colnames(matNew)<-colnames(listaStep1[[3]])
dev.off()
coreDay<-c();namesDay<-c();r2Day<-c()
coreNight<-c();namesNight<-c();r2Night<-c()
for(i in 1:length(peakEpidermis2)){
	if(coreG[i]!="ARNTL" & coreG[i]!="DBP"){
		if(0<=peakEpidermis2[i] & peakEpidermis2[i]<pi ){
			coreDay<-c(coreDay,peakEpidermis2[i])
			r2Day<-c(r2Day,listaStep1[[22]][i,1])
			namesDay<-c(namesDay,coreG[i])
		}else{
			coreNight<-c(coreNight,peakEpidermis2[i])
			r2Night<-c(r2Night,listaStep1[[22]][i,1])
			namesNight<-c(namesNight,coreG[i])
		}
	}
}
return(list(oNew,escNew,matNew,peakEpidermis2,listaStep1[[22]][,1],parCore,coreDay,r2Day,namesDay,coreNight,r2Night,namesNight,reverse))
}




basicPreOderMod_cores<-function(listaStep1,tejidoN,coreG){
peaksEpidermis2<-listaStep1[[25]]


peakRefEpidermis2<-peaksEpidermis2[match("ARNTL",coreG)]
peakDBPEpidermis2<-peaksEpidermis2[match("DBP",coreG)]
escTEpidermis<-order((listaStep1[[2]]-peakRefEpidermis2+pi)%%(2*pi))
esc2<-((listaStep1[[2]]-peakRefEpidermis2+pi)%%(2*pi))[escTEpidermis]
#x11()
png(filename=paste0("12CorePre_",tejidoN,".png"))
par(mfrow=c(compPerfSq(1,length(coreG)),compPerfSq(1,length(coreG))))
par(mar=c(2,2,2,2))
peakEpidermis2<-c()
a=TRUE
parCore<-matrix(0,length(coreG),5)
for(i in 1:length(coreG)){
	dat_i<-listaStep1[[3]][match(coreG[i],rownames(listaStep1[[3]])),]
	plot(0,0,type="n",main=paste0(tejidoN," ",coreG[i]," R2=",round(listaStep1[[22]][i,1],2)),xaxt="n",yaxt="n",ylab="Exp",xlab="Time",ylim=c(-1,1),xlim=c(0,2*pi))
		newAlphaEpidermis<-(listaStep1[[19]][[i]][[2]]-peakRefEpidermis2+pi)%%(2*pi)
		if((peakDBPEpidermis2-peakRefEpidermis2+pi)%%(2*pi)<pi){
		#if(!a){
			lines(esc2,dat_i[escTEpidermis],col=rainbowColor(coreG)[i])
			points(esc2,dat_i[escTEpidermis],col=rainbowColor(coreG)[i])
			lines(seq(0,2*pi,length.out=length(dat_i)),fitMob(listaStep1[[19]][[i]][[5]],listaStep1[[19]][[i]][[6]],newAlphaEpidermis,
				listaStep1[[19]][[i]][[3]],listaStep1[[19]][[i]][[4]],seq(0,2*pi,length.out=length(dat_i))),col=1)
			peakEpidermis2[i]<-compUU(newAlphaEpidermis,listaStep1[[19]][[i]][[3]],listaStep1[[19]][[i]][[4]])
			parCore[i,]<-c(listaStep1[[19]][[i]][[5]],listaStep1[[19]][[i]][[6]],
				newAlphaEpidermis,listaStep1[[19]][[i]][[3]],listaStep1[[19]][[i]][[4]])
		}else{
			lines(2*pi-rev(esc2),rev(dat_i[escTEpidermis]),col=rainbowColor(coreG)[i])
			points(2*pi-rev(esc2),rev(dat_i[escTEpidermis]),col=rainbowColor(coreG)[i])
			dd<-outF2(rev(dat_i),c(listaStep1[[19]][[i]][[5]],listaStep1[[19]][[i]][[6]],
				2*pi-newAlphaEpidermis,2*pi-listaStep1[[19]][[i]][[3]],listaStep1[[19]][[i]][[4]]),2*pi-rev(esc2))
			lines(2*pi-rev(esc2),dd,col=1)
			peakEpidermis2[i]<-compUU(2*pi-newAlphaEpidermis,2*pi-listaStep1[[19]][[i]][[3]],listaStep1[[19]][[i]][[4]])
			parCore[i,]<-c(listaStep1[[19]][[i]][[5]],listaStep1[[19]][[i]][[6]],
				2*pi-newAlphaEpidermis,2*pi-listaStep1[[19]][[i]][[3]],listaStep1[[19]][[i]][[4]])
		}
}
matNew<-matrix(0,nrow(listaStep1[[3]]),ncol(listaStep1[[3]]))
if(peakEpidermis2[10]<pi){#if((peakDBPEpidermis2-peakRefEpidermis2+pi)%%(2*pi)<pi){#cambio 01/11/21
	oNew<-escTEpidermis
	escNew<-esc2
	matNew<-listaStep1[[3]][,escTEpidermis]
}else{
	oNew<-rev(escTEpidermis)
	escNew<-2*pi-rev(esc2)
	for(i in 1:nrow(listaStep1[[3]])){
		matNew[i,]<-rev(listaStep1[[3]][i,escTEpidermis])
	}
}
rownames(matNew)<-rownames(listaStep1[[3]])
colnames(matNew)<-colnames(listaStep1[[3]])
dev.off()
coreDay<-c();namesDay<-c();r2Day<-c()
coreNight<-c();namesNight<-c();r2Night<-c()
for(i in 1:length(peakEpidermis2)){
	if(coreG[i]!="ARNTL" & coreG[i]!="DBP"){
		if(0<=peakEpidermis2[i] & peakEpidermis2[i]<pi ){
			coreDay<-c(coreDay,peakEpidermis2[i])
			r2Day<-c(r2Day,listaStep1[[22]][i,1])
			namesDay<-c(namesDay,coreG[i])
		}else{
			coreNight<-c(coreNight,peakEpidermis2[i])
			r2Night<-c(r2Night,listaStep1[[22]][i,1])
			namesNight<-c(namesNight,coreG[i])
		}
	}
}
return(list(oNew,escNew,matNew,peakEpidermis2,listaStep1[[22]][,1],parCore,coreDay,r2Day,namesDay,coreNight,r2Night,namesNight))
}

#salida1<-list(oNew,escNew,matNew,peakEpidermis2,listaStep1[[22]][,1],parCore,coreDay,r2Day,namesDay,coreNight,r2Night,namesNight)
#segundo orden
#tejidoN<-tissues40[kk]
#objPre<-salida1

#tejidoN="Epidermis"
#objPre=preOrdRefSkin
#coreG=coreG2

#objPre<-preOrdRefG2
basicOderMod_cores<-function(tejidoN,objPre,coreG){
	peaksPre<-objPre[[4]]
	aviso1<-FALSE;aviso2<-FALSE;aviso3<-FALSE;cond31<-FALSE;cond32<-FALSE
	peakD<-0
	p6am<-1-cos(peaksPre[match("DBP",coreG)])
	p6pm<-1-cos(peaksPre[match("DBP",coreG)]-pi)
	for(i in 1:length(coreG)){
		if(i!=match("ARNTL",coreG) & i!=match("DBP",coreG) & peaksPre[i]>0 & peaksPre[i]<pi)peakD<-peakD+1
	}
	mitad<-floor(length(coreG)-2)/2
	if(p6am<=0.1 )aviso1<-TRUE# | peakD<mitad)aviso1<-TRUE
	if(peaksPre[match("DBP",coreG)]>pi/2 & p6pm>0.1)aviso2<-TRUE
	if(p6pm<=0.1 )aviso3<-TRUE
	cond1<-(aviso1 & peakD<=mitad &
	!(peaksPre[match("DBP",coreG)]<peaksPre[match("CRY1",coreG)] & peaksPre[match("CRY1",coreG)]<peaksPre[match("ARNTL",coreG)]) )|
	(aviso1 & peakD<mitad &
	(peaksPre[match("DBP",coreG)]<peaksPre[match("CRY1",coreG)] & peaksPre[match("CRY1",coreG)]<peaksPre[match("ARNTL",coreG)]) )
	cond2<-aviso2 & peakD<mitad &
		!(peaksPre[match("DBP",coreG)]<peaksPre[match("ARNTL",coreG)] & peaksPre[match("ARNTL",coreG)]<peaksPre[match("CRY1",coreG)])
	if(peaksPre[match("TEF",coreG)]>pi & (peaksPre[match("CRY1",coreG)]>peaksPre[match("TEF",coreG)] | peaksPre[match("CRY1",coreG)]<pi))cond31<-TRUE
	if(peaksPre[match("TEF",coreG)]<pi & peaksPre[match("CRY1",coreG)]>peaksPre[match("TEF",coreG)] & peaksPre[match("CRY1",coreG)]<pi)cond32<-TRUE
	cond3<- aviso3 & !cond31 & !cond32 
	cond4<-!aviso1 & !aviso2 &  (1-cos(peaksPre[match("CRY1",coreG)]-pi))<0.1 & peakD<mitad
	cambiaOri<-FALSE
	peaksPreNew<-c()

	if(cond1 | cond2 |cond3 | cond4){
		peaksPreNew<-(2*pi-peaksPre)%%(2*pi)
		print(tejidoN)
		oNewNew<-rev(objPre[[1]])
		escNewNew<-2*pi-rev(objPre[[2]])
		matNewNew<-matrix(0,nrow(objPre[[3]]),ncol(objPre[[3]]))
		pars<-matrix(0,length(coreG),5)
		for( i in 1:length(coreG)){
			pars[i,]<-c(objPre[[6]][i,1],objPre[[6]][i,2],2*pi-objPre[[6]][i,3],2*pi-objPre[[6]][i,4],objPre[[6]][i,5])
		}
		for(i in 1:nrow(objPre[[3]])){
			matNewNew[i,]<-rev(objPre[[3]][i,])
		}
		rownames(matNewNew)<-rownames(objPre[[3]])
		indNewNew<-rev(1:length(oNewNew))
		cambiaOri<-TRUE

	}else{
		peaksPreNew<-peaksPre
		oNewNew<-objPre[[1]]
		escNewNew<-objPre[[2]]
		matNewNew<-objPre[[3]]
		rownames(matNewNew)<-rownames(objPre[[3]])
		indNewNew<-1:length(oNewNew)
		pars<-matrix(0,length(coreG),5)
		for( i in 1:length(coreG)){
			pars[i,]<-c(objPre[[6]][i,1],objPre[[6]][i,2],objPre[[6]][i,3],objPre[[6]][i,4],objPre[[6]][i,5])
		}
	}

	#x11()
if(!cambiaOri)png(filename=paste0("12CorePre_",tejidoN,".png"))
if(cambiaOri)png(filename=paste0("!!12CorePre_",tejidoN,".png"))
par(mfrow=c(compPerfSq(1,length(coreG)),compPerfSq(1,length(coreG))))
par(mar=c(2,2,2,2))
peakEpidermis2<-c()
for(i in 1:length(coreG)){
	dat_i<-matNewNew[match(coreG[i],rownames(matNewNew)),]
	plot(0,0,type="n",main=paste0(tejidoN," ",coreG[i]," R2=",round(objPre[[5]][i],2)),xaxt="n",yaxt="n",ylab="Exp",xlab="Time",ylim=c(-1.2,1.2),xlim=c(0,2*pi))
	lines(escNewNew,dat_i,col=rainbowColor(coreG)[i])
	points(escNewNew,dat_i,col=rainbowColor(coreG)[i])
	ss<-seq(0,2*pi,length.out=100)#length(dat_i))
	lines(ss,fitMob(pars[i,1],pars[i,2],pars[i,3],pars[i,4],pars[i,5],seq(0,2*pi,length.out=100)),col=1)#length(dat_i)#cambio
	axis(1,at=c(ss[1],ss[ceiling(length(ss)/2)]),labels=c("6AM","6PM"))
}
dev.off()

if(!cambiaOri)png(filename=paste0("12PeaksPre_",tejidoN,".png"))
if(cambiaOri)png(filename=paste0("!!12PeaksPre_",tejidoN,".png"))
plantillaCircular_cores(tejidoN,peaksPreNew,coreG,objPre[[5]])
dev.off()


	return(list(oNewNew,escNewNew,matNewNew,peaksPreNew,objPre[[5]],pars,cambiaOri,indNewNew))
}

#i=1
#	load(file=paste("F:/gTEX/database/store/giveMatIniRefG",tejChar[i],".RData",sep=""))
#	load(file=paste("F:/gTEX/database/store/preOrdRefG2",tejChar[i],".RData",sep=""))
#	load(file=paste("F:/gTEX/database/store/basicOrdRefG2",tejChar[i],".RData",sep=""))


#tejidoN<-tejChar[i]
#objPre<-preOrdRefG2
#coreG<-coreG
basicOderMod2_cores<-function(tejidoN,objPre,coreG){
	peaksPre<-objPre[[4]]
	aviso1<-FALSE;aviso2<-FALSE;aviso3<-FALSE;cond31<-FALSE;cond32<-FALSE;cond33<-FALSE
	peakD<-0
	p6am<-1-cos(peaksPre[match("DBP",coreG)])
	p6pm<-1-cos(peaksPre[match("DBP",coreG)]-pi)
	for(i in 1:length(coreG)){
		if(i!=match("ARNTL",coreG) & i!=match("DBP",coreG) & peaksPre[i]>0 & peaksPre[i]<pi)peakD<-peakD+1
	}
	mitad<-floor(length(coreG)-2)/2
	if(p6am<=0.1 )aviso1<-TRUE# | peakD<mitad)aviso1<-TRUE
	if(peaksPre[match("DBP",coreG)]>pi/2 & p6pm>0.1)aviso2<-TRUE
	if(p6pm<=0.1 )aviso3<-TRUE
	cond1<-(aviso1 & peakD<mitad &
	!(peaksPre[match("DBP",coreG)]<peaksPre[match("CRY1",coreG)] & peaksPre[match("CRY1",coreG)]<peaksPre[match("ARNTL",coreG)]) )#|
	#(aviso1 & peakD<mitad &
	#(peaksPre[match("DBP",coreG)]<peaksPre[match("CRY1",coreG)] & peaksPre[match("CRY1",coreG)]<peaksPre[match("ARNTL",coreG)]) )
	cond2<-aviso2 & peakD<mitad &
		!(peaksPre[match("DBP",coreG)]<peaksPre[match("ARNTL",coreG)] & peaksPre[match("ARNTL",coreG)]<peaksPre[match("CRY1",coreG)])
	if(peaksPre[match("TEF",coreG)]>pi | peaksPre[match("CRY1",coreG)]>pi | (peaksPre[match("CRY1",coreG)]<pi & peaksPre[match("CRY1",coreG)]<peaksPre[match("TEF",coreG)]))cond33<-TRUE
	cond3<-(1-cos(peaksPre[match("CRY1",coreG)]-pi))<0.1 & cond33 & !objPre[[13]]
	cond4<-!aviso1 & !aviso2 &  (1-cos(peaksPre[match("CRY1",coreG)]-pi))<0.1 & peakD<mitad
	cambiaOri<-FALSE
	peaksPreNew<-c()

	if(cond1 | cond2 |cond3 | cond4){
		peaksPreNew<-(2*pi-peaksPre)%%(2*pi)
		print(tejidoN)
		oNewNew<-rev(objPre[[1]])
		escNewNew<-2*pi-rev(objPre[[2]])
		matNewNew<-matrix(0,nrow(objPre[[3]]),ncol(objPre[[3]]))
		pars<-matrix(0,length(coreG),5)
		for( i in 1:length(coreG)){
			pars[i,]<-c(objPre[[6]][i,1],objPre[[6]][i,2],2*pi-objPre[[6]][i,3],2*pi-objPre[[6]][i,4],objPre[[6]][i,5])
		}
		for(i in 1:nrow(objPre[[3]])){
			matNewNew[i,]<-rev(objPre[[3]][i,])
		}
		rownames(matNewNew)<-rownames(objPre[[3]])
		indNewNew<-rev(1:length(oNewNew))
		cambiaOri<-TRUE

	}else{
		peaksPreNew<-peaksPre
		oNewNew<-objPre[[1]]
		escNewNew<-objPre[[2]]
		matNewNew<-objPre[[3]]
		rownames(matNewNew)<-rownames(objPre[[3]])
		indNewNew<-1:length(oNewNew)
		pars<-matrix(0,length(coreG),5)
		for( i in 1:length(coreG)){
			pars[i,]<-c(objPre[[6]][i,1],objPre[[6]][i,2],objPre[[6]][i,3],objPre[[6]][i,4],objPre[[6]][i,5])
		}
	}

	#x11()
if(!cambiaOri)png(filename=paste0("12CorePre_",tejidoN,".png"))
if(cambiaOri)png(filename=paste0("!!12CorePre_",tejidoN,".png"))
par(mfrow=c(compPerfSq(1,length(coreG)),compPerfSq(1,length(coreG))))
par(mar=c(2,2,2,2))
peakEpidermis2<-c()
for(i in 1:length(coreG)){
	dat_i<-matNewNew[match(coreG[i],rownames(matNewNew)),]
	plot(0,0,type="n",main=paste0(tejidoN," ",coreG[i]," R2=",round(objPre[[5]][i],2)),xaxt="n",yaxt="n",ylab="Exp",xlab="Time",ylim=c(-1,1),xlim=c(0,2*pi))
	lines(escNewNew,dat_i,col=rainbowColor(coreG)[i])
	points(escNewNew,dat_i,col=rainbowColor(coreG)[i])
	ss<-seq(0,2*pi,length.out=length(dat_i))
	lines(ss,fitMob(pars[i,1],pars[i,2],pars[i,3],pars[i,4],pars[i,5],seq(0,2*pi,length.out=length(dat_i))),col=1)
	axis(1,at=c(ss[1],ss[ceiling(length(ss)/2)]),labels=c("6AM","6PM"))
}
dev.off()

if(!cambiaOri)png(filename=paste0("12PeaksPre_",tejidoN,".png"))
if(cambiaOri)png(filename=paste0("!!12PeaksPre_",tejidoN,".png"))
plantillaCircular_cores(tejidoN,peaksPreNew,coreG,objPre[[5]])
dev.off()


	return(list(oNewNew,escNewNew,matNewNew,peaksPreNew,objPre[[5]],pars,cambiaOri,indNewNew))
}





basicOder_cores<-function(tejidoN,objPre,coreG){
	peaksPre<-objPre[[4]]
	aviso<-FALSE
	peakD<-0
	p6am<-1-cos(peaksPre[match("DBP",coreG)])
	p6pm<-1-cos(peaksPre[match("DBP",coreG)]-pi)
	for(i in 1:length(coreG)){
		if(i!=match("ARNTL",coreG) & i!=match("DBP",coreG) & peaksPre[i]>0 & peaksPre[i]<pi)peakD<-peakD+1
	}
	mitad<-floor(length(coreG)-2)/2
	if(p6am<=0.1 | p6pm<=0.1 | peakD<mitad)aviso<-TRUE
	cambiaOri<-FALSE
	peaksPreNew<-c()
	if(aviso & !(peaksPre[match("DBP",coreG)]<peaksPre[match("CRY1",coreG)] & peaksPre[match("CRY1",coreG)]<peaksPre[match("ARNTL",coreG)])){
		peaksPreNew<-(2*pi-peaksPre)%%(2*pi)
		print(tejidoN)
		oNewNew<-rev(objPre[[1]])
		escNewNew<-2*pi-rev(objPre[[2]])
		matNewNew<-matrix(0,nrow(objPre[[3]]),ncol(objPre[[3]]))
		pars<-matrix(0,length(coreG),5)
		for( i in 1:length(coreG)){
			pars[i,]<-c(objPre[[6]][i,1],objPre[[6]][i,2],2*pi-objPre[[6]][i,3],2*pi-objPre[[6]][i,4],objPre[[6]][i,5])
		}
		for(i in 1:nrow(objPre[[3]])){
			matNewNew[i,]<-rev(objPre[[3]][i,])
		}
		rownames(matNewNew)<-rownames(objPre[[3]])
		indNewNew<-rev(1:length(oNewNew))
		cambiaOri<-TRUE

	}else{
		peaksPreNew<-peaksPre
		oNewNew<-objPre[[1]]
		escNewNew<-objPre[[2]]
		matNewNew<-objPre[[3]]
		rownames(matNewNew)<-rownames(objPre[[3]])
		indNewNew<-1:length(oNewNew)
		pars<-matrix(0,length(coreG),5)
		for( i in 1:length(coreG)){
			pars[i,]<-c(objPre[[6]][i,1],objPre[[6]][i,2],objPre[[6]][i,3],objPre[[6]][i,4],objPre[[6]][i,5])
		}
	}

	#x11()
if(!cambiaOri)png(filename=paste0("12CorePre_",tejidoN,".png"))
if(cambiaOri)png(filename=paste0("!!12CorePre_",tejidoN,".png"))
par(mfrow=c(compPerfSq(1,length(coreG)),compPerfSq(1,length(coreG))))
par(mar=c(2,2,2,2))
peakEpidermis2<-c()
for(i in 1:length(coreG)){
	dat_i<-matNewNew[match(coreG[i],rownames(matNewNew)),]
	plot(0,0,type="n",main=paste0(tejidoN," ",coreG[i]," R2=",round(objPre[[5]][i],2)),xaxt="n",yaxt="n",ylab="Exp",xlab="Time",ylim=c(-1,1),xlim=c(0,2*pi))
	lines(escNewNew,dat_i,col=rainbowColor(coreG)[i])
	points(escNewNew,dat_i,col=rainbowColor(coreG)[i])
	ss<-seq(0,2*pi,length.out=length(dat_i))
	lines(ss,fitMob(pars[i,1],pars[i,2],pars[i,3],pars[i,4],pars[i,5],seq(0,2*pi,length.out=length(dat_i))),col=1)
	axis(1,at=c(ss[1],ss[ceiling(length(ss)/2)]),labels=c("6AM","6PM"))
}
dev.off()

if(!cambiaOri)png(filename=paste0("12PeaksPre_",tejidoN,".png"))
if(cambiaOri)png(filename=paste0("!!12PeaksPre_",tejidoN,".png"))
plantillaCircular_cores(tejidoN,peaksPreNew,coreG,objPre[[5]])
dev.off()


	return(list(oNewNew,escNewNew,matNewNew,peaksPreNew,objPre[[5]],pars,cambiaOri,indNewNew))
}
#salida2<-list(oNewNew,escNewNew,matNewNew,peaksPre,objPre[[5]],objPre[[6]])

rainbowColor<-function(coreG){
	colorea<-c()
	are<-match(coreG,c("PER1","PER2","PER3","CRY1","CRY2","ARNTL","CLOCK","NR1D1","RORA","DBP","TEF","STAT3"))
	for(i in 1:length(coreG)){
		colorea[i]<-rainbow(12)[are[i]]
	}
	return(colorea)
}

cambiaPeaks<-function(peaks){
	peaks2<-c()
	for(i in 1:length(peaks)){
		if(0<peaks[i] & peaks[i]<=pi/2)peaks2[i]<-(3*pi/2-(peaks[i]-0))%%(2*pi)
		if(pi/2<peaks[i] & peaks[i]<=pi)peaks2[i]<-(pi-(peaks[i]-pi/2))%%(2*pi)
		if(pi<peaks[i] & peaks[i]<=3*pi/2)peaks2[i]<-(pi/2-(peaks[i]-pi))%%(2*pi)
		if(3*pi/2<peaks[i] & peaks[i]<=2*pi)peaks2[i]<-(2*pi-(peaks[i]-3*pi/2))%%(2*pi)
	}
	return(peaks2)
}

#titulo="Epidermis Pre"
#peaks=preOrdRefSkin[[4]]
#namesPeaks=coreG2
#r2P=preOrdRefSkin[[5]]
plantillaCircular_cores<-function(titulo,peaks,namesPeaks,r2P){
	peaks2<-c()
	for(i in 1:length(peaks)){
		if(0<peaks[i] & peaks[i]<=pi/2)peaks2[i]<-(3*pi/2-(peaks[i]-0))%%(2*pi)
		if(pi/2<peaks[i] & peaks[i]<=pi)peaks2[i]<-(pi-(peaks[i]-pi/2))%%(2*pi)
		if(pi<peaks[i] & peaks[i]<=3*pi/2)peaks2[i]<-(pi/2-(peaks[i]-pi))%%(2*pi)
		if(3*pi/2<peaks[i] & peaks[i]<=2*pi)peaks2[i]<-(2*pi-(peaks[i]-3*pi/2))%%(2*pi)
	}
	require(plotrix)
	#x11()
	par(mfrow=c(1,1))
	par(mar=c(0.05,0.05,0.05,0.05),mgp=c(-0.5, 0.5, 0))
	plot.circular(as.circular(0),type="n",axes=FALSE,cex.main=0.5)
	title(titulo,cex=0.25)
	points(as.circular(pi/2-0.3),pch=17,cex=1.5)

	for(i in 1:length(namesPeaks)){
		points(as.circular(peaks2[i]),col=rainbowColor(namesPeaks)[i])	
	}
	for(i in 1:length(namesPeaks)){
		axis.circular(at=peaks2[i],labels=namesPeaks[i],cex=0.9,col=rainbowColor(namesPeaks)[i])
	}
	
	axis.circular(at=peaks2,labels=round(r2P,2),tcl.text = -0.095)
	axis.circular(at=c(pi/2,3*pi/2,0,pi),labels=c("6PM","6AM","Midnight","Noon"),tcl.text = -0.055,cex=0.9)#,col=4)
	axis.circular(at=c(pi/2,3*pi/2),labels=c(expression(pi),expression(paste("2",pi,"/0"))),tcl.text = 0.075)
}

plantillaCircular_cores_screen<-function(titulo,peaks,namesPeaks,r2P){
	peaks2<-c()
	for(i in 1:length(peaks)){
		if(0<peaks[i] & peaks[i]<=pi/2)peaks2[i]<-(3*pi/2-(peaks[i]-0))%%(2*pi)
		if(pi/2<peaks[i] & peaks[i]<=pi)peaks2[i]<-(pi-(peaks[i]-pi/2))%%(2*pi)
		if(pi<peaks[i] & peaks[i]<=3*pi/2)peaks2[i]<-(pi/2-(peaks[i]-pi))%%(2*pi)
		if(3*pi/2<peaks[i] & peaks[i]<=2*pi)peaks2[i]<-(2*pi-(peaks[i]-3*pi/2))%%(2*pi)
	}
	require(plotrix)
	#x11()
	#par(mfrow=c(1,1))
	par(mar=c(0.05,0.05,0.05,0.05),mgp=c(-0.5, 0.5, 0))
	plot.circular(as.circular(0),type="n",axes=FALSE,cex.main=0.5)
	title(titulo,cex=0.25)
	points(as.circular(pi/2-0.3),pch=17,cex=1.5)

	for(i in 1:length(namesPeaks)){
		points(as.circular(peaks2[i]),col=rainbowColor(namesPeaks)[i])	
	}
	for(i in 1:length(namesPeaks)){
		axis.circular(at=peaks2[i],labels=namesPeaks[i],cex=0.9,col=rainbowColor(namesPeaks)[i])
	}
	
	axis.circular(at=peaks2,labels=round(r2P,2),tcl.text = -0.095)
	axis.circular(at=c(pi/2,3*pi/2,0,pi),labels=c("6PM","6AM","Midnight","Noon"),tcl.text = -0.055,cex=0.9)#,col=4)
	axis.circular(at=c(pi/2,3*pi/2),labels=c(expression(pi),expression(paste("2",pi,"/0"))),tcl.text = 0.075)
}


cambioOri<-function(peaksOri){
	peaksNew<-c()
	for(i in 1:length(peaksOri)){
		# primer cuadrante -> segundo cuadrante
		if(peaksOri[i]>=0 & peaksOri[i]<pi/2)peaksNew[i]<-pi-peaksOri[i]
		# segundo cuadrante -> primer cuadrante
		if(peaksOri[i]>=pi/2 & peaksOri[i]<pi)peaksNew[i]<-pi-peaksOri[i]
		# tercer cuadrante -> cuarto cuadrante
		if(peaksOri[i]>=pi & peaksOri[i]<3*pi/2)peaksNew[i]<-3*pi/2+(peaksOri[i]-pi)
		# cuarto cuadrante -> tercer cuadrante
		if(peaksOri[i]>=3*pi/2 & peaksOri[i]<2*pi)peaksNew[i]<-pi+(2*pi-peaksOri[i])
	}
	return(peaksNew)
}

compPerfSq<-function(i,f){
	t<-length(i:f)
	j<-0
	sq<-0
	while(sq<t){
		j<-j+1
		sq<-j*j
	}
	return(j)
}


plotTop<-function(plotIni,plotFin,repi,esci,mFini,paramFi,nameTissue,r2){
	#x11()
	side<-compPerfSq(plotIni,plotFin)
	par(mfrow=c(side,side))
	#par(mfrow=c(7,8))
	par(mar=c(1,1,1,1))
	esci<-esci[repi,]
	mFini<-mFini[[repi]]
	paramFi<-paramFi[[repi]]
	r2Ord<-order(r2[-c(1:12)],decreasing=TRUE)+12
	ii<-1
	for(j in plotIni:plotFin){
		if(j<=12){
			plot(esci,mFini[j,],type="b",main=paste(nameTissue," k=",repi," ",rownames(mFini)[j]," R2=",round(r2[j],3),sep=""),
				xaxt="n",yaxt="n")
			tes<-seq(esci[1],esci[length(esci)],length.out=100)
			lines(tes,outF2(mFini[j,],paramFi[j,1:5],tes),col=4)
			#if(rownames(mFini)[j]=="PER1")abline(v=pi,col=3)
		}else{
			plot(esci,mFini[r2Ord[ii],],type="b",main=paste(nameTissue," k=",repi," ",rownames(mFini)[r2Ord[ii]]," R2=",round(r2[r2Ord[ii]],3),sep=""),
				xaxt="n",yaxt="n")
			tes<-seq(esci[1],esci[length(esci)],length.out=100)
			lines(tes,outF2(mFini[r2Ord[ii],],paramFi[r2Ord[ii],1:5],tes),col=4)
			#if(rownames(mFini)[r2Ord[ii]]=="PER1")abline(v=pi,col=3)
			ii<-ii+1
		}
	}
}


function1Local<-function(v){
  v2<-c(v,v)
  i2<-c(1:length(v),1:length(v))
  candL<-c()
  candU<-c()
  for(i in 2:(length(v)+1)){
    #maximos locales
    if(v2[i-1]<=v2[i] & v2[i]>=v2[i+1]){
      candiU<-i2[i]
      candU<-unique(c(candU,candiU))
    }
    #minimos locales
    if(v2[i-1]>=v2[i] & v2[i]<=v2[i+1]){
      candiL<-i2[i]
      candL<-unique(c(candL,candiL))
    }
  }
  
  
  
  mseFin<-9999
  for( i in 1:length(candL) ){
    for (j in 1:length(candU) ){
      #print(c(i,j))
      if(candL[i]<candU[j]){
        indexPavaLU<-candL[i]:candU[j]
        if(length(indexPavaLU)>2){
          pavaLU<-pava(v[(candL[i]+1):(candU[j]-1)])
          pavaLU<-c(v[candL[i]],pavaLU,v[candU[j]])
        }else{
          pavaLU<-c(v[candL[i]],v[candU[j]])
        }
        if(candU[j]!=length(v) & candL[i]!=1){
          indexPavaUL1<-(candU[j]+1):length(v)
          indexPavaUL2<-1:(candL[i]-1)
          pavaUL<-pava(v2[(candU[j]+1):(length(v)+candL[i]-1)],decreasing=TRUE)
          pavaAux<-c(pavaUL[(length(indexPavaUL1)+1):length(pavaUL)],pavaLU,pavaUL[1:length(indexPavaUL1)])
          if((length(v)-length(pavaAux))!=0){
            #print("#######################################################error")
          }
          mseAux<-sum((v-pavaAux)^2)/length(v)
          if(pavaLU[2]>=v[candL[i]] & pavaLU[length(pavaLU)-1]<=v[candU[j]] & 
             pavaUL[1]<v[candU[j]] & pavaUL[length(pavaUL)]>v[candL[i]] & mseAux<mseFin){
            mseFin<-mseAux
            pavaFin<-pavaAux
            Lopt<-candL[i]
            Uopt<-candU[j]
          }
        }
        if(candU[j]==length(v) & candL[i]!=1){
          pavaUL<-pava(v[1:(candL[i]-1)],decreasing=TRUE)
          pavaAux<-c(pavaUL,pavaLU)
          if((length(v)-length(pavaAux))!=0){
            #print("#######################################################error")
          }
          mseAux<-sum((v-pavaAux)^2)/length(v)
          if(pavaLU[2]>v[candL[i]] & pavaLU[length(pavaLU)-1]<v[candU[j]] & 
             pavaUL[1]<v[candU[j]] & pavaUL[length(pavaUL)]>v[candL[i]] & mseAux<mseFin){
            mseFin<-mseAux
            pavaFin<-pavaAux
            Lopt<-candL[i]
            Uopt<-candU[j]
          }
        }
        if(candL[i]==1 & candU[j]!=length(v)){
          pavaUL<-pava(v[(candU[j]+1):(length(v))],decreasing=TRUE)
          pavaAux<-c(pavaLU,pavaUL)
          if((length(v)-length(pavaAux))!=0){
            #print("#######################################################error")
          }
          mseAux<-sum((v-pavaAux)^2)/length(v)
          if(pavaLU[2]>v[candL[i]] & pavaLU[length(pavaLU)-1]<v[candU[j]] & 
             pavaUL[1]<v[candU[j]] & pavaUL[length(pavaUL)]>v[candL[i]] & mseAux<mseFin){
            mseFin<-mseAux
            pavaFin<-pavaAux
            Lopt<-candL[i]
            Uopt<-candU[j]
          }
        }
        if(candL[i]==1 & candU[j]==length(v)){
          pavaUL<-c()
          pavaAux<-pavaLU
          if((length(v)-length(pavaAux))!=0){
            #print("#######################################################error")
          }
          mseAux<-sum((v-pavaAux)^2)/length(v)
          if(pavaLU[2]>v[candL[i]] & pavaLU[length(pavaLU)-1]<v[candU[j]] & mseAux<mseFin){
            mseFin<-mseAux
            pavaFin<-pavaAux
            Lopt<-candL[i]
            Uopt<-candU[j]
          }
        }
      }else{#candU<=candL#revisar
        #agnado un if oara cuando son inguales los candidatos el else es lo que habia para candU<candL
        if(candL[i]==candU[j]){
		indexPavaLU<-c(candL[i]:length(v),1:candU[j])#new
          if(length(indexPavaLU)>2){
            pavaLU<-pava(v2[(candL[i]+1):(length(v)+candU[j]-1)])
            pavaLU<-c(v[candL[i]],pavaLU,v[candU[j]])
          }else{
            pavaLU<-c(v[candL[i]],v[candU[j]])
          }
          pavaAux<-rep(mean(v),length(v))
          mseAux<-sum((v-pavaAux)^2)/length(v)
          if(pavaLU[2]>=v[candL[i]] & pavaLU[length(pavaLU)-1]<=v[candU[j]] &  mseAux<mseFin){#ojo lo he camabiado
            mseFin<-mseAux
            pavaFin<-pavaAux
            Lopt<-candL[i]
            Uopt<-candU[j]
          }else{#agnado es to porque sigue dando error, hay que revisar
			mseFin<-mseAux
            	pavaFin<-pavaAux
            	Lopt<-candL[i]
            	Uopt<-candU[j]
		}
	    
        }else{
          indexPavaLU<-c(candL[i]:length(v),1:candU[j])
          if(length(indexPavaLU)>2){
            pavaLU<-pava(v2[(candL[i]+1):(length(v)+candU[j]-1)])
            pavaLU<-c(v[candL[i]],pavaLU,v[candU[j]])
          }else{
            pavaLU<-c(v[candL[i]],v[candU[j]])
          }
          if(length(pavaLU)==length(v)){
            pavaAux<-c(pavaLU[(length(candL[i]:length(v))+1):length(pavaLU)],pavaLU[1:length(candL[i]:length(v))])
            if((length(v)-length(pavaAux))!=0){
              #print("#######################################################error")
            }
            mseAux<-sum((v-pavaAux)^2)/length(v)
            if( pavaLU[2]>v[candL[i]] & pavaLU[length(pavaLU)-1]<v[candU[j]]  & mseAux<mseFin){
              mseFin<-mseAux
              pavaFin<-pavaAux
              Lopt<-candL[i]
              Uopt<-candU[j]
            }
          }else{
            if(candU[j]!=1 & candL[i]!=length(v)){
              #if(length(pavaLU)!=length(v)){
              pavaUL<-pava(v[(candU[j]+1):(candL[i]-1)],decreasing=TRUE)
              pavaAux<-c(pavaLU[(length(candL[i]:length(v))+1):length(pavaLU)],pavaUL,pavaLU[1:length(candL[i]:length(v))])
              if((length(v)-length(pavaAux))!=0){
               # print("#######################################################error")
              }
              mseAux<-sum((v-pavaAux)^2)/length(v)
              if(pavaLU[2]>=v[candL[i]] & pavaLU[length(pavaLU)-1]<=v[candU[j]] & 
                 pavaUL[1]<v[candU[j]] & pavaUL[length(pavaUL)]>v[candL[i]] & mseAux<mseFin){
                mseFin<-mseAux
                pavaFin<-pavaAux
                Lopt<-candL[i]
                Uopt<-candU[j]
              }
              #}else{
              #pavaUL<-c()
              #\tpavaAux<-pavaLU[c((length(candL[i]:length(v))+1):(length(v)),1:(length(candL[i]:length(v))))]
              #\tmseAux<-sum((v-pavaAux)^2)/length(v)
              #\tif(pavaLU[2]>v[candL[i]] & pavaLU[length(pavaLU)-1]<v[candU[j]]  & mseAux<mseFin){
              #\t\tmseFin<-mseAux
              #\t\tpavaFin<-pavaAux
              #\t\tLopt<-candL[i]
              #\t\tUopt<-candU[j]
              #\t}
              #}
            }
            if(candU[j]==1 & candL[i]!=length(v)){
              pavaUL<-pava(v[2:(candL[i]-1)],decreasing=TRUE)
              pavaAux<-c(pavaLU[length(pavaLU)],pavaUL,pavaLU[-length(pavaLU)])
              if((length(v)-length(pavaAux))!=0){
                #print("#######################################################error")
              }
              mseAux<-sum((v-pavaAux)^2)/length(v)
              if(pavaLU[2]>v[candL[i]] & pavaLU[length(pavaLU)-1]<v[candU[j]] & 
                 pavaUL[1]<v[candU[j]] & pavaUL[length(pavaUL)]>v[candL[i]] & mseAux<mseFin){
                mseFin<-mseAux
                pavaFin<-pavaAux
                Lopt<-candL[i]
                Uopt<-candU[j]
              }
            }
            if(candL[i]==length(v) & candU[j]!=1){
              pavaUL<-pava(v[(candU[j]+1):(length(v)-1)],decreasing=TRUE)
              pavaAux<-c(pavaLU[2:length(pavaLU)],pavaUL,pavaLU[1])
              if((length(v)-length(pavaAux))!=0){
                #print("#######################################################error")
              }
              mseAux<-sum((v-pavaAux)^2)/length(v)
              if(pavaLU[2]>v[candL[i]] & pavaLU[length(pavaLU)-1]<v[candU[j]] & 
                 pavaUL[1]<v[candU[j]] & pavaUL[length(pavaUL)]>v[candL[i]] & mseAux<mseFin){
                mseFin<-mseAux
                pavaFin<-pavaAux
                Lopt<-candL[i]
                Uopt<-candU[j]
              }
            }
            if(candL[i]==length(v) & candU[j]==1){
              pavaUL<-pava(v[(candU[j]+1):(candL[i]-1)],decreasing=TRUE)
              pavaAux<-c(pavaLU[2],pavaUL,pavaLU[1])
              if((length(v)-length(pavaAux))!=0){
                #print("#######################################################error")
              }
              mseAux<-sum((v-pavaAux)^2)/length(v)
              if(pavaLU[2]>v[candL[i]] & pavaLU[length(pavaLU)-1]<v[candU[j]] & mseAux<mseFin){
                mseFin<-mseAux
                pavaFin<-pavaAux
                Lopt<-candL[i]
                Uopt<-candU[j]
              }
            }
          }
        }
        
      }
    }
  }
  return(list(pavaFin,mseFin,candL,candU,Lopt,Uopt))
}


















