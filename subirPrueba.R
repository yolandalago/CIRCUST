

P<-"PER1"
M<-mLiverNorm[1:2,]
o<-orderCPCALiver8
p<-escalaPhiLiver8
cyc<-cycGenes[1:2]


x11()
par(mfrow=c(1,1))
plot(p,M[match(P,cyc),o],type="p")
lines(p,fittingRef[[1]],col=4)
x11()
plot((pN3-pN3[1])%%(2*pi),vN3,type="b")
lines((pN3-pN3[1])%%(2*pi),fit3N[[1]],col=4)
reajustarPlotsParam_v2<-function(P,M,o,p,cyc){
	muevoR<-FALSE;muevoL<-FALSE
	mEstPar3N<-matrix(0,nrow(M),10)
	#order de referencia ORIN
	muevoR<-FALSE;muevoL<-FALSE
	#gen de referencia para saber sii sube primero y luego baja o no, Siempre orden ORIN
	#Fijandonos en el ajuste para metrico 
	fittingRef<-fitFMM_Par2(M[match(P,cyc),o],p)
	al<-fittingRef[[2]]
	be<-fittingRef[[3]]
	om<-fittingRef[[4]]
	LN<-(((al+2*atan2(1/om*sin((pi-be)/2),cos((pi-be)/2))))%%(2*pi))/(2*pi-(2*pi)/ncol(M))*ncol(M)+1#compL(al=all,be=bee,om=omm,tod01Kidny=length(o[[1]]))
	UN<-(((al+2*atan2(1/om*sin(-be/2),cos(-be/2))))%%(2*pi))/(2*pi-(2*pi)/ncol(M))*ncol(M)+1
	
	#time Point medio 
	medio<-length(o)/2+0.5
	if((length(fittingRef[[1]])%%2)==0){
		a<-p[medio-0.5]/(2*pi-(2*pi)/ncol(M))*ncol(M)+1
		b<-p[medio+0.5]/(2*pi-(2*pi)/ncol(M))*ncol(M)+1
		if(((medio-a)^2)<((medio-b)^2)){
			medio<-medio-0.5
		}else{
			medio<-medio+0.5
		}
	}
	if(UN<medio){
		muevoR<-TRUE
		cuantos<-trunc((medio-UN))
		inicio<-(medio-UN)-trunc((medio-UN))
	}else{
		muevoL<-TRUE
		cuantos<-trunc(UN-medio)
		inicio<-(UN-medio)-trunc(UN-medio)
	}

	
	vN3<-c()
	oN3<-c()
	vN<-M[match(P,cyc),o]
	vN2<-c(vN,vN)
	oN<-o
	oN2<-c(oN,oN)
	vN3<-vN
	oN3<-oN
	pN2<-c(p,p)
	if(muevoR){
		if(cuantos>0){
			vNRemp<-vN2[length(vN2):(length(vN2)-cuantos+1)]
			oNRemp<-oN2[length(vN2):(length(vN2)-cuantos+1)]
			vN3[1:cuantos]<-rev(vNRemp)
			vN3[(cuantos+1):length(vN)]<-vN2[1:(length(vN)-cuantos)]
			oN3[1:cuantos]<-rev(oNRemp)
			oN3[(cuantos+1):length(vN)]<-oN2[1:(length(oN)-cuantos)]
			x<-which(oN3[1]==oN2,arr.ind=TRUE)[1]
			pN3<-pN2[x:(x+length(oN)-1)]
		}
		fit3N<-fitFMM_Par2(vN3,(pN3-pN3[1])%%(2*pi))
	}
	if(muevoL){
		if(cuantos>0){
			vN3[(length(vN)-cuantos+1):length(vN)]<-vN2[1:cuantos]
			vN3[1:(length(vN)-cuantos)]<-vN2[(cuantos+1):length(vN)]
			oN3[(length(oN)-cuantos+1):length(oN)]<-oN2[1:cuantos]
			oN3[1:(length(oN)-cuantos)]<-oN2[(cuantos+1):length(oN)]
			x<-which(oN3[1]==oN2,arr.ind=TRUE)[2]
			pN3<-oN2[(x+length(oN)+1):x]
		}
		fit3N<-fitFMM_Par(vN3,pN3)
	}
	fitAllN<-matrix(0,nrow(M),ncol(M))
	mFitP3N<-matrix(0,nrow(M),ncol(M))
	for(i in 1:nrow(M)){
		vi<-M[i,oN3]
		if( i!=match(P,cyc))fitAuxN<-fitFMM_Par2(vi,pN3)
		if(i==match(P,cyc))fitAuxN<-fit3N
		fitAllN<-fitAuxN[[1]]
		al<-fitAuxN[[2]];be<-fitAuxN[[3]];om<-fitAuxN[[4]];M<-fitAuxN[[5]];A<-fitAuxN[[6]]
		peak<-compU(al,be,om,vi);trough<-compL(al,be,om,vi);peakRel<-peak/length(vi);troughRel<-trough/length(vi)
		sigmaPar<-sum((vi-fitAllN)^2)/(length(vi)-5)
		mEstPar3N[i,]<-c(M,A,al,be,om,peak,trough,peakRel,troughRel,sigmaPar)
		colnames(mEstPar3N)<-c("M","A","al","be","om","peak","trough","peakRel","troughRel","sigmaPar")
		mFitP3N[i,]<-fitAllN
	}
	return(list(oN3,pN3,mEstPar3N,mFitP3N))
}
