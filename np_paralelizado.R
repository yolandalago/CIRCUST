library(doParallel)

CutBySize <- function(m, block.size, nb = ceiling(m/block.size)) {
  if(nb > m)
    nb <- m
  int <-  m/nb
  upper <- round(1:nb * int)
  lower <- c(1, upper[-nb]+1)
  size <- c(upper[1], diff(upper))
  return(cbind(lower,upper,size))
}

#Funci?n para el c?lculo del mejor estimador Up-Down-Up
function1Local_par<-function(v){
  
  
  #Buscamos los m?ximos y m?nimos locales que existen
  extremos <- extremosLocales(v)
  candL <- extremos$candL
  candU <- extremos$candU
  
  #N?mero de cores a usar en la paralelizaci?n de la b?squeda
  ncores <- min(detectCores()-1,length(candL))
  if(is.na(maxCores) == FALSE){
    ncores <- min(ncores,maxCores)
  }
  programacion<-CutBySize(length(candL),nb = ncores)
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  
  #Divisi?n del trabajo en los procesos hijos
  mejores <- foreach(pr = iter(programacion,by='row'), .export = c("busquedaMejorC")) %dopar% {
    
    dyn.load("pava.dll")
    dyn.load("busquedaMejor.dll")

    candLInt <- candL[pr[1]:pr[2]]

    mejor <- busquedaMejorC(v,candLInt,candU)
    return(mejor)
  }

  stopCluster(cl)
  
  #Busamos el mejor de entre los mejores
  mseFin <- Inf
  for(elemento in mejores){
    if(!is.null(elemento) && elemento$mse < mseFin){
      mseFin<-elemento$mse
      pavaFin<-elemento$pava
      Lopt<-elemento$Lopt
      Uopt<-elemento$Uopt
    }
  }
  
  registerDoSEQ()
  return(list(pava=pavaFin,mse=mseFin,candL=candL,candU=candU,Lopt=Lopt,Uopt=Uopt))
}
