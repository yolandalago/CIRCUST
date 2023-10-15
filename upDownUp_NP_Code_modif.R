source("upDownUp_NP_Code_Paralelizado.R")
source("upDownUp_NP_Code_NoParalelizado.R")
source("NucleoComun.R")

leerTamRentableVersionParalela <- function(){
  tryCatch({
    load("tamRentableVersionParalela.RData")
    return(tamRentableVersionParalela)
  }, warning = function (w){
    return(NA)
  }, error = function(e){
    return(NA)
  })
}

leernRentableVersionMultiple <- function(){
  tryCatch({
    load("tamRentableVersionMultiple.RData")
    return(nRentableVersionMultiple)
  }, warning = function (w){
    return(NA)
  }, error = function(e){
    return(NA)
  })
}

leerPeriodoRentableVersionMultiple <- function(){
  tryCatch({
    load("tamRentableVersionMultiple.RData")
    return(periodoRentableVersionMultiple)
  }, warning = function (w){
    return(NA)
  }, error = function(e){
    return(NA)
  })
}

leerMaxCores <- function(){
  tryCatch({
    load("maxCores.RData")
    return(maxCores)
  }, warning = function (w){
    return(NA)
  }, error = function(e){
    return(NA)
  })
}

tamRentableVersionParalela <- leerTamRentableVersionParalela()
nRentableVersionMultiple <- leernRentableVersionMultiple()
periodoRentableVersionMultiple <- leerPeriodoRentableVersionMultiple()
maxCores <- leerMaxCores()

#Funci?n para el c?lculo del mejor estimador Up-Down-Up
function1Local_modif<-function(v, parallel = "auto"){
  
  if(parallel == "si"){
    return(function1Local_par(v))
    
  } else if (parallel == "no"){
    return(function1Local_sec(v))
    
  } else if (parallel == "auto") {
    if( !is.na(tamRentableVersionParalela) && length(v) > tamRentableVersionParalela){
      return(function1Local_par(v))
    } else {
      return(function1Local_sec(v))
    }
    
  } else {
    stop("Value of parallel not allowed\n")
  }
  
}


#Ejecuci?n en paralelo siempre
function1Local_multiple <- function(datos){
  
  n <- nrow(datos)
  
  if( !is.na(nRentableVersionMultiple) && n < nRentableVersionMultiple){
    aux <- n >= as.numeric(names(periodoRentableVersionMultiple))
    perAux <- periodoRentableVersionMultiple[aux]
    tamRentable <- perAux[length(perAux)]
    if(all(aux) == FALSE || ncol(datos) < tamRentable){
      cat("La version multiple podria no ser adecuada para tan pocos datos")
      flush.console()
    }
  }
  ncores <- min(detectCores()-1,n)
  if(is.na(maxCores) == FALSE){
    ncores <- min(ncores,maxCores)
  }
  programacion <- CutBySize(n,nb=ncores)
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  
  resultados <- foreach(pr = iter(programacion,by='row'), 
                        .inorder=T, .combine = c, 
                        .export = c("function1Local_sec","extremosLocales","busquedaMejorC")) %dopar% {
    
      dyn.load("busquedaMejor.dll")
      datosInt <- datos[pr[1]:pr[2],]
      resultadosInt <- list()
    if(is.null(nrow(datosInt))){
      resultadosInt[[1]] <- function1Local_sec(datosInt)
    } else {
      for(i in 1:nrow(datosInt)){
        resultadosInt[[i]] <- function1Local_sec(datosInt[i,])
      }
    }
   return(resultadosInt)
    
  }
  
  stopCluster(cl)
  registerDoSEQ()
  return(resultados)
  
}


