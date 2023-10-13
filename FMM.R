# Funcón para representar de forma sencilla un FMM
# Recibe como argumentos:
#  - M --> Valor para M
#  - A --> valor  para A
#  - a --> valor para alpha
#  - b --> valor para beta
#  - w --> valor para omega
#  - from --> instante temporal de inicio. Por defecto es 0.
#  - to --> instante temporal de fin Por defecto es 2*pi.
#  - timePoints --> También es posible definir unos instantes de tiempo prefijados
#  - plot --> en caso afirmativo, se dispondrá por la salida de gráficos estándar de una
#             representación gráfica del FMM
#  - outvalues --> en caso afirmativo, se devolverá una lista con los instantes temporales 
#                  y los valores producidos.
#  - length.out --> longitud de la simulación
#  - sigmaNoise --> en caso afirmativo, se añade un ruido normal con s.d. especificada.
# En lugar de establecer valores individuales para M, A, alpha, beta y omega,
# se puede establecer un array. En este caso se genera una combinación lineal de FMMs.
# Se añaden los que faltan por replicación
FMM <- function(M,A,a,b,w,from=0,to=2*pi,length.out=100,timePoints=seq(from,to,length=length.out),
                plot=T,outvalues=F,sigmaNoise=0){
  
  narg <- max(length(M),length(A),length(a),length(b),length(w))
  M <- rep(M,length.out=narg)
  A <- rep(A,length.out=narg)
  a <- rep(a,length.out=narg)
  b <- rep(b,length.out=narg)
  w <- rep(w,length.out=narg)
  
  t <- timePoints
  
  phi <- list()
  for(i in 1:narg){
    phi[[i]] <- b[i]+2*atan(w[i]*tan((t-a[i])/2))
  }
  
  ym <- list()
  for(i in 1:narg){
    ym[[i]] <- M[i]+A[i]*cos(phi[[i]])
  }
  
  y <- rep(0,length(t))
  for(i in 1:narg){
    y <- y + ym[[i]]
  }
  
  if (sigmaNoise > 0) y <- y + rnorm(length.out,0,sigmaNoise)
  
  if(plot) plot(t,y,type="l",lwd=2,col=2,xlab="tiempo",ylab="Respuesta",
       main=paste("FMM:   M=",M,"   A=",A,"   a=",a,"   b=",b,"   w=",w,sep=""))
  
  if(outvalues) return(list(t=t,y=y))
}


predict_FMM <- function(M,A,alpha,beta,omega,timePoints){
  return(FMM(M,a,alpha,beta,omega,timePoints = timePoints,plot=F,outvalues = T)$y)
}


# Función para ajustar un modelo FMM sobre unos datos.
# Recibe como entrada los datos y los instantes de tiempo donde se observan:
#  vData --> valores para el ajuste
#  timePoints --> instantes de tiempo en el que se produce el ajuste.
#                  Por defecto están equiespaciados
#  lengthAlphaGrid --> Longitud del GRID utilizado para alpha. Debe indicarse con un número entero.
#                     No es necesario indicarlo si se va a utilizar un GRID personalizado mediante
#                     el argumento alphaGrid. Por defecto vale 24.
#  lengthOmegaGrid --> Longitud del GRID utilizado para omega Debe indicarse con un número entero.
#                     No es necesario indicarlo si se va a utilizar un GRID personalizado mediante
#                     el argumento omegaGrid. Por defecto vale 24.
#  alphaGrid --> Valores de alpha utilizados para la estimación inicial
#                Si se conoce el valor aproximado de alpha, es recomendable indicarlo aquí
#                Por defecto, equiespaciados entre 0 y 2pi, utilizando 24 elementos
#  omegaMax --> valor máximo que puede tomar el parámetro omega. En ningún caso se devuelve
#               un omega estimado mayor del especificado.
#               Si se conoce el valor aproximado de omega, es recomendable indicar aquí una cota
#  omegaGrid --> Valores de omega utilizados para la estimación inicial
#                Por defecto toma valores entre 0.0001 y 1 en escala logarítmica
#  numReps --> número de veces que se repite el paso 1 + 2. La salida del paso 2 se toma como una nueva entrada
#              en torno a la que construir un nuevo GRID. Se repite la estimación un total de numReps.
#              Por defecto vale 3, es decir, se hace tres veces.
# REESCRITA CON RESPECTO A LA ORIGINAL
# PROBABLEMENTE HAYA QUE CAMBIARLE EL NOMBRE
fitFMM_Par<-function(vData, timePoints = seq(0,2*pi,length.out = length(vData)),
                     lengthAlphaGrid = 48, lengthOmegaGrid = 24,
                     alphaGrid = seq(0,2*pi,length.out = lengthAlphaGrid), omegaMax = 1,
                     omegaGrid = exp(seq(log(0.0001),log(omegaMax),length.out=lengthOmegaGrid)),
                     numReps = 3){
  
  n <- length(vData)
  
  ## Step 1: Calcular los estimadores iniciales de M, A, alpha, beta y omega
  #Para ello se fija alpha y omega y se calculan los demás utilizando un modelo Cosinor
  grid <- expand.grid(alphaGrid,omegaGrid)
  
  #Para evitar un doble bucle, ineficientes en R, se utiliza la función apply con una función auxiliar
  #Esta función auxiliar se encarga de la estimación de beta, M, A, y también devuelve como añadido el RSS
  #para facilitar la búsqueda del mejor a posteriori
  step1 <- t(apply(grid,1,step1FMM, vData=vData, timePoints=timePoints))
  colnames(step1) <- c("M","A","alpha","beta","omega","RSS")
  
  #De todos los ajustados, buscamos el que tiene un menor RSS pero que cumpla ciertas
  #condiciones de estabilidad. Para esto se usa la función bestStep1.
  bestPar <- bestStep1(vData,step1)
  #Los estimadores iniciales para M, A, alpha, beta y omega ya se han encontrado
  
  #Step-2: Optimización de Nelder-Mead para el RSS
  nelderMead <- optim(par = bestPar[1:5], fn = step2FMM, vData = vData, timePoints = timePoints, omegaMax = omegaMax)
  parFinal <- nelderMead$par
  SSE <- nelderMead$value*n
  
  #Llevar alpha y beta al intervalo [0,2pi]
  parFinal[3] <- parFinal[3]%%(2*pi)
  parFinal[4] <- parFinal[4]%%(2*pi)
  
  #Se repite el GRID el número de veces especificado (por defecto NO se aplica esta opción)
  numReps <- numReps - 1
  while(numReps > 0){
    
    #Nuevo GRID para alpha, de la misma longitud del anterior. Se asegura que esté entre 0 y 2pi.
    nAlphaGrid <- length(alphaGrid)
    amplitudeAlphaGrid <- 1.5*mean(diff(alphaGrid))
    alphaGrid <- seq(parFinal[3]-amplitudeAlphaGrid,parFinal[3]+amplitudeAlphaGrid,length.out = nAlphaGrid)
    alphaGrid <- alphaGrid%%(2*pi)
    
    #Nuevo GRID para omega, de la misma longitud del anterior. Se asegura que esté entre 0 y omegaMax
    nOmegaGrid <- length(omegaGrid)
    amplitudeOmegaGrid <- 1.5*mean(diff(omegaGrid))
    omegaGrid <- seq(max(parFinal[5]-amplitudeOmegaGrid,0),
                     min(omegaMax,parFinal[5]+amplitudeOmegaGrid),
                     length.out = nOmegaGrid)
    
    
    ## Step 1
    grid <- as.matrix(expand.grid(alphaGrid,omegaGrid))
    step1 <- t(apply(grid,1,step1FMM, vData=vData, timePoints=timePoints))
    colnames(step1) <- c("M","A","alpha","beta","omega","RSS")
    antBestPar <- bestPar
    bestPar <- bestStep1(vData,step1)
    
    #Es posible que no haya ninguno que satisfaga las condiciones
    if(is.null(bestPar)){
      bestPar <- antBestPar
      numReps <- 0
      warning("Es posible que el FMM no sea adecuado en este caso")
    }
    
    
    ## Step 2
    nelderMead <- optim(par = bestPar[1:5], fn = step2FMM, vData = vData, timePoints = timePoints, omegaMax = omegaMax)
    parFinal <- nelderMead$par
    
    #Llevar alpha y beta al intervalo [0,2pi]
    parFinal[3] <- parFinal[3]%%(2*pi)
    parFinal[4] <- parFinal[4]%%(2*pi)
    
    
    numReps <- numReps - 1
  }
  
  names(parFinal) <- c("M","A","alpha","beta","omega")
  
  #Devuelve, además de los parámetros, un ajuste en los intervalos de tiempo especificados 
  #y el SSE (sum squared error).
  adjMob <- parFinal["M"] + parFinal["A"]*cos(parFinal["beta"] + 2*atan(parFinal["omega"]*tan((timePoints-parFinal["alpha"])/2)))
  SSE <- sum((adjMob-vData)^2)
  
  outMobius<-list(FMM=adjMob,alpha=parFinal[3],beta=parFinal[4],omega=parFinal[5],M=parFinal["M"],
                  A=parFinal["A"], SSE = SSE)
  return(outMobius)
}



# Función para ajustar un modelo FMM múltiple (backfitting) sobre unos datos.
# Recibe como entrada los datos y los instantes de tiempo donde se observan:
#  vData --> valores para el ajuste
#  timePoints --> instantes de tiempo en el que se produce el ajuste.
#                  Por defecto están equiespaciados
#  nback --> número de componentes a ajustar mediante backfitting
#  maxiter --> número máximo de iteraciones a realizar. Por defecto vale lo mismo que nback.
#  stopFunction --> parámetro opcional. Debe especificarse una función que recibe tres argumentos:
#                   los valores reales (vData), el valor ajustado por el FMM en la iteración actual,
#                   y el valor ajustado por el FMM en la anterior iteración. Debe devolver TRUE si se desea
#                   detener el backfitting, o FALSE si debe continuar porque aún no se ha alcanzado la convergencia.
#                   Por defecto siempre devuelve FALSE: stopFunction = function(vData,pred,prevPred){return(FALSE)}
#  objectFMM --> un objeto FMM, como máximo de nback componentes, sobre el que continuar realizando iteraciones.
#                Esta opción permite continuar un ajuste a partir de un objetoFMM anterior
#                ahorrando en tiempo y recursos al no comenzar desde el principio.
#                Si el número de componentes no es el mismo, los demás parámetros se rellenan con cero.
#  lengthAlphaGrid --> Longitud del GRID utilizado para alpha. Debe indicarse con un número entero.
#                     No es necesario indicarlo si se va a utilizar un GRID personalizado mediante
#                     el argumento alphaGrid. Por defecto vale 24.
#  lengthOmegaGrid --> Longitud del GRID utilizado para omega Debe indicarse con un número entero.
#                     No es necesario indicarlo si se va a utilizar un GRID personalizado mediante
#                     el argumento omegaGrid. Por defecto vale 24.
#  alphaGrid --> Valores de alpha utilizados para la estimación inicial
#                Si se conoce el valor aproximado de alpha, es recomendable indicarlo aquí
#                Por defecto, equiespaciados entre 0 y 2pi, utilizando 24 elementos
#  omegaMax --> valor máximo que puede tomar el parámetro omega. En ningún caso se devuelve
#               un omega estimado mayor del especificado.
#               Si se conoce el valor aproximado de omega, es recomendable indicar aquí una cota
#  omegaGrid --> Valores de omega utilizados para la estimación inicial
#                Por defecto toma valores entre 0.0001 y 1 en escala logarítmica
#  numReps --> número de veces que se repite el paso 1 + 2. La salida del paso 2 se toma como una nueva entrada
#              en torno a la que construir un nuevo GRID. Se repite la estimación un total de numReps.
#              Por defecto vale 3, es decir, se hace tres veces.
# REESCRITA CON RESPECTO A LA ORIGINAL
# PROBABLEMENTE HAYA QUE CAMBIARLE EL NOMBRE
fitFMM_back<-function(vData, timePoints = seq(0,2*pi,length.out = length(vData)), nback, maxiter=nback,
                      stopFunction = function(vData,pred,prevPred){return(FALSE)}, objectFMM = NULL,
                     lengthAlphaGrid = 48, lengthOmegaGrid = 24,
                     alphaGrid = seq(0,2*pi,length.out = lengthAlphaGrid), omegaMax = 1,
                     omegaGrid = exp(seq(log(0.0001),log(omegaMax),length.out=lengthOmegaGrid)),
                     numReps = 3, showProgress = T){
  
  n <- length(vData)
  
  if(showProgress){
    marcasTotales <- 50
    granoInforme <- 2
    cat("|")
    for(m in 1:marcasTotales) cat("-")
    cat("|\n")
    cat("|")
    porcentajeCompletado <- 0.00001
  }
  
  #Inicialización de las estructuras para almacenar las componentes.
  predichosComponente <- list()
  ajusteComponente <- list()
  
  #En caso de comenzar desde el principio (no se aporta un objectFMM)
  if(is.null(objectFMM)){
    for(i in 1:nback){
      predichosComponente[[i]] <- rep(0,n)
    }
    prevAdjMob <- NULL
    
  #En caso de que se aporte un objectFMM anterior para seguir con el ajuste
  } else {
    prevAdjMob <- objectFMM$FMM
    nbackAnterior <- length(objectFMM$alpha)
    if(nbackAnterior > nback){
      stop("Impossible to reduce dimensions from input objectFMM")
    }
    for(i in 1:nback){
      if(i <= nbackAnterior){
        predichosComponente[[i]] <- objectFMM$M/nbackAnterior + objectFMM$A[i]*cos(objectFMM$beta[i] + 
                                      2*atan(objectFMM$omega[i]*tan((timePoints-objectFMM$alpha[i])/2)))
      } else {
        predichosComponente[[i]] <- rep(0,n)
      }
    }
  }
  
  
  #Backfitting hasta un máximo de maxiter
  for(i in 1:maxiter){
    
    #Backfitting por componentes
    for(j in 1:nback){
      
      #Obtenemos la componente j a ajustar como la diferencia de los datos originales con respecto
      #al ajuste de todas las demás componentes
      vDataAjuste <- vData
      for(k in 1:nback){
        if(j != k){
          vDataAjuste <- vDataAjuste - predichosComponente[[k]]
        }
      }
      
      #Realizamos el ajuste de la componente j con las opciones especificadas
      ajusteComponente[[j]] <- fitFMM_Par(vDataAjuste,timePoints = timePoints, lengthAlphaGrid = lengthAlphaGrid,
                                          lengthOmegaGrid = lengthOmegaGrid, alphaGrid = alphaGrid, omegaMax = omegaMax,
                                          omegaGrid = omegaGrid, numReps = numReps)
      predichosComponente[[j]] <- ajusteComponente[[j]]$FMM
      
      #Esto se utiliza para mostrar progreso visible
      if(showProgress){
        porcentajeAntes <- porcentajeCompletado
        porcentajeCompletado <- porcentajeCompletado + 100/(nback*maxiter)
        numMarcas <- sum((seq(ceiling(porcentajeAntes),floor(porcentajeCompletado),by=1)%%granoInforme == 0))
        if (numMarcas > 0) for(m in 1:numMarcas) cat("=")
      }
      
    }
    
    #Comprobación de la condición de parada
    #Calculamos los valores predichos por el modelo como la suma de todas las componentes
    adjMob <- rep(0,n)
    for(j in 1:nback){
      adjMob <- adjMob + predichosComponente[[j]]
    }
    if(!is.null(prevAdjMob)){
      if(stopFunction(vData,adjMob,prevAdjMob)){
        break
      }
    }
    prevAdjMob <- adjMob
    
    
  }
  
  #Esto se utiliza para mostrar progreso visible
  if(showProgress){
    if(porcentajeCompletado < 100){
      porcentajeAntes <- porcentajeCompletado
      porcentajeCompletado <- 100
      numMarcas <- sum((seq(ceiling(porcentajeAntes),floor(porcentajeCompletado),by=1)%%granoInforme == 0))
      for(m in 1:numMarcas) cat("=")
    }
    cat("|\n")
    if(i == maxiter){
      cat("Stopped by reaching maximum iterations\n")
    } else {
      cat("Stopped by the stopFunction\n")
    }
  }

  #El parámetro M es común a todas las componentes (indica el desplazamiento)
  #Debe ser así para evitar singularidades
  #Se toma el parámetro M como la suma de todos los M de cada componente
  #Los parámetros A, alpha, beta y omega son únicos de cada componente
  alpha <- rep(0,nback)
  beta <- rep(0,nback)
  omega <- rep(0,nback)
  for(j in 1:nback){
    alpha[j] <- ajusteComponente[[j]]$alpha
    beta[j] <- ajusteComponente[[j]]$beta
    omega[j] <- ajusteComponente[[j]]$omega
  }
  
  #Se recalculan A y M con una regresion lineal
  phi <- list()
  for(j in 1:nback){
    phi[[j]] <- cos(beta[j] + 2*atan(omega[j]*tan((timePoints-alpha[j])/2)))
  }
  M <- matrix(unlist(phi),ncol=nback)
  regresion <- lm(vData ~ M)
  M <- coefficients(regresion)[1]
  A <- coefficients(regresion)[-1]
  
  #Calculamos los valores predichos por el modelo
  adjMob <- predict(regresion)

  #Calculo del SSE
  SSE <- sum((adjMob-vData)^2)

  
  #Se devuelve una salida similar al caso sin backfitting
  outMobius<-list(FMM=adjMob,alpha=alpha,beta=beta,omega=omega,M=M,A=A,SSE = SSE)
  return(outMobius)

}


#dropComponents_FMM <- function(objectFMM,numOfComponents=0,components=NULL,revise=T){
  
  
#  FMM(20,2,c(1.5,3.4,5),c(0.2,2.3,4.7),c(0.1,0.2,0.05))
#  res <- FMM(20,2,c(1.5,3.4,5),c(0.2,2.3,4.7),c(0.1,0.2,0.05),outvalues = T,plot=F,sigmaNoise = 0.1)
#  fit <- fitFMM_back(res$y,res$t,nback = 4)
#  plot(res$t,res$y,xlab="Tiempo",ylab="Respuesta")
#  points(res$t,fit$FMM,col=2,type="l",lwd=2)
#  plot_FMM(fit,res$y,res$t,components = T)
  
#  nback <- length(objectFMM$alpha)
#  if(is.null(components)){
#    components <- seq(len)
#  }
  
#}


# Función auxiliar que calcula el porcentaje de variabilidad recogido
# Recibe como entradas:
#   vData --> datos originales
#   pred --> datos predichos por el modelo
# FUNCTION OCULTA
PV <- function(vData,pred){
  meanVData <- mean(vData)
  return(1 - sum((vData-pred)^2)/sum((vData-meanVData)^2))
}

# Función posible para utilizar como criterio de parada en el backfitting.
# Devuelve TRUE (parada) si se cumple al menos una de estas dos condiciones:
# cond1 -> ((PV(i) - PV(i-1)) > 0.0001) & ((PV(i)-PV(i-1)) < 0.025*(1-PV(i-1)))
# cond2 -> (PV(i) - PV(i-1)) <= 0.0001
# Donde PV(i) es el porcentaje de variabilidad recogido en la iteración i
# Recibe como entradas:
#   vData --> datos originales
#   pred --> valores predichos en la iteración actual
#   prevPred --> valores predichos en la iteración anterior
stopFunction1 <- function(vData,pred,prevPred){
  PV_pred <- PV(vData,pred)
  PV_prevPred <- PV(vData,prevPred)
  cond1 <- ((PV_pred-PV_prevPred) > 0.0001) & ((PV_pred-PV_prevPred) < 0.025*(1-PV_prevPred))
  cond2 <- (PV_pred-PV_prevPred) <= 0.0001
  if(cond1 | cond2){
    return(TRUE)
  } else {
    return(FALSE)
  }
}


# Función auxiliar para el primer paso de la estimación FMM
# Recibe como entradas:
#   alphaOmega --> un vector que contiene dos parámetros: alpha y omega
#   vData --> valores para el ajuste
#   timePoints --> instantes de tiempo en el que se produce el ajuste.
#La salida producida es un vector de seis componentes: M, A, alpha, beta, omega y RSS
# FUNCION OCULTA
step1FMM <- function(alphaOmega, vData, timePoints) {
  
  
  alpha <- alphaOmega[1]
  omega <- alphaOmega[2]
  
  parteMobius <- 2*atan(omega*tan((timePoints-alpha)/2))
  t_star <- alpha + parteMobius
  
  #Modelo cosinor suponiendo fijos alpha y omega
  xx<-cos(t_star)
  zz<-sin(t_star)
  fit<-lm((vData)~xx+zz)
  M <-fit$coefficients[1] #intercept
  bb<-fit$coefficients[2] #Coeficiente del coseno
  gg<-fit$coefficients[3] #Coeficiente del seno
  phiEst<-atan2(-gg,bb) #acrophase (phi)
  A<-sqrt(bb^2+gg^2) #Amplitud de la onda
  beta <- (phiEst+alpha)%%(2*pi)
  
  dataReg<-cos(beta+parteMobius) #Resultado de Mobius sin amplitud e intercept
  
  adj0<-M+A*dataReg #Mobius Reg
  RSS<-sum((vData-adj0)^2)/length(timePoints) #suma de cuadrados de residuos (medio)
  
  devolver <- c(M,A,alpha,beta,omega,RSS)
  return(devolver)
  
}




#Función para uso interno que permite encontrar el parámetro con menor RSE. Toma como entradas:
#  vData --> el vector que contiene los datos originales
#  step1 --> la salida de t(apply(grid,1,step1FMM, vData=vData, timePoints=timePoints)),
#            que es un data.frame que tiene por columnas estimaciones de M, A, alpha, beta, omega, RSS
#  Devuelve la fila de step1 que se haya elegido.
# FUNCION OCULTA
bestStep1 <- function(vData,step1){
  
  #De todos los ajustados, buscamos el que tiene un menor RSS pero que cumpla ciertas
  #condiciones de estabilidad. 
  ordenMinRSS <- order(step1[,"RSS"])
  maxVData <- max(vData)
  minVData <- min(vData)
  n <- length(vData)
  
  #Empezamos desde el que tiene un RSS menor y vamos buscando el primero que cumpla las condiciones
  condicionContinuar <- TRUE
  i <- 1
  while(condicionContinuar){
    
    #Recuperar los parámetros y el RSS
    M <- step1[ordenMinRSS[i],"M"]
    A <- step1[ordenMinRSS[i],"A"]
    alpha <- step1[ordenMinRSS[i],"alpha"]
    beta <- step1[ordenMinRSS[i],"beta"]
    omega <- step1[ordenMinRSS[i],"omega"]
    sigma <- sqrt(step1[ordenMinRSS[i],"RSS"]*n/(n-5))
    
    #Restricciones de estabilidad
    maxi <- M + A
    mini <- M - A
    rest1 <- maxi <= maxVData+0.1*(range(vData)[2]-range(vData)[1])#1.96*sigma
    rest2 <- mini >= minVData-0.1*(range(vData)[2]-range(vData)[1])#1.96*sigma
    
    #Si se cumple, parar aquí. Hay que comprobar que todos son distintos de NA, porque puede ser un ajuste extremo
    if(is.na(rest1)) rest1 <- FALSE
    if(is.na(rest2)) rest2 <- FALSE
    
    if(rest1 & rest2){
      condicionContinuar <- FALSE
    } else {
      i <- i+1
    }
    
    if(i > nrow(step1)){
      return(NULL)
    }
  }
  
  return(step1[ordenMinRSS[i],])
  
}

# Función que realiza el segundo paso del FMM. En este caso, es una implementación
# de la función objetivo a minimizar con Nelder-Mead. Calcula el RSS dados los parámetros
# y los datos. Los argumentos que recibe son:
#   param --> un vector que contiene, en orden, M, A, alpha, beta, omega
#   vData --> los datos para ser ajustados
#   timePoints --> instantes de tiempo donde se ajustan los datos
# FUNCION OCULTA
step2FMM <- function(param, vData, timePoints, omegaMax){
  
  n <- length(timePoints)
  
  #Ajuste FMM
  ffMob<-param[1] + param[2] * cos(param[4]+2*atan2(param[5]*sin((timePoints-param[3])/2),cos((timePoints-param[3])/2)))
  
  #La media de la suma de residuos al cuadrado
  RSS<-sum((ffMob - vData)^2)/n
  sigma <- sqrt(RSS*n/(n-5))
  
  #En caso de que se cumpla la condición de la amplitud, se devuelve el RSS.
  #En caso contrario, se devuelve infinito.
  maxi<-param[1]+param[2]
  mini<-param[1]-param[2]
  rest1 <- maxi <= max(vData)+0.1*(range(vData)[2]-range(vData)[1])#1.96*sigma
  rest2 <- mini >= min(vData)-0.1*(range(vData)[2]-range(vData)[1])#1.96*sigma
  
  #Otras condiciones de integridad que deben cumplirse
  rest3 <- param[2] > 0 #A > 0
  rest4 <- param[5] > 0 #omega > 0
  rest5 <- param[5] <= omegaMax #omega <= omegaMax
  if(rest1 & rest2 & rest3 & rest4 & rest5){
    return(RSS)
  }else{
    return(Inf)
  }
  
}


#Dibuja el FMM. Si components = F, dibuja las componentes centradas sin los datos.
plot_FMM <- function(object_FMM,vData,timePoints = seq(0,2*pi,length.out = length(vData)),components=F){
  
  if(components){
    
    nc <- length(object_FMM$alpha)
    predicted <- list()
    
    for(i in 1:nc){
      predicted[[i]] <- object_FMM$A[i]*cos(object_FMM$beta[i] + 2*atan(object_FMM$omega[i]*tan((timePoints-object_FMM$alpha[i])/2)))
    }
    
    mini <- min(vData)
    maxi <- max(vData)
    for(i in 1:nc){
      #predicted[[i]] <- predicted[[i]] - median(predicted[[i]])
      predicted[[i]] <- predicted[[i]] - predicted[[i]][1] + vData[1]
      mini <- min(mini,min(predicted[[i]]))
      maxi <- max(maxi,max(predicted[[i]]))
    }
    
    plot(timePoints,vData,ylim=c(mini,maxi),xlab="Time",ylab="Response",main="Components FMM",type="n")
    for(i in 1:nc){
      points(timePoints,predicted[[i]],type="l",lwd=2,col=i+1)
    }
    
  } else {
    
    plot(timePoints,vData,xlab="Time",ylab="Response",main="Adjusted FMM")
    points(timePoints,object_FMM$FMM,type="l",col=2,lwd=2)
    
  }
  
}

