
#Función para el cálculo de la regresión isotónica simple
pavaC <- function(y, w = NULL, decreasing = F){
  
  if (any(is.na(y))) 
    stop("missing values not allowed")
  if ((!is.null(w)) && (length(w) != length(y))){
    stop("length of weights and points differ")
  }
  
  if(is.null(w)){
    w <- rep(1,length(y))
  }
  
  resultado <- rep(0, length(y))
  n<-length(y)
  
  if(decreasing){
    return(.C("pavaCAdriDecreasing",nin = as.integer(n),y = as.double(y),w = as.double(w),
              resultado = as.double(resultado))$resultado)
  } else {
    return(.C("pavaCAdri",nin = as.integer(n),y = as.double(y),w = as.double(w),
              resultado = as.double(resultado))$resultado)
  }
}

#Determina si el elemento está en el vector
contiene <- function(elemento, vector){
  for(e in vector){
    if(elemento == e) return (TRUE)
  }
  return (FALSE)
}

#Encontrar los mínimos y máximos locales
extremosLocales <- function(v){
  candL<-c()
  candU<-c()
  for(i in 1:length(v)){
    if(i == length(v)){
      iant <- i-1
      ipos <- 1
    } else if (i == 1){
      iant <- length(v)
      ipos <- i+1
    } else {
      iant <- i-1
      ipos <- i+1
    }
    #minimos locales
    if(v[iant]<=v[i] & v[i]>=v[ipos]){
      candU<-c(candU,i)
    }
    #maximos locales
    if(v[iant]>=v[i] & v[i]<=v[ipos]){
      candL<-c(candL,i)
    }
  }
  return(list(candL=candL,candU=candU))
}

busquedaMejorC <- function(v, candL, candU){
  
  if(length(v) == 0){
    stop("Empty data")
  }
  
  if(length(candL)==0 || length(candU)==0){
    return(NULL)
  }
  
  sol <- .Call("busquedaMejor",as.double(v),as.integer(candL-1),as.integer(candU-1))
  
  if(sol$mseFin < 0){
    return(NULL)
  } else {
    devolver <- list(pava = sol$pavaFin, mse = sol$mseFin, candL = candL, candU = candU, Lopt = sol$Lopt, Uopt = sol$Uopt)
    return(devolver)
  }
  
}

busquedaMejor <- function(v,candL,candU){
  
  #Vector que contine los datos concatenados dos veces; y el índice correspondiente
  v2<-c(v,v)
  mseFin <- Inf
  
  #cat("Lista de candidatos a mínimo:\n")
  #print(candL)
  #cat("\n")
  #cat("Lista de candidatos a máximo:\n")
  #print(candU)
  #cat("\n")
  
  #Para cada combinación de máximo local y mínimo local se busca el ajuste
  for(indL in candL){
    
    #cat("\nProbando los candidatos a U válidos con indL = ",indL,"\n",sep="")
    
    cond <- (candU >= indL)
    ordenComprobacion <- c(candU[cond],candU[!cond])
    
    nCandUValidos <- 0
    candUValidos <- c()
    pavasCrecientes <- list()
    
    #Determinamos todos los posibles candidatos a U, dado indL
    for(indU in ordenComprobacion){
      
      #cat("Comprobando si indU = ",indU," es válido... \n",sep="")
      
      ######################################
      #Primer caso: L se produce antes que U
      ######################################
      if(indL<indU){
        
        #Índices de las partes donde se va a ejecutar el pava creciente, y cálculo correspondiente
        indexPavaLU<-indL:indU
        if(length(indexPavaLU)>2){
          pavaLU<-pavaC(v[(indL+1):(indU-1)])
          pavaLU<-c(v[indL],pavaLU,v[indU])
        }else{
          pavaLU<-c(v[indL],v[indU])
        }
        
        #Si el ajuste no es el adecuado a la derecha del mínimo local, no hace falta seguir probando
        #con más U's, puesto que el inicio del PAVA no se verá modificado
        if(pavaLU[2]<=v[indL]){
          #cat("\nNo se cumple la condición por la derecha de la parte creciente\n")
          break
        }
        
        #Para saber si indU es un candidato válido, se debe cumplir la condición de que
        #el máximo local es mayor que el inmediato anterior. En este caso, se guarda este
        #candidato y el pava creciente usado
        if(pavaLU[length(pavaLU)-1]<v[indU]){
          
          #cat("(sí lo es)\n")
          nCandUValidos <- nCandUValidos + 1
          candUValidos[nCandUValidos] <- indU
          pavasCrecientes[[nCandUValidos]] <- pavaLU
          
        }
     
      ########################################
      #Segundo caso: L y U son los mismos
      ########################################   
      } else if(indL == indU){
        
        #En este caso el pava siempre será valido
        #cat("(sí lo es)\n")
        nCandUValidos <- nCandUValidos + 1
        candUValidos[nCandUValidos] <- indU
        pavasCrecientes[[nCandUValidos]] <- rep(mean(v),length(v))
        
      ########################################
      #Tercer caso: L se produce después que U
      ########################################  
      }else {
        
        #Cálculo del PAVA creciente
        indexPavaLU<-c(indL:length(v),1:indU)
        if(length(indexPavaLU)>2){
          pavaLU<-pavaC(v2[(indL+1):(length(v)+indU-1)])
          pavaLU<-c(v[indL],pavaLU,v[indU])
        }else{
          pavaLU<-c(v[indL],v[indU])
        }
        
        #Si el ajuste no es el adecuado a la derecha del mínimo local, no hace falta seguir probando
        #con más U's, puesto que el inicio del PAVA no se verá modificado
        if(pavaLU[2]<=v[indL]){
          #cat("\nNo se cumple la condición por la derecha de la parte creciente\n")
          break
        }
        
        #Para saber si indU es un candidato válido, se debe cumplir la condición de que
        #el máximo local es mayor que el inmediato anterior. En este caso, se guarda este
        #candidato y el pava creciente usado
        if(pavaLU[length(pavaLU)-1]<v[indU]){
            
          #cat("(sí lo es)\n")
          nCandUValidos <- nCandUValidos + 1
          candUValidos[nCandUValidos] <- indU
          pavasCrecientes[[nCandUValidos]] <- pavaLU
            
        }
      
      }
      
    }
  

    ##########################################################################################
    #Ahora hay que calcular los pavas decrecientes de los candidatos guardados anteriormente
    ##########################################################################################
    
    #cat("\n\nCálculo de los PAVAS decrecientes\n")
    
    nV <- 0
    for (indU in rev(candUValidos)){
      
      nV <- nV + 1
      
      #cat("Probando con indL = ",indL," y indU = ",indU,"\n",sep="")
      
      ######################################
      #Primer caso: L se produce antes que U
      ######################################
      if(indL<indU){
        
        #Recuperamos el pava correspondiente
        pavaLU <- pavasCrecientes[[nCandUValidos - nV + 1]]
        
        #Ahora hay que calcular el PAVA en la parte decreciente.
          
        #Si el pava creciente no era desde el primer hasta el último punto,
        #hay que construir la parte de pavaUL en dos partes
        if(indU!=length(v) && indL!=1){
          
          #Ajuste del pava decreciente
          indexPavaUL1<-(indU+1):length(v)
          pavaUL<-pavaC(v2[(indU+1):(length(v)+indL-1)],decreasing=TRUE)
          
          #Si no se cumple la condición izquierda a indL, entonces ningún U a la izquierda
          #del actual será válido
          if(pavaUL[length(pavaUL)] < v[indL]){
            #cat("No se cumple la condición por la izquierda de indL\n")
            break
          }
          
          #Si es válido, se calcula el MSE y se toma como mejor solución si es el mínimo hasta ahora
          if(pavaUL[1]<=v[indU]){
            
            #Ajuste completo
            pavaAux<-c(pavaUL[(length(indexPavaUL1)+1):length(pavaUL)],pavaLU,pavaUL[1:length(indexPavaUL1)])
            mseAux<-sum((v-pavaAux)^2)/length(v)
            
            if(mseAux < mseFin){
              mseFin<-mseAux
              pavaFin<-pavaAux
              Lopt<-indL
              Uopt<-indU
            }
          }
          
          # En este caso se considera el PAVA en el caso donde U es el último punto
        } else if(indU==length(v) & indL!=1){
          
          #Calcular el PAVA decreciente
          pavaUL<-pavaC(v[1:(indL-1)],decreasing=TRUE)
          
          #Si no se cumple la condición izquierda a indL, entonces ningún U a la izquierda
          #del actual será válido
          if(pavaUL[length(pavaUL)] < v[indL]){
            #cat("No se cumple la condición por la izquierda de indL\n")
            break
          }
          
          #Si es válido, se calcula el MSE y se toma como mejor solución si es el mínimo hasta ahora
          if(pavaUL[1]<=v[indU]){
            
            pavaAux<-c(pavaUL,pavaLU)
            mseAux<-sum((v-pavaAux)^2)/length(v)
            
            if(mseAux < mseFin){
              mseFin<-mseAux
              pavaFin<-pavaAux
              Lopt<-indL
              Uopt<-indU
            }
          }
          
          #En este caso se considera el PAVA donde L es el primer punto 
        } else if(indL==1 & indU!=length(v)){
          
          pavaUL<-pavaC(v[(indU+1):(length(v))],decreasing=TRUE)
          
          #Si no se cumple la condición izquierda a indL, entonces ningún U a la izquierda
          #del actual será válido
          if(pavaUL[length(pavaUL)] < v[indL]){
            #cat("No se cumple la condición por la izquierda de indL\n")
            break
          }
          
          #Si es válido, se calcula el MSE y se toma como mejor solución si es el mínimo hasta ahora
          if(pavaUL[1]<=v[indU]){
            
            pavaAux<-c(pavaLU,pavaUL)
            mseAux<-sum((v-pavaAux)^2)/length(v)
            
            if(mseAux < mseFin){
              mseFin<-mseAux
              pavaFin<-pavaAux
              Lopt<-indL
              Uopt<-indU
            }
          }
          
          #En este caso se considera el PAVA cuando L es el mínimo y U es el máximo
        } else {
          
          # En este caso no hay que hacer pava decreciente
          
          pavaAux<-pavaLU
          mseAux<-sum((v-pavaAux)^2)/length(v)
          
          if(mseAux < mseFin){
            mseFin<-mseAux
            pavaFin<-pavaAux
            Lopt<-indL
            Uopt<-indU
          }
          
        }  

      ######################################
      #Segundo caso: L es igual a U
      ######################################     
      }else if(indL == indU){
        
        #El estimador es constante
        pavaAux<-pavasCrecientes[[nCandUValidos - nV + 1]]
        mseAux<-sum((v-pavaAux)^2)/length(v)
        
        if( mseAux<mseFin){          
          mseFin<-mseAux
          pavaFin<-pavaAux
          Lopt<-indL
          Uopt<-indU
        }
        
      ########################################
      #Tercer caso: L se produce después que U
      ########################################  
      }else {
        
        #Recuperamos el pava creciente
        pavaLU <- pavasCrecientes[[nCandUValidos - nV + 1]]
        
        
        #Ahora hay que calcular el PAVA en la parte decreciente.
          
        #Si el pava creciente abarca todos los puntos, entonces sólo es necesario el creciente
        if(length(pavaLU)==length(v)){
            
            pavaAux<-c(pavaLU[(length(indL:length(v))+1):length(pavaLU)],pavaLU[1:length(indL:length(v))])
            mseAux<-sum((v-pavaAux)^2)/length(v)
            
            if(mseAux < mseFin){
              mseFin<-mseAux
              pavaFin<-pavaAux
              Lopt<-indL
              Uopt<-indU
            }
          
        #En caso contrario, hay que hacer un pava decreciente
        }else{
          
          #Si el trozo decreciente no empieza en el primer índice, ni termina en el último
          if(indU!=1 & indL!=length(v)){
            
            pavaUL<-pavaC(v[(indU+1):(indL-1)],decreasing=TRUE)
            
            #Si no se cumple la condición izquierda a indL, entonces ningún U a la izquierda
            #del actual será válido
            if(pavaUL[length(pavaUL)] < v[indL]){
              #cat("No se cumple la condición por la izquierda de indL\n")
              break
            }
            
            #Si es válido, se guarda si es el mejor hasta ahora
            if(pavaUL[1]<=v[indU]){
              
              pavaAux<-c(pavaLU[(length(indL:length(v))+1):length(pavaLU)],pavaUL,pavaLU[1:length(indL:length(v))])
              mseAux<-sum((v-pavaAux)^2)/length(v)
              
              if(mseAux < mseFin){
                mseFin<-mseAux
                pavaFin<-pavaAux
                Lopt<-indL
                Uopt<-indU
              }
            }
            
            #Si el trozo decreciente sí empieza en el primer índice, pero no acaba en el último
          } else if(indU==1 & indL!=length(v)){
            
            pavaUL<-pavaC(v[2:(indL-1)],decreasing=TRUE)
            
            #Si no se cumple la condición izquierda a indL, entonces ningún U a la izquierda
            #del actual será válido
            if(pavaUL[length(pavaUL)] < v[indL]){
              #cat("No se cumple la condición por la izquierda de indL\n")
              break
            }
            
            if(pavaUL[1]<=v[indU]){
              
              pavaAux<-c(pavaLU[length(pavaLU)],pavaUL,pavaLU[-length(pavaLU)])
              mseAux<-sum((v-pavaAux)^2)/length(v)
              
              if(mseAux < mseFin){
                mseFin<-mseAux
                pavaFin<-pavaAux
                Lopt<-indL
                Uopt<-indU
              }
            }
            
            #Si el trozo decreciente no empieza en el primer índice pero termina en el último
          } else if(indL==length(v) & indU!=1){
            
            pavaUL<-pavaC(v[(indU+1):(length(v)-1)],decreasing=TRUE)
            
            #Si no se cumple la condición izquierda a indL, entonces ningún U a la izquierda
            #del actual será válido
            if(pavaUL[length(pavaUL)] < v[indL]){
              #cat("No se cumple la condición por la izquierda de indL\n")
              break
            }
            
            if(pavaUL[1]<=v[indU]){
              
              pavaAux<-c(pavaLU[2:length(pavaLU)],pavaUL,pavaLU[1])
              mseAux<-sum((v-pavaAux)^2)/length(v)
              
              if(mseAux<mseFin){
                mseFin<-mseAux
                pavaFin<-pavaAux
                Lopt<-indL
                Uopt<-indU
              }
            }
          }
          
          #Si el trozo decreciente abarca todos los puntos
          if(indL==length(v) & indU==1){
            
            pavaUL<-pavaC(v[(indU+1):(indL-1)],decreasing=TRUE)
            
            #Si no se cumple la condición izquierda a indL, entonces ningún U a la izquierda
            #del actual será válido
            if(pavaUL[length(pavaUL)] < v[indL]){
              #cat("No se cumple la condición por la izquierda de indL\n")
              break
            }
            
            if(pavaUL[1]<=v[indU]){ 
              
              pavaAux<-c(pavaLU[2],pavaUL,pavaLU[1])
              mseAux<-sum((v-pavaAux)^2)/length(v)
              
              if(mseAux < mseFin){
                mseFin<-mseAux
                pavaFin<-pavaAux
                Lopt<-indL
                Uopt<-indU
              }
            }
          }
        }
      }
    }
    
    if(mseFin == 0){
      break
    }
    
  }
  
  if(is.infinite(mseFin)){
    return (NULL)
  } else {
    return(list(pava=pavaFin,mse=mseFin,candL=candL,candU=candU,Lopt=Lopt,Uopt=Uopt))
  }
  
}
