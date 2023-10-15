
#Carga din?mica de c?digo en C
dyn.load("pava.dll")
dyn.load("busquedaMejor.dll")

#Funci?n para el c?lculo del mejor estimador Up-Down-Up
function1Local_sec<-function(v){
  
  #Buscamos los m?ximos y m?nimos locales que existen
  extremos <- extremosLocales(v)
  candL <- extremos$candL
  candU <- extremos$candU
  
  mejor <- busquedaMejorC(v,candL,candU)
  
  return(mejor)
}
