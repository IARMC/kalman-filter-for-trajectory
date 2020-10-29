#####     Configuracion Inicial     #####
rm(list = ls())
options( digits = 15 )

#####     Conexion a la base de datos     #####
library(RPostgreSQL)
dbdriver <- "PostgreSQL";
host <- 'localhost';
port <- '5432';
dbname <- 'CPPP';
user <- 'postgres';
pass <- 'toor';
drv <- dbDriver(dbdriver)
con <- dbConnect(drv, host = host, port = port, dbname = dbname, user = user, pass = pass)

#####     Obtencion de la informacion desde la base de datos     #####
dataset <- "brasil"
data <- dbGetQuery(con, "SELECT longitude, latitude, unixtime, file_name FROM public.brasil WHERE  longitude BETWEEN  '-37.095' and '-37.04' and latitude BETWEEN  '-10.99' and '-10.89' and unixtime > '13/09/14' and unixtime <= '29/07/15' ")
dbDisconnect(con)

#####     Preprocesamiento     #####
sapply(data, function(x) sum(is.na(x)))
data <- data[!is.na(data$longitude),]
data <- data[!is.na(data$latitude),]
data <- data[!is.na(data$file_name),]
data <- data[!is.na(data$unixtime),]

trajectory_after_kalman <- list()
binnacle_for_latitude <- data.frame(Medicion = 0,MedicionK = 0,Ruido_Proc = 0,Ruido_Med = 0,GainK = 0,Cov = 0,trayectoria = 0)
binnacle_for_longitude <- data.frame(Medicion = 0,MedicionK = 0,Ruido_Proc = 0,Ruido_Med = 0,GainK = 0,Cov = 0,trayectoria = 0)

#####     Division de los datos por trayectorias     #####
array_name <- data.frame(file_name=unique(data$file_name)) 

todas <- data.frame(longitude = c(1:nrow(data)), latitude = 0, file_name = 0, speed = 0, unixtime = 0)
todas$longitude <- data$longitude
todas$latitude  <- data$latitude
todas$file_name <- data$file_name
todas$unixtime  <- data$unixtime
todas$speed     <- data$speed

lista_trayectorias_sinprocesar <- list ()

for(count in 1:length(unique(todas$file_name))){
  tray_corresp <- subset(todas, todas$file_name == array_name$file_name[count])
  lista_trayectorias_sinprocesar[[count]] <-  tray_corresp
}

#####     Funciones de tratameinto de las bitacoras     #####
add_binnacle_longitude <- function(newBinnacle){
  binnacle_for_longitude <<-  rbind(binnacle_for_longitude, newBinnacle)
}

add_binnacle_latitude <- function(newBinnacle){
  binnacle_for_latitude <<-  rbind(binnacle_for_latitude, newBinnacle)
}

#####     Funciones Kalman     #####
get_dimentional_kalman <- function(dimention){
  #input: class dimention => vector
  tamanio_dimention <- length(dimention)
  
  binnacle <- data.frame(Medicion = c(1:tamanio_dimention), 
                         MedicionK = 0, 
                         Ruido_Proc = 0, 
                         Ruido_Med = 0, 
                         GainK = 0, 
                         Cov = 0, 
                         trayectoria = 0)
  
  Q <- 0  #Process noise covariance - Ruido del proceso     X
  R <- 0  #Measurement noise covariance - Ruido del sensor  X
  X <- 0  #Value - Valor filtrado
  P <- 0  #Estimation error covariance - Error estimado
  K <- 0  #Kalman gain - Ganancia de Kalman [0~1]
  M <- 0  #Measurement MEDICION
  
  #Valores iniciales
  Q <- 0.5 #0.125   0.33   0.5
  R <- 0.125
  P <- 1 #No importante, se ajusta durante el proceso
  
  #Calibrar valores
  for (i in 1:tamanio_dimention) {
    Medicion <- dimention[i]
    M <- Medicion
    P <- P + Q
    K <- P / (P + R)
    X <- X + K * (M - X)
    P <- (1 - K) * P
  }
  
  for (i in 1:tamanio_dimention) {
    Medicion <- dimention[i]
    binnacle[i, 'Medicion'] <- Medicion
    M <- Medicion
    
    binnacle[i, 'Ruido_Proc'] <- Q
    binnacle[i, 'Ruido_Med'] <- R
    
    #Prediccion
    P <- P + Q
    
    #Correccion
    K <- P / (P + R)  #[0~1]
    binnacle[i, 'GainK'] <- K
    
    X <- X + K * (M - X)
    if(i>2){
      binnacle[i, 'MedicionK'] <- X
    }else{
      binnacle[i, 'MedicionK'] <- Medicion
    }
    
    P <- (1 - K) * P
    binnacle[i, 'Cov'] <- P
  }
  
  return(binnacle)
}

create_dataframe_kalman <- function(longitude, latitude, file_name, unixtime){
  kalman_dataframe <- data.frame(longitude, latitude, file_name, unixtime)
  return(kalman_dataframe)
}

apply_kalman <- function(trayectoria) {
  #Aplicar kalman a la latitud
  binnacle_latitude <- get_dimentional_kalman(trayectoria$latitude)
  binnacle_latitude$trayectoria <- trayectoria$file_name
  add_binnacle_latitude(binnacle_latitude)
  
  #Aplicar kalman a la longitud
  binnacle_longitude <- get_dimentional_kalman(trayectoria$longitude)
  binnacle_longitude$trayectoria <- trayectoria$file_name
  add_binnacle_longitude(binnacle_longitude)
  
  new_trajectory <- create_dataframe_kalman(longitude=binnacle_longitude$MedicionK, 
                                            latitude=binnacle_latitude$MedicionK, 
                                            file_name=trayectoria$file_name, 
                                            unixtime=trayectoria$unixtime)
  
  return(new_trajectory)
}

#####     Funciones de almacenamiento en disco duro de los resultados     #####
#Establecer directorio
establecer_directorio_kalman <- function(dataset){
  simbol <- "\\"
  root <- "C:"
  setwd(paste("C:",simbol,sep=""))
  
  folder <- "Filtro-Kalman"
  dir.create(folder);
  setwd(paste(getwd(),simbol,folder,sep=""))
  
  algoritmo <- paste("Kalman_", Sys.Date(),sep = "")
  dir.create(algoritmo);
  setwd(paste(getwd(),simbol,algoritmo,sep=""))
  
  dir.create(dataset)
  setwd(paste(getwd(),simbol,dataset,sep=""))
  
  return(getwd())
}

#Almacenamiento de las bitacoras


#Almacenamiento de los resultados


#Almacenameinto de los graficos
guardar_trayectoria_individual <- function(idtr, tr_original, tr_kalman){
  nombre <- paste("trayectoria_",idtr,".png",sep = "")
  png(filename=nombre, width = 1024, height = 720, units = 'px', pointsize = 20 ,res = NA)
  plot(x=tr_original$longitude, y=tr_original$latitude,type="l",col="red", xlab="Longitud", ylab="Latitud", main=paste("Results Kalman Filter_TR",idtr,sep = ""))
  points(x=tr_kalman$longitude, y=tr_kalman$latitude,type="l",col="blue")
  #https://www.rdocumentation.org/packages/graphics/versions/3.6.2/topics/legend
  legend("bottomleft",legend=c("Original","Filtro de Kalman"),col=c("red","blue"), pch=1,bty="n",ncol=1,cex=1,pt.cex=1)
  #Aplicar zoom al grafico
  #zoom::zm()
  dev.off()
}

guardar_graficos <- function(original, kalman){
  if(length(original)==length(kalman)){
    print("Guardando graficos...")
    for (k in 1:length(original)) {
      guardar_trayectoria_individual(idtr=k, tr_original = original[[k]], tr_kalman = kalman[[k]])
    }
    print("Graficos guardados.")
  }else{
    print("Error en la cantidad de trayectorias")
  }
}

#####     Recorrido de trayectorias para aplicar Kalman     #####
cant_tr <- length(unique(lista_trayectorias_sinprocesar))

for (count in 1:cant_tr){
  print(paste("Aplicando Filtro de Kalman. Trayectoria ",count,"...", sep = ""))
  
  trayectoria_individual <- lista_trayectorias_sinprocesar[[count]]
  cant_points <- nrow(trayectoria_individual)
  
  if(cant_points <= 2){
    tr_kalman <- trayectoria_individual
  }else{
    tr_kalman <- apply_kalman(trayectoria_individual)
  }
  
  trajectory_after_kalman[[count]] <- tr_kalman
}

#####     Almacenar Datos      #####
establecer_directorio_kalman(dataset)

guardar_graficos(lista_trayectorias_sinprocesar, trajectory_after_kalman)

binnacle_prueba()
#TODO: UNIR BITACORAS
#TODO: GUARDAR BITACORAS


binnacle_for_latitude

