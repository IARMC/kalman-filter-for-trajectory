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
vk <- c()
vlat <- c()
vlon <- c()

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

#####     Funciones Kalman     #####
prueba_funcion <- function(pkg) { #TODO:QUIT
  #ambito local de la funcion
  valor <- 5
  print(valor)
}

get_dimentional_kalman <- function(dimention){
  #input: class dimention => vector
  #dimention <- trayectoria$latitude #TODO: QUIT
  tamanio_dimention <- length(dimention)
  
  binnacle <- data.frame(Medicion = c(1:tamanio_dimention), 
                         MedicionK = 0, 
                         Ruido_Proc = 0, 
                         Ruido_Med = 0, 
                         GainK = 0, 
                         Cov = 0, 
                         trayectoria = 0)
  
  results_from_kalman <- c()
  
  Q <- 0  #Process noise covariance - Ruido del proceso
  R <- 0  #Measurement noise covariance - Ruido del sensor
  X <- 0  #Value - Valor filtrado
  P <- 0  #Estimation error covariance - Error estimado
  K <- 0  #Kalman gain - Ganancia de Kalman [0~1]
  M <- 0  #Measurement MEDICION
  
  #valores iniciales
  Q <- 0.5 #0.125   0.33   0.5
  R <- 0.125
  P <- 1 #No importante, se ajusta durante el proceso
  results_from_kalman <- c()
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
    #print(K)
    X <- X + K * (M - X)
    binnacle[i, 'MedicionK'] <- X
    #print(X)
    P <- (1 - K) * P
    binnacle[i, 'Cov'] <- P
    #print(P)
    results_from_kalman <- c(results_from_kalman, X)
  }
  
  #TODO: MED
  
  binnacle$MedicionK <- results_from_kalman
  #par(lwd=0.5)
  plot(binnacle[,'Medicion'],type="o",col="red", ylab="",main="Latitude Output")
  
  return(binnacle)
}

create_dataframe_kalman <- function(longitude, latitude, file_name, unixtime){
  kalman_dataframe <- data.frame(longitude, latitude, filename, unixtime)
  return(kalman_dataframe)
}

apply_kalman <- function(trayectoria) {
  print(trayectoria) #TODO: QUIT
  
  trayectoria <- trayectoria_individual #TODO: QUIT
  
  binnacle_latitude <- get_dimentional_kalman(trayectoria$latitude)
  binnacle_latitude$trayectoria <- trayectoria$file_name
  binnacle_latitude$MedicionK
  
  binnacle_longitude <- get_dimentional_kalman(trayectoria$longitude)
  binnacle_longitude$trayectoria <- trayectoria$file_name
  
  #TODO: UNIR BITACORAS
  #TODO: CREAR GRAFICOS POR TRAYECTORIA
  
  new_trajectory <- create_dataframe_kalman(binnacle_longitude$MedicionK, binnacle_latitude$MedicionK, trayectoria$file_name, trayectoria$unixtime)
  
  #return cbind entre kalman_longitude, kalman_latitude, file_name, unixtime
  return(new_trajectory)
}


#####     Recorrido de trayectorias     #####
cant_tr <- length(unique(lista_trayectorias_sinprocesar))
count <- 1 #TODO: QUIT
for (count in 1:cant_tr){
  trayectoria_individual <- lista_trayectorias_sinprocesar[[count]]
  cant_points <- nrow(trayectoria_individual)
  
  if(cant_points <= 2){
    trajectory_after_kalman[[count]] <- trayectoria_individual
    
  }else{
    
    tr_kalman <- apply_kalman(trayectoria_individual)
    trajectory_after_kalman[[count]] <- tr_kalman
    
  }
}
#TODO: GUARDAR BITACORAS




