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
print("Conexion establecida")

#####     Obtencion de la informacion desde la base de datos     #####
dataset <- "brasil"
print("Obtencion de los datos")
data <- dbGetQuery(con, "SELECT longitude, latitude, unixtime, file_name FROM public.brasil WHERE  longitude BETWEEN  '-37.095' and '-37.04' and latitude BETWEEN  '-10.99' and '-10.89' and unixtime > '13/09/14' and unixtime <= '29/07/15' ")
dbDisconnect(con)

#####     Preprocesamiento     #####
sapply(data, function(x) sum(is.na(x)))
data <- data[!is.na(data$longitude),]
data <- data[!is.na(data$latitude),]
data <- data[!is.na(data$file_name),]
data <- data[!is.na(data$unixtime),]
print("Dataset Limpiado")

trajectory_after_kalman <- list()
binnacle_for_latitude <- data.frame(Medicion = 0,MedicionK = 0,Ruido_Proc = 0,Ruido_Med = 0,GainK = 0,Cov = 0,trayectoria = 0)
binnacle_for_longitude <- data.frame(Medicion = 0,MedicionK = 0,Ruido_Proc = 0,Ruido_Med = 0,GainK = 0,Cov = 0,trayectoria = 0)
print("Varaibles de almacenameinto declaradas")

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
print("Division por trayectorias finalizada")

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
    binnacle[i, 'GainK'] <- K
    X <- X + K * (M - X)
    P <- (1 - K) * P
    binnacle[i, 'Cov'] <- P
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
    
    X <- X + K * (M - X)
    if(i>2){
      binnacle[i, 'MedicionK'] <- X
    }else{
      binnacle[i, 'MedicionK'] <- Medicion
    }
    
    P <- (1 - K) * P
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
  
  print("Directorio Establecido")
  return(getwd())
}

#Almacenamiento de las bitacoras
save_binnacles <- function(){
  write.table(binnacle_for_latitude, file ="binnacle-latitude.csv" , sep = ";", row.names = FALSE, col.names = TRUE)
  write.table(binnacle_for_longitude, file ="binnacle-longitude.csv" , sep = ";", row.names = FALSE, col.names = TRUE)
  print("Bitacoras guardadas")
}

#Almacenamiento de los resultados
save_trajectories <- function(list){
  write.table(do.call(rbind,list), file ="trajectories-after-kalman.csv" , sep = ";", row.names = FALSE, col.names = TRUE)
  print("Trayectorias guardadas")
}

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
  
  print(paste("Se han guardado ", length(kalman), " graficos", sep=""))
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

print(paste("Finalizado. Filtro de Kalman aplicado a ",cant_tr,"Trayectorias.", sep = ""))

#####     Almacenar Datos      #####
establecer_directorio_kalman(dataset)

guardar_graficos(lista_trayectorias_sinprocesar, trajectory_after_kalman)

save_binnacles()

save_trajectories(trajectory_after_kalman)



#####     Comprobacion con libreria     #####
library(dlm)

kalman_libreria <- list()

nileBuild <- function(par) {
  dlmModPoly(order = 1, dV = exp(par[1]), dW = exp(par[2]))
  #dlmModPoly(order = 1, dV = 0.5, dW = 0.125)
}

cant_tr <- length(unique(lista_trayectorias_sinprocesar))

for (count in 1:cant_tr){
  count <- 1
  print(paste("Aplicando DLM Kalman. Trayectoria ",count,"...", sep = ""))
  tr_i <- lista_trayectorias_sinprocesar[[count]]
  cant_points <- nrow(tr_i)
  frame_tray <- data.frame(x = lista_trayectorias_sinprocesar[[count]]$longitude, 
                           y = lista_trayectorias_sinprocesar[[count]]$latitude, 
                           t = lista_trayectorias_sinprocesar[[count]]$unixtime,
                           g = 0)
  if(cant_points <= 6){
    new_frame_kalman <- frame_tray
    print("Menos de 6")
  }else{
    matrix1 <- frame_tray
    
    nileMLE1 <- dlmMLE(matrix1[,1],parm = c(1, 1),nileBuild)
    nileMLE2 <- dlmMLE(matrix1[,2],parm = c(1, 1),nileBuild)
    
    nileMod1 <- nileBuild(nileMLE1$par)
    nileMod2 <- nileBuild(nileMLE2$par)
    
    nileFilt1 <- dlmFilter(matrix1[,1], nileMod1)
    nileFilt2 <- dlmFilter(matrix1[,2], nileMod2)
    
    nileSmooth1 <- dlmSmooth(nileFilt1)
    nileSmooth2 <- dlmSmooth(nileFilt2)
    
    tamano_kalman <- length(nileSmooth1$s)
    new_frame_kalman <- data.frame(x = nileSmooth1$s[2:tamano_kalman], 
                                   y = nileSmooth2$s[2:(tamano_kalman)], 
                                   t = frame_tray[,3], 
                                   g = frame_tray[,4])
  }
  
  kalman_libreria[[count]] <- new_frame_kalman
}

#####     Comparativa Original - Kalman - DLM     #####
sel_tr <- 1
#kalman_libreria[[5]]
nombre <- paste("Comparativa trayectoria ",sel_tr,".png",sep = "")
#png(filename = nombre, width = 1024, height = 720, units = 'px', pointsize = 20 ,res = NA)
plot(x = lista_trayectorias_sinprocesar[[sel_tr]]$longitude, 
     y = lista_trayectorias_sinprocesar[[sel_tr]]$latitude,
     type = "l", col = "red", xlab = "Longitud", ylab = "Latitud", 
     main = paste("Results Kalman Filter_TR",sel_tr,sep = ""))
points(x = trajectory_after_kalman[[sel_tr]]$longitude, 
       y = trajectory_after_kalman[[sel_tr]]$latitude,
       type = "l",col = "blue")
points(x = kalman_libreria[[sel_tr]]$x, 
       y = kalman_libreria[[sel_tr]]$y, 
       type = "l", col = "green")
#https://www.rdocumentation.org/packages/graphics/versions/3.6.2/topics/legend
legend("bottomleft",legend = c("Original","Filtro de Kalman","Libreria"), 
       col = c("red","blue","green"), pch = 1, bty = "n", 
       ncol = 1, cex = 1, pt.cex = 1)
#Aplicar zoom al grafico
#zoom::zm()
#dev.off()
#dev.list()



#####     Graficacion en Mapa     #####
#packages <- c("RPostgreSQL","rlang", "kmlShape","brotli","dlm", "ggplot2", "caret", "class", "mapview",  "compare", "pracma" , "stringr","SpatialTools","matlib", "dplyr","chron","lubridate","zoom","RgoogleMaps","ggmap") #librerias 
#packages <- c("RPostgreSQL","rlang", "ggplot2", "caret", "class", "mapview",  "compare", "pracma", "leaflet" , "stringr","SpatialTools","matlib", "dplyr","chron","lubridate","zoom","RgoogleMaps","ggmap") #librerias
#install.packages("leaflet")
#install.packages("mapview")

packages <- c("leaflet","mapview","webshot")
ipak <- function(pkg) { #cargar e instalar paquetes de manera dinamica
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
ipak(packages)

n_tr <- 42
m <- leaflet()
m <- addTiles(m)
m <- addCircles(map = m, lng = lista_trayectorias_sinprocesar[[n_tr]]$longitude, lat =  lista_trayectorias_sinprocesar[[n_tr]]$latitude, radius = 0.05, weight = 5, color = "red")
m <- addCircles(map = m, lng = trajectory_after_kalman[[n_tr]]$longitude, lat =  trajectory_after_kalman[[n_tr]]$latitude, radius = 0.05, weight = 5, color = "blue")
#m <- addCircles(map = m, lng = trajectory_after_kalman[[n_tr]]$longitude, lat =  trajectory_after_kalman[[n_tr]]$latitude, radius = 0.05, weight = 5, color = "yellow")
show(m)


install_phantomjs(version = "2.1.1",baseURL = "https://github.com/wch/webshot/releases/download/v0.3.1/")
#webshot::install_phantomjs()
print(getwd())
#mapshot(m, url = "map.html", file = "Mapa localidad y trayectorias.png")
mapshot(m, url = paste0(getwd(), "/map.html"), file = paste0(getwd(), "/map.png"))