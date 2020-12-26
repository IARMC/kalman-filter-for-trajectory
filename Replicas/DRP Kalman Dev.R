#####     Configuracion Inicial     #####
rm(list = ls())
options( digits = 15 )

#####     Funciones Generales     #####
#Importacion de librerias
ipak <- function(pkg) { #cargar e instalar paquetes de manera dinamica
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

#####     Funciones de Tratamiento de Bitacoras del Filtro de Kalman     #####
add_binnacle_longitude <- function(newBinnacle){
  binnacle_for_longitude <<-  rbind(binnacle_for_longitude, newBinnacle)
}

add_binnacle_latitude <- function(newBinnacle){
  binnacle_for_latitude <<-  rbind(binnacle_for_latitude, newBinnacle)
}

#####     Funciones de Tratamiento del Filtro de Kalman     #####
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

#####     Funciones de Almacenamiento de Resultados del Filtro de Kalman     #####
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
  nombre <- paste("trayectoria_", idtr, ".png", sep = "")
  png(filename = nombre, width = 1024, height = 720, units = 'px', pointsize = 20 ,res = NA)
  plot(x = tr_original$longitude, 
       y = tr_original$latitude,
       type = "l",col = "red", 
       xlab = "Longitud", ylab = "Latitud", 
       main = paste("Results Kalman Filter_TR",idtr,sep = ""))
  points(x = tr_kalman$longitude, 
         y = tr_kalman$latitude, 
         type = "l",col = "blue")
  #https://www.rdocumentation.org/packages/graphics/versions/3.6.2/topics/legend
  legend("bottomleft", legend = c("Original","Filtro de Kalman"), 
         col = c("red","blue"), pch = 1, bty = "n", 
         ncol = 1, cex = 1, pt.cex = 1)
  #Aplicar zoom al grafico
  #zoom::zm()
  dev.off()
}

guardar_graficos <- function(original, kalman){
  if(length(original)==length(kalman)){
    print("Guardando graficos...")
    for (k in 1:length(original)) {
      guardar_trayectoria_individual(idtr = k, tr_original = original[[k]], tr_kalman = kalman[[k]])
    }
    print("Graficos guardados.")
  }else{
    print("Error en la cantidad de trayectorias")
  }
  
  print(paste("Se han guardado ", length(kalman), " graficos", sep=""))
}

#Guardado de las trayectorias en los mapas
guardar_mapa_individual <- function(idtr, tr_original, tr_kalman){
  nombre <- paste("Mapas\\","Trayectoria_", idtr, sep = "")
  mapa <- leaflet()
  mapa <- addTiles(mapa)
  mapa <- addCircles(map = mapa, lng = tr_original[idtr]$longitude, lat = tr_original[idtr]$latitude, radius = 0.05, weight = 5, color = "red")
  mapa <- addCircles(map = mapa, lng = tr_kalman[idtr]$longitude, lat = tr_kalman[idtr]$latitude, radius = 0.05, weight = 5, color = "blue")
  #print(getwd())
  mapshot(mapa, url = paste0(getwd(), nombre, ".html"), file = paste0(getwd(), nombre, ".png"))
}

guardar_en_mapa <- function(original, kalman){
  if(length(original)==length(kalman)){
    print("Guardando mapas...")
    for (k in 1:length(original)) {
      guardar_mapa_individual(idtr = k, tr_original = original[[k]], tr_kalman = kalman[[k]])
    }
    print("Mapas guardados.")
  }else{
    print("Error en la cantidad de trayectorias.")
  }
  
  print(paste("Se han guardado ", length(kalman), " mapas", sep=""))
}

#####     Funciones de Tratamiento de Algoritmos de Simplificacion     #####
crear_directorio<- function(folders, algorithm, dataset){
  
  simbol<-"\\";
  root<- "D:";  	##directorio raiz
  setwd(paste("D:",simbol,sep="")); ##C:
  
  dir.create(folders);
  setwd(paste(getwd(),simbol,folders,sep=""));#C:folders
  
  dir.create(algorithm);
  setwd(paste(getwd(),simbol,algorithm,sep=""));#C:folders/algorithm
  
  dir.create(dataset);							## modificacion del directorio por base a experimentar
  setwd(paste(getwd(),simbol,dataset,sep="")); #C:folders/algorithm/base
  
  return(getwd())
}

loc_dist <-function(end, start) {
  # """ Spatial distance between two points (end-start)
  
  #   Args:
  #       start (:obj:`Point`)
  #       end (:obj:`Point`)
  #   Returns:
  #       float, distance in m
  #  """
  return (distance_f(end,start))
}

time.difference<- function(self, previous){
  #""" Calcultes the time difference against another point
  
  #     Args:
  #          previous (:obj:`Point`): Point before
  #      Returns:
  #          Time difference in seconds
  #      """
  return(abs(self[,c('unixtime')]-previous[,c('unixtime')]) )
}

time_dist<- function(end, start){
  #""" Temporal distance between two points (end-start)
  
  #  Args:
  #      start (:obj:`Point`)
  #      end (:obj:`Point`)
  # Returns:
  #      float, time difference in seconds
  #  """
  return(time.difference(end,start))
}

distance_f<-function(p_a, p_b){
  #  p_a=points[i,]
  #  p_b= point
  #""" Euclidean distance, between two points
  
  #  Args:
  #      p_a (:obj:`Point`)
  #      p_b (:obj:`Point`)
  #  Returns:
  #      float: distance, in degrees
  #  """
  return (sqrt((p_a$latitude-p_b$latitude)^2 + (p_a$longitude-p_b$longitude)^2) )
}

td_tr<-function(points, dist_threshold){
  #points = lista_trayectorias[[j]]; dist_threshold = dist_threshold[[j]];
  #""" Top-Down Time-Ratio Trajectory Compression Algorithm
  #Detailed in https://www.itc.nl/library/Papers_2003/peer_ref_conf/meratnia_new.pdf
  #Args:
  #   points (:obj:`list` of :obj:`Point`): trajectory or part of it
  #   dist_threshold (float): max distance error, in meters
  #Returns:
  #  :obj:`list` of :obj:`Point`, compressed trajectory
  #"""
  #http://tinyurl.com/rbd7n2r
  if (nrow(points)<=2) {
    return(points)
  }else {
    max_dist_threshold = 0
    found_index = 0
    delta_e = as.numeric(time_dist(points[nrow(points),],points[1,]), units = "secs") * I_3600 #segundos
    d_lat = points[nrow(points),][,'latitude'] - points[1,][,'latitude']
    d_lon = points[nrow(points),][,'longitude'] - points[1,][,'longitude']
    
    for (i in 2:(nrow(points)-1)) {
      delta_i= as.numeric(time_dist(points[i,],points[1,]), units = "secs")  * I_3600
      di_de = delta_i / delta_e
      if(is.nan(di_de) == FALSE){ 
        point = data.frame(longitude = points[1,]$longitude + d_lon * di_de,
                           latitude = points[1,]$latitude + d_lat * di_de,
                           unixtime=0)
        dist = loc_dist(points[i,], point)#loc_dist(points[1,], points[nrow(points),])
        if (dist > max_dist_threshold) {
          max_dist_threshold = dist
          found_index = i
        }
      }
    }
  }  
  if ( max_dist_threshold > dist_threshold ) {
    one = td_tr(points =  points[1:found_index,], dist_threshold = dist_threshold)
    two = td_tr(points =  points[found_index:nrow(points),], dist_threshold = dist_threshold)
    two=two[2:nrow(two),]
    one = rbind(one,two)
    return(one)
  }else{
    return(rbind(points[1,],points[nrow(points),]))
  }
}

point_line_distance<- function(point, start, end){
  # """ Distance from a point to a line, formed by two points
  # 
  #   Args:
  #       point (:obj:`Point`)
  #       start (:obj:`Point`): line point
  #       end (:obj:`Point`): line point
  #   Returns:
  #       float: distance to line, in degrees
  #   """
  if (isTRUE(compare(start,end))== TRUE ) {
    return(distance_f(point,start))
  }else{
    un_dist<-abs(((end$latitude-start$latitude) * (start$longitude-point$longitude))
                 - ((start$latitude-point$latitude) * (end$longitude-start$longitude)))
    n_dist<- sqrt((end$latitude-start$latitude)^2 + (end$longitude-start$longitude)^2)
    if (n_dist==0) {
      return(0)
    }else{
      return(un_dist/n_dist)
    }
  }
}

drp<-function(points, dist_threshold){
  #""" Douglas ramer peucker
  #
  #   Based on https://en.wikipedia.org/wiki/Ramer%E2%80%93Douglas%E2%80%93Peucker_algorithm
  #
  #   Args:
  #      points (:obj:`list` of :obj:`Point`)
  #     dist_threshold (float): drp threshold
  # Returns:
  #    :obj:`list` of :obj:`Point`
  #"""
  #http://tinyurl.com/rbd7n2r
  dmax = 0
  index = 0
  if (nrow(points) >1) {
    for (i in 2:nrow(points)-1) {
      dist<-point_line_distance(points[i,],points[1,],points[nrow(points),])
      
      if (dist > dmax) {
        index = i
        dmax = dist
      }
    }
    if (dmax > dist_threshold) {
      one=drp(points[1:index+1,],dist_threshold)[-nrow(points),]
      two=drp(points[index:nrow(points),],dist_threshold)
      
      return(rbind(one,two))
      #return(rbind(drp(points[1:index+1,],epsilon)[-nrow(points),],
      #       drp(points[index:nrow(points),],epsilon)))
    }else{
      return(rbind(points[1,],points[nrow(points),]))
    }
  }else{return(points)}
}

graficar_inicial <- function(lista_trayectoria, trayectorias_informacion){
  par(lwd=2) 
  plot(0 ,ylim=c(min_latitud , max_latitud) , xlim= c(men_longitud , max_longitud ), xlab = 'Longitud', ylab = 'Latitud',main=dataset)
  for(a in 1:length(lista_trayectoria)){
    color <- trayectorias_informacion[a,'Centroide']  			## color para pintar, cuando tengan un mismo centroide quedan del mismo color
    if(color == 0 || color == 'None'){ 									## al incializar pinta trayectorias de todos los colores
      color <- 'black'
    }
    tray_aux <- lista_trayectoria[[a]]
    points(tray_aux[,1:2], pch=46, col=color, cex=2.5)# color 64
  }
}

preprocesar_informacion <-  function (lista_trayectorias){
  cant_trayectorias <- length(lista_trayectorias)
  trayectorias_informacion<-data.frame( No_Trayectoria = c(1:cant_trayectorias),Cant_puntos_originales = 0,Cant_puntos_simplificados = 0, Centroide = 0, silhouette = 0,Velocidad_Prom = 0, Men_Latitud =0 ,Long_Corresp1 = 0 ,May_Latitud  =0 ,Long_Corresp2=0 ,Men_Longitud=0, Lat_Corresp1 = 0, May_Longitud=0 , Lat_Corresp2 = 0 , Dist_ancho=0, Dist_alto = 0,tiempo_recorrido=0)
  
  for(i in 1:cant_trayectorias){
    tray_actual <- lista_trayectorias[[i]]
    cant_puntos <- nrow(tray_actual)
    
    ##ordeno de menor a mayor   y 
    tray_order1 <- tray_actual[order(tray_actual$latitude, decreasing=FALSE),]
    trayectorias_informacion$Men_Latitud[i] <- tray_order1$latitude[1]
    trayectorias_informacion$Long_Corresp1[i] <- tray_order1$longitude[1]
    ## la mayor es la ultima, pues esta ordenado de menor a mayor  
    trayectorias_informacion$May_Latitud[i] <- tray_order1$latitude[cant_puntos]
    trayectorias_informacion$Long_Corresp2[i] <- tray_order1$longitude[cant_puntos]
    
    ##distancia alto
    trayectorias_informacion$Dist_alto[i] <- dist(rbind(tray_order1[1,c('longitude','latitude')], tray_order1[cant_puntos,c('longitude','latitude')]), method = "euclidean")
    
    ##ordeno de menor a mayor la longitud  x
    tray_order2 <- tray_actual[order(tray_actual$longitude,decreasing=FALSE),]
    trayectorias_informacion$Men_Longitud[i] <- tray_order2$longitude[1]
    trayectorias_informacion$Lat_Corresp1[i] <- tray_order2$latitude[1]
    ##ordeno de mayor a menor la longitud  x
    trayectorias_informacion$May_Longitud[i] <- tray_order2$longitude[cant_puntos]
    trayectorias_informacion$Lat_Corresp2[i] <- tray_order2$latitude[cant_puntos]
    
    ##distancia ancho
    trayectorias_informacion$Dist_ancho[i] <- dist(rbind(tray_order2[1,c('longitude','latitude')], tray_order2[cant_puntos,c('longitude','latitude')]), method = "euclidean")
    ##almaceno la cant de puntos de la tray
    trayectorias_informacion$Cant_puntos_originales[i] <- cant_puntos
    
    #trayectorias_informacion$Velocidad_Prom[i] <- mean(tray_actual$unixtime)
    
    #trayectorias_informacion$tiempo_recorrido[i] <- as.numeric(time_dist(tray_actual[1,],tray_actual[cant_puntos,]), units = "secs") 
  }
  return (trayectorias_informacion)
}

calcular_metricas <- function(){
  datos <-data.frame ( Grupos = c(1:(length(unique(trayectorias_informacion$Centroide))-1)), Cantidad_Trayectorias = 0, silhouette = 0) 
  for(i in 1:(length(unique(trayectorias_informacion$Centroide))-1)){
    grupo_informacion <- subset(trayectorias_informacion , Centroide  == i )
    datos[i,'Cantidad_Trayectorias'] <- nrow(grupo_informacion)
    datos[i,'silhouette'] <- mean(grupo_informacion$silhouette)
  }
  print("Datos de los grupos formados:")
  print(datos)
  return(datos)
}

calcular_angulo <-  function (coord1 , coord2){
  ## angulo en grados
  anguloGrados <- (atan2(coord2[1,2] - coord1[1,2] , coord2[1,1] -coord1[1,1]) * 180) / pi
  return (anguloGrados)
}

#####     MAIN - Filtro de Kalman     #####
packages <- c("RPostgreSQL", "rlang", "ggplot2", "caret", "class", 
              "mapview",  "compare", "pracma" , "stringr", "SpatialTools",
              "matlib", "dplyr", "chron", "lubridate", "zoom",
              "RgoogleMaps", "ggmap","leaflet","webshot")
ipak(packages)

#####     Datos de la Conexion a la base de datos     #####
dbdriver <- "PostgreSQL";
host <- 'localhost';
port <- '5432';
dbname <- 'CPPP2';
user <- 'postgres';
pass <- 'toor';
drv <- dbDriver(dbdriver)
con <- dbConnect(drv, host = host, port = port, dbname = dbname, user = user, pass = pass)
print("Conexion establecida")

#####     Obtencion de la informacion desde la base de datos     #####
print("Obtencion de los datos")
#DataSets Disponibles:
# 1. Brasil     14096
# 2. Beijing    62138
# 3. Guayaquil  ~DATASET: TDTR
# 4. Guayaquil  1460
# 5. Quito      1460
# 6. California 914684
# 7. San Francisco
# 8. Circular total
# 9. Media circular
# 10. Circular lineal
# 11. Mitaad circular
selected_dataset <- 4

if (selected_dataset == 1) { 
  dataset <- "Brasil"
  time_in_timestamp = 'FALSE'
  data <- dbGetQuery(con, "SELECT longitud as longitude, latitud as latitude, fecha as unixtime, id as file_name FROM brasil")
  
}else if (selected_dataset == 2) {
  dataset <- "Beijing"
  time_in_timestamp = 'FALSE'
  data <- dbGetQuery(con, "SELECT longitud as longitude, latitud as latitude, fecha as unixtime, id as file_name FROM beijing")
  
}else if (selected_dataset == 3) {
  dataset <- "Guayaquil"
  time_in_timestamp = 'FALSE'
  data <- dbGetQuery(con, "SELECT longitud as longitude, latitud as latitude, fecha as unixtime, id as file_name FROM guayaquil")
  
}else if (selected_dataset == 4) {
  dataset <- "Guayaquil"
  time_in_timestamp = 'FALSE'
  data <- dbGetQuery(con, "SELECT longitud as longitude, latitud as latitude, fecha as unixtime, id as file_name FROM guayaquilmin")
  
}else if (selected_dataset == 5) { 
  dataset <- "Quito"
  time_in_timestamp = 'FALSE'
  data <- dbGetQuery(con, "SELECT longitud as longitude, latitud as latitude, fecha  as unixtime, id as file_name FROM quitomin")
  
}else if (selected_dataset == 6) { 
  dataset <- "California"
  time_in_timestamp = 'TRUE'
  data <- dbGetQuery(con, "SELECT longitud as longitude, latitud as latitude, fecha  as unixtime, id as file_name 	FROM california")
  
}else if (selected_dataset == 7){ 
  dataset <- "Sanfrancisco"
  time_in_timestamp = 'TRUE'
  data <- dbGetQuery(con, "SELECT longitud as longitude, latitud as latitude, fecha as unixtime, id as file_name FROM sanfrancisco")
  
}else if (selected_dataset == 8) { 
  dataset <- "Circulartotal"
  time_in_timestamp = 'FALSE'
  data <- dbGetQuery(con, "SELECT longitud as longitude, latitud as latitude, fecha as unixtime, id as file_name FROM circulartotal")
  
}else if (selected_dataset == 9) {
  dataset <- "Mediacircular"
  time_in_timestamp = 'FALSE'
  data <- dbGetQuery(con, "SELECT longitud as longitude, latitud as latitude, fecha as unixtime, id as file_name FROM mediacircular")
  
}else if (selected_dataset == 10) {
  dataset <- "Circularlineal"
  time_in_timestamp = 'FALSE'
  data <- dbGetQuery(con, "SELECT longitud as longitude, latitud as latitude, fecha as unixtime, id as file_name FROM circularlineal")
  
}else if (selected_dataset == 11) {
  dataset <- "Mitadcircular"
  time_in_timestamp = 'FALSE'
  data <- dbGetQuery(con, "SELECT longitud as longitude, latitud as latitude, fecha as unixtime, id as file_name FROM mitadcircular")
  
}

dbDisconnect(con)

#####     Preprocesamiento para Kalman     #####
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
array_name <- data.frame(file_name = unique(data$file_name)) 

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

#####     Recorrido de la lista de trayectorias y aplicación de Kalman     #####
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

#####     Conteo Coincidencias     #####

coincidencia_per_dim <- function(dimention){
  cant_points <- nrow(dimention)
  coincide <- 0
  no_coincide <- 0
  for (count in 1:cant_points){
    if(dimention$Medicion[count] == dimention$MedicionK[count]){
      coincide <- coincide + 1
    }else{
      no_coincide <- no_coincide + 1
    }
  }
  return (c(coincide = coincide, no_coincide = no_coincide, total = cant_points))
}

contar_coinc <- function(latitude, longitude){
  c_latitud <- coincidencia_per_dim(latitude)
  c_longitud <- coincidencia_per_dim(longitude)
  
  return (t(data.frame(Latitud = c_latitud,Longitud = c_longitud)))
}

coincidencias <- contar_coinc(latitude = binnacle_for_latitude,longitude = binnacle_for_longitude)



cant_points <- nrow(binnacle_for_latitude)
dimention <- binnacle_for_latitude

coincide_lat <- 0
coincide_long <- 0
no_coincide_lat <- 0
no_coincide_long <- 0

for (count in 1:cant_points){
  if(dimention$Medicion[count] == dimention$MedicionK[count]){
    #print("Coincide")
    coincide <- coincide + 1
  }else{
    #print("No coincide")
    no_coincide <- no_coincide + 1
  }
}

coincidencias <- c(coincide = coincide, no_coincide = no_coincide)
#####     RMSE     #####

#####     Guardado/Almacenamiento de Resultados     #####
establecer_directorio_kalman(dataset)

guardar_graficos(original = lista_trayectorias_sinprocesar, kalman = trajectory_after_kalman)

#install_phantomjs(version = "2.1.1",baseURL = "https://github.com/wch/webshot/releases/download/v0.3.1/")
#webshot::install_phantomjs()
guardar_en_mapa(original = lista_trayectorias_sinprocesar, kalman = trajectory_after_kalman)

save_binnacles()

save_trajectories(trajectory_after_kalman)

#####     MAIN - Algoritmo de Simplificacion     #####
print("Cuantas veces quieres ejecutar el algortimo?")
cant_ejecuciones <- length(lista_trayectorias)

print("Ingresa el Epsilon para cada experimento\n")
for (i in 1:as.numeric(cant_ejecuciones) ) {
  dist_threshold[i]<-Epsilon
}
#####     Variables de aplicación     #####
I_3600 = 1 / 3600.0
Epsilon<- 0.00001

result<-list()
time_proc<- list()
dist_threshold<-list()
lista_trayectorias <- list ()
lista_trayectorias_sinprocesar <- list ()
#####     Aplicacion del Algoritmo  de simplificacion     ######
##  Adaptar a seleccion
for (j in 1:as.numeric(cant_ejecuciones)) {
  print(paste("RDP Trayectoria ",j ," de ", cant_trayectorias, sep = ""))
  ini_iter<- as.POSIXct(Sys.time()) 
  result[[j]]<-  drp(points = kalman_libreria[[j]], dist_threshold = dist_threshold[[j]])
  end_iter<- as.POSIXct(Sys.time())
  time_proc[[j]]<-c(ini_iter,end_iter)
}

#Remover duplicados
for (k in 1:as.numeric(cant_ejecuciones)) {
  result[[k]]<- result[[k]][!duplicated(result[[k]]),]
}

##  Actualizacion de infomracion de trayectorias
for (p in  1:as.numeric(cant_ejecuciones)) {
  aux <- nrow(result[[p]])
  trayectorias_informacion$Cant_puntos_simplificados[p] <- aux
}

#####     Almacenamiento de los resultados del algoritmo de simplificación     #####
write.table(trayectorias_informacion, file ="trayectoria informacion.csv" , sep = ";", row.names = FALSE, col.names = TRUE)
write.table(do.call(rbind,lista_trayectorias), file ="lista_trayectorias.csv" , sep = ";", row.names = FALSE, col.names = TRUE)
write.table(do.call(rbind,result), file ="lista_trayectorias_simplificada.csv" , sep = ";", row.names = FALSE, col.names = TRUE)

#####     Graficacion de cada trayectoria (original y simplificada)     ######
for (count in trayectorias_informacion$No_Trayectoria) {
  original<- lista_trayectorias[[count]]
  simplificada<-result[[count]]
  png(paste("trayectoria",count,".png",sep = ""), width = 1024, height = 720, units = 'px', pointsize = 20 ,res = NA)
  plot(original[,1:2], type="b", cex=1.2,col=1, main=paste("trayectoria_",count,sep = ""))
  points(simplificada[,1:2], type="b", cex=1.2,col=2)
  legend("bottomright",legend=c("Original","Simplificada"),col=c(1,2), pch=1,bty="n",ncol=1,cex=1,pt.cex=1)
  dev.off()
}

#####     Margen de Error     #####
##  Tabla de tabulacion
list_size<- list()
for (a in 1:cant_ejecuciones) {
  list_size[[a]]<-nrow(result[[a]]) 
}

##  Razon de compresion 
values_z<-data.frame(valor_z=c(1.28,1.44,1.65,1.96,2.58),valor_confianza=c("80%","85%","90%","95%","99%"))

print("Indique el valor de confianza:")
print(values_z)
opcion <- 1.96

rc<- list()
for (count in 1:as.numeric(cant_ejecuciones)) {
  rc[[count]] <- (1 - trayectorias_informacion$Cant_puntos_simplificados[count]/ trayectorias_informacion$Cant_puntos_originales[count]) *100
  #rc[[count]]<-(1 - (list_size[[count]]/nrow(todas)))*100
}

## TODO: Funcion para mergen de error

#####     Resultados de los experimentos     #####

detalles<- data.frame(n_experimento= 1:as.numeric(cant_ejecuciones), trayectoria_inicial= trayectorias_informacion$Cant_puntos_originales,
                      epsilon= format(do.call(rbind,dist_threshold), scientific=F), trayectoria_final= trayectorias_informacion$Cant_puntos_simplificados,
                      tiempo_proceso_Segundos= paste(do.call(rbind,time_proc)[,2]-do.call(rbind,time_proc)[,1],sep = ""),margen_Error=do.call(rbind,margen_e),razon_compresion=paste(do.call(rbind,rc),sep = ""))

write.table(detalles, file = "Resultados de los experimentos.csv", sep = ";", row.names = FALSE, col.names = TRUE)

png(paste("grafica_original.png",sep = ""), width = 1024, height = 720, units = 'px', pointsize = 20 ,res = NA)
graficar_inicial(lista_trayectoria = lista_trayectorias,trayectorias_informacion = trayectorias_informacion)
dev.off()

png(paste("grafica_simplificada.png",sep = ""), width = 1024, height = 720, units = 'px', pointsize = 20 ,res = NA)
graficar_inicial(lista_trayectoria = result,trayectorias_informacion = trayectorias_informacion)
dev.off()