rm(list = ls())
options( digits = 15 ) ## para cambiar el numero de decimales. Por defecto viene 7.

#Librerias
ipak <- function(pkg) { #cargar e instalar paquetes de manera dinamica
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
#Funciones
crear_directorio<- function(folders, algorithm, dataset){
  
  simbol<-"\\";
  root<- "C:";  	##directorio raiz
  setwd(paste("C:",simbol,sep="")); ##C:
  
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







packages <- c("dlm","RPostgreSQL","rlang", "ggplot2", "caret", "class", "mapview",  "compare", "pracma" , "stringr","SpatialTools","matlib", "dplyr","chron","lubridate","zoom","RgoogleMaps","ggmap") #librerias 
ipak(packages)  #llama las librerias que se requeriran para el correcto funcionamiento del algortimo

#Datos Globales
carpeta<- "PRUEBA TD-TR"								                                        #Directorio experimentos
algoritmo<- "TDTR"						                              #Algorithm
algoritmo <- paste(algoritmo, Sys.Date(),sep = "")
#dataset<- "Sanfrancisco"                                                      #Base de datos
I_3600 = 1 / 3600.0
Epsilon<- 0.001






result<-list()                                                        #listas de datos que contendra los resultados por iteracion del algortimo
time_proc<- list()                                                          #tiempo del proceso por cada iteracion
dist_threshold<-list()                                                      #lista de las distancias limite / umbral para cada experimento / iteracion 
lista_trayectorias <- list () 						                                  #lista que almacena todas las trayectorias
lista_trayectorias_sinprocesar <- list ()


# Configuracion Data Source ####
dbdriver <- "PostgreSQL";
host <- 'localhost';
port='5432';
dbname='CPPP2';
user='postgres';
pass='toor';



##CREATE DATA SOURCE
print("CARGANDO DATOS DE LA BD")
drv <- dbDriver(dbdriver)
con <- dbConnect(drv, host = host, port = port, dbname = dbname, user = user, pass = pass)


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
selected_dataset <- 2

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

array_name <- data.frame(file_name=unique(data$file_name)) 
dbDisconnect(con)

#Creacion del directorio
crear_directorio(folders = carpeta, algorithm = algoritmo, dataset = dataset)
#Crea el directorio donde se alojaran los resultadso del algortimo


todas <- data.frame(longitude = c(1:nrow(data)), latitude = 0, file_name = 0, speed = 0, unixtime = 0) ## creo un data frame para almacenar todas las trayectorias
todas$longitude <- data$longitude
todas$latitude  <- data$latitude
todas$file_name <- data$file_name
todas$unixtime  <- if(time_in_timestamp == 'TRUE'){ as.POSIXct(data[,'unixtime']/1000, origin="1970-01-01")}else{data[,'unixtime']}

pto_medio <- matrix(0, nrow=1, ncol=2)
min_latitud <- min(todas[,'latitude']); max_latitud <- max(todas[,'latitude']); 
men_longitud <- min(todas[,'longitude']); max_longitud <- max(todas[,'longitude']);
pto_medio[2] <- (max_latitud + min_latitud)/2   ## valor de x para el centro de la ventana
pto_medio[1] <- (max_longitud + men_longitud)/2## valor de y para el centro de la ventana


#Grafico el mapa del Dataset a evaluar
#plot(data[,c('longitude','latitude')],xlab="Longitude",ylab="Latitude",main=dataset,pch=46, cex = 1.5)

#Separarlo en listas de tray
for(count in 1:length(unique(todas$file_name))){
  tray_corresp <-subset(todas, todas$file_name == array_name$file_name[count])		##divido por trayectorias
  lista_trayectorias_sinprocesar[[count]] <-  tray_corresp
}
lista_trayectorias <- lista_trayectorias_sinprocesar
cant_trayectorias <- length(lista_trayectorias)
print('Cantidad de trayectorias')
print(cant_trayectorias)

#Prepocesar informacion
trayectorias_informacion <-preprocesar_informacion(lista_trayectorias)





#frame_tray <- data.frame(x = lista_trayectorias[[1]]$longitude, y = lista_trayectorias[[1]]$latitude, t=lista_trayectorias[[1]]$unixtime,g=0)
#matrix1 <- frame_tray 

#####    KALMAN

kalman_libreria <- list()

nileBuild <- function(par) {
  dlmModPoly(order = 1, dV = exp(par[1]), dW = exp(par[2]))
}

cant_tr <- length(unique(lista_trayectorias_sinprocesar))

for (count in 1:cant_tr){
  #count <- 1
  print(paste("Aplicando DLM Kalman. Trayectoria ",count," de ", cant_tr, ep = ""))
  tr_i <- lista_trayectorias_sinprocesar[[count]]
  cant_points <- nrow(tr_i)
  frame_tray <- data.frame(longitude = lista_trayectorias_sinprocesar[[count]]$longitude, 
                           latitude = lista_trayectorias_sinprocesar[[count]]$latitude, 
                           unixtime = lista_trayectorias_sinprocesar[[count]]$unixtime,
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
    new_frame_kalman <- data.frame(longitude = nileSmooth1$s[2:tamano_kalman], 
                                   latitude = nileSmooth2$s[2:(tamano_kalman)], 
                                   unixtime = frame_tray[,3], 
                                   g = frame_tray[,4])
  }
  
  kalman_libreria[[count]] <- new_frame_kalman
}



print("Cuantas veces quieres ejecutar el algortimo?")
cant_ejecuciones <- length(lista_trayectorias)

print("Ingresa el Epsilon para cada experimento\n")
for (i in 1:as.numeric(cant_ejecuciones) ) {
  dist_threshold[i]<-Epsilon
}

#Aplicacion del Algoritmo  de simplificacion
for (j in 1:as.numeric(cant_ejecuciones)) {
  print(paste("TD-TR Trayectoria ",j ," de ", cant_trayectorias, sep = ""))
  ini_iter<- as.POSIXct(Sys.time()) 
  result[[j]]<-  td_tr(points = kalman_libreria[[j]], dist_threshold = dist_threshold[[j]])
  end_iter<- as.POSIXct(Sys.time())
  time_proc[[j]]<-c(ini_iter,end_iter)
}


#lista_trayectorias[[42]]
#kalman_libreria[[42]]


#Remover duplicados
for (k in 1:as.numeric(cant_ejecuciones)) {
  result[[k]]<- result[[k]][!duplicated(result[[k]]),]
}



#Actualizacion de la trayectoria informacion
for (p in  1:as.numeric(cant_ejecuciones)) {
  aux <- nrow(result[[p]])
  trayectorias_informacion$Cant_puntos_simplificados[p] <- aux
}

write.table(trayectorias_informacion, file ="trayectoria informacion.csv" , sep = ";", row.names = FALSE, col.names = TRUE)

write.table(do.call(rbind,lista_trayectorias), file ="lista_trayectorias.csv" , sep = ";", row.names = FALSE, col.names = TRUE)

write.table(do.call(rbind,result), file ="lista_trayectorias_simplificada.csv" , sep = ";", row.names = FALSE, col.names = TRUE)



#Graficas cada trauectoria (original y simplificada)

#for (count in trayectorias_informacion$No_Trayectoria) {
#  original<- lista_trayectorias[[count]]
#  simplificada<-result[[count]]
#  png(paste("trayectoria",count,".png",sep = ""), width = 1024, height = 720, units = 'px', pointsize = 20 ,res = NA)
#  plot(original[,1:2], type="b", cex=1.2,col=1, main=paste("trayectoria_",count,sep = ""))
#  points(simplificada[,1:2], type="b", cex=1.2,col=2)
#  legend("bottomright",legend=c("Original","Simplificada"),col=c(1,2), pch=1,bty="n",ncol=1,cex=1,pt.cex=1)
#  dev.off()
#}



#tabla de tabulacion
list_size<- list()
for (a in 1:cant_ejecuciones) {
  list_size[[a]]<-nrow(result[[a]]) 
}

##Razon de compresion 
values_z<-data.frame(valor_z=c(1.28,1.44,1.65,1.96,2.58),valor_confianza=c("80%","85%","90%","95%","99%"))

print("Indique el valor de confianza:")
print(values_z)
opcion<-1.96

rc<- list()
for (count in 1:as.numeric(cant_ejecuciones)) {
  rc[[count]] <- (1 - trayectorias_informacion$Cant_puntos_simplificados[count]/ trayectorias_informacion$Cant_puntos_originales[count]) *100
  #rc[[count]]<-(1 - (list_size[[count]]/nrow(todas)))*100
}


##Margen de error
margen_e<- list()
for (count in 1:as.numeric(cant_ejecuciones)){
  dist_or <- lista_trayectorias[[count]]
  dist_si <- result[[count]]
  
  dist_sed <- 0
  if (nrow(dist_or)<=2) {
    return(dist_or)
  }else{
    delta_e = as.numeric(time_dist(dist_or[nrow(dist_or),],dist_or[1,]), units = "secs") * I_3600 #segundos
    d_lat = dist_or[nrow(dist_or),][,'latitude'] - dist_or[1,][,'latitude']
    d_lon = dist_or[nrow(dist_or),][,'longitude'] - dist_or[1,][,'longitude']
    
    
    delta_e = as.numeric(time_dist(dist_or[nrow(dist_si),],dist_si[1,]), units = "secs") * I_3600 #segundos
    d_lat = dist_si[nrow(dist_si),][,'latitude'] - dist_si[1,][,'latitude']
    d_lon = dist_si[nrow(dist_si),][,'longitude'] - dist_si[1,][,'longitude']
    
    for (i in 1:(nrow(dist_or)-1)) {
      
      delta_i= as.numeric(time_dist(dist_or[i,],dist_or[1,]), units = "secs")  * I_3600
      di_de = delta_i / delta_e
      prim = data.frame(longitude = dist_or[1,]$longitude + d_lon * di_de,
                        latitude = dist_or[1,]$latitude + d_lat * di_de,
                        unixtime=0)
      dist_sed <- dist_sed + loc_dist(dist_or[i,], prim)
    }
  }
  P_dist_euclidean_or <- dist_sed/nrow(dist_or)
  
  dist_sed <- 0
  if (nrow(dist_si)<=2) {
    return(dist_si)
  }else{
    for (i in 1:(nrow(dist_si)-1)) {
      
      delta_i= as.numeric(time_dist(dist_si[i,],dist_si[1,]), units = "secs")  * I_3600
      di_de = delta_i / delta_e
      prim = data.frame(longitude = dist_si[1,]$longitude + d_lon * di_de,
                        latitude = dist_si[1,]$latitude + d_lat * di_de,
                        unixtime=0)
      dist_sed <- dist_sed + loc_dist(dist_si[i,], prim)
    }
    P_dist_euclidean_sim <- dist_sed/nrow(dist_si)
    
  }
  margen_e<- list()
  margen_e[[count]]<- abs(P_dist_euclidean_sim - P_dist_euclidean_or)
}


detalles<- data.frame(n_experimento= 1:as.numeric(cant_ejecuciones), 
                      trayectoria_inicial= trayectorias_informacion$Cant_puntos_originales,
                      epsilon= as.vector(unlist(dist_threshold)), 
                      trayectoria_final= trayectorias_informacion$Cant_puntos_simplificados,
                      tiempo_proceso_Segundos= do.call(rbind,time_proc)[,2]-do.call(rbind,time_proc)[,1],
                      margen_Error=do.call(rbind,margen_e),
                      #razon_compresion=paste(do.call(rbind,rc),sep = ""))
                      razon_compresion=do.call(rbind,rc))

#convert list to vector
#as.vector(unlist(myList))

write.table(detalles, file = "Resultados de los experimentos.csv", sep = ";", row.names = FALSE, col.names = TRUE)

#png(paste("grafica_original.png",sep = ""), width = 1024, height = 720, units = 'px', pointsize = 20 ,res = NA)
#graficar_inicial(lista_trayectoria = lista_trayectorias,trayectorias_informacion = trayectorias_informacion)
#dev.off()

#png(paste("grafica_simplificada.png",sep = ""), width = 1024, height = 720, units = 'px', pointsize = 20 ,res = NA)
#graficar_inicial(lista_trayectoria = result,trayectorias_informacion = trayectorias_informacion)
#dev.off()


#contador de puntos finales
contador <- 0

for (i in 1:length(lista_trayectorias)) {
  contador <- contador + nrow(result[[i]])
}
total_puntos <- c(Cant_Puntos = contador)

medias <- c(sapply(detalles,mean))
#medias <- summary(detalles)
#print(detalles)
#sumatoria de las columnas del dataset detalles
sumatorias <- sapply(detalles,sum)

write.table(total_puntos, file = "Metrica total puntos.txt", sep = ";", row.names = TRUE, col.names = TRUE)
write.table(medias, file = "Metricas Medias.txt", sep = ";", row.names = TRUE, col.names = TRUE)
write.table(sumatorias, file = "Metricas Sumatorias.txt", sep = ";", row.names = TRUE, col.names = TRUE)
#epsilon, 
#summary(detalles)
#info <- sapply(detalles,mean)
#mean(detalles$razon_compresion)
#mean(detalles[,"razon_compresion"])
#detalles$razon_compresion