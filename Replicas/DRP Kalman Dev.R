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


#####     Funciones de Tratamiento del Filtro de Kalman     #####


#####     Funciones de Almacenamiento de Resultados del Filtro de Kalman     #####


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

#####     ABC     #####