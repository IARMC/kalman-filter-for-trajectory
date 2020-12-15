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
