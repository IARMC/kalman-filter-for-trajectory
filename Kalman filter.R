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
sapply(dataset, function(x) sum(is.na(x)))
dataset <- dataset[!is.na(dataset$longitude),]
dataset <- dataset[!is.na(dataset$latitude),]
dataset <- dataset[!is.na(dataset$file_name),]
dataset <- dataset[!is.na(dataset$unixtime),]

vk <- c()
vlat <- c()
vlon <- c()

#####     Division de los datos por trayectorias     #####
array_name <- data.frame(file_name=unique(data$file_name)) 