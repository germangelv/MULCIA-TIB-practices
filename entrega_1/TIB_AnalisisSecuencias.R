## Carga de bibliotecas adicionales
## library("biblioteca") o a través del menú Packages

## Variables y asignación 
## variable <- valor

## Funciones
## función(datos, ..., arg1 = valor1, ...)

## Vectores: Todos los elementos de la misma naturaleza
## * Creación:
##   - c(dato1, dato2, ...)
##   - inicio:fin siendo inicio y fin dos números cualesquiera
## * Selección:
##   - Por los índices de las posiciones: vector[vector.índices]
##   - Por exclusión de las posiciones: vector[-vector.índices]

## Matrices
## Acceso:
## - Por el índice de la columna: matriz[, índice.col]

## Listas: Cada elemento puede tener una naturaleza diferente.
## Acceso:
## - Por el índice del elemento: lista[[índice]]
## - Por el nombre del elemento (no siempre disponible): lista$nombre

## Gráficos: n.color, n.línea y n.punto son un entero indicando el
## color, el tipo de línea y el tipo de punto, respectivamente (por
## defecto 1 en ambos casos, negro, continua respectivamente).
## * En general: main, xlab, ylab (título y etiquetas de los ejes x e
##   y).
## * De líneas:
##   plot(coord.x, coord.y, type = "l", col = n.color, lty = n.línea)
## * Histograma (de densidad):
##   hist(datos, freq = FALSE)
## * Diagrama de caja y bigote
##   boxplot(datos1, datos2, ..., col = vector.colores,
##           names = vector.nombres)

## Añadidos a un gráfico:
## * Líneas verticales a distancias vector.x del origen
##   abline(v = vector.x, col = vector.color, lty = vector.línea)
## * Líneas horizontales a distancias vector.y del origen
##   abline(h = vector.y, col = vector.color, lty = vector.línea)
## * Líneas
##   lines(x.coord, y.coord, col = n.color, lty = n.línea)
## * Puntos
##   points(x.coord, y.coord, col = n.color, pch = n.punto)
## * Leyenda para un gráfico de líneas
##   legend(x.pos, y.pos, vector.nombres, col = vector.colores, 
##          lty = vector.líneas)

## plot(1:25, pch = 1:25, col = rep(1:8, length.out = 25))
## abline(h = 1:6, col = 1:6, lty = 1:6)

######################################################################
## Haemophilus influenzae (NC_000907) es una bacteria gram-negativa
## que fue erróneamente considerada la causa de la gripe hasta que en
## 1933 se descubrió la causa viral de la misma.

## El genoma de H. influenzae fue el primer genoma en secuenciarse de
## un organismo de vida libre por el grupo de Craig Venter en 1995.
######################################################################

######################################################################
## (1) Descargar su genoma de la base de datos de NCBI (como un
## fichero fasta) y, utilizando la función read.fasta de la biblioteca
## seqinr, guardar dicha información como hi.
library("seqinr")
hi <- read.fasta("D:\\Proyecto\\MULCIA-TIB-practices\\entrega_1\\NC_000907.fasta", seqtype="AA")
######################################################################

######################################################################
## (2) Utilizar la función getSequence para obtener, del primer
## elemento (y único) de hi (un objeto de la clase SeqFastadna
## anterior), la secuencia de ADN. Guardar dicha información como
## genomaHi
genomaHi <- getSequence(hi)[[1]]
######################################################################

######################################################################
## La siguiente función calcula (para una secuencia, seq) el valor
## local de una función, func, utilizando una ventada de longitud
## window.length y un desplazamiento offset.

local.composition <- function (seq, func, window.length, offset) {
  
  ## Paso 1: Inicialización de variables antes del bucle:
  ## low.extrems ← serie aritmética de razón offset con inicio en 1
  ## y con (tamaño de seq) - window.length + 1 como tope superior
  ## results ← matriz vacía (numérico)
  low.extrems <- seq(1, length(seq) - window.length + 1, 
                     by = offset)
  results <- numeric()
  
  ## Paso 2: Para cada uno de los elementos de low.extrems hace:
  for (i in low.extrems) {
    
    ## Paso 2.1: Añadir a results una nueva fila con los valores
    ## de func sobre el trozo de seq entre i e (i + window.length
    ## - 1)
    results <- rbind(results, func(seq[i:(i+window.length-1)]))
  }
  
  ## Paso 3: Devolver una lista con los valores de results y
  ## low.extrems como los elementos results y positions
  ## (respectivamente).
  list(results = results, positions = low.extrems)
}
######################################################################

######################################################################
## (3) Calcular en contenido en GC, global y local utilizando una
## ventana de longitud 20000 y un desplazamiento de 2000. Representar, 
## utilizando un gráfico de líneas, los resultados obtenidos.
## Nota: Utilizar la función GC de la biblioteca seqinr.
GC(c(genomaHi))
longitud <- 20000
desplazamiento <- 2000
local.composition(seq = genomaHi, func = GC, window.length = longitud, offset = desplazamiento)
analisis <- local.composition(seq = genomaHi, func = GC, window.length = 20000, offset = 2000)
lapply(analisis[1][sapply(analisis[1], is.numeric)], max)
convertframe <- data.frame(X=c(analisis[2])  , Y=c(analisis[1]), stringsAsFactors = FALSE)
plot(convertframe, 
     type = "l",
     pch=".",
     col="blue")
######################################################################

######################################################################
## El contenido en GC de un genoma, la proporción de los nucleótidos g
## y c, es una característica muy específica de cada organismo en
## concreto.

## La transferencia horizontal de genes hace referencia a la
## transmisión de genes entre organismos vivos (comúnmente de
## diferentes especies) sin que uno sea el progenitor del otro.

## Debido a que el contenido en GC es una de las características más
## específicas de cada especie, la identificación de zonas en el
## genoma cuyo contenido en GC diverge significativamente del
## contenido GC global puede indicar la existencia de un evento de
## transferencia horizontal de genes.

## ¿Se observa en el gráfico anterior una desviación significativa?
## Si, se observa una desviacion significativa cuando es analizada la cadena superando el gen
## 1500000 por encima del 47.075% dada que al emdia es de 38.15029%

## ¿Puede esto indicar transferencia horizontal de genes?
## Si, la zona donde diverge a un 46% es un indicador de tranbsferencia horizontal de genes
######################################################################

######################################################################
## (4) Identificar los puntos en los que se observa la mayor
## desviación
## Se puede ver un punto de maxima desviacion entre las posiciones 1568001 a 1572001 en aproximadamente 0.4707500
######################################################################

######################################################################
## (5) Comprobar, visualizando los datos y con un test de
## Shapiro-Wilk, si el contenido en GC local en los distintos tramos
## sigue una distribución normal.
lapply(analisis[1][sapply(analisis[1], is.numeric)], shapiro.test)
## Dado que 2.2e-16 es inferior a 0.05 los datos de la muestra no provienen de una distribuci�n normal
######################################################################

######################################################################
## (6) Determinar, utilizando un test adecuado, si la diferencia entre
## los tramos es significativa.
# install.packages("tidyverse")
# install.packages("rstatix")
# install.packages("ggpubr")
library(tidyverse)
library(rstatix)
library(ggpubr)
media <- convertframe %>% get_summary_stats(results, type = "median_iqr")
media
stat <- convertframe %>% select(results) %>% rstatix::wilcox_test(results, mu = 0.3815)
stat
# .y.     group1 group2         n statistic     p
#  * <chr>   <chr>  <chr>      <int>     <dbl> <dbl>
#  1 results 1      null model   906   184574. 0.011

## Se calcul� una prueba Wilcoxon con signo para evaluar si la media GC local era diferente a la media GC y
## dado que p es menor a 0.05 pero no menor a 0.01 o 0.001 podemos decir que tenemos evidencias estadisticas
## significativas para aceptar la hipotesis y considerar que existen datos similares a la media

library(coin) # wilcox.test wilcoxsign_test
# Preparo datos para funcion
rfirst <- convertframe %>% select(results)
rfirst = rfirst[-1,]
rend<- convertframe %>% select(results)
rend = rend[-nrow(rend),]
datos <- data.frame(antes = rend, despues = rfirst)

# calculo la media de variacion entre segmentos
r1 <- 0
for (i1 in 1:length(rend)){
  r1 = r1+ c(abs(rfirst[1]-rend[1]))
}
r1 = r1/length(rend)

stat2 <- wilcox.test(x = rend, y = rfirst, alternative = "two.sided", mu = r1, paired = TRUE)
stat2
# data:  rend and rfirst
# V = 175633, p-value = 0.000348
# alternative hypothesis: true location shift is not equal to 0.00045
## Se calcul� una prueba Wilcoxon con signo para evaluar cada cadena con la subsecuentes y 
## dado que p es menor a 0.05, 0.01 y 0.001 podemos decir que tenemos evidencias estadisticas
## significativas para descartar la hipotesis y aceptar que los datos son diferentes a la media de variacion
## de la cadena
# require(coin)
# wilcoxsign_test(antes ~ despues, data = datos, distribution = "exact", mu = 0.3815)
######################################################################

######################################################################
## (7) Visualizar, utilizando diagramas de caja y bigote, la
## significancia del resultado obtenido en el apartado anterior.

## - Un asterisco p-valor < 0.05
## - Dos asteriscos p-valor < 0.01
## - Tres asteriscos p-valor < 0.001
bxp <- ggboxplot(
  convertframe$results, width = 0.5, add = c("mean", "jitter"), 
  ylab = "Resultados GC", xlab = FALSE
)
bxp
######################################################################

######################################################################
## Tarea:
## Realizar un análisis parecido para el genoma de la bacteria
## Methanocaldococcus jannaschii (NC_000909). Se trata de una archaea
## metanógena termofílica que habita ventanas hidrotermales creciendo
## usando como fuente de energía dióxido de carbono e hidrógeno y
## produciendo metano como producto secundario de su metabolismo.

## El grupo de Craig Venter fue el primero en secuenciar su genoma en
## 1996 que constituyó el primer genoma de archaea en secuenciarse
## completamente. La secuenciación de su genoma produjo evidencias
## claves para la existencia de los tres dominios de la vida (Archaea,
## Bacteria y Eukarya).

## Ref: "Compositional Biases of Bacterial Genomes and Evolutionary
## Implications", S. Karlin, J. Mrázek, A.M. Campbell, Journal of
## Bacteriology, June 1997, p. 3899–3913


######################################################################

## Analisis NC_000909
hi <- read.fasta("D:\\Proyecto\\MULCIA-TIB-practices\\entrega_1\\NC_000909.fasta", seqtype="AA")
genomaHi <- getSequence(hi)[[1]]
GC(c(genomaHi)) # 0.3142
longitud <- 20000
analisis <- local.composition(seq = genomaHi, func = GC, window.length = longitud, offset = desplazamiento)
lapply(analisis[1][sapply(analisis[1], is.numeric)], max) # 0.4093
convertframe <- data.frame(X=c(analisis[2])  , Y=c(analisis[1]), stringsAsFactors = FALSE)
plot(convertframe, 
     type = "l",
     pch=".",
     col="blue")
## Se pueden ver dos maximos en la secuencia en la posicion 144001 a 148001 en aproximadamente 0.4093000 y
## en la posicion 632001 a 636001 en aproximadamente 0.4073500
lapply(analisis[1][sapply(analisis[1], is.numeric)], shapiro.test)
## Dado que 2.2e-16 es inferior a 0.05 los datos de la muestra no provienen de una distribuci�n normal

media <- convertframe %>% get_summary_stats(results, type = "median_iqr")
media #??? 0.311

stat <- convertframe %>% select(results) %>% rstatix::wilcox_test(results, mu = 0.3815)
stat

# Preparo datos para funcion
rfirst <- convertframe %>% select(results)
rfirst = rfirst[-1,]
rend<- convertframe %>% select(results)
rend = rend[-nrow(rend),]
datos <- data.frame(antes = rend, despues = rfirst)

# Calculo la media de variacion entre segmentos
r1 <- 0
for (i1 in 1:length(rend)){
  r1 = r1+ c(abs(rfirst[1]-rend[1]))
}
r1 = r1/length(rend)
r1 # 0.0041

stat2 <- wilcox.test(x = rend, y = rfirst, alternative = "two.sided", mu = r1, paired = TRUE)
stat2
# data:  rend and rfirst
# V = 20774, p-value < 2.2e-16
# alternative hypothesis: true location shift is not equal to 0.0041
## Se calcul� una prueba Wilcoxon con signo para evaluar cada cadena con la subsecuentes y 
## dado que p es menor a 0.05, 0.01 y 0.001 podemos decir que tenemos evidencias estadisticas
## significativas para descartar la hipotesis y aceptar que los datos son diferentes a la media de variacion

# Grafico
bxp <- ggboxplot(
  convertframe$results, width = 0.5, add = c("mean", "jitter"), 
  ylab = "Resultados GC", xlab = FALSE
)
bxp

