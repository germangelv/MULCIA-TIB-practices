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
  ## low.extrems â† serie aritmética de razón offset con inicio en 1
  ## y con (tamaño de seq) - window.length + 1 como tope superior
  ## results â† matriz vacía (numérico)
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
# Calculo la media para despues
media <- GC(c(genomaHi))
media # 0.3815029
# Preparo para local.composition
longitud <- 20000
desplazamiento <- 2000
# Calculo contenido local y global
analisis <- local.composition(seq = genomaHi, func = GC, window.length = longitud, offset = desplazamiento)
# Verifco o convierto a numerico
lapply(analisis[1][sapply(analisis[1], is.numeric)], max)
# Lo paso a dataframe para poder manipular
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
## 1500000 por encima del 47.075% dada que la emedia es de 38.15029%

## ¿Puede esto indicar transferencia horizontal de genes?
## Si, la zona donde diverge a un 46% es un indicador de transferencia horizontal de genes
######################################################################

######################################################################
## (4) Identificar los puntos en los que se observa la mayor
## desviación
## Se puede ver un punto de maxima desviacion entre las posiciones 1568001 a 1572001 en aproximadamente 0.47075%
######################################################################

######################################################################
## (5) Comprobar, visualizando los datos y con un test de
## Shapiro-Wilk, si el contenido en GC local en los distintos tramos
## sigue una distribución normal.
lapply(analisis[1][sapply(analisis[1], is.numeric)], shapiro.test)
## Dado que 2.2e-16 es inferior a 0.05 los datos de la muestra no provienen de una distribución normal
######################################################################

######################################################################
## (6) Determinar, utilizando un test adecuado, si la diferencia entre
## los tramos es significativa.
# Utilice 1.2 porque de esa manera podia facilmente separar y unir la secuencia
fragmento <-convertframe [ convertframe$results > media*1.2, ]
resto <- convertframe [ convertframe$results <= media*1.2, ]
library(stats)
stat2 <- wilcox.test(fragmento$results, resto$results)
stat2
## Se calculo una prueba Wilcoxon entre el fragmento de interes y el resto de 
## la secuencia y hay mucha evidencia en contra de la hipotesis nula de que ambos
## fragmentos puedan ser similares por lo que tenemos evidencia estadistica para
## considerar una divergencia y por lo tanto una posible transferencia horizontal de genes
######################################################################

######################################################################
## (7) Visualizar, utilizando diagramas de caja y bigote, la
## significancia del resultado obtenido en el apartado anterior.

## - Un asterisco p-valor < 0.05
## - Dos asteriscos p-valor < 0.01
## - Tres asteriscos p-valor < 0.001
boxplot(fragmento$results, convertframe$results,
        main = "Comparación entre fragmento divergente y secuencia",
        at = c(1,3),
        names = c("Fragmento","Secuencia"),
        col = c("orange", "yellow"),
        pch=8,
        sub = "W = 3608, p-value = 0.0005532",
        vertical = TRUE)
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
## Bacteriology, June 1997, p. 3899â€“3913
######################################################################

## Analisis NC_000909
hi <- read.fasta("D:\\Proyecto\\MULCIA-TIB-practices\\entrega_1\\NC_000909.fasta", seqtype="AA")
genomaHi <- getSequence(hi)[[1]]
# Obtengo media
media <- GC(c(genomaHi))
media # 0.3142
# Preparo para local.composition
longitud <- 20000
desplazamiento <- 2000
# Calculo contenido local y global
analisis <- local.composition(seq = genomaHi, func = GC, window.length = longitud, offset = desplazamiento)
# Verifco o convierto a numerico
lapply(analisis[1][sapply(analisis[1], is.numeric)], max) # 0.4093
# Lo paso a dataframe para poder manipular
convertframe <- data.frame(X=c(analisis[2])  , Y=c(analisis[1]), stringsAsFactors = FALSE)
plot(convertframe, 
     type = "l",
     pch=".",
     col="blue")
## Se pueden ver dos maximos en la secuencia en la posicion 144001 a 148001 en aproximadamente 0.4093000 y
## en la posicion 632001 a 636001 en aproximadamente 0.4073500

# Analizo distribución
lapply(analisis[1][sapply(analisis[1], is.numeric)], shapiro.test) # W = 0.74377, p-value < 2.2e-16
## Dado que 2.2e-16 es inferior a 0.05 los datos de la muestra no provienen de una distribución normal

#Separo divergencias
fragmento1 <-convertframe [ convertframe$positions < 500000, ]
fragmento1_analyzed <-fragmento1 [ fragmento1$results > media*1.2, ]
fragmento1_analyzed
fragmento2 <-convertframe [ convertframe$positions > 500000, ]
fragmento2_analyzed <-fragmento2 [ fragmento2$results > media*1.2, ]
fragmento2_analyzed
resto <- convertframe [ convertframe$results <= media*1.2, ]

library(stats)
stat1 <- wilcox.test(fragmento1_analyzed$results, resto$results)
stat1 # W = 6456, p-value = 1.11e-06
stat2 <- wilcox.test(fragmento2_analyzed$results, resto$results)
stat2 # W = 6456, p-value = 1.11e-06
## Se calculo una prueba Wilcoxon entre el fragmento de interes y el resto de 
## la secuencia y hay mucha evidencia en contra de la hipotesis nula de que ambos
## fragmentos puedan ser similares por lo que tenemos evidencia estadistica para
## considerar una divergencia y por lo tanto una posible transferencia horizontal de genes
## en ambos fragmentos

boxplot(fragmento1_analyzed$results, fragmento2_analyzed$results, convertframe$results,
        main = "Comparación entre fragmentos divergentes y secuencia",
        at = c(1,2,4),
        names = c("Fragmento 1","Fragmento 2","Secuencia"),
        col = c("orange", "purple", "yellow"),
        pch=8,
        sub = "W = 6456, p-value = 1.11e-06",
        vertical = TRUE)
