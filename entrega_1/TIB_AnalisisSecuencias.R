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



######################################################################

######################################################################
## (2) Utilizar la función getSequence para obtener, del primer
## elemento (y único) de hi (un objeto de la clase SeqFastadna
## anterior), la secuencia de ADN. Guardar dicha información como
## genomaHi



######################################################################

######################################################################
## La siguiente función calcula (para una secuencia, seq) el valor
## local de una función, func, utilizando una ventada de longitud
## window.length y un desplazamiento offset.


##si tengo secuencia ACCAAAAAAAACA
##puedo ver el contenido de aaaaes en el camino con u offset
##ws=5, offset =2
## eso el contenido es loq ue hace esta funcion
## a secuencia sq
## la funcion de analisis
## la ventana de offset
## el offset de dicha ventana ala nalizar al seq



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

library("seqinr")
GC(c("G","T","A","C"))
## GC nos dice la proporcion de nucleotidos de la cadena
## la R la considero como una G o C


## Aplicamos GC con ventana 20k y desplazamiento (1k creo)no se cuanto y graficamos un plot
## el % de G o Cs y el grafico sera de 2Millones porque es largo

## a
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
## ¿Puede esto indicar transferencia horizontal de genes?
######################################################################

######################################################################
## (4) Identificar los puntos en los que se observa la mayor
## desviación



######################################################################

######################################################################
## (5) Comprobar, visualizando los datos y con un test de
## Shapiro-Wilk, si el contenido en GC local en los distintos tramos
## sigue una distribución normal.

r <- shapira.test()
## menor a 95% 
rip.value


######################################################################

######################################################################
## (6) Determinar, utilizando un test adecuado, si la diferencia entre
## los tramos es significativa.

r <- wilcox.test(fragmento, resto)

rip.value

######################################################################

######################################################################
## (7) Visualizar, utilizando diagramas de caja y bigote, la
## significancia del resultado obtenido en el apartado anterior.

## - Un asterisco p-valor < 0.05
## - Dos asteriscos p-valor < 0.01
## - Tres asteriscos p-valor < 0.001



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

