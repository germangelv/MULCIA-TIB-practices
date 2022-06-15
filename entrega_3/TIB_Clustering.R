## Carga de bibliotecas adicionales
## library("biblioteca") o a través del menú packages

## Algunas funciones de interés

## daisy(dataframe) calcula todas las distancias entre observaciones
## del conjunto de datos.
## cutree(árbol, k) corta los datos provinientes del árbol en k
## grupos distintos.
## sapply(vector.x, function(y) ...) aplica la función function(y)
## sobre todos los elementos de x, siendo y el elemento en cuestión.

######################################################################
## El citocromo c es una proteína que funciona como transportador
## mitocondrial entre dos complejos respiratorios.
######################################################################

######################################################################
## (1) Cargar los datos de amino.acid.sequence.1972 (incluídos en la
## biblioteca cluster.datasets). Una vez cargados, visualizar las
## primeras entradas de la tabla para ver cómo están estructurados los
## datos.
######################################################################
install.packages("cluster.datasets")
library(cluster)
library(cluster.datasets)
data(amino.acid.sequence.1972)
aminoacid <- amino.acid.sequence.1972
head(aminoacid)
#        species p.1 p.2 p.3 p.4 p.5 p.6 p.7 p.8 p.9 p.10 p.11 p.12 p.13 p.14 p.15 p.16 p.17
# ...
#    p.18 p.19 p.20 p.21 p.22 p.23 p.24 p.25 p.26 p.27 p.28 p.29 p.30 p.31 p.32 p.33 p.34 p.35
# ...
#    p.36
######################################################################
## (2) Calcular la tabla de distancias entre los elementos de la
## tabla. Si se intenta hacer directamente, R devolverá una tabla
## llena de elementos NA. Esto es debido a que los valores tienen un
## valor alfanumérico. Para ello, previo al cálculo debemos
## transformar los valores de manera adecuada para que puedan ser
## tratados. Además, si uno se fija, se puede observar como en los
## datos hay alguna columna que no se precisa para la ejecución de
## cualquiera de los algoritmos estudiados.
######################################################################
# Normalizo datos
aminoacid_numeric <- sapply(aminoacid, as.numeric)
# Aplico calculo de distancias
distances <- daisy(aminoacid_numeric)
distances
# Dissimilarities :
#           1        2        3        4        5        6        7        8        9       10
# 2  2.027588                                                                                 
# 3  5.459446 5.823802                                                                        
# 4  5.364492 5.734884 1.013794                                                               
# 5  3.511885 4.055175 3.362374 3.205897                                                      
# 6  3.511885 4.055175 3.655285 3.511885 1.433721                                             
# 7  3.793269 4.301163 5.823802 5.734884 3.205897 3.205897                                    
# 8  3.926406 4.419025 4.055175 3.926406 2.266912 2.266912 2.266912                           
# 9  3.362374 3.926406 3.793269 3.926406 2.266912 2.266912 3.362374 2.867442                  
# 10 6.491447 6.800735 4.755114 4.645787 4.861984 5.823802 7.096556 6.249444 5.734884         
# 11 4.533824 4.966555 7.096556 7.023769 4.533824 4.055175 4.301163 5.267827 4.645787 7.787098
# 12 5.267827 5.644565 4.055175 3.926406 4.179979 3.655285 5.997685 4.966555 4.301163 5.169354
# 13 5.169354 5.552777 7.096556 7.023769 4.533824 4.755114 4.533824 5.267827 4.645787 6.950220
# 14 5.911383 6.249444 5.644565 5.552777 5.169354 5.364492 5.911383 5.267827 5.644565 5.997685
# 15 6.411795 6.411795 6.800735 6.724747 5.911383 6.411795 7.310571 7.096556 5.644565 7.096556
# 16 6.875884 7.168604 6.491447 6.411795 5.911383 6.411795 7.023769 6.800735 5.997685 7.653975
# 17 8.236099 7.982620 6.491447 6.411795 7.023769 7.023769 8.110350 7.239936 7.239936 7.787098
#          11       12       13       14       15       16
# 2                                                       
# 3                                                       
# 4                                                       
# 5                                                       
# 6                                                       
# 7                                                       
# 8                                                       
# 9                                                       
# 10                                                      
# 11                                                      
# 12 5.997685                                             
# 13 4.966555 6.331140                                    
# 14 7.310571 5.823802 8.110350                           
# 15 7.023769 7.653975 6.875884 8.602325                  
# 16 7.586538 8.173467 7.720823 7.982620 5.364492         
# 17 8.602325 7.096556 9.067647 7.310571 8.110350 7.720823
# 
# Metric :  euclidean 
# Number of objects : 17
# Warning message:
# In daisy(aminoacid_numeric) :
#   binary variable(s) 2, 3, 6, 7, 8, 9, 10, 12, 14, 15, 17, 21, 22, 23, 24, 25, 26, 27, 31, 32 # # treated as interval scaled
######################################################################
## (3) Para comenzar el análisis de estos datos, comenzaremos usando
## un algoritmo de clústering jerárquico, mediante la función hclust
## de la biblioteca estándar stats. Generar una gráfica que represente
## la jerarquía de clústeres del conjunto de datos. Haga que las
## etiquetas de cada rama sea la especie que representa.
######################################################################
# Aplico clustering
cls <- hclust(distances)
plot(cls)
######################################################################
## (4) Mirando la gráfica anterior se puede observar, en cierta
## manera, cómo están estructurados los datos de la tabla. Viendo la
## jerarquía generada se pueden diferenciar n grupos distintos, con
## los que podemos comenzar a trabajar. Vamos a «cortar» el árbol
## en el número de clústeres correspondientes, y a mostrar una tabla
## cuántos elementos hay en cada uno de los grupos.
######################################################################
data.frame(cutree(cls, 5))
######################################################################
## (5) De hecho, se pueden visualizar estos datos para distintas
## agrupaciones de una vez. Recuerda la función sapply para aplicar
## la misma función a cada uno de los elementos de una lista.
######################################################################
sapply(cls, cutree)
######################################################################
## (6) Ahora, intenta visualizar los individuos de cada uno de los
## grupos. Puedes hacerlo de uno en uno o puedes ver varios grupos
## usando la función sapply.
######################################################################

######################################################################
## (7) Prueba a visualizar la misma salida con distinto número de
## clústeres.
######################################################################
plot(cls)

######################################################################
## (8) Una vez tenemos el número de clústeres deseado para representar
## los grupos de individuos, pasamos a usar clústering no jerárquico.
## En este caso, usaremos la función pam de la biblioteca cluster.
## Una vez usado, muestre en una tabla los resultados obtenidos
## comparados con los grupos obtenidos en el ejercicio anterior.
## ¿Qué se puede concluir de esto? ¿Esperabas este resultado?
######################################################################
# Eleji 5 porque eran 5 la separabilidad pero despues de ver el grafico en (9) deje 3
# No esperaba ese resultado
pam.res <- pam(distances, 3)
######################################################################
## (9) Por último, muestre un gráfico de silueta del último clústering
## realizado. En el gráfico aparece la anchura media de la silueta,
## que indica la «naturalidad» de la estructura encontrada.
## De 0,71 a 1, se ha encontrado una estructura fuerte.
## De 0,51 a 0,7, se encuentra una estructura razonable.
## De 0,26 a 0,5, la estructura es débil y podría ser artificial.
## Si es menor que 0,25, no se ha encontrado una estructura
## sustancial.
## ¿Por qué crees que has obtenido los datos que has obtenido?
######################################################################
# Ploteo pam
plot(pam.res)
# Observo como queda
head(pam.res$clustering)
# La anchura media de la silueta es 0.27. No encuentro una estructura 
# sustancial por la naturaleza de los datos
