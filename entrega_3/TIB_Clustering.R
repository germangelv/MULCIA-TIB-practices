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



######################################################################
## (3) Para comenzar el análisis de estos datos, comenzaremos usando
## un algoritmo de clústering jerárquico, mediante la función hclust
## de la biblioteca estándar stats. Generar una gráfica que represente
## la jerarquía de clústeres del conjunto de datos. Haga que las
## etiquetas de cada rama sea la especie que representa.
######################################################################



######################################################################
## (4) Mirando la gráfica anterior se puede observar, en cierta
## manera, cómo están estructurados los datos de la tabla. Viendo la
## jerarquía generada se pueden diferenciar n grupos distintos, con
## los que podemos comenzar a trabajar. Vamos a «cortar» el árbol
## en el número de clústeres correspondientes, y a mostrar una tabla
## cuántos elementos hay en cada uno de los grupos.
######################################################################



######################################################################
## (5) De hecho, se pueden visualizar estos datos para distintas
## agrupaciones de una vez. Recuerda la función sapply para aplicar
## la misma función a cada uno de los elementos de una lista.
######################################################################



######################################################################
## (6) Ahora, intenta visualizar los individuos de cada uno de los
## grupos. Puedes hacerlo de uno en uno o puedes ver varios grupos
## usando la función sapply.
######################################################################



######################################################################
## (7) Prueba a visualizar la misma salida con distinto número de
## clústeres.
######################################################################



######################################################################
## (8) Una vez tenemos el número de clústeres deseado para representar
## los grupos de individuos, pasamos a usar clústering no jerárquico.
## En este caso, usaremos la función pam de la biblioteca cluster.
## Una vez usado, muestre en una tabla los resultados obtenidos
## comparados con los grupos obtenidos en el ejercicio anterior.
## ¿Qué se puede concluir de esto? ¿Esperabas este resultado?
######################################################################



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


