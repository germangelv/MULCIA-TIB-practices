######################################################################
## (1) Descargar el genoma del bacterio fago lambda, identificado con
## NC_001416.1 de la base de datos NCBI, y cargar la función
## local.composition definida en la práctica anterior (se encuentra en
## el fichero HMM.RData).
######################################################################

# Cargamos RData:
load("HMM.RData")

# Descargamos y leemos el fichero sequence.fasta que contiene el genoma 
# del bacterio faga lambda

library(seqinr)
lambda <- read.fasta("NC_001416.1.fasta",seqtype="DNA")

lambda.sequences <- getSequence(lambda)[[1]]

######################################################################
## (2) Crear un gráfico para comparar el contenido global y local en
## GC del genoma del bacterio fago lambda usando una ventana de
## longitud 5000 y un desplazamiento de 500.
######################################################################

# Aplicamos la función local.coposition a la secuencia usando una ventana
# de longitud 5000 y un desplazamiento de 500:
lambdaGC.local <- local.composition(lambda.sequences, GC, 5000,500)

# Veamos el contenido glocal aplicando directamente la función GC:
labmdaGC.global <- GC(lambda.sequences)

# Gráficamente:

plot(lambdaGC.local[["positions"]],lambdaGC.local[["results"]],
     type="l",pch=0,col="red",cex=1.5,lwd=1,
     xlab="Genome Posicion",ylab="GC result", cex.lab=1.5,
     cex.axis=1.25,
     main="Local GC result in lambda fago")
lines(lambdaGC.local[["positions"]], rep(labmdaGC.global, length(lambdaGC.local[["positions"]])),
      type="l",pch=1,col="blue",cex=1.5,lwd=2)
legend("topright",legend=c("Local G C","Global G C"),col=c("red","blue"),
       lwd=2,cex=0.75)

# Se observan diferencias entre la zona alta y la zona baja de la gráfica
######################################################################
## (3) Escribir una función para el entrenamiento de un HMM (es decir,
## que estime las probabilidades iniciales, matriz de transición y
## matriz de emisión. Como argumentos debe recibir una secuencia de
## símbolos y otra (de la misma longitud) con los correspondientes
## estados.
######################################################################

prueba = hidro.initial.hmm

estimate.hmm <- function (sec.simb, sec.est) {
  ## Paso 1: Definir las siguientes variables
  alf.e <- sort(unique(sec.est))
  m <- length(alf.e)
  alf.x <- sort(unique(sec.simb))
  n <- length(alf.x)
  
  ## Paso 2: Calcular el número de veces que aparece cada estado,
  cuentaEstados <- count(sec.est, 1, alphabet = alf.e)
  
  ## Paso 3: Definir pi como el vector con las probabilidades
  #pi <- cuentaEstados/ length(sec.est)
  
  ## iniciales de la variable E.
  prueba$pi[1] <-cuentaEstados[1]/ length(sec.est)
  prueba$pi[2]<-cuentaEstados[1]/ length(sec.est)
  pi <- prueba$pi
  
  ## Probabilidades de transición
  ## Paso 4: Calcular el número de veces que aparece cada posible
  ## pareja de estados en la secuencia de estados.
  cuentaEstados2 <- count(sec.est, 2, alphabet = alf.e)
  
  ## Paso 5: Definir A como la matriz mXm adecuada que contiene
  ## las frecuencias anteriores (nombrar adecuadamente sus filas y
  ## columnas).
  A = matrix(cuentaEstados2, nrow = m, byrow = T, 
             dimnames = list(alf.e, alf.e))
  
  ## Paso 6: Normalizar cada fila de la matriz A.
  for (i in 1:m) {
    A[i,] <- A[i,] / sum(A[i,])
  }
  
  ## Probabilidades de emisión
  ## Paso 7: Definir B como una matriz mXn (nombrar adecuadamente
  ## sus filas y columnas).
  B = matrix(nrow = m, ncol = n, 
             dimnames = list(alf.e, alf.x))
  
  ## Paso 8: Para cada elemento e de alf.e
  ##   Paso 8.1: Almacenar en la fila e de B el número de veces que
  ##   aparece cada símbolo, alf.x, en las posiciones de la
  ##   secuencia de símbolos que se corresponde con las apariciones
  ##   de e en la secuencia de estados.
  for (i in alf.e) {
    B[i,] <- count(sec.simb[sec.est == i], 1, alphabet = alf.x)
    # Paso 8.2: Normalizar la fila e de la matriz B
    B[i,] <- B[i,] / sum(B[i,])
  }         
  
  ## Devolver una lista con pi, A y B.
  list(pi = pi, A = A, B = B)
}

dna.seq <- c("a", "c", "t", "a", "a", "t", "g", "t", "c", "t", "g", "a", 
             "c", "t", "a", "t", "c", "a", "g", "c", "g", "g", "c", "a", 
             "g", "c", "a", "g", "c", "t", "g", "a", "c", "c", "g", "t")


hidden.seq <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
                1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2,
                2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2)

estimate.hmm(dna.seq, hidden.seq)
## $pi
##   1   2 
## 0.5 0.5 
## $A
##           1          2
## 1 0.9444444 0.05555556
## 2 0.0000000 1.00000000
## $B
##           a         c         g         t
## 1 0.3333333 0.2222222 0.1111111 0.3333333
## 2 0.1666667 0.3333333 0.3888889 0.1111111

######################################################################
## La siguiente función estima la sucesión de estados más probable
## dada una secuencia símbolos y el hmm correspondiente.

viterbiLog <- function (obs, mom) {
  l <- length(obs)
  estados <- names(mom$pi)
  pr <- matrix(nrow = l-1, ncol = length(mom$pi),
               dimnames = list(2:l, estados))
  B <- log(mom$B)
  A <- log(mom$A)
  pi <- log(mom$pi)
  nu <- B[,obs[1]] + pi 
  for (k in 2:l) {
    aux <- nu + A
    nu <- B[,obs[k]] + apply(aux, 2, max)
    pr[as.character(k),] <- apply(aux, 2, which.max)
  }
  q <- c()
  q[l] <- estados[which.max(nu)]
  for (k in (l-1):1) {
    q[k] <- estados[pr[as.character(k+1),q[k+1]]]
  }
  q
}

######################################################################
## (4) Identificar utilizando una cadena de Markov oculta las regiones
## con alto y bajo contenido en GC. Utiliza lambda.initial.hmm como
## modelo (se encuentra definido en el archivo HMM.RData).
######################################################################

secuenciaMasLogica <- viterbiLog (lambda.sequences, lambda.initial.hmm)
table(secuenciaMasLogica)
# Tenemos dos estados, el estado bajo contenido en GC, con 16089 elementos 
# y el otro estado h2 de alto contenido en GC

######################################################################
## (5) Representa en un panel dos gráficos: uno con el gráfico
## anterior con el contenido local en GC y otro gráfico que represente
## la estimación de la secuencia oculta. Comenta estos gráficos.
######################################################################

par(mfrow=c(2,1), oma = c(0,0,0,2))
plot(lambdaGC.local[["positions"]],lambdaGC.local[["results"]],
     type="l",pch=0,col="red",cex=1.5,lwd=1,
     xlab="Genome Posicion",ylab="GC result", cex.lab=1.5,
     cex.axis=1.25,
     main="Local GC result in lambda fago")
lines(lambdaGC.local[["positions"]],rep(labmdaGC.global,length(lambdaGC.local[["positions"]])),
      type="l",pch=1,col="blue",cex=1.5,lwd=2)

# Definimos el vector x como una secuencia de la misma longitud que secuencia mas logica
x <- seq(1,length(secuenciaMasLogica))

# Calculamos el vector y que estará formado por -1 y 1, correspondiente al
# estado h1 y h2 respectivamente
y <- c()
for (i in x){
  if (secuenciaMasLogica[i] == "h1") y <- c(y,-1) else y <- c(y,1)
}

plot(x,y,
     type="l",pch=0,col="green",cex=1.5,lwd=1,
     xlab="Genome Posicion",ylab="State", cex.lab=1.5,
     cex.axis=1.25,
     main="Estimation of the hidden sequence")
par(xpd=NA)
legend(x=45000, y=50,legend=c("Local GC","Global GC","EHS"),col=c("red","blue","green"),
       lwd=2,cex=0.75)

# Se observa como las bajadas y subidas son similares en ambas gráficas.
# Cuando el contenido local en GC es alto, la estimación de la secuencia
# oculta es 1, lo que equivale, según la asignación que hemos hecho al
# estado h2. Sin embargo, cuando hay bajadas del contenido en GC en la
# primera gráfica vemos que el estado de la secuencia oculta obtenida
# de la cadena de Markov es -1, es decir, estado h1.
# Por lo tanto, podemos ver una clara relación entre el contenido local en GC
# y la estimación de la cadena oculta de estados.
dev.off()

######################################################################
## Tarea:

## En general se desconocen exactamente los parámetros pi, A y B que
## definen un HMM. Sin embargo, dada una secuencia de símbolos y una
## aproximación inicial pi0, A0 y B0 se puede seguir un proceso
## iterativo de esperanza-maximización o algoritmo EM
## (Expectation-Maximization Algorithm, EM Algorithm) para obtener un
## modelo más verosímil.

## Este tipo de algoritmos se organizan en dos fases diferentes:

## Entrada: Una secuencia de símbolos observable y un primer HMM
## asociado a la misma.

## * Fase 1, Paso M o Paso de máxima verosimilitud:
##   Dadas una secuencia de símbolos observable y un HMM calcular la
##   secuencia de estados oculta de máxima verosimilitud.

## * Fase 2, Paso E o Paso de esperanza:
##   Dadas una secuencia de símbolos observable y una secuencia oculta
##   estimar el HMM más probable.

## Estas dos fases se repiten un número de iteraciones produciendo en
## cada vuelta del bucle una mejor aproximación del modelo

## Salida: La estimación del HMM de la última iteración y la secuencia
## oculta más probable.

## (a) Implementar el algoritmo EM descrito. Debe devolver el modelo
## HMM estimado y la secuencia de estados oculta más probable

algoritmoEM <- function(seq.simb, HMM_0, k){
  hmm.est <- HMM_0
  for (i in 1:k) {
    # Guardamos la estimación para hacer la condición de parada
    hmm.est.ant <- hmm.est
    
    # Cálculo secuencia máxima verosimilitud
    seq.est.ma <- viterbiLog(seq.simb,  hmm.est)
    
    # Recalculamos las matrices de probabilidad
    hmm.est <- estimate.hmm(seq.simb, seq.est.ma)
    
    # Parada si la estimación obtenida es igual a la que ya teníamos
    if (hmm.est.ant$pi == hmm.est$pi && hmm.est.ant$A == hmm.est$A
        && hmm.est.ant$B == hmm.est$B) break
  }
  # Devolvemos el HMM estimado y la secuencia de estados oculta más probable:
  list(hmm.est = hmm.est, seq.est.ma = seq.est.ma,k=k)
}

# ESTO ESTA MAL
# seq.simb = olfac_14C36.sequences
# hmm.est = hidro.initial.hmm
# mom = hmm.est
# obs=seq.simb

## (b) Mejorar el modelo proporcionado en el apartado 4 y comparar los
## resultados del apartado 5 con los que se obtengan con el nuevo
## modelo.

## (c) Estudiar la hidrofobicidad local y global de proteínas
## transmembranales.

## Las proteínas transmembranales son aquellas que se encuentran
## insertadas en las membranas atravesando una o varias veces la
## bicapa lipídica. De esta forma en las secuencias de aminoácidos que
## forman estas proteínas se pueden distinguir diferentes
## segmentos. Por una parte, los dominios insertados en la bicapa
## lipídica son altamente hidrofóbicos. Mientras que las regiones
## citosólicas y extracitosólicas son altamente hidrofílicas. Las
## siguientes funciones nos permitirán detectar de forma aproximada
## los dominios hidrofílicos e hidrofóbicos de proteínas
## transmembranales.


## (c.1) Descargar de la web de Uniprot las secuencias de la proteínas
## transmembranales correspondientes a los receptores olfativos del
## Homo sapiens: olfactory receptor 14C36, olfactory receptor 2J3 y
## olfactory receptor 1F1. Lee estas secuencias en formato fasta en R
## siguiendo la convención de usar letras mayúsculas para representar
## aminoácidos.

olfac_14C36 <- read.fasta("14C36.fasta",seqtype="AA")
olfac_14C36.sequences <- getSequence(olfac_14C36)[[1]]

olfac_2J3 <- read.fasta("2J3.fasta",seqtype="AA")
olfac_2J3.sequences <- getSequence(olfac_2J3)[[1]]

olfac_1F1 <- read.fasta("1F1.fasta",seqtype="AA")
olfac_1F1.sequences <- getSequence(olfac_1F1)[[1]]

## El vector hydrophobic.values (definido en el archivo HMM.RData)
## contiene los valores hidrofóbicos de cada aminoácido.

hydrophobic.values

## (c.2) Implementar una función, hydrophobicity, que reciba como
## entrada una secuencia de aminoácidos y devuelva su hidrofobicidad
## global. Esta hidrofobicidad se calcula como la suma de la
## hidrofibicidad de cada aminoácido dividida por la longitud de la
## secuencia de aminoácidos.

alfabet <- sort(unique(olfac_14C36.sequences))

hydrophobicity <- function (seq.ami) {
  # Calculamos el número de veces que aparece cada aminoácido:
  cuenta <- count(seq.ami, 1, alphabet = alfabet)
  
  sum(hydrophobic.values * cuenta) / length(seq.ami)
}

## (c.3) Calcular la hidrofobicidad local y global de los receptores
## olfativos del apartado (c.1). Usa una ventana de longitud 20 con un
## desplazamiento de 1 para el calculo de la hidrofobicidad
## local. Representar en tres gráficos de líneas diferentes la
## hidrofobicidad global y local de cada proteína. ¿Cuántas regiones
## hidrofóbicas se distinguen en cada proteína? ¿Cuántas veces
## atraviesan estas proteínas la capa lipídica?

# Global:
hidro.global.14C36 <- hydrophobicity(olfac_14C36.sequences)
hidro.global.2J3 <- hydrophobicity(olfac_2J3.sequences)
hidro.global.1F1 <- hydrophobicity(olfac_1F1.sequences)

# Local:
hidro.local.14C36 <- local.composition(olfac_14C36.sequences, hydrophobicity, 20,1)
hidro.local.2J3 <- local.composition(olfac_2J3.sequences, hydrophobicity, 20,1)
hidro.local.1F1 <- local.composition(olfac_1F1.sequences, hydrophobicity, 20,1)

# Representacion:
plot(hidro.local.14C36[["positions"]],hidro.local.14C36[["results"]],
     type="l",pch=0,col="red",cex=1.5,lwd=1,
     xlab="Position",ylab="Hydrophobicity", cex.lab=1.5,
     cex.axis=1.25,
     main="Local hydrophobicity in olfactory receptor 14C36")
lines(hidro.local.14C36[["positions"]],rep(hidro.global.14C36,length(hidro.local.14C36[["positions"]])),
      type="l",pch=1,col="blue",cex=1.5,lwd=2)
legend("topright",legend=c("Local hydrophobicity","Global hydrophobicity"),
       col=c("red","blue"),lwd=2,cex=0.75)

# Las proteinas se encuentran atravesando la bicapa lipidica, asi podemos
# diferenciar diferentes segmentos en las secuencias de aminoacidos. En
# el receptor 14C36 vemos 6 regiones altamente hidrofobicas y 6 cruces de la capa lipidica


plot(hidro.local.2J3[["positions"]],hidro.local.2J3[["results"]],
     type="l",pch=0,col="red",cex=1.5,lwd=1,
     xlab="Position",ylab="Hydrophobicity", cex.lab=1.5,
     cex.axis=1.25,
     main="Local hydrophobicity in olfactory receptor 2J3")
lines(hidro.local.2J3[["positions"]],rep(hidro.global.2J3,length(hidro.local.2J3[["positions"]])),
      type="l",pch=1,col="blue",cex=1.5,lwd=2)

legend("topright",legend=c("Local hydrophobicity","Global hydrophobicity"),
       col=c("red","blue"),lwd=2,cex=0.75)

# En el receptor 2J3 encontramos con 4 regiones altamente hidrofobicas 
# y que atraviesa 5 veces la capa lipídica.


plot(hidro.local.1F1[["positions"]],hidro.local.1F1[["results"]],
     type="l",pch=0,col="red",cex=1.5,lwd=1,
     xlab="Position",ylab="Hydrophobicity", cex.lab=1.5,
     cex.axis=1.25,
     main="Local hydrophobicity in olfactory receptor 1F1")
lines(hidro.local.1F1[["positions"]],rep(hidro.global.1F1,length(hidro.local.1F1[["positions"]])),
      type="l",pch=1,col="blue",cex=1.5,lwd=2)

legend("topright",legend=c("Local hydrophobicity","Global hydrophobicity"),
       col=c("red","blue"),lwd=2,cex=0.75)

# En el receptor 1F1 tenemos 6 regiones altamente hidrofobicas


## (c.4) Usa la función del apartado (a) para realizar un aprendizaje
## no supervisado con los tres receptores olfativos tomando como
## modelo inicial hidro.initial.hmm y 10 iteraciones.

olfac_14C36.em <- algoritmoEM(olfac_14C36.sequences,hidro.initial.hmm, 10)
olfac_2J3.em <- algoritmoEM(olfac_2J3.sequences,hidro.initial.hmm, 10)
olfac_1F1.em <- algoritmoEM(olfac_1F1.sequences,hidro.initial.hmm, 10)

## (c.5) Muestra gráficamente los resultados obtenidos. Comenta los
## resultados. ¿Qué propones para mejorar la eficiencia de la cadena
## de Markov oculta utilizada?


# olfac_14C36.em
par(mfrow=c(2,1), oma = c(0,0,0,2))
plot(hidro.local.14C36[["positions"]],hidro.local.14C36[["results"]],
     type="l",pch=0,col="red",cex=1.5,lwd=1,
     xlab="Position",ylab="Hydrophobicity", cex.lab=1.5,
     cex.axis=1.25,
     main="Local hydrophobicity in olfactory receptor 14C36")
lines(hidro.local.14C36[["positions"]],rep(hidro.global.14C36,length(hidro.local.14C36[["positions"]])),
      type="l",pch=1,col="blue",cex=1.5,lwd=2)

legend("topright",legend=c("Local hydrophobicity","Global hydrophobicity"),
       col=c("red","blue"),lwd=2,cex=0.75)


# Definimos el vector x como una secuencia de la misma longitud que secuencia mas logica
x1 <- seq(1,length(olfac_14C36.em$seq.est.ma))

# Calculamos el vector y que estará formado por -1 y 1, correspondiente al
# estado h1 y h2 respectivamente
y1 <- c()
for (i in x1){
  if (olfac_14C36.em$seq.est.ma[i] == 1) y1 <- c(y1,-1) else y1 <- c(y1,1)
}

plot(x1,y1,
     type="l",pch=0,col="green",cex=1.5,lwd=1,
     xlab="Genome Posicion",ylab="State", cex.lab=1.5,
     cex.axis=1.25,
     main="Estimation of the hidden sequence")
par(xpd=NA)


# olfac_1F1.em
plot(hidro.local.1F1[["positions"]],hidro.local.1F1[["results"]],
     type="l",pch=0,col="red",cex=1.5,lwd=1,
     xlab="Position",ylab="Hydrophobicity", cex.lab=1.5,
     cex.axis=1.25,
     main="Local hydrophobicity in olfactory receptor 1F1")
lines(hidro.local.1F1[["positions"]],rep(hidro.global.1F1,length(hidro.local.1F1[["positions"]])),
      type="l",pch=1,col="blue",cex=1.5,lwd=2)

legend("topright",legend=c("Local hydrophobicity","Global hydrophobicity"),
       col=c("red","blue"),lwd=2,cex=0.75)


# Definimos el vector x como una secuencia de la misma longitud que secuencia mas logica
x2 <- seq(1,length(olfac_1F1.em$seq.est.ma))

# Calculamos el vector y que estará formado por -1 y 1, correspondiente al
# estado h1 y h2 respectivamente
y2 <- c()
for (i in x1){
  if (olfac_1F1.em$seq.est.ma[i] == 1) y2 <- c(y2,-1) else y2 <- c(y2,1)
}

plot(x2,y2,
     type="l",pch=0,col="green",cex=1.5,lwd=1,
     xlab="Genome Posicion",ylab="State", cex.lab=1.5,
     cex.axis=1.25,
     main="Estimation of the hidden sequence")
par(xpd=NA)


# olfac_2J3.em
plot(hidro.local.2J3[["positions"]],hidro.local.2J3[["results"]],
     type="l",pch=0,col="red",cex=1.5,lwd=1,
     xlab="Position",ylab="Hydrophobicity", cex.lab=1.5,
     cex.axis=1.25,
     main="Local hydrophobicity in olfactory receptor 2J3")
lines(hidro.local.2J3[["positions"]],rep(hidro.global.2J3,length(hidro.local.2J3[["positions"]])),
      type="l",pch=1,col="blue",cex=1.5,lwd=2)

legend("topright",legend=c("Local hydrophobicity","Global hydrophobicity"),
       col=c("red","blue"),lwd=2,cex=0.75)


# Definimos el vector x como una secuencia de la misma longitud que secuencia mas logica
x3 <- seq(1,length(olfac_2J3.em$seq.est.ma))

# Calculamos el vector y que estará formado por -1 y 1, correspondiente al
# estado h1 y h2 respectivamente
y3 <- c()
for (i in x1){
  if (olfac_2J3.em$seq.est.ma[i] == 1) y3 <- c(y3,-1) else y3 <- c(y3,1)
}

plot(x3,y3,
     type="l",pch=0,col="green",cex=1.5,lwd=1,
     xlab="Genome Posicion",ylab="State", cex.lab=1.5,
     cex.axis=1.25,
     main="Estimation of the hidden sequence")
par(xpd=NA)
dev.off()

# Para mejorar la estimacion propongo aumentar k iteraciones
