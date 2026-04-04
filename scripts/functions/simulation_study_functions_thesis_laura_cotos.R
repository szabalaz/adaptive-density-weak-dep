# Selector ventana de Silverman
h.Silverman<-function(x, K) {
  n<-length(x)
  sd<-sqrt(var(x))
  ri<-diff(quantile(x, probs=c(0.25, 0.75)))/(qnorm(0.75)-qnorm(0.25))
  # el rango intercuartilico muestral tambien se puede obtener con IQR (x)
  sigma.hat<-min(sd,ri)
  
  if (K=="Gauss"){
    mu_2<- 1
    Rk<-1/(2*sqrt(pi))
  } else if (K=="Epa"){
    mu_2<-1/5
    Rk<-3/5
  }
  
  mu_2<-mu_2^2
  h.Silverman<-((8*sqrt(pi)*Rk)/(3*mu_2*n))^(1/5)*sigma.hat
  return (h.Silverman)
}

#Selector de ventana de Sheather y Jones
h.SJ<-function (x,K) {
  n<-length(x)
  if (K=="Gauss"){
    K6<-function (x) {exp((-x^2)/2)*(x^6-15*x^4+45*x^2-15)/(sqrt(2*pi))}
    K4<-function (x) {exp((-x^2)/2)*(x^4-6*x^2+3)/(sqrt(2*pi))}
    mu_2<-1
    Rk<-1/(2*sqrt(pi))
  } else if (K=="Epa"){
    K6<-function (x) {0}
    K4<-function (x) {0} # la derivada segunda ya es una constante
    mu_2<-1/5
    Rk<-3/5 #mu_2 y Rk sacados de la tabla B.2 de Kernel Smoothing
  }
  
  #1º estimamos fi8 utilizando la regla del pulgar
  sd<-sd(x) #sqrt(var(x))
  ri<-diff(quantile(x, probs=c(0.25, 0.75)))/(qnorm(0.75)-qnorm(0.25))
  sigma.hat<-min(sd,ri)
  fi8<-105/(32*sqrt(pi)*sigma.hat^9)
  
  #ahora calculamos la ventana óptima de fi6 (que es la que depende de fi8)
  g1<-((2*K6(0))/(-mu_2*fi8*n))^(1/9)
  
  #ahora puedo calcular la estimación de fi6 con la ventana óptima g1
  a<-outer(x,x,"-")/g1 #matriz nxn con (Xi-Xj)/g1
  Keval<-matrix(K6(a), nc=n, nrow=n) #evalúo la matriz en la derivada correspondiente de K
  fi6<-(1/(n^2*g1^7))*sum(Keval)
  
  #calculamos ahora la ventana optima g2 de fi4
  g2<-((2*K4(0))/(-mu_2*fi6*n))^(1/7)
  
  #calculamos fi4
  a<-outer(x,x,"-")/g2 #matriz nxn con (Xi-Xj)/g2
  Keval2<-matrix(K4(a), nc=n,nrow=n) #evalúo la matriz en la derivada correspondiente de K
  fi4<-(1/(n^2*g2^5))*sum(Keval2)
  
  h.SJ<-(Rk/((mu_2)^2*fi4*n))^(1/5)
  return(h.SJ)
}

#Estimador tipo núcleo
f.kernel<-function (x1,x0,h,K) {
  n<-length(x1)
  a<-outer(x0,x1, FUN="-")/h # me da una matriz dxn con d la longitud de x0
  if (K=="Gauss"){Kf=K_gauss} else if (K=="Epa"){Kf=K_epa}
  Keval<-matrix(Kf(a), nc=n, nrow=length(x0)) # funcion nucleo de a
  yy<-rowSums(Keval)/(n*h) # promedio los nucleos por filas
  return(data.frame('xrej'=x0,'y'=yy))
}

K_gauss <- function(x) {1/sqrt(2*pi)*exp(-x^2/2)} #con -infty<x<infty

# Estudio de simulación
library(nor1mix)
n<-50 # n=50,100,200,500
M<-1000 #numero de muestras con cada tamaño muestral
n.rej<-100 #longitud de la rejilla

#creamos unas matrices para guardar los datos de las ventanas y los errores
mat.err <- matrix(NA, ncol=4, nrow=M)
colnames(mat.err) <- c("Teorica", "Silverman", "SJ", "CV")

mat.e1 <- matrix(NA, ncol=4, nrow=1)
colnames(mat.e1) <- c("Teorica", "Silverman", "SJ", "CV")

mat.e2 <- matrix(NA, ncol=4, nrow=1)
colnames(mat.e2) <- c("Teorica", "Silverman", "SJ", "CV")

mat.h <- matrix(NA, ncol=4, nrow=M)
colnames(mat.h) <- c("Teorica", "Silverman", "SJ", "CV")

mat.dif <- matrix(NA, ncol=3, nrow=M)
colnames(mat.dif) <- c("Teorica-Silverman", "Teorica-SJ", "Teorica-CV")

mat.e3 <- matrix(NA, ncol=3, nrow=1)
colnames(mat.e3) <- c("Silverman", "SJ", "CV")

mat.e4 <- matrix(NA, ncol=3, nrow=1)
colnames(mat.e4) <- c("Silverman", "SJ", "CV")

#hacemos una función con la ventana teórica dada por el ECMIA ya que esta
#solo va a depender del tamaño muestral y del modelo teórico
h.teorica<- function (n, modelo) {
  #estoy suponiendo que K es el núcleo gaussiano
  Rk<-1/(2*sqrt(pi))
  mu2 <-1
  
  if (modelo=="MW.nm1"){ #N(0,1)
    Rf<-3/(8*sqrt(pi)) #sacado de Wand & Jones
    h <- (Rk/(n*mu2^2*Rf))^0.2
  } else {
    if (modelo=="MW.nm4") {
      f4 <- function (x) {(2/(3*sqrt(2*pi)) *exp(-x^2/2) *exp(-50*x^2) * (-10^3+10^5*x^2))^2}
      Rf <- integrate(f4,-Inf, Inf)
      h <- (Rk/(n*mu2^2*Rf$value))^0.2
    } else {
      if (modelo=="MW.nm8") {
        f8<-function(x) {(0.75* 1/sqrt(2*pi) *exp(-x^2/2)*(-1+x^2) + 0.25*1/sqrt(2*pi) * exp(-9/2*(x-1.5)^2) * (-27+243*(x-1.5)^2))^2}
        Rf <- integrate(f8,-Inf, Inf)
        h <- (Rk/(n*mu2^2*Rf$value))^0.2
      } else{
        if (modelo=="MW.nm10") {
          f10 <- function (x) {
            mu<-c(-1,-0.5,0,0.5,1)
            sig <- 0.1
            res <- 0.5 * (x^2-1) * 1/sqrt(2*pi) * exp(-x^2/2) +
              0.1/sig^3 * (((x-mu[1])/sig)^2-1) * 1/sqrt(2*pi) * exp(-0.5*((x-mu[1])/sig)^2) +
              0.1/sig^3 * (((x-mu[2])/sig)^2-1) * 1/sqrt(2*pi) * exp(-0.5*((x-mu[2])/sig)^2) +
              0.1/sig^3 * (((x-mu[3])/sig)^2-1) * 1/sqrt(2*pi) * exp(-0.5*((x-mu[3])/sig)^2) +
              0.1/sig^3 * (((x-mu[4])/sig)^2-1) * 1/sqrt(2*pi) * exp(-0.5*((x-mu[4])/sig)^2) +
              0.1/sig^3 * (((x-mu[5])/sig)^2-1) * 1/sqrt(2*pi) * exp(-0.5*((x-mu[5])/sig)^2)
            return (res^2)
          }
          Rf <- integrate(f10,-Inf, Inf)
          h<- (Rk/(n*mu2^2*Rf$value))^0.2
        }
      }
    }
  }
  return (h)
}

############################ MODELO 1
#calculamos entonces la densidad teórica:
h.teor<-h.teorica (n, modelo="MW.nm1")
rejilla <- seq(qnorMix(0.01, MW.nm1), qnorMix(0.99, MW.nm1), len=n.rej)
d.teor <- dnorMix(rejilla, obj=nor1mix::MW.nm1) #densidad

set.seed (12345) #fijamos la semilla

for (m in 1:M) {
  #generamos una muestra
  x <- nor1mix::rnorMix(n, obj=nor1mix::MW.nm1)
  
  # Aplicamos cada uno de los selectores:
  mat.h[m, 1] <- h.teor
  mat.h[m, 2] <- h.Silverman(x, K="Gauss")
  mat.h[m, 3] <- h.SJ(x, K="Gauss")
  mat.h[m, 4] <- bw.bcv(x)
  
  # Calculamos la estimacion con cada uno de los selectores
  d.teor.est <- f.kernel(x,x0=rejilla, h=mat.h[m, 1], K="Gauss")
  d.silverman <- f.kernel(x,x0=rejilla, h=mat.h[m, 2], K="Gauss")
  d.sj <- f.kernel(x,x0=rejilla, h=mat.h[m, 3], K="Gauss")
  d.cv <- f.kernel(x,x0=rejilla, h=mat.h[m, 4], K="Gauss")
  
  # Calculamos medida de error ECI (de cada muestra)
  mat.err[m, 1] <- mean((d.teor.est[,2]-d.teor)^2)
  mat.err[m, 2] <- mean((d.silverman[,2]-d.teor)^2)
  mat.err[m, 3] <- mean((d.sj[,2]-d.teor)^2)
  mat.err[m, 4] <- mean((d.cv[,2]-d.teor)^2)
  
  # Calculamos el error de las ventanas (de cada muestra)
  mat.dif[m, 1] <- mat.h[m,1]-mat.h[m,2]
  mat.dif[m, 2] <- mat.h[m,1]-mat.h[m,3]
  mat.dif[m, 3] <- mat.h[m,1]-mat.h[m,4]
}

# Calculamos la media (ECMI) y sd de los errores (de las M muestras)
mat.e1[1] <- mean(mat.err[,1])
mat.e1[2] <- mean(mat.err[,2])
mat.e1[3] <- mean(mat.err[,3])
mat.e1[4] <- mean(mat.err[,4])

mat.e2[1] <- sd(mat.err[,1])
mat.e2[2] <- sd(mat.err[,2])
mat.e2[3] <- sd(mat.err[,3])
mat.e2[4] <- sd(mat.err[,4])

# Calculamos media y sd de los errores (diferencias de ventanas)
mat.e3[1] <- mean(mat.dif[,1]) # media de h.silverman-h.teorica
mat.e3[2] <- mean(mat.dif[,2]) # media de h.sj-h.teorica
mat.e3[3] <- mean(mat.dif[,3]) 

mat.e4[1] <- sd(mat.dif[,1])
mat.e4[2] <- sd(mat.dif[,2])
mat.e4[3] <- sd(mat.dif[,3])

#boxplot de los errores (ECI)
boxplot(mat.err*100)

######################## MODELO 4
#calculamos entonces la densidad teórica:
h.teor<-h.teorica(n, modelo="MW.nm4")
rejilla <- seq(qnorMix(0.01, MW.nm4), qnorMix(0.99, MW.nm4), len=n.rej)
d.teor <- dnorMix(rejilla, obj=nor1mix::MW.nm4)

set.seed (12345) #fijamos la semilla

for (m in 1:M) {
  # generamos una muestra
  x <- nor1mix::rnorMix(n, obj=nor1mix::MW.nm4)
  
  # Aplicamos cada uno de los selectores:
  mat.h[m, 1] <- h.teor
  mat.h[m, 2] <- h.Silverman(x, K="Gauss")
  mat.h[m, 3] <- h.SJ(x, K="Gauss")
  mat.h[m, 4] <- bw.bcv(x)
  
  # Calculamos la estimacion con cada uno de los selectores
  d.teor.est <- f.kernel(x,x0=rejilla, h=mat.h[m, 1], K="Gauss")
  d.silverman <- f.kernel(x,x0=rejilla, h=mat.h[m, 2], K="Gauss")
  d.sj <- f.kernel(x,x0=rejilla, h=mat.h[m, 3], K="Gauss")
  d.cv <- f.kernel(x,x0=rejilla, h=mat.h[m, 4], K="Gauss")
  
  # Calculamos medida de error ECI (de cada muestra)
  mat.err[m, 1] <- mean((d.teor.est[,2]-d.teor)^2)
  mat.err[m, 2] <- mean((d.silverman[,2]-d.teor)^2)
  mat.err[m, 3] <- mean((d.sj[,2]-d.teor)^2)
  mat.err[m, 4] <- mean((d.cv[,2]-d.teor)^2)
  
  # Calculamos tambien el error de las ventanas (de cada muestra)
  mat.dif[m, 1] <- mat.h[m,1]-mat.h[m,2]
  mat.dif[m, 2] <- mat.h[m,1]-mat.h[m,3]
  mat.dif[m, 3] <- mat.h[m,1]-mat.h[m,4]
}

# Calculamos la media (ECMI) y sd de los errores (de las M muestras)
mat.e1[1] <- mean(mat.err[, 1])
mat.e1[2] <- mean(mat.err[, 2])
mat.e1[3] <- mean(mat.err[, 3])
mat.e1[4] <- mean(mat.err[, 4])

mat.e2[1] <- sd(mat.err[, 1])
mat.e2[2] <- sd(mat.err[, 2])
mat.e2[3] <- sd(mat.err[, 3])
mat.e2[4] <- sd(mat.err[, 4])

# Calculamos media y sd de los errores (diferencias de ventanas)
mat.e3[1] <- mean(mat.dif[,1]) # media de h.silverman-h.teorica
mat.e3[2] <- mean(mat.dif[,2]) # media de h.sj-h.teorica
mat.e3[3] <- mean(mat.dif[,3])

mat.e4[1] <- sd(mat.dif[,1])
mat.e4[2] <- sd(mat.dif[,2])
mat.e4[3] <- sd(mat.dif[,3])

#boxplot de los errores (ECI)
boxplot(mat.err*100)

######################## MODELO 8
#calculamos entonces la densidad teórica:
h.teor<-h.teorica(n,modelo="MW.nm8")
rejilla <- seq(qnorMix(0.01, MW.nm8), qnorMix(0.99, MW.nm8), len=n.rej)
d.teor <- dnorMix(rejilla, obj=nor1mix::MW.nm8)

set.seed (12345) #fijamos la semilla

for (m in 1:M) {
  # generamos una muestra
  x <- nor1mix::rnorMix(n, obj=nor1mix::MW.nm8)
  
  # Aplicamos cada uno de los selectores:
  mat.h[m, 1] <- h.teor
  mat.h[m, 2] <- h.Silverman(x, K="Gauss")
  mat.h[m, 3] <- h.SJ(x, K="Gauss")
  mat.h[m, 4] <- bw.bcv(x)
  
  # Calculamos la estimacion con cada uno de los selectores
  d.teor.est <- f.kernel(x,x0=rejilla, h=mat.h[m, 1], K="Gauss")
  d.silverman <- f.kernel(x,x0=rejilla, h=mat.h[m, 2], K="Gauss")
  d.sj <- f.kernel(x,x0=rejilla, h=mat.h[m, 3], K="Gauss")
  d.cv <- f.kernel(x,x0=rejilla, h=mat.h[m, 4], K="Gauss")
  
  # Calculamos medida de error ECI (de cada muestra)
  mat.err[m, 1] <- mean((d.teor.est[,2]-d.teor)^2)
  mat.err[m, 2] <- mean((d.silverman[,2]-d.teor)^2)
  mat.err[m, 3] <- mean((d.sj[,2]-d.teor)^2)
  mat.err[m, 4] <- mean((d.cv[,2]-d.teor)^2)
  
  # Calculamos tambien el error de las ventanas (de cada muestra)
  mat.dif[m, 1] <- mat.h[m,1]-mat.h[m,2]
  mat.dif[m, 2] <- mat.h[m,1]-mat.h[m,3]
  mat.dif[m, 3] <- mat.h[m,1]-mat.h[m,4]
}

# Calculamos la media (ECMI) y sd de los errores (de las M muestras)
mat.e1[1] <- mean(mat.err[,1])
mat.e1[2] <- mean(mat.err[,2])
mat.e1[3] <- mean(mat.err[,3])
mat.e1[4] <- mean(mat.err[,4])

mat.e2[1] <- sd(mat.err[,1])
mat.e2[2] <- sd(mat.err[,2])
mat.e2[3] <- sd(mat.err[,3])
mat.e2[4] <- sd(mat.err[,4])

# Calculamos media y sd de los errores (diferencias de ventanas)
mat.e3[1] <- mean(mat.dif[,1]) # media de h.silverman-h.teorica
mat.e3[2] <- mean(mat.dif[,2]) # media de h.sj-h.teorica
mat.e3[3] <- mean(mat.dif[,3])

mat.e4[1] <- sd(mat.dif[,1])
mat.e4[2] <- sd(mat.dif[,2])
mat.e4[3] <- sd(mat.dif[,3])

#boxplot de los errores (ECI)
boxplot(mat.err*100)

######################## MODELO 10
#calculamos la densidad teórica:
h.teor<-h.teorica(n, modelo="MW.nm10")
rejilla <- seq(qnorMix(0.01,MW.nm10), qnorMix(0.99, MW.nm10),len=n.rej)
d.teor <- dnorMix(rejilla, obj=nor1mix::MW.nm10)

set.seed (12345) #fijamos la semilla

for (m in 1:M){
  # generamos una muestra
  x <- nor1mix::rnorMix(n, obj=nor1mix::MW.nm10)
  
  # Aplicamos cada uno de los selectores:
  mat.h[m, 1] <- h.teor
  mat.h[m, 2] <- h.Silverman(x, K="Gauss")
  mat.h[m, 3] <- h.SJ(x, K="Gauss")
  mat.h[m, 4] <- bw.bcv(x)
  
  # Calculamos la estimacion con cada uno de los selectores
  d.teor.est <- f.kernel(x,x0=rejilla, h=mat.h[m, 1], K="Gauss")
  d.silverman <- f.kernel(x,x0=rejilla, h=mat.h[m, 2], K="Gauss")
  d.sj <- f.kernel(x,x0=rejilla, h=mat.h[m, 3], K="Gauss")
  d.cv <- f.kernel(x,x0=rejilla, h=mat.h[m, 4], K="Gauss")
  
  # Calculamos medida de error ECI (de cada muestra)
  mat.err[m, 1] <- mean((d.teor.est[,2]-d.teor)^2)
  mat.err[m, 2] <- mean((d.silverman[,2]-d.teor)^2)
  mat.err[m, 3] <- mean((d.sj[,2]-d.teor)^2)
  mat.err[m, 4] <- mean((d.cv[,2]-d.teor)^2)
  
  # Calculamos tambien el error de las ventanas (de cada muestra)
  mat.dif[m, 1] <- mat.h[m,1]-mat.h[m,2]
  mat.dif[m, 2] <- mat.h[m,1]-mat.h[m,3]
  mat.dif[m, 3] <- mat.h[m,1]-mat.h[m,4]
}

# Calculamos la media (ECMI) y sd de los errores (de las M muestras)
mat.e1[1] <- mean(mat.err[,1])
mat.e1[2] <- mean(mat.err[,2])
mat.e1[3] <- mean(mat.err[,3])
mat.e1[4] <- mean(mat.err[,4])

mat.e2[1] <- sd(mat.err[,1])
mat.e2[2] <- sd(mat.err[,2])
mat.e2[3] <- sd(mat.err[,3])
mat.e2[4] <- sd(mat.err[,4])

# Calculamos media y sd de los errores (diferencias de ventanas)
mat.e3[1] <- mean(mat.dif[, 1]) # media de h.silverman-h.teorica
mat.e3[2] <- mean(mat.dif[, 2]) # media de h.sj-h.teorica
mat.e3[3] <- mean(mat.dif[, 3])

mat.e4[1] <- sd(mat.dif[, 1])
mat.e4[2] <- sd(mat.dif[, 2])
mat.e4[3] <- sd(mat.dif[, 3])

#boxplot de los errores (ECI)
boxplot(mat.err*100)