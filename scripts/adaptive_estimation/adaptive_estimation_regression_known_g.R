# ############################################################
# ############################################################
### Programa : Estimacion adaptativa de la funcion de regresion,
### Usando nucleo gaussiano, con g conocida.
# ############################################################
# ############################################################

n <- 2200 # Longitud de la muestra que se generara.
sig1 <- 0.5 # Desviacion Estandar del Error de la regresion.

# ##########################################################
# ##########################################################
# ########### INICIO CALIBRACION #######################
# ##########################################################
# ##########################################################

# ##########################################################
### Generacion de numeros aleatorios de distribucion N(0,1)
### truncada por 1_{(-a,a)}(x) dependientes usando datos
### AR(1) con phi=0.75 y el metodo de la transformada
### inversa.
# ###########################################################

Z = arima.sim(list(ar = 0.75), n = n) # Generacion de AR(1) con
# phi=0.75 y sigma=1

U = pnorm(Z, mean = 0, sd = sqrt(1/(1-0.75^2))) 
# U=F(Z), donde Z~N(0, 1/(1-0.75^2))
# y F es su funcion de distribucion.
# U~U(0,1) y es dependiente.

a = 2; # Soporte para la normal truncada, dependiente.
p = pnorm(a) - pnorm(-a); # Area bajo la normal (0,1) truncada
# en a, tal densidad se denota por g
# y ademas g=G'.

X = qnorm(p*U+pnorm(-a), mean = 0, sd = 1) 
# X=G^(-1)(U)=(G^(-1)oF)(Z).
# Normal truncada dependiente.

plot(seq(1, n, 1), X, xlab="", ylab="") # Grafico nube de puntos
# de datos generados,
# variable explicativa.

acf(X, lag=100) # Grafico de la funcion de correlacion de la
# variable explicativa

T1 = X
Y0t = 0.7*T1+2*exp(-10*T1^2) # Funcion de regresion a estimar
# r(x)=0.7x+2exp(-10x^2).

Yt = Y0t+rnorm(n, 0, sig1) # Variable respuesta

ne <- n-200 # Numero de datos con los que se construye el
# estimador

X = T1[1:(n-200)] # Datos de la variable explicativa con los que
# se construye el estimador

Y0 = Y0t[1:(n-200)] # Y0=r(X), donde X es normal truncada
# dependiente.

Y = Yt[1:(n-200)] # Modelo Y=r(X)+ e donde e~N(0,1), datos de
# la variable respuesta con los que se
# construye el estimador

# ##########################################################
# ####### Estimacion desviacion estandar de e_j ##########
# ####### donde e_j es el error de la regresion ##########
# ##########################################################
library(KernSmooth)
h <- dpill(X,Y) # Metodo plug-in de Ruppert, Sheather y
# Wand (1995), esta ventana obtenida por el
# metodo del pulgar la sustituire en el
# estimador de N-W, para asi tener una
# estimacion de r en cada punto de X en R^2000.

# ###########################################################
# ############## MODULO 1 m ########################
# ###########################################################
Estimmm <- function(x, X, Y, h){
  # Estimador de la funcion m con kernel gaussiano, ventana
  # h, muestra aleatoria de la variable explicativa X,
  # muestra aleatoria de la variable respuesta Y, y variable
  # independiente x.
  # m(x)=sum_{i=1}^n(1/(nh*sqrt(2*pi))) Y[i] exp(-(((x-
  # X[i])/h)^2)/2).
  n <- length(X);
  mm = 0;
  for(i in 1:n){
    mm = mm+(1/(n*h*sqrt(2*pi)))*Y[i]*exp(-(((x-X[i])/h)^2)/2);
  }
  return(mm)
}

# ############################################################
# ################ MODULO 1 g ########################
# ############################################################
Estimg <- function(x, X, h){
  # Estimador de la densidad g con kernel gausiano, ventana h
  # para la muestra aleatoria X envaluada en x.
  # g(x)=sum_{i=1}^n(1/(nh sqrt(2*pi))) exp(-(((x-X[i])/h)^2)/2).
  n <- length(X);
  g = 0;
  for(i in 1:n){
    g = g+(1/(n*h*sqrt(2*pi)))*exp(-(((x-X[i])/h)^2)/2);
  }
  return(g)
}

# ############################################################
# ############# Estimador mpul, gpul y rpul #################
# ############################################################
mpul <- 0; gpul <- 0;
for(i in 1:ne){
  mpul[i] <- Estimmm(X[i],X, Y, h); # Numerador del estimador de
  # Nadaraya-W
  
  gpul[i] <- Estimg(X[i],X, h); # Denominador del estimador de
  # Nadaraya-W
}
rpul <- mpul/gpul # Estimador de Nadaraya-W, en la ventana h,
# estimada por el metodo del pulgar.
# #############################################################

sigqchu <- (1/(ne-1))*sum((Y-rpul)^2) # Estimacion Var(e_j)
sigchu <- sqrt(sigqchu) # Estimacion desviacion estandar de e_j
sig <- sigchu # Desviacion estandar de e_j (Estimacion
# del error de la regresion)

# ############################################################
# ####### Fin Estimacion desviacion estandar e_j ############
# ############################################################

t1 = T1[((n-100)+1):n]
t2 = Yt[((n-100)+1):n]
ii = order(t1)
x = t1[ii]
Yct <- t2[ii]

# ############################################################
# #### Formula de la densidad g, normal (0,1) truncada ########
# ############################################################
mu <- 0; sigm <- 1; # Extremo del soporte es a=2.
dTnorm <- function(x, mu, sigm, a){
  # Densidad N(mu, sigm) truncada de -a hasta a.
  p <- pnorm(mu+a, mu, sigm, lower.tail=TRUE, log.p=FALSE) -
    pnorm(mu-a, mu, sigm, lower.tail=TRUE, log.p=FALSE);
  if(mu-a<=x && x<=mu+a){
    g <- dnorm(x, mu, sigm, log=FALSE)/p;
  }else{g=0;}
  return(g)
}
# ###########################################################

# ###########################################################
# ################ MODULO 1 ########################
# ###########################################################
EstimR <- function(x, X, Y, h){
  # Estimador de la regresion con kernel gausiano, ventana h
  # para la muestra aleatoria (X,Y) envaluada en x.
  # r(x)=sum_{i=1}^n(1/(nh sqrt(2*pi))) Y[i] exp(-(((x-
  # X[i])/h)^2)/2)g^{-1}(X[i]).
  n <- length(X);
  r = 0;
  for(i in 1:n){
    r = r+(1/(n*h*sqrt(2*pi)))*Y[i]*exp(-(((x-X[i])/h)^2)/2) *dTnorm(X[i], mu, sigm, a)^(-1);
  }
  return(r)
}
# ###########################################################

# ###########################################################
# ####################### MODULO 2 ######################
# ###########################################################
EstimRSS <- function(x, X, Y, h, h1){
  # Estimador sobresuavisado con nucleo gausiano, muestra
  # aleatoria simple (X,Y), ventanas h y h1, evaluado en x.
  n <- length(X);
  rss = 0
  for(i in 1:n){
    rss = rss+(1/(n*sqrt(2*pi*(h^2+h1^2))))*Y[i]*exp(-((x-X[i])^2)/(2*(h^2+h1^2)))*dTnorm(X[i], mu, sigm, a)^(-1);
  }
  return(rss)
}
# ###########################################################

# ############################################################
# ####################### MODULO 3 #####################
# ############################################################
# Esta funcion determina el valor absoluto de la diferencia
# entre el estimador de la regresion sobresuavizado rhh1 y rh1
# en un valor x del soporte.
VAdifEst <- function(x, X, Y, h, h1){
  VA <- abs(EstimRSS(x, X, Y, h, h1) - EstimR(x, X, Y, h1))
  return(VA)
}
# ############################################################

ne <- length(X); # Longitud de la muestra de entrenamiento
I <- 0
xt <- 0
Yc <- 0
for(i in 1:100){
  if(-1<=x[i] & x[i]<=1){
    I <- I+1;
    xt[I] <- x[i];
    Yc[I] <- Yct[i]
  }
}
x <- xt

s <- length(x) # Numero de elementos en la particion del soporte
# de la densidad de X

M <- floor(log(ne))*(2/3) # Definicion de M, donde M es la parte
# entera n/log(n),

H <- exp(-seq(0, M, 0.1)) # Familia de ventanas, donde H={h_i}_0^M,
# con M definido en la linea anterior

m <- length(H) # Numero de ventanas.

# ########################################################
# ########################################################
#### Rutina para almacenar los datos
#### usando la funcion VAdifEst
#### en una lista de nombre l que
#### contiene s matrices donde cada
#### matriz es de orden (mxm)
#### s=length(x) y m=length(H).
#### El almacenamiento en cada matriz
#### se realiza por filas.
# ########################################################

l <- list() # Nombre de la lista que almacenara las matrices
# donde cada una es de orden mxm.

# Mi comentario: Esta parte demora algo (fuera bueno medir cuanto y porque tarda)
for(k in 1:s){
  l0 <- 0 # Variable auxiliar que almacenara temporalmente
  # cada matriz de datos.
  l00 <- 0 # Variable auxiliar que almacena temporalmente
  # cada fila de la matriz l0.
  for(j in 1:m){
    l00 <- VAdifEst(x[k],X, Y, H, H[j])
    if(j==1){l0<-l00}else{l0<-cbind(l0, l00)}
  }
  l[[k]] <- l0
}
# ###########################################################
# ###########################################################

# ###########################################################
# ####### Cambio Local de rsup y ginf IMPORTANTE ###
# ####### Ahora rsup y ginf son vectores de dimension s ###
# ####### Cada rsup[i] y ginf[i] corresponde con x[i], ###
# ####### para i=1:s ###
# ###########################################################
d <- rep(0, s) #d distancias para el calculo de Int(rchu-Y)^2
d[51] <- ((x[51]+x[52])/2)-(-1)
d[s] <- 1-((x[s-1]+x[s])/2)
for(i in 2:(s-1)){
  d[i] <- ((x[i]+x[i+1])/2)-((x[i-1]+x[i])/2)
}

e <- rep(0, s+1)
se <- length(e)
e[51] <- (-1)
e[se] <- 1

for(i in 2:(se-1)){
  e[i] <- (x[i-1]+x[i])/2
}

rsup <- rep(0, s)
ginf <- rep(0, s)

for(i in 1:s){
  rsup[i] <- max(abs(Y*(e[i]<X)*(X<e[i+1]))) # Maximo Local de
  # la funcion de
  # regresion.
  W1 <- 0
  W1 <- X*(e[i]<X)*(X<e[i+1])
  W2 <- 0
  I1 <- 1
  for(j in 1:ne){
    if(W1[j]!=0){W2[I1]<-X[j]; I1<-I1+1}
  }
  valorg <- rep(0, length(W2))
  for(j in 1:length(W2)){
    valorg[j] <- dTnorm(W2[j], 0, 1, a)
  }
  ginf[i] <- min(valorg) #Minimo de la funcion de densidad de X.
}
# ############### FIN IMPORTANTE, CAMBIO LOCAL ##############
# ############### FIN IMPORTANTE, CAMBIO LOCAL ##############

del <- (log(ne))^(-1/5)
N2 <- 1/sqrt(2*sqrt(pi)) #Norma 2 del nucleo gaussiano.

# ###########################################################
# ##### Cambio Local de A3 IMPORTANTE ###
# ##### Ahora A3 es un vector de dimension s ###
# ##### Cada A3[i] corresponde con x[i], para cada i=1:s. ###
# ###########################################################
A3 <- rep(0, s)
for(i in 1:s){
  A3[i] <- (rsup[i]^2+sig^2)*ginf[i]^(-1)*N2^2
}
# ############## FIN IMPORTANTE, CAMBIO LOCAL ###############
# ############## FIN IMPORTANTE, CAMBIO LOCAL ###############

GAM1 <- 0.00000005
GAM21 <- 0.05
dd <- (GAM21-GAM1)/20
GAM <- seq(GAM1, GAM21, dd) # Valores de la variable gamma,
# usados para calibrar el metodo
Norma2Cuadrado <- rep(0, length(GAM));

for(q in 1:length(GAM)){
  print(q)
  gam <- GAM[q] # Parametro gamma > 2. (valor que he cambiado)
  B4 <- sqrt(2*gam*A3)*(2) # El 2 que multiplica despues de la
  # raiz cuadrada es 2 = \|K\|_1 +1,
  # donde \|K\|_1=1
  # "IMPORTANTE ahora B4 es una variable LOCAL"
  # " Es claro que B4 es un vector de dimension s".
  # Cada B4[i] corresponde a x[i], para cada i=1:s
  # IMPORTANTE en esta parte no hubo cambios.
  
  # ########################################################
  # ##### Calculo del estimador V(x,h), para cada x en ##
  # ##### la malla de longitud s y h \in H. ##
  # ##### Es decir V es una matriz de orden sxm. ##
  # ########################################################
  # ##### IMPORTANTE aca se hacen cambios locales, ##
  # ##### Ahora V tambien depende de x, ademas de ##
  # ##### depender de H (como antes) ##
  # ##### sera una matriz de oden sxm, donde cada fila ##
  # ##### corresponde a un valor de x y cada columna a ##
  # ##### los valores de H. V[i,j] donde i=1:s y j=1:m ##
  # ########################################################
  V <- matrix(0, nrow=s, ncol=m)
  for(k in 1:m){
    V[,k] <- B4*(1+del)*sqrt(log(ne))/sqrt(ne*H[k]) # Aparece ne, pues es la longitud
    # de la muestra de entrenamiento
  }
  # ######### Fin del cambio local de V ##################
  # #######################################################
  
  # #######################################################
  # #### Calculo del Maximo de cada fila ##############
  # #### de la matriz {l[k]-V}_{+}, ##############
  # #### para cada k de 1:s. ##############
  # #### Lo cual se guarda en la lista ##############
  # #### A={A[k]}_{k=1;s} con A[k] en R^m ##############
  # #### Se genera la lista AV, compuesta ##############
  # #### de s vectores m dimensionales ##############
  # #### AV[[k]]=A[[k]]+V en R^m para ##############
  # #### k = 1, ..., s. ##############
  # #### Se geneta el vector hopt en R^s ##############
  # #### donde hopt[k]=argmin(AV[[k]]+V) ##############
  # #######################################################
  A <- list() # Es una lista de s vectores, donde cada
  # vector A[k] en R^m, corresponde a un
  # valor x[k] en de la malla, y A[k][i]
  # es el valor en H[i]
  
  AV <- list() # Es una lista de vectores, donde cada
  # vector AV[k] coresponde a la suma de
  # vectores A[k]+V[k,].
  
  arg <- 0 # Esta Variable sera un vector donde
  # cada componente es
  # arg[k]=argmin_{h \in H}(A[k]+V[k,])
  
  hopt <- 0 # Variable que almacena las ventanas optimas
  # para cada x[k] en la malla.
  
  for(k in 1:s){
    A0 <- 0 # Esta variable auxiliar, sera una matriz que
    # almacenara temporalmente los datos contenidos
    # en cada matriz A[[k]] de la lista A
    
    for(i in 1:m){
      A0[i] <- max(l[[k]][i,]-V[k,])
      if(A0[i]<0){A0[i]<-0}
    }
    A[[k]] <- A0 # Se almacenan los datos en la lista A
    AV[[k]] <- A0+V[k,] # Se almacenan los datos en la lista AV
    arg[k] <- which.min(AV[[k]]) # Se almacena la posicion de la
    # ventana optima para cada x[k]
    # de la malla de x.
  }
  
  hopt <- H[arg] # Vector de ventanas optimas,
  # es un vector en R^s
  
  # ########################################################
  # ################## Estimador rchu #####################
  # ########################################################
  rchu <- 0;
  TMSE0 <- 0; # Variable auxiliar que almacenara temporalmente
  # los terminos del MSE de la estimacion.
  for(i in 1:s){
    rchu[i] <- EstimR(x[i],X, Y, hopt[i]);
    TMSE0[i] <- (rchu[i]-Yc[i])^2;
  }
  # #########################################################
  
  # #########################################################
  # ##### Grafico del estimador rchu para cada GAM[q], ######
  # ##### con q=1:\length(GAM) ######
  # #########################################################
  plot(x, rchu, ylim=c(-1.5, 2.5))
  lines(x, rchu, lwd=2, col="red")
  points(X, Y0, col="black") # Curva Y0=r(X)
  
  r <- EstimR(x, X, Y, 1/6); # Estimador de nucleo en una ventana
  # fija MODULO 1
  
  lines(x, r, lwd=2, col='blue')
  points(x, Yc, lwd=2, col='green')
  # ##########################################################
  
  Norma2Cuadrado[q] <- sum(d*TMSE0)
}

plot(Norma2Cuadrado, main="Calibración Método GL", xlab="", ylab="")
lines(Norma2Cuadrado)

i0 <- which.min(Norma2Cuadrado)
i0
GAM[i0]
# ###########################################################
# ###########################################################
# ################### Final Calibracion ################
# ###########################################################
# ###########################################################

# ###########################################################
# ###########################################################
# ############### INICIO AJUSTE ##################
# ###########################################################
# ###########################################################
a <- 2
a1 <- 1 # Limite superior del soporte de la densidad de X
x1 <- seq(-a1, a1, 0.1) # Soporte de la densidad de X
s1 <- length(x1) # Numero de elementos en la particion del
# soporte de la densidad de X

# ###########################################################
# ###########################################################
# ######### Rutina para almacenar los datos ##############
# ######### usando la funcion VAdifEst ##############
# ######### en una lista de nombre l que ##############
# ######### contiene s matrices donde cada ##############
# ######### matriz es de orden (mxm) ##############
# ######### s=length(x) y m=length(H). ##############
# ######### El alamcenamiento en cada matriz ##############
# ######### se realiza por filas. ##############
# ###########################################################
l <- 0 # Nombre de la lista que almacenara s matrices
# donde cada una es de orden mxm.

# Mi comentario: Esta parte también demora
for(k in 1:s1){
  l0 <- 0 # Variable auxiliar que almacenara temporalmente
  # cada matriz de datos.
  l00 <- 0 # Variable auxiliar que almacena temporalmente
  # cada fila de la matriz l0.
  
  for(j in 1:m){
    l00 <- VAdifEst(x1[k],X, Y, H, H[j])
    if(j==1){l0<-l00}else{l0<-cbind(l0, l00)}
  }
  if(k==1){l<-list(l0)}else{l[[k]]<-l0}
}

# ###########################################################
# ###########################################################
# #### Cambio Local de rsup y ginf IMPORTANTE ######
# #### Ahora rsup y ginf son vectores de dimension s ######
# #### Cada rsup[i] y ginf[i] corresponde con x[i], ######
# #### para i=1:s ######
# ###########################################################
# ###########################################################
rsup <- rep(0, s1)
ginf <- rep(0, s1)

for(i in 1:s1){
  rsup[i] <- max(abs(Y*((x1[i]-0.05)<X)*(X<(x1[i]+0.05)))) # Maximo Local de
  # la funcion de regresion. #####
  
  W1 <- 0
  W1 <- X*((x1[i]-0.05)<X)*(X<(x1[i]+0.05))
  W2 <- 0
  I1 <- 1
  for(j in 1:ne){
    if(W1[j]!=0){W2[I1]<-X[j]; I1<-I1+1}
  }
  valorg <- rep(0, length(W2))
  for(j in 1:length(W2)){
    valorg[j] <- dTnorm(W2[j], 0, 1, a)
  }
  ginf[i] <- min(valorg) #Minimo de la funcion de densidad de X.
}
# ############## FIN IMPORTANTE, CAMBIO LOCAL ##############
# ############## FIN IMPORTANTE, CAMBIO LOCAL ##############
# ##########################################################

gam <- GAM[i0] # Parametro gamma
del <- (log(ne))^(-1/5)
N2 <- 1/sqrt(2*sqrt(pi)) # Norma 2 del nucleo gaussiano.

# ##########################################################
# ######### Cambio Local de A3 IMPORTANTE ######
# ######### Ahora A3 es un vector de dimension s ######
# ######### Cada A3[i] corresponde con x[i], ######
# ######### para cada i=1:s. ######
# ##########################################################
A3 <- rep(0, s1)
for(i in 1:s1){
  A3[i] <- (rsup[i]^2+sig^2)*ginf[i]^(-1)*N2^2
}
# ############### FIN IMPORTANTE, CAMBIO LOCAL #############
# ############### FIN IMPORTANTE, CAMBIO LOCAL #############
# ##########################################################

B4 <- sqrt(2*gam*A3)*(2) # El 2 que multiplica despues de
# la raiz cuadrada es
# 2 = \|K\|_1 +1, donde \|K\|_1=1

# ##########################################################
# ##### Calculo del estimador V(h), para cada h \in H. #####
# ##### IMPORTANTE aca se hacen cambios locales, #####
# ##### Ahora V tambien depende de x, ademas de h #####
# ##### V sera una matriz de oden sxm, donde cada #####
# ##### fila corresponde a un valor de x y cada #####
# ##### columna a los valores de H. #####
# ##### V[i,j] donde i=1:s y j=1:m #####
# ##########################################################
V <- matrix(0, nrow=s1, ncol=m)
for(k in 1:m){
  V[,k] <- B4*(1+del)*sqrt(log(ne))/sqrt(ne*H[k])
}
# ##########################################################
# ############## Fin del cambio local de V. ################
# ##########################################################

# ##########################################################
# #### Calculo del Maximo de cada fila #################
# #### de la matriz {l[k]-V}_{+}, #################
# #### para cada k de 1:s. #################
# #### Lo cual se guarda en la lista #################
# #### A={A[k]}_{k=1;s} con A[k] en R^m #################
# #### Se genera la lista AV, compuesta #################
# #### de s vectores m dimencionales #################
# #### AV[[k]]=A[[k]]+V en R^m para #################
# #### k = 1, ..., s. #################
# #### Se geneta el vector hopt en R^s #################
# #### donde hopt[k]=argmin(AV[[k]]+V) #################
# ##########################################################
A <- list() # Es una lista de s vectores, donde cada
# vector A[k] en R^m, corresponde a un valor
# x[k] de la malla, y A[k][i] es el
# valor en H[i]

AV <- list() # Es una lista de s, ectores, donde cada vector
# AV[k] coresponde a la suma de vectores A[k]+V.

arg <- 0 # Esta Variable sera un vector donde cada
# componente es arg[k]=argmin_{h \in H}(A[k]+V)

hopt <- 0 # Variable que almacena las ventanas optimas
# para cada x[k] en la malla.

for(k in 1:s1){
  A0 <- 0 # Esta variable auxiliar, sera una matriz que
  # almacenara temporalmente los datos contenidos
  # en cada matris A[[k]] de la lista A
  
  for(i in 1:m){
    A0[i] <- max(l[[k]][i,]-V[k,])
    if(A0[i]<0){A0[i]<-0}
  }
  # Se genera la lista A y AV
  A[[k]] <- A0; # Se guardan resultados en la lista A
  AV[[k]] <- A0+V[k,] # Se guardan resultados en la lista AV
  arg[k] <- which.min(AV[[k]]) # Se almacena la ventana optima
  # para cada x[k] de la malla de x
}
hopt <- H[arg]
# ###########################################################

# ############################################################
# ###################### Estimador rchu #####################
# ############################################################
rchu1 <- 0;
for(i in 1:s1){
  rchu1[i] <- EstimR(x1[i],X, Y, hopt[i]);
}
# ############################################################

# ############################################################
# ############################################################
# #################### FINALIZO EL AJUSTE ####################
# ############################################################
# ############################################################

# ##############################################################
# ##################### GRAFICOS FINALES ######################
# ##############################################################
plot(x1, rchu1, ylim=c(-2, 3), type="n", main="Estimacion Adaptativa Método GL para r", xlab="", ylab="")
points(X, Y, col="gray")
lines(x1, rchu1, lwd=2, lty="dotted")

#### Funcion de regresion a estimar r(x)=0.7x+2exp(-10x^2)
#### evaluada en la malla x1 del intervalo [-1,1]
r=0.7*x1+2*exp(-10*x1^2)
lines(x1, r, lwd=2, col='black')
