# ########################################################### 
# ###### Estimador de r con densidad g desconocida. ####### 
# ###### Calibración de gchu usando                 ####### 
# ###### ||gchu-HHg||_2^2                           ####### 
# ###### donde gchu es la estimacion de la          ####### 
# ###### funcion de densidad y HHg un histograma    ####### 
# ###### V(h) distinto al articulo                  ####### 
# ###### de Bertin y Klutchnikoff                   ####### 
# ###########################################################

n <- 2200 # Numero de datos a generar, para la muestra
sig1 <- 0.5 # Desviacion Estandar del Error de la regresion.

# ########################################################### 
# ############### CALIBRACION Y ESTIMACION g ################ 
# ########################################################### 
for (HH1 in 1:1) {
  print(HH1) 
  # #########################################################
  # ######## Generacion de numeros aleatorios       ######### 
  # ######## de distribucion N(0,1) truncada por    ######### 
  # ######## 1_{-a,a}(x) dependientes usando        ######### 
  # ######## datos AR(1) con phi=0.75 y el          ######### 
  # ######## metodo de la transformada inversa.     ######### 
  # ######################################################### 
  Z = arima.sim(list(ar=0.75), n=n) # Generacion de AR(1) con
  # phi=0.75 y sigma=1
  
  U = pnorm(Z, mean=0, sd=sqrt(1/(1-0.75^2)))
  # U=F(Z), donde Z~N(0, 1/(1-0.75^2))
  # y F es su funcion de distribucion.
  # U~U(0,1) y es dependiente.
  
  a = 2; # Soporte para la normal truncada, dependiente.
  p = pnorm(a) - pnorm(-a); # Area bajo la normal (0,1) truncada
  # en a, tal densidad se denota por g
  # y ademas g=G'.
  
  X = qnorm(p*U+pnorm(-a), mean=0, sd=1) # X=G^(-1)(U)=(G^(-1)oF)(Z).
  # Normal truncada dependiente.
  
  plot(seq(1, n, 1), X, xlab="", ylab="") # Diagrama de dispersion de los
  # datos generados {X_i}_{i=1}^n
  
  acf(X, lag=100) # Grafico de la funcion de correlacion
  # de los primeros 100 datos X_i
  
  T1 = X
  
  ### Funcion de regresion a estimar r(x)=0.7x+2exp(-10x^2).
  Y0t = 0.7*T1+2*exp(-10*T1^2)
  
  Yt = Y0t+rnorm(n, 0, sig1) # Generacion de la variable
  # respuesta del modelo de regresion
  
  ne <- n-200 # Numero de datos para estimar
  
  X = T1[1:(n-200)] # Datos de la variable explicativa para estimar
  
  Y0 = Y0t[1:(n-200)] # Y0=r(X), donde X es normal truncada dependiente.
  
  Y = Yt[1:(n-200)] # Modelo Y=r(X)+e, donde e~N(0,1)
  # o datos de la variable respuesta para estimar
  
  # ######################################################### 
  # ######################################################### 
  # ######## Estimacion desviacion estandar e_j     ######### 
  # ######################################################### 
  # #########################################################
  
  library(KernSmooth) 
  h <- dpill(X,Y) # Metodo plug-in de Ruppert, Sheather y Wand (1995)
  # Esta ventana obtenida por el metodo del pulgar,
  # la sustituire en el estimador de N-W para asi
  # tener una pre-estimacion de r en cada punto de
  # X en R^2000. Lo que permitira estimar e_j
  
  # ##########################################################
  
  # ######################################################## 
  # ############## MODULO 1 m          ################### 
  # ############## Estimador de m      ################### 
  # ######################################################## 
  Estimmm <- function(x, X, Y, h) {
    # Estimador de m con kernel gausiano, ventana h para la
    # muestra aleatoria X envaluada en x.
    # m(x)=sum_{i=1}^n (1/(nh*sqrt(2*pi))) Y[i] exp(-(((x-X[i])/h)^2)/2).
    n <- length(X);
    mm = 0;
    for(i in 1:n) {
      mm = mm + (1/(n*h*sqrt(2*pi))) * Y[i] * exp(-(((x-X[i])/h)^2)/2);
    }
    return(mm)
  }
  
  # ######################################################### 
  # ############# MODULO 1 g           ############## 
  # ############# Estimador de g       ############## 
  # ######################################################### 
  Estimg <- function(x, X, h) {
    # Estimador de g con kernel gausiano, ventana h para la
    # muestra aleatoria X envaluada en x.
    # g(x)=sum_{i=1}^n (1/(nh sqrt(2*pi))) exp(-(((x-X[i])/h)^2)/2).
    n <- length(X);
    g = 0;
    for(i in 1:n) {
      g = g + (1/(n*h*sqrt(2*pi))) * exp(-(((x-X[i])/h)^2)/2);
    }
    return(g)
  }
  
  # ######################################################### 
  # ############# Estimador mpul, gpul y rpul ############## 
  # ######################################################### 
  mpul <- 0; gpul <- 0; 
  for(i in 1:ne) {
    mpul[i] <- Estimmm(X[i], X, Y, h); 
    gpul[i] <- Estimg(X[i], X, h);
  }
  rpul <- mpul/gpul
  # ########################################################
  
  sigqchu <- (1/(ne-1))*sum((Y-rpul)^2) # Estimacion Var(e_j) 
  sigchu <- sqrt(sigqchu) # Estimacion desviacion estandar de e_j 
  sig <- sigchu # Desviacion estandar de e_j
  
  # ######################################################### 
  # ######################################################### 
  # ###### Fin Estimacion desviacion estandar e_j    ######## 
  # ######################################################### 
  # #########################################################
  
  t1 = T1[((n-100)+1):n] # Ultimos 100 datos de la variable explicativa para calibrar.
  t2 = Yt[((n-100)+1):n] # Ultimos 100 datos de la variable respuesta para calibrar.
  
  ii = order(t1) # Comando para ordenar
  
  x = t1[ii] # Ultimos 100 datos de la variable explicativa para calibrar, "ORDENADOS"
  Yct <- t2[ii] # Ultimos 100 datos de la variable respuesta para calibrar, "ORDENADOS"
  
  # ########################################################## 
  ### Formula de la densidad g, normal (0,1) truncada ##### 
  # ########################################################## 
  mu <- 0; sigm <- 1; # Extremo del soporte es a=2. 
  dTnorm <- function(x, mu, sigm, a) {
    # Densidad N(mu, sigm) truncada de -a hasta a.
    p <- pnorm(mu+a, mu, sigm, lower.tail=TRUE, log.p=FALSE) - 
      pnorm(mu-a, mu, sigm, lower.tail=TRUE, log.p=FALSE);
    if(mu-a<=x && x<=mu+a) {
      g <- dnorm(x, mu, sigm, log=FALSE)/p;
    } else {g=0;}
    return(g)
  }
  #**ACA PARA QUE CORRA**#
  # ###########################################################
  
  # ########################################################### 
  # ################# MODULO 2 g                #################### 
  # ########### Estimador de g sobresuavizado   ########## 
  # ########################################################### 
  EstimgSS <- function(x, X, h, h1) {
    # Estimador sobresuavizado con nucleo gausiano, muestra
    # aleatoria simple X, ventanas h y h1.
    n <- length(X);
    gss = 0
    for(i in 1:n) {
      gss = gss + (1/(n*sqrt(2*pi*(h^2+h1^2)))) * exp(-((x-X[i])^2)/(2*(h^2+h1^2)));
    }
    return(gss)
  }
  #**ACA PARA QUE CORRA**#
  # ###########################################################
  
  # ########################################################### 
  # ##################### MODULO 3 g ###################### 
  # ########################################################### 
  # Esta funcion determina el valor absoluto de la diferencia 
  # entre el estimador de la regresion sobresuavizado rhh1 y rh1 
  # en un valor x del soporte. 
  VAdifEstg <- function(x, X, h, h1) {
    VA <- abs(EstimgSS(x, X, h, h1) - Estimg(x, X, h1))
    return(VA)
  }
  # ########################################################### 
  # ###########################################################
  
  histg <- hist(X, freq=FALSE) # Calculo y grafico del histograma
  # usado para calibrar la estimacion de la densidad
  
  s <- length(x) # Numero de elementos en la particion del
  # soporte de la densidad de X
  
  go <- 0 # Variable que almacena la funcion de densidad original
  for(i in 1:s) {
    go[i] <- dTnorm(x[i], mu, sigm, a)
  }
  lines(x, go, lwd=2, col='black') # Al histograma se le agrega la
  # curva de la densidad original
  
  #eh <- histg[[4]]; EstHist <- histg[[5]];  # SEGUN HAY QUE COMENTAR ESTO
  eh <- histg$breaks
  EstHist <- histg$density
  nh <- length(EstHist) 
  gsup <- max(EstHist) # Se determina el valor maximo del histohrama
  
  # ######################################################### 
  # ################ MODULO 4 Histograma      ############### 
  # ################ Funcion Histograma       ############### 
  # ######################################################### 
  Histograma <- function(x, nh, EstHist, eh) {
    Hg <- 0
    for(i in 1:(nh-1)) {
      if(eh[i]<=x && x<eh[i+1]) {Hg <- EstHist[i]}
    }
    if(eh[nh]<=x && x<=eh[nh+1]) {Hg <- EstHist[nh]}
    return(Hg)
  }
  # ######################################################### 
  # #########################################################
  
  ne <- length(X); # Numero de datos en la muestra de
  # estimacion (Para construir el Estimador)
  
  # a1<-1 # Limite superior del soporte de la densidad de X 
  # x<- seq(-a1, a1, 0.1) # Soporte de la densidad de X
  
  M <- floor(log(ne))*(2/3) # Definicion de M, donde M es la parte entera n/log(n),
  
  H <- exp(-seq(0, M, 0.1)) # Familia de ventanas, donde
  # H={h_i}_0^M, con M definido en Linea anterior.
  
  m <- length(H) # Numero de ventanas.
  
  # ########################################################## 
  # ########################################################## 
  #### Rutina para almacenar los datos        ################## 
  #### usando la funcion VAdifEst             ################## 
  #### en una lista de nombre l que           ################## 
  #### contiene s matrices donde cada         ################## 
  #### matriz es de orden (mxm)               ################## 
  #### s=length(x) y m=length(H).             ################## 
  #### El alamcenamiento en cada matriz       ################## 
  #### se realiza por filas.                  ################## 
  # ########################################################## 
  l <- list() # Nombre de la lista que almacenara s matrices
  # donde cada una es de orden mxm.
  
  for(k in 1:s) {
    l0 <- 0 # Variable auxiliar que almacena temporalmente
    # cada matriz de datos.
    l00 <- 0 # Variable auxiliar que almacena temporalmente
    # cada fila de la matriz l0.
    
    for(j in 1:m) {
      l00 <- VAdifEstg(x[k], X, H, H[j])
      if(j==1) {l0 <- l00} else {l0 <- cbind(l0, l00)}
    }
    l[[k]] <- l0
  }
  # ######################################################### 
  # #########################################################
  
  # ######################################################### 
  #### Cambio Local de rsup y ginf IMPORTANTE #### 
  #### Ahora rsup y ginf son vectores de dimension s #### 
  #### Cada rsup[i] y ginf[i] corresponde con x[i], #### 
  #### para i=1:s #### 
  # ######################################################### 
  d <- rep(0, s) #d distancias para el calculo de Int(rchu-Y)^2 
  d[4] <- ((x[4]+x[6])/2)-(-2) 
  d[s] <- 2-((x[s-1]+x[s])/2) 
  for(i in 2:(s-1)) {
    d[i] <- ((x[i]+x[i+1])/2)-((x[i-1]+x[i])/2)
  }
  
  delg <- (log(ne))^(-1/2);
  N2 <- 1/sqrt(2*sqrt(pi)); #Norma 2 del nucleo gaussiano. 
  N1 <- 1; #Norma 1 del nucleo gasussiano
  
  GAM1 <- 0.00000005 # Extremo inferior de la malla para calibrar 
  GAM21 <- 0.08 # Extremo superior de la malla para calibrar 
  dd <- (GAM21-GAM1)/20 # Numero de elementos de la malla 
  GAM <- seq(GAM1, GAM21, dd) # Malla para calibrar el estimador
  # de la funcion de densidad g
  
  Norma2gchuHHg <- rep(0, length(GAM)); # Variable que alamacena
  # la funcion Error(GAM)
  
  GCHU <- matrix(0, nrow=s, ncol=length(GAM)) # Matriz que almacena
  # las estimaciones gchu para cada valor GAM[i]
  
  # ######################################################### 
  # ######################################################### 
  # ############ INICIO CALIBRACION ############### 
  # ######################################################### 
  # ######################################################### 
  for(q in 1:length(GAM)) {
    print(q) 
    gam <- GAM[q] # Parametro gamma > 2. (valor que he cambiado)
    
    # ####################################################### 
    ### Calculo del estimador V(x,h), para cada x en ### 
    ### la malla de longitud s y h \in H.            ### 
    ### Es decir V es una matriz de orden sxm        ### 
    # ####################################################### 
    ### IMPORTANTE aca se hacen cambios locales,     ### 
    ### Ahora V tambien depende de x, ademas de      ### 
    ### depender de H(como antes)                    ### 
    ### V sera una matriz de oden sxm, donde cada    ### 
    ### fila corresponde a un valor de x y cada      ### 
    ### columna a los valores de H.                  ### 
    ### V[i,j] donde i=1:s y j=1:m                   ### 
    # ####################################################### 
    V <- matrix(0, nrow=s, ncol=m) 
    for(i in 1:s) {
      for(k in 1:m) {
        V[i,k] <- sqrt(2*gsup*gam)*N2*(1+N1)*(1+delg) * ((log(ne))/(ne*H[k]))^(1/2)
      }
    }
    #**ACA PARA QUE CORRA**#
    # ######################################################## 
    ### Fin del cambio local de V. ### 
    # ########################################################
    
    # ######################################################## 
    # ######################################################## 
    # #### Calculo del Maximo de cada fila ################# 
    # #### de la matriz {l[k]-V}_{+},      ################# 
    # #### para cada k de 1:s.             ################# 
    # #### Lo cual se guarda en la lista   ################# 
    # #### A={A[k]}_{k=1;s} con A[k] en R^m. ################# 
    # #### Se genera la lista AV, compuesta ################# 
    # #### de s vectores m dimencionales    ################# 
    # #### AV[[k]]=A[[k]]+V en R^m para     ################# 
    # #### k = 1, ..., s.                   ################# 
    # #### Se geneta el vector hopt en R^s  ################# 
    # #### donde hopt[k]=argmin(AV[[k]]+V)  ################# 
    # ######################################################## 
    # ######################################################## 
    A <- list() # Es una lista de s vectores, donde cada
    # vector A[k] en R^m, corresponde a un valor
    # x[k] en de la malla, y A[k][i] es el valor en H[i]
    
    AV <- list() # Es una lista de s, vectores, donde cada
    # vector AV[k] coresponde a la suma de vectores A[k]+V[k,].
    
    arg <- 0 # Esta Variable sera un vector donde cada
    # componente es arg[k]=argmin_{h \in H}(A[k]+V[k,])
    
    hoptg <- 0 # Variable que almacena las ventanas optimas
    # para cada x[k] en la malla.
    
    ### IMPORTANTE, En la parte que sigue se hicieron #### 
    ### CAMBIOS LOCALES, correspondiente a la variable V. ### 
    for(k in 1:s) {
      A0 <- 0 # Esta variable auxiliar, sera una matriz que
      # almacenara temporalmente los datos contenidos
      # en cada matriz A[[k]] de la lista A
      
      for(i in 1:m) {
        A0[i] <- max(l[[k]][i,] - V[k,])
        if(A0[i]<0) {A0[i] <- 0}
      }
      A[[k]] <- A0 # Se almacenan los datos en la lista A 
      AV[[k]] <- A0+V[k,] # Se almacenan los datos en la lista AV
      
      arg[k] <- which.min(AV[[k]]) # Se almacena la posicion de la ventana optima
      # para cada x[k] de la malla de x.
    }
    # FIN CAMBIOS LOCALES, correspondientes a la variable V. # 
    # FIN CAMBIOS LOCALES, correspondientes a la variable V. #
    
    hoptg <- H[arg] # Vector de ventanas optimas, es un vector en R^s
    
    # ######################################################### 
    # ######################################################### 
    # ################ Estimador gchu ################## 
    # ######################################################### 
    gchu <- 0; HHg <- 0; TNorma2 <- 0 
    for(i in 1:s) {
      gchu[i] <- Estimg(x[i], X, hoptg[i]); 
      HHg[i] <- Histograma(x[i], nh, EstHist, eh); 
      TNorma2[i] <- d[i]*(gchu[i]-HHg[i])^2;
    }
    
    GCHU[,q] <- gchu; # Matriz de estimadores de gchu para
    # cada gam cada columna tiene un estimado de gchu
    
    # ######################################################### 
    # ######################################################### 
    # ##### Grafico del estimador gchu para cada       ###### 
    # ##### GAM[q], con q=1:\length(GAM)               ###### 
    # ######################################################### 
    plot(x, gchu, ylim=c(0, 0.6)) 
    lines(x, gchu, lwd=2, col="red")
    lines(x, go, lwd=2, col='black') 
    lines(x, HHg, lwd=2, col='blue') 
    # ########################################################
    
    Norma2gchuHHg[q] <- sum(TNorma2)
  }
  
  # ########################################################### 
  # ######### Grafico de la funcion Error(\gamma) ########### 
  # ######### Que determina el criterio para      ########### 
  # ######### calibrar la funcion de densidad g   ########### 
  # ###########################################################
  plot(Norma2gchuHHg, main="Calibración Método GL para g", xlab="", ylab="") 
  lines(Norma2gchuHHg)
  
  #**ACA PARA QUE CORRA**#
  # ###########################################################
  
  i0g <- which.min(Norma2gchuHHg) # Posicion del valor de gamma
  # que calibra al estimador gchu de g.
  
  # ########################################################### 
  # ########################################################### 
  # ############# Final Calibracion     ########### 
  # ############# Del estimador de g    ########### 
  # ########################################################### 
  # ###########################################################
  
  # ########################################################### 
  # ########################################################### 
  # ########### INICIO AJUSTE           ########### 
  # ########### CALIBRADO DE g          ########### 
  # ########################################################### 
  # ########################################################### 
  a <- 2 
  x1 <- seq(-a, a, 0.1); # Soporte para la estimacion. 
  s1 <- length(x1) # Numero de elementos en la particion del soporte de la densidad de X
  
  # ########################################################### 
  # ########################################################### 
  #### Rutina para almacenar los datos        ################### 
  #### usando la funcion VAdifEst             ################### 
  #### en una lista de nombre l que           ################### 
  #### contiene s matrices donde cada         ################### 
  #### matriz es de orden (mxm)               ################### 
  #### s=length(x) y m=length(H).             ################### 
  #### El alamcenamiento en cada matriz       ################### 
  #### se realiza por filas.                  ################### 
  # ########################################################### 
  l <- 0 # Nombre de la lista que almacenara s matrices
  # donde cada una es de orden mxm.
  
  for(k in 1:s1) {
    l0 <- 0 # Variable auxiliar que almacenara temporalmente
    # cada matriz de datos.
    l00 <- 0 # Variable auxiliar que almacena temporalmente
    # cada fila de la matriz l0.
    
    for(j in 1:m) {
      l00 <- VAdifEstg(x1[k], X, H, H[j])
      if(j==1) {l0 <- l00} else {l0 <- cbind(l0, l00)}
    }
    if(k==1) {l <- list(l0)} else {l[[k]] <- l0}
  }
  # ########################################################### 
  # ###########################################################
  
  gam <- GAM[i0g] # Parametro gamma > 2. (valor que he cambiado) 
  Q <- i0g # Posicion en la malla GAM, donde se calibro gchu,
  # es decir la columna
  # Q de GCHU (GCHU[,Q]) tiene la estimacion calibrada en una
  # malla de longitud 100.
  
  # ########################################################## 
  #### Calculo del estimador V(h), para cada h \in H. ##### 
  #### IMPORTANTE aca se hacen cambios locales,       ##### 
  #### Ahora V tambien depende de x, ademas de        ##### 
  #### depender de H(como antes)                      ##### 
  #### V sera una matriz de oden sxm, donde cada      ##### 
  #### fila corresponde a un valor de x y cada        ##### 
  #### columna a los valores de H.                    ##### 
  #### V[i,j] donde i=1:s y j=1:m                     ##### 
  # ########################################################## 
  V <- matrix(0, nrow=s1, ncol=m) 
  for(i in 1:s1) {
    for(k in 1:m) {
      V[i,k] <- sqrt(2*gsup*gam)*N2*(1+N1)*(1+delg) * ((log(ne))/(ne*H[k]))^(1/2)
    }
  }
  # ########################################################## 
  # ######### Fin del cambio local de V. #################### 
  # ##########################################################
  
  # ########################################################## 
  # #### Calculo del Maximo de cada fila ################# 
  # #### de la matriz {l[k]-V}_{+},      ################# 
  # #### para cada k de 1:s.             ################# 
  # #### Lo cual se guarda en la lista   ################# 
  # #### A={A[k]}_{k=1;s} con A[k] en R^m. ################# 
  # #### Se genera la lista AV, compuesta ################# 
  # #### de s vectores m dimencionales    ################# 
  # #### AV[[k]]=A[[k]]+V en R^m para     ################# 
  # #### k = 1, ..., s.                   ################# 
  # #### Se geneta el vector hopt en R^s  ################# 
  # #### donde hopt[k]=argmin(AV[[k]]+V)  ################# 
  # ########################################################## 
  A <- list() # Es una lista de s vectores, donde cada
  # vector A[k] en R^m, corresponde a un valor
  # x[k] de la malla, y A[k][i] es el valor en H[i]
  
  AV <- list() # Es una lista de s, ectores, donde cada
  # vector AV[k] coresponde a la suma de vectores A[k]+V.
  
  arg <- 0 # Esta Variable sera un vector donde cada
  # componente es arg[k]=argmin_{h \in H}(A[k]+V)
  
  hoptg <- 0 # Variable que almacena las ventanas optimas
  # para cada x[k] en la malla.
  
  for(k in 1:s1) {
    A0 <- 0 # Esta variable auxiliar, sera una matriz que
    # almacenara temporalmente los datos contenidos
    # en cada matris A[[k]] de la lista A
    
    for(i in 1:m) {
      A0[i] <- max(l[[k]][i,] - V[k,])
      if(A0[i]<0) {A0[i] <- 0}
    }
    # Se genera la lista A y AV 
    A[[k]] <- A0; # Se guardan resultados en la lista A 
    AV[[k]] <- A0+V[k,] # Se guardan resultados en la lista AV 
    arg[k] <- which.min(AV[[k]]) # Se almacena la ventana optima
    # para cada x[k] de la malla de x.
  }
  hoptg <- H[arg]
  
  # ########################################################## 
  # ########################################################## 
  # ###################### Estimador gchu ################### 
  # ########################################################## 
  # ########################################################## 
  gchu <- 0; 
  for(i in 1:s1) {
    gchu[i] <- Estimg(x1[i], X, hoptg[i]); # Coloque hopt^{-1} 
  }
  # ##########################################################
  
  # ########################################################## 
  # ############## FINALIZO EL AJUSTE             ############## 
  # ############## DE LA FUNCION DE DENSIDAD g    ############## 
  # ########################################################## 
  # ##########################################################
  
  # ############################################################ 
  # ################ GRAFICOS FINALES                ############# 
  # ################ PARA LA FUNCION DE DENSIDAD g   ############# 
  # ############################################################ 
  plot(x1, gchu, ylim=c(0, 0.6), type="n", main="Estimacion Adaptativa Método GL para g", xlab="", ylab="") 
  lines(x1, gchu, lwd=2, lty="dotted")
  
  go <- 0 
  HHg <- 0; 
  for(i in 1:s1) {
    go[i] <- dTnorm(x1[i], mu, sigm, a) 
    HHg[i] <- Histograma(x1[i], nh, EstHist, eh)
  }
  lines(x1, go, lwd=2, col='black') 
  # lines(x1, HHg, lwd=2, col='blue')
  
  #**ACA PARA QUE CORRA**#
  # ############################################################# 
  # ################ FIN GRAFICOS FINALES          ################ 
  # ################ DE LA FUNCION DE DENSIDAD g   ################ 
  # #############################################################
  
  # ############################################################# 
  # ################ CALIBRACION y ESTIMACION m    ############### 
  # ############################################################# 
  for(HH2 in 1:1) {
    print(HH2) 
    # ########################################################### 
    # ##################### MODULO 2 mss         ################## 
    # ########################################################### 
    EstimmmSS <- function(x, X, Y, h, h1) {
      # Estimador sobresuavisado con nucleo gausiano, muestra
      # aleatoria simple X, ventanas h y h1.
      n <- length(X);
      mmss = 0
      for(i in 1:n) {
        mmss = mmss + (1/(n*sqrt(2*pi*(h^2+h1^2)))) * Y[i] * exp(-((x-X[i])^2)/(2*(h^2+h1^2)));
      }
      return(mmss)
    }
    #**ACA PARA QUE CORRA**#
    # ########################################################### 
    # ###########################################################
    
    # ########################################################### 
    # ####################### MODULO 3 m         ################### 
    # ########################################################### 
    # Esta funcion determina el valor absoluto de la diferencia 
    # entre el estimador de la regresion sobresuavizado rhh1 y rh1 
    # en un valor x del soporte. 
    VAdifEstmm <- function(x, X, Y, h, h1) {
      VA <- abs(EstimmmSS(x, X, Y, h, h1) - Estimmm(x, X, Y, h1))
      return(VA)
    }
    # ######################################################### 
    # #########################################################
    
    # ######################################################### 
    # ######################################################### 
    #### Rutina para almacenar los datos        ################## 
    #### usando la funcion VAdifEst             ################## 
    #### en una lista de nombre l que           ################## 
    #### contiene s matrices donde cada         ################## 
    #### matriz es de orden (mxm)               ################## 
    #### s=length(x) y m=length(H).             ################## 
    #### El alamcenamiento en cada matriz       ################## 
    #### se realiza por filas.                  ################## 
    # ######################################################### 
    l <- list() # Nombre de la lista que almacenara s matrices
    # donde cada una es de orden mxm.
    
    for(k in 1:s) {
      l0 <- 0 # Variable auxiliar que almacenara temporalmente
      # cada matriz de datos.
      l00 <- 0 # Variable auxiliar que almacena temporalmente
      # cada fila de la matriz l0.
      
      for(j in 1:m) {
        l00 <- VAdifEstmm(x[k], X, Y, H, H[j])
        if(j==1) {l0 <- l00} else {l0 <- cbind(l0, l00)}
      }
      l[[k]] <- l0
    }
    # ########################################################## 
    # ##########################################################
    
    # ########################################################## 
    ## Cambio Local de rsup y ginf IMPORTANTE                 ## 
    ## Ahora rsup y ginf son vectores de dimension s          ## 
    ## Cada rsup[i] y ginf[i] corresponde con x[i], para i=1:s## 
    # ########################################################## 
    e <- rep(0, s+1) 
    se <- length(e) 
    e[4] <- (-2) 
    e[se] <- 2
    
    for(i in 2:(se-1)) {
      e[i] <- (x[i-1]+x[i])/2
    }
    
    rsup <- rep(0, s) 
    gchusup <- max(gchu)
    
    for(i in 1:s) {
      rsup[i] <- max(abs(Y*(e[i]<X)*(X<e[i+1]))) # Maximo Local de la funcion de regresion
    }
    # ########################################################### 
    # ############## FIN IMPORTANTE, CAMBIO LOCAL ############### 
    # ############## FIN IMPORTANTE, CAMBIO LOCAL ############### 
    # ###########################################################
    
    delm <- (log(ne))^(-1/5)
    N2 <- 1/sqrt(2*sqrt(pi)) #Norma 2 del nucleo gaussiano.
    
    # ########################################################### 
    #### Cambio Local de A3 IMPORTANTE                       #### 
    #### Ahora A3 es un vector de dimension s                #### 
    #### Cada A3[i] corresponde con x[i], para cada i=1:s    #### 
    # ########################################################### 
    A3 <- rep(0, s) 
    for(i in 1:s) {
      A3[i] <- (rsup[i]^2+sig^2)*gchusup*N2^2
    }
    
    # ########################################################## 
    # ############# FIN IMPORTANTE, CAMBIO LOCAL ############## 
    # ############# FIN IMPORTANTE, CAMBIO LOCAL ############## 
    # ##########################################################
    
    Norma2mchuYHHg <- rep(0, length(GAM)); # Variable que almacenara los valores de la funcion Error_m(gamma)
    
    # ########################################################## 
    # ############ CALIBRACION m                 ############# 
    # ########################################################## 
    for(q in 1:length(GAM)) {
      gam <- GAM[q] # Parametro gamma > 2. (valor que he cambiado)
      
      # ######################################################## 
      # ##### Calculo del estimador V(h), para cada h \in H. ### 
      # ##### IMPORTANTE aca se hacen cambios locales,       ### 
      # ##### Ahora V tambien depende de x, ademas de        ### 
      # ##### depender de H(como antes)                      ### 
      # ##### V sera una matriz de oden sxm, donde           ### 
      # ##### cada fila corresponde                          ### 
      # ##### a un valor de x y cada columna a los valores de H. ### 
      # ##### V[i,j] donde i=1:s y j=1:m                     ### 
      # ######################################################## 
      V <- matrix(0, nrow=s, ncol=m) 
      for(k in 1:m) {
        V[,k] <- sqrt(2*gam*A3)*(N1+1)*(1+delm) * sqrt(log(ne))/sqrt(ne*H[k])
      }
      #**ACA PARA QUE CORRA**#
      # ######################################################## 
      # ############# Fin del cambio local de V. ############### 
      # ########################################################
      
      # ######################################################## 
      # #### Calculo del Maximo de cada fila       ############## 
      # #### de la matriz {l[k]-V}_{+},            ############## 
      # #### para cada k de 1:s.                   ############## 
      # #### Lo cual se guarda en la lista         ############## 
      # #### A={A[k]}_{k=1;s} con A[k] en R^m      ############## 
      # #### Se genera la lista AV, compuesta      ############## 
      # #### de s vectores m dimencionales         ############## 
      # #### AV[[k]]=A[[k]]+V en R^m para          ############## 
      # #### k = 1, ..., s.                        ############## 
      # #### Se geneta el vector hopt en R^s       ############## 
      # #### donde hopt[k]=argmin(AV[[k]]+V)       ############## 
      # ######################################################## 
      A <- 0 # Es una lista de s vectores, donde cada
      # vector A[k] en R^m, corresponde a un valor x[k]
      # de la malla, y A[k][i] es el valor en H[i]
      
      AV <- 0 # Es una lista de s, vectores, donde cada vector
      # AV[k] coresponde a la suma de vectores A[k]+V.
      
      arg <- 0 # Esta Variable sera un vector donde cada
      # componente es arg[k]=argmin_{h \in H}(A[k]+V)
      
      hoptm <- 0 # Variable que almacena las ventanas optimas
      # para cada x[k] en la malla.
      
      # ######################################################## 
      ### IMPORTANTE, En la parte que sigue se hicieron      ### 
      ### CAMBIOS LOCALES, correspondiente a la variable V.  ### 
      # ######################################################## 
      for(k in 1:s) {
        A0 <- 0 # Esta variable auxiliar, sera una matriz que
        # almacenara temporalmente los datos contenidos
        # en cada matriz A[[k]] de la lista A
        
        for(i in 1:m) {
          A0[i] <- max(l[[k]][i,] - V[k,])
          if(A0[i]<0) {A0[i] <- 0}
        }
        
        if(k==1) {A <- list(A0); AV <- list(A0+V[k,])} # Se genera la lista A y AV
        else {A[[k]] <- A0; AV[[k]] <- A0+V[k,]} # Se generan las listas A y AV
        
        arg[k] <- which.min(AV[[k]]) # Se almacena la ventana
        # optima para cada x[k] de la malla de x.
      }
      # ######################################################### 
      # FIN CAMBIOS LOCALES, correspondientes a la variable V. # 
      # FIN CAMBIOS LOCALES, correspondientes a la variable V. # 
      # ######################################################### 
      hoptm <- H[arg]
      
      # ######################################################### 
      # ###################### Estimador mchu ################## 
      # ######################################################### 
      mchu <- 0; 
      TNorma2mchuYHHg <- 0; # Variable auxiliar que alamcenara
      # temporalmente los terminos del MSE de la estimacion.
      
      for(i in 1:s) {
        mchu[i] <- Estimmm(x[i], X, Y, hoptm[i]); 
        TNorma2mchuYHHg[i] <- d[i]*(mchu[i] - Yct[i]*GCHU[i,Q])^2;
      }
      # ########################################################
      
      # ######################################################## 
      # ##### Grafico del estimador mchu para cada GAM[q], ##### 
      # ##### con q=1:\length(GAM)                         ##### 
      # ########################################################
      plot(x, mchu, ylim=c(-0.5, 1)) 
      lines(x, mchu, lwd=2, col="red") 
      points(x, Yct*GCHU[,Q], col="yellow") 
      # #########################################################
      
      Norma2mchuYHHg[q] <- sum(TNorma2mchuYHHg);
    }
    
    # ########################################################### 
    # #### Grafico de la funcion Error_m(gamma) que permite ##### 
    # #### saber donde se calibra el estimador mchu de la   ##### 
    # #### funcion m.                                       ##### 
    # ########################################################### 
    plot(Norma2mchuYHHg, main="Calibración Método GL para m", xlab="", ylab="") 
    lines(Norma2mchuYHHg)
    
    #**ACA PARA QUE CORRA**#
    # ###########################################################
    
    i0m <- which.min(Norma2mchuYHHg) # En esta posicion de la malla
    # GAM se calibra mchu
    
    # ########################################################### 
    # ########################################################### 
    # ############ Final Calibracion m              ############# 
    # ########################################################### 
    # ###########################################################
    
    # ########################################################### 
    # ########################################################### 
    # ############# INICIO AJUSTE m                 ############## 
    # ########################################################### 
    # ###########################################################
    
    # ########################################################## 
    # ########################################################## 
    #### Rutina para almacenar los datos         ############### 
    #### usando la funcion VAdifEst              ############### 
    #### en una lista de nombre l que            ############### 
    #### contiene s matrices donde cada          ############### 
    #### matriz es de orden (mxm)                ############### 
    #### s=length(x) y m=length(H).              ############### 
    #### El alamcenamiento en cada matriz        ############### 
    #### se realiza por filas.                   ############### 
    # ##########################################################
    
    l <- 0 # Nombre de la lista que almacenara s matrices
    # donde cada una es de orden mxm.
    
    for(k in 1:s1) {
      l0 <- 0 # Variable auxiliar que almacenara temporalmente cada matriz de datos.
      l00 <- 0 # Variable auxiliar que almacena temporalmente cada fila de la matriz l0.
      
      for(j in 1:m) {
        l00 <- VAdifEstmm(x1[k], X, Y, H, H[j])
        if(j==1) {l0 <- l00} else {l0 <- cbind(l0, l00)}
      }
      if(k==1) {l <- list(l0)} else {l[[k]] <- l0}
    }
    # ######################################################### 
    # #########################################################
    
    gam <- GAM[i0m] # Parametro gamma > 2. (valor que he cambiado)
    
    # ######################################################### 
    #### Calculo del estimador V(h), para cada h \in H.    #### 
    #### IMPORTANTE aca se hacen cambios locales,          #### 
    #### Ahora V tambien depende de x, ademas de           #### 
    #### depender de H(como antes)                         #### 
    #### V sera una matriz de oden sxm, donde              #### 
    #### cada fila corresponde a un valor de x             #### 
    #### y cada columna a los valores de H.                #### 
    #### V[i,j] donde i=1:s y j=1:m                        #### 
    # ######################################################### 
    e <- rep(0, s1+1) 
    se <- length(e) 
    e[4] <- (-2) 
    e[se] <- 2
    
    for(i in 2:(se-1)) {
      e[i] <- (x1[i-1]+x1[i])/2
    }
    
    rsup <- rep(0, s1) 
    gchusup <- max(gchu)
    
    for(i in 1:s1) {
      rsup[i] <- max(abs(Y*(e[i]<X)*(X<e[i+1]))) # Maximo Local de la funcion de regresion
    }
    
    A3 <- rep(0, s1) 
    for(i in 1:s1) {
      A3[i] <- (rsup[i]^2+sig^2)*gchusup*N2^2
    }
    
    V <- matrix(0, nrow=s1, ncol=m) 
    for(k in 1:m) {
      V[,k] <- sqrt(2*gam*A3)*(N1+1)*(1+delm) * sqrt(log(ne))/sqrt(ne*H[k])
    }
    # ########################################################## 
    # ######### Fin del cambio local de V.          ############ 
    # ##########################################################
    
    # ########################################################## 
    # #### Calculo del Maximo de cada fila        ################# 
    # #### de la matriz {l[k]-V}_{+},             ################# 
    # #### para cada k de 1:s.                    ################# 
    # #### Lo cual se guarda en la lista          ################# 
    # #### A={A[k]}_{k=1;s} con A[k] en R^m       ################# 
    # #### Se genera la lista AV, compuesta       ################# 
    # #### de s vectores m dimencionales          ################# 
    # #### AV[[k]]=A[[k]]+V en R^m para           ################# 
    # #### k = 1, ..., s.                         ################# 
    # #### Se geneta el vector hopt en R^s        ################# 
    # #### donde hopt[k]=argmin(AV[[k]]+V)        ################# 
    # ########################################################## 
    A <- list() # Es una lista de s vectores, donde cada
    # vector A[k] en R^m, corresponde a un valor x[k]
    # de la malla, y A[k][i] es el valor en H[i]
    
    AV <- list() # Es una lista de s, ectores, donde cada vector
    # AV[k] coresponde a la suma de vectores A[k]+V.
    
    arg <- 0 # Esta Variable sera un vector donde cada
    # componente es arg[k]=argmin_{h \in H}(A[k]+V)
    
    hoptm <- 0 # Variable que almacena las ventanas optimas
    # para cada x[k] en la malla.
    
    for(k in 1:s1) {
      A0 <- 0 # Esta variable auxiliar, sera una matriz que
      # almacenara temporalmente los datos contenidos
      # en cada matris A[[k]] de la lista A
      
      for(i in 1:m) {
        A0[i] <- max(l[[k]][i,] - V[k,])
        if(A0[i]<0) {A0[i] <- 0}
      }
      # Se genera la lista A y AV 
      A[[k]] <- A0; # Se guardan resultados en la lista A 
      AV[[k]] <- A0+V[k,] # Se guardan resultados en la lista AV 
      arg[k] <- which.min(AV[[k]]) # Se almacena la ventana optima
      # para cada x[k] de la malla de x
    }
    
    hoptm <- H[arg]
    # ##########################################################
    
    # ########################################################## 
    # ###################### Estimador mchu  ################### 
    # ########################################################## 
    mchu <- 0; 
    for(i in 1:s1) {
      mchu[i] <- Estimmm(x1[i], X, Y, hoptm[i]);
    }
    # ##########################################################
    
    # ########################################################## 
    # ############# FINALIZO EL AJUSTE m           ############## 
    # ##########################################################
  }
  
  # ############################################################ 
  # ###### Grafico de mchu estimador de m               ######## 
  # ###### Junto a la nube de puntos de entrenamiento   ######## 
  # ###### y la funcion m, en el intervalo [-2,2]       ######## 
  # ############################################################ 
  plot(x1, mchu, ylim=c(-0.5, 1), type="n", main="Estimacion Adaptativa Método GL para m", xlab="", ylab="") 
  lines(x1, mchu, lwd=2, lty="dotted") # Curva del estimador mchu 
  points(x, Yct*GCHU[,Q], col="yellow") # Nube de puntos de entrenamiento
  
  r = 0.7*x1+2*exp(-10*x1^2) # Funcion de regresion a estimar
  # r(x)=0.7x+2exp(-10x^2) evaluada en la malla x1 del intervalo [-2,2]
  
  lines(x1, r*go, lwd=2, col='black') # Grafico de m=r*go evaluado en la malla x1 del intervalo [-2,2]
  
  #**ACA PARA QUE CORRA**#
  # ############################################################
  
  # ############################################################ 
  ### GRAFICO de rchu=mchu/gchu en el intervalo [-2,2]       ### 
  # ############################################################ 
  rchu <- mchu/gchu; # Estimador rchu de la funcion de regrion r
  
  plot(x1, rchu, ylim=c(-2.5, 3), type="n", main="Estimacion Adaptativa Método GL para r", xlab="", ylab="") 
  points(X, Y, col="yellow") # Nube de puntos de estimacion 
  lines(x1, rchu, lwd=2, lty="dotted") # Curva del estimador mchu de la regresion r
  
  lines(x1, r, lwd=2, col='black') # Curva de la funcion de regresion r
  
  #**ACA PARA QUE CORRA**#
  # ############################################################
  
  # ############################################################ 
  # ####### Grafico de rchu en el intervalo [-1,1]     ######### 
  # ############################################################ 
  j2 <- seq(11, 31, 1) # Indices para que la malla del soporte este en [-1,1]
  
  plot(x1[j2], rchu[j2], ylim=c(-2, 3), type="n", main="Estimacion Adaptativa Método GL para r", xlab="", ylab="") 
  points(X, Y, col="yellow") # Nube de puntos de estimacion 
  lines(x1[j2], rchu[j2], lwd=2, lty="dotted") # Curva del estimador mchu en [-1,1]
  
  r = 0.7*x1[j2]+2*exp(-10*x1[j2]^2) # Funcion de regresion a estimar
  # r(x)=0.7x+2exp(-10x^2) evaluada en la malla x1[j2] del intervalo [-1,1]
  
  lines(x1[j2], r, lwd=2, col='black') # Curva de la funcion de regresion r, en el intervalo [-1,1]
  
  #**ACA PARA QUE CORRA**#
  # ############################################################
}
