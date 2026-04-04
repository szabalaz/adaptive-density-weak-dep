# =================================================================================
# EFECTO DEL PARÁMETRO VENTANA (h) CON VENTANA FIJA (h) EN UNA MIXTURA DE NORMALES
# CON LA VERDADERA DENSIDAD SUPERPUESTA
# =================================================================================

# 1. Fijar semilla para reproducibilidad
set.seed(2023)

# 2. Generar datos simulados de una mixtura bimodal
#    50% de N(-2, 1) y 50% de N(2, 1)
n <- 250
componentes <- sample(1:2, prob = c(0.5, 0.5), size = n, replace = TRUE)
medias <- c(-2, 2)
desviaciones <- c(1, 1)
datos <- rnorm(n, mean = medias[componentes], sd = desviaciones[componentes])

# 3. Calcular la VERDADERA curva de densidad teórica
#    Creamos una secuencia fina de valores de x para trazar la curva
x_seq <- seq(-6, 6, length.out = 500)
verdadera_densidad <- 0.5 * dnorm(x_seq, mean = -2, sd = 1) + 
  0.5 * dnorm(x_seq, mean = 2, sd = 1)

# 4. Definir los valores de h 
h_valores <- c(0.1, 0.3, 0.7)

# 5. Configurar el lienzo para 1 fila y 3 columnas
par(mfrow = c(1, 3), mar = c(5, 4, 4, 1) + 0.1)

# 6. Generar los tres gráficos
for (h in h_valores) {
  # Calcular la estimación por núcleo (Gaussiano por defecto)
  estimacion <- density(datos, bw = h, kernel = "gaussian")
  
  # Crear el gráfico con la densidad ESTIMADA
  plot(estimacion, 
       main = paste("Ancho de ventana h =", h),
       xlab = "Valores de x", 
       ylab = "Densidad", 
       col = "darkblue", 
       lwd = 2,
       ylim = c(0, 0.25)) # Eje Y fijo para poder comparar
  
  # SUPERPONER la densidad VERDADERA
  lines(x_seq, verdadera_densidad, col = "red", lwd = 2)
  
  # Añadir marcas de los datos reales en el eje X ("rug")
  rug(datos, col = "gray40", ticksize = 0.03)
  
  # Añadir leyenda
  legend("topright", 
         legend = c("Estimada", "Verdadera"), 
         col = c("darkblue", "red"), 
         lwd = 2, 
         bty = "n", # Sin caja alrededor de la leyenda
         cex = 0.9)
}

# Restaurar los parámetros gráficos a su valor original
par(mfrow = c(1, 1))

# =====================================================================================================
# COMPARACIÓN DE MÉTODOS CLÁSICOS PARA SELECCIONAR EL PARÁMETRO VENTANA (h) EN UNA MIXTURA DE NORMALES
# CON LA VERDADERA DENSIDAD SUPERPUESTA
# =====================================================================================================
# Seleccionando h por cada método
h1 <- bw.SJ(datos) # Sheather y Jones (1991)
h2 <- bw.nrd0(datos) # Silverman(1986)
h3 <- bw.nrd(datos)  # Scott(1992)
h4 <- bw.ucv(datos)  # CV insesgada
h5 <- bw.bcv(datos)  # CV sesgada
# h_valores_mc <- c(bw.SJ(datos), # Sheather y Jones (1991)
#                   bw.nrd0(datos), # Silverman(1986)
#                   bw.nrd(datos),  # Scott(1992)
#                   bw.ucv(datos),  # CV insesgada
#                   bw.bcv(datos))  # CV sesgada
# Estimando la densidad para cada h
npden1 <- density(datos, bw = h1) 
npden2 <- density(datos, bw = h2)
npden3 <- density(datos, bw = h3) 
npden4 <- density(datos, bw = h4) 
npden5 <- density(datos, bw = h5) 

# Crear el gráfico con la densidad ESTIMADA
hist(datos, freq = FALSE, main = " Comparación de métodos clásicos para seleccionar h",
     xlab = "Valores de x", ylab = "Densidad",
     border = "darkgray")

lines(x_seq, verdadera_densidad, col = "red", lwd = 2)
lines(npden1, lwd = 2, col = "blue")
lines(npden2, lwd = 2, col = "lightblue")
lines(npden3, lwd = 2, col = "green")
lines(npden4, lwd = 2, col = "darkblue")
lines(npden5, lwd = 2, col = "darkgreen")
#rug(datos, col = "darkgray")

# Añadir leyenda
legend("topright", 
       legend = c("Verdadera","SJ", "Silverman", "Scott", "UCV", "BCV"), 
       col = c("red","blue", "lightblue", "green","darkblue", "darkgreen"), 
       lwd = 2, 
       bty = "n", # Sin caja alrededor de la leyenda
       cex = 0.9)

# =====================================================================
# INFLUENCIA MARGINAL DE LA FUNCIÓN NÚCLEO (KERNEL)
# CON ANCHO DE VENTANA (h) FIJO
# =====================================================================
# 1. Fijar semilla para reproducibilidad
set.seed(123)

# 2. Generar datos de la mixtura: 50% N(-2, 1) y 50% N(2, 1)
n <- 300
pesos <- c(0.5, 0.5)
medias <- c(-2, 2)
desviaciones <- c(1, 1)

componentes <- sample(1:2, prob = pesos, size = n, replace = TRUE)
datos <- rnorm(n, mean = medias[componentes], sd = desviaciones[componentes])

# Función de la verdadera densidad teórica
verdadera_densidad <- function(x) {
  pesos[1] * dnorm(x, mean = medias[1], sd = desviaciones[1]) +
    pesos[2] * dnorm(x, mean = medias[2], sd = desviaciones[2])
}

# 3. Fijar el ancho de ventana (h)
h_fijo <- 0.6 
h_fijos <- c(0.6, 0.8, 3)

# 4. Definir los tipos de núcleo a comparar
kernels <- c("gaussian", "epanechnikov", "rectangular", "biweight")
nombres_kernels <- c("Gaussiano", "Epanechnikov", "Rectangular", "Biweight")

# 5. Configurar el lienzo para 2x2
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1) + 0.1)

# 6. Calcular límites de los ejes (ylim hasta 0.30 da espacio perfecto arriba)
eje_x <- c(-6, 6)
eje_y <- c(0, 0.30) 

# 7.1 Bucle para generar los gráficos con h = 0.6
for (i in 1:length(kernels)) {
  estimacion <- density(datos, bw = h_fijos[1], kernel = kernels[i])
  
  plot(estimacion, 
       main = paste("Núcleo:", nombres_kernels[i], "| h =", h_fijos[1]),
       xlab = "x", 
       ylab = "Densidad", 
       col = "darkred", 
       lwd = 2,
       xlim = eje_x,
       ylim = eje_y,
       las = 1) 
  
  curve(verdadera_densidad(x), col = "blue", lwd = 2, add = TRUE)
  
  rug(datos, col = "gray40", ticksize = 0.05)
  
  # Leyenda en el centro superior, donde está el "valle" entre las dos campanas
  legend("top", legend = c("Estimada", "Teórica"), 
         col = c("darkred", "blue"), lwd = 2, 
         bty = "n", cex = 0.8, horiz = TRUE)
}

# 7.2 Bucle para generar los gráficos con h = 0.8
for (i in 1:length(kernels)) {
  estimacion <- density(datos, bw = h_fijos[2], kernel = kernels[i])
  
  plot(estimacion, 
       main = paste("Núcleo:", nombres_kernels[i], "| h =", h_fijos[2]),
       xlab = "x", 
       ylab = "Densidad", 
       col = "darkred", 
       lwd = 2,
       xlim = eje_x,
       ylim = eje_y,
       las = 1) 
  
  curve(verdadera_densidad(x), col = "blue", lwd = 2, add = TRUE)
  
  rug(datos, col = "gray40", ticksize = 0.05)
  
  # Leyenda en el centro superior, donde está el "valle" entre las dos campanas
  legend("top", legend = c("Estimada", "Teórica"), 
         col = c("darkred", "blue"), lwd = 2, 
         bty = "n", cex = 0.8, horiz = TRUE)
}

# 7.3 Bucle para generar los gráficos con h = 3
for (i in 1:length(kernels)) {
  estimacion <- density(datos, bw = h_fijos[3], kernel = kernels[i])
  
  plot(estimacion, 
       main = paste("Núcleo:", nombres_kernels[i], "| h =", h_fijos[3]),
       xlab = "x", 
       ylab = "Densidad", 
       col = "darkred", 
       lwd = 2,
       xlim = eje_x,
       ylim = eje_y,
       las = 1) 
  
  curve(verdadera_densidad(x), col = "blue", lwd = 2, add = TRUE)
  
  rug(datos, col = "gray40", ticksize = 0.05)
  
  # Leyenda en el centro superior, donde está el "valle" entre las dos campanas
  legend("top", legend = c("Estimada", "Teórica"), 
         col = c("darkred", "blue"), lwd = 2, 
         bty = "n", cex = 0.8, horiz = TRUE)
}

# Restaurar lienzo
par(mfrow = c(1, 1))


# # 1. Fijar semilla para reproducibilidad
# set.seed(123)
# 
# # 2. Generar datos simulados de la mixtura (3/4 N(0,1) y 1/4 N(1.5, 1/3))
# n <- 300
# pesos <- c(0.75, 0.25)
# medias <- c(0, 1.5)
# desviaciones <- c(1, 1/3)
# 
# componentes <- sample(1:2, prob = pesos, size = n, replace = TRUE)
# datos <- rnorm(n, mean = medias[componentes], sd = desviaciones[componentes])
# 
# # Función de la verdadera densidad teórica
# verdadera_densidad <- function(x) {
#   pesos[1] * dnorm(x, mean = medias[1], sd = desviaciones[1]) +
#     pesos[2] * dnorm(x, mean = medias[2], sd = desviaciones[2])
# }
# 
# # 3. Fijar el ancho de ventana (h)
# h_fijo <- 0.6 
# 
# # 4. Definir los tipos de núcleo a comparar
# kernels <- c("gaussian", "epanechnikov", "rectangular", "biweight")
# nombres_kernels <- c("Gaussiano", "Epanechnikov", "Rectangular", "Biweight")
# 
# # 5. Configurar el lienzo para 2x2
# par(mfrow = c(2, 2), mar = c(4, 4, 3, 1) + 0.1)
# 
# # 6. LÍMITES FIJOS (La clave para que no se superponga nada)
# eje_x <- c(-4, 4)   # Espacio de sobra a los lados
# eje_y <- c(0, 0.45) # Techo alto garantizado (las curvas llegan hasta ~0.30)
# 
# # 7. Bucle para generar los gráficos
# for (i in 1:length(kernels)) {
#   # Calcular la estimación
#   estimacion <- density(datos, bw = h_fijo, kernel = kernels[i])
#   
#   # Crear el gráfico
#   plot(estimacion, 
#        main = paste("Núcleo:", nombres_kernels[i], "| h =", h_fijo),
#        xlab = "x", 
#        ylab = "Densidad", 
#        col = "darkred", 
#        lwd = 2,
#        xlim = eje_x,
#        ylim = eje_y,
#        las = 1) # Pone los números del eje Y en horizontal (más estético)
#   
#   # Añadir la curva teórica
#   curve(verdadera_densidad(x), col = "blue", lwd = 2, lty = 2, add = TRUE)
#   
#   # Marcas de los datos reales en la base
#   rug(datos, col = "gray40", ticksize = 0.05)
#   
#   # Leyenda arriba a la derecha, pequeña y sin marco
#   legend("topright", legend = c("Estimada", "Teórica"), 
#          col = c("darkred", "blue"), lwd = 2, lty = c(1, 2), 
#          bty = "n", cex = 0.8)
# }
# 
# # Restaurar lienzo a su estado original
# par(mfrow = c(1, 1))
# 
# #####################################
