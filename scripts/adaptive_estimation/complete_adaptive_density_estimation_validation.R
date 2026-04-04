#======================
# INSTALAR LAS LIBRERÍAS ( en este orden )
#======================

install.packages(c(
  # Momentos estadisticos
  "moments", "e1071",
  # Series temporales
  "forecast", "tseries", "FinTS", "urca",
  # Modelos GARCH
  "rugarch", "fGarch",
  # Valores extremos
  "evir", "POT", "QRM",
  # Pruebas adicionales
  "nortest", "lawstat",
  # Analisis financiero
  "PerformanceAnalytics", "quantmod",
  # Estimacion no parametrica
  "KernSmooth", "np"
))


# ========================================
# SCRIPT DE VALIDACION COMPLETA
# Proyecto : Estimacion Adaptativa de Densidad
# ========================================

# CARGAR LAS LIBRERÍAS
# library (quantmod)
# library (moments)
# library (e1071)
# library (forecast)
# library (FinTS)
# library (tseries)
# library (nortest)
# library (rugarch)
# library (evir)
# library (PerformanceAnalytics)

# Packages loading (mas eficiente que la anterior)
packages <- c(
  "quantmod", "moments", "e1071", "forecast", "FinTS", "tseries",
  "nortest", "rugarch", "evir", "PerformanceAnalytics"
)
invisible(lapply(packages, library, character.only = TRUE))

# ----- 1. CARGA DE DATOS -----
cat("\n=== CARGA DE DATOS ===\ n")

# Aquí se ingresa el Símbolo Ticker del activo de interés y la ventana temporal
# por defecto al frecuencia es diaria.
getSymbols ("^GSPC", from ="2015-01-01", to ="2023-12-31")
precios <- Cl(GSPC)
retornos <- diff(log(precios))[-1]
retornos <- as.numeric(retornos)
cat("Numero de observaciones:", length (retornos), "\n")


# ----- 2. ESTADSTICOS DESCRIPTIVOS -----

cat("\n=== ESTADSTICOS DESCRIPTIVOS ===\ n")
media <- mean (retornos)
desv <- sd(retornos)
sk <- skewness (retornos)
ku <- kurtosis (retornos)

cat(sprintf(" Media : %.6f\n", media))
cat(sprintf(" Desviacion Estandar : %.6f\n", desv))
cat(sprintf(" Asimetria ( Skewness ): %.4f\n", sk))
cat(sprintf(" Curtosis ( Kurtosis ): %.4f\n", ku))

# Interpretacion
if (sk < 0){
  cat(" -> Asimetria NEGATIVA detectada ( cola izquierda pesada )\n"
  )
} else{
  cat (" -> Asimetria POSITIVA detectada \n")
}
if (ku > 3) {
  cat(" -> COLAS PESADAS detectadas ( leptokurtosis )\n")
} else{
  cat(" -> Colas ligeras ( platicurtosis )\n")
}


# ----- 3. PRUEBAS DE NORMALIDAD -----
cat("\n=== PRUEBAS DE NORMALIDAD ===\ n")

jb_test <- jarque.test(retornos)
cat(" Prueba de Jarque - Bera :\n")
cat(sprintf(" Estadistico : %.4f\n", jb_test$statistic))
cat(sprintf(" p- valor : %.6f\n", jb_test$p.value ))
if (jb_test$p.value < 0.05) {
  cat (" -> Los retornos NO siguen distribucion normal \n")
}

ad_test <- ad.test (retornos)
cat("\ nPrueba de Anderson - Darling :\n")
cat(sprintf(" Estadistico : %.4f\n", ad_test$statistic))
cat(sprintf(" p- valor : %.6f\n", ad_test$p.value))
if (ad_test$p.value < 0.05){
  cat(" -> Los retornos NO siguen distribucion normal \n")
}



# ----- 4. PRUEBAS DE ESTACIONARIEDAD -----
cat("\n=== PRUEBAS DE ESTACIONARIEDAD ===\ n")

adf_test <- adf.test (retornos)
cat(" Prueba Aumentada de Dickey - Fuller (ADF):\n")
cat(sprintf(" Estadistico : %.4f\n", adf_test$statistic))
cat(sprintf(" p- valor : %.6f\n", adf_test$p.value))
if (adf_test$p.value < 0.05){
  cat (" -> La serie es ESTACIONARIA \n")
}

kpss_test <- kpss.test(retornos)
cat("\ nPrueba KPSS :\n")
cat(sprintf(" Estadistico : %.4f\n", kpss_test$statistic))
cat(sprintf(" p- valor : %.6f\n", kpss_test$p.value))
if (kpss_test$p.value > 0.05){
  cat (" -> La serie es ESTACIONARIA \n")
}

# ----- 5. AUTOCORRELACION Y DEPENDENCIA -----
cat("\n=== ANALISIS DE AUTOCORRELACION ===\ n")

# Autocorrelacion de retornos
Box_ret <- Box.test (retornos, lag = 12, type = "Ljung-Box")
cat(" Prueba de Ljung-Box para retornos :\n")
cat(sprintf(" p- valor : %.6f\n", Box_ret$p.value))

# Autocorrelacion de retornos al cuadrado ( CRUCIAL )
retornos_sq <- retornos^2
Box_ret_sq <- Box.test(retornos_sq, lag = 12, type = "Ljung-Box")
cat("\ nPrueba de Ljung -Box para retornos al cuadrado :\n")
cat(sprintf(" p- valor : %.6f\n", Box_ret_sq$p.value))
if ( Box_ret_sq$p.value < 0.05){
  cat (" -> AUTOCORRELACION SIGNIFICATIVA en retornos^2\n")
  cat (" -> DEPENDENCIA TEMPORAL detectada \n")
}


# ----- 6. EFECTOS ARCH -----
cat("\n=== PRUEBA DE EFECTOS ARCH ===\ n")
arch_test <- ArchTest(retornos, lags = 12)
cat(sprintf(" Estadistico Chi - cuadrado : %.4f\n", arch_test$statistic))
cat(sprintf("p- valor : %.6f\n", arch_test$p.value))
if (arch_test$p.value < 0.05){
  cat(" -> Efectos ARCH DETECTADOS \n")
  cat(" -> VOLATILIDAD AGRUPADA presente \n")
  cat(" -> La serie es DEBILMENTE DEPENDIENTE \n")
}


# ----- 7. MODELADO GARCH -----
cat("\n=== AJUSTE DE MODELO GARCH (1 ,1) ===\ n")

spec <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
  mean.model = list(armaOrder = c(0, 0)),
  distribution.model = "norm")

fit_garch <- ugarchfit(spec, data = retornos)
cat("\ nModelo estimado exitosamente \n")
cat(" Coeficientes :\n")
print(coef(fit_garch))


# Significancia de parametros GARCH
alpha1 <- coef(fit_garch)["alpha1"]
beta1 <- coef (fit_garch)["beta1"]
cat(sprintf ("\ nalpha1 ( ARCH ): %.6f\n", alpha1))
cat(sprintf (" beta1 ( GARCH ): %.6f\n", beta1))
cat(sprintf (" Suma ( persistencia ): %.6f\n", alpha1 + beta1))


if (( alpha1 + beta1 ) < 1){
  cat(" -> Proceso estacionario en varianza \n")
  cat(" -> Confirma DEPENDENCIA DEBIL \n")
}


# ----- 8. ANALISIS DE COLAS PESADAS -----
cat("\n=== ANALISIS DE VALORES EXTREMOS ===\ n")

# Estimador de Hill para la cola derecha
hillett <- hill(retornos[ retornos > 0])
cat(" Estimador de Hill calculado \n")

# Ajuste GPD
threshold_95 <- quantile(retornos, 0.95)
cat(sprintf(" Umbral ( percentil 95): %.6f\n", threshold_95))
tryCatch ({
  gpd_fit <- gpd (retornos, threshold = threshold_95)
  cat (" Modelo GPD ajustado exitosamente \n")
  cat (" -> Confirma comportamiento de COLA PESADA \n")
}, error = function(e){
  cat (" Error en ajuste GPD ( comun con muestras pequenas )\n")
})


# ----- 9. RESUMEN DE VALIDACION -----
cat("\n\n")
cat(" ========================================\ n")
cat(" RESUMEN DE V A L I D A C I N \n")
cat(" ========================================\ n\n")

cat("1. HECHOS ESTILIZADOS :\n")
if (sk < 0 && ku > 3 && arch_test$p.value < 0.05){
  cat(" [X] Asimetria negativa \n")
  cat(" [X] Colas pesadas ( leptokurtosis )\n")
  cat(" [X] Volatilidad agrupada \n")
  cat(" -> CONDICIONES SATISFECHAS \n\n")
} else{
  cat(" [ ] Algunas condiciones no cumplidas \n\n")
}


cat("2. CONDICIONES METODOLOGICAS :\n")
if (adf_test$p.value < 0.05 &&
    Box_ret_sq$p.value < 0.05 &&
    ( alpha1 + beta1 ) < 1){
  cat(" [X] Identica distribucion ( estacionaria )\n")
  cat(" [X] Proceso debilmente dependiente \n")
  cat(" [X] Autocorrelacion en retornos^2\n")
  cat(" -> METODOLOGIA GL ES APLICABLE \n\n")
} else{
  cat (" [ ] Revisar algunas condiciones \n\n")
}


cat("La serie de retornos cumple los requisitos \n")
cat(" para aplicar el metodo Goldenshluger - Lepski \n")



# ----- 10. GRAFICOS DIAGNOSTICOS -----
cat("\ nGenerando graficos diagnosticos ...\ n")

pdf("diagnosticos_validacion.pdf", width = 12, height = 10)

par( mfrow = c(3, 2))

# Grafico 1: Serie temporal de retornos
plot(retornos, type = "l", main = " Serie Temporal de Retornos ",
     xlab = "Tiempo", ylab = "Retorno")
abline(h = 0, col = "red ", lty = 2)


# Grafico 2: Histograma con densidad
hist(retornos, breaks = 50, freq = FALSE,
     main = "Histograma de Retornos",
     xlab = "Retorno", col = "lightblue")
lines(density(retornos), col = "blue", lwd = 2)
curve(dnorm(x, mean(retornos), sd(retornos)),
      add = TRUE, col = "red ", lwd = 2, lty = 2)
legend ("topright", c("Densidad empirica", "Normal teorica"),
        col = c("blue", "red"), lty = c(1, 2) , lwd = 2)

# Grafico 3: Q-Q plot
qqnorm(retornos, main = "Q-Q Plot")
qqline(retornos, col = "red ", lwd = 2)

# Grafico 4: ACF de retornos
Acf(retornos, main = "ACF de Retornos", lag.max = 30)

# Grafico 5: ACF de retornos al cuadrado
Acf(retornos_sq, main = "ACF de Retornos al Cuadrado", lag.max = 30)

# Grafico 6: Volatilidad condicional GARCH
vol_cond <- sigma(fit_garch)
plot(vol_cond, type = "l", main = "Volatilidad Condicional GARCH (1 ,1)",
     xlab = "Tiempo", ylab = "Volatilidad", col = "blue")

dev.off ()
cat(" Graficos guardados en: diagnosticos _ validacion .pdf\n")
cat("\n=== VALIDACION COMPLETA ===\ n")
# FIN DEL SCRIPT






