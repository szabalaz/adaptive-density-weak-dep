# Comparación de tasas de convergencia
# Autor: Sebastián Zabala
# Fecha: 21/04/2026

library(latex2exp)

# Parámetro de tamaño muestral
n_max <- 10000
n <- seq_len(n_max)

# Tasas de convergencia
tasa_parametrica <- n^(-1)
tasa_np_optima    <- n^(-2/5)
tasa_np_glm       <- (log(n) / n)^(2/5)
tasa_np_alt       <- (n / log(n))^(-2/5)

# Configuración gráfica
op <- par(mfrow = c(2, 2), mar = c(4, 4, 3, 1) + 0.1)
on.exit(par(op), add = TRUE)

plot(n, tasa_parametrica, type = "l", col = "blue", lwd = 2,
     main = TeX(r"(Evolución de la tasa paramétrica $n^{-1}$)"),
     xlab = "n", ylab = TeX(r"($n^{-1}$)"), log = "x")

plot(n, tasa_np_optima, type = "l", col = "red", lwd = 2,
     main = TeX(r"(Evolución de la tasa no paramétrica $n^{-2/5}$)"),
     xlab = "n", ylab = TeX(r"($n^{-2/5}$)"), log = "x")

plot(n, tasa_np_glm, type = "l", col = "darkgreen", lwd = 2,
     main = TeX(r"(Tasa no paramétrica $(\log n / n)^{2/5}$)"),
     xlab = "n", ylab = TeX(r"($(log n / n)^{2/5}$)"), log = "x")

plot(n, tasa_np_alt, type = "l", col = "purple", lwd = 2,
     main = TeX(r"(Tasa no paramétrica $(n / \log n)^{-2/5}$)"),
     xlab = "n", ylab = TeX(r"($(n / \log n)^{-2/5}$)"), log = "x")