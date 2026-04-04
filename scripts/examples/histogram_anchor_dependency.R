# ==============================================================================
# DEPENDENCIA DEL HISTOGRAMA AL ANCLA
# ==============================================================================
library(ggplot2)
library(patchwork)

# 1. Simulación de datos
set.seed(2023) 
n <- 40        
datos <- data.frame(x = rnorm(n, mean = 0, sd = 1.2))

# Parámetros globales
h <- 1.0       # Ancho de ventana
y_limit <- 16  # Límite fijo eje Y

# 2. Función de graficado minimalista
plot_clean <- function(data, ancla, label_text) {
  ggplot(data, aes(x = x)) +
    geom_histogram(
      binwidth = h,
      boundary = ancla,
      fill = "white",    
      color = "black",   
      linewidth = 0.6    # CORRECCIÓN: 'linewidth' en lugar de 'size'
    ) +
    # Marcas de los datos en el eje X
    geom_rug(sides = "b", alpha = 0.6, linewidth = 0.5) +
    # Límites fijos para permitir comparación
    coord_cartesian(ylim = c(0, y_limit), xlim = c(-4, 4)) +
    # Etiqueta discreta del ancla
    annotate("text", x = -3.8, y = y_limit * 0.9, 
             label = label_text, hjust = 0, size = 3.5, fontface = "italic") +
    labs(x = NULL, y = NULL) + 
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey92"),
      axis.text = element_text(size = 9, color = "black")
    )
}

# 3. Generar los 4 gráficos
p1 <- plot_clean(datos, ancla = 0.15,        label_text = "(a) Origen: 0.15") + labs(y = "Frecuencia")
p2 <- plot_clean(datos, ancla = 0.25,     label_text = "(b) Origen: 0.25")
p3 <- plot_clean(datos, ancla = 0.50,     label_text = "(c) Origen: 0.50") + labs(y = "Frecuencia", x = "x")
p4 <- plot_clean(datos, ancla = 0.75,     label_text = "(d) Origen: 0.75") + labs(x = "x")

# 4. Composición final (2x2)
grafico_final <- (p1 + p2) / (p3 + p4)

# 5. Guardar
ggsave("histograma_ancla.png", plot = grafico_final, width = 7, height = 5, dpi = 300)

# ==============================================================================
# SCRIPT: DEPENDENCIA DEL HISTOGRAMA AL ANCHO DE VENTANA (h)
# ==============================================================================


# 1. Simulación de datos
set.seed(2023) 
n <- 40        # Aumentamos ligeramente n para que el efecto del h pequeño sea visible
datos <- data.frame(x = rnorm(n, mean = 0, sd = 1.2))

# Parámetro fijo
ancla_fija <- 0  # El origen se mantiene constante en todos los gráficos

# 2. Función de graficado
plot_bandwidth <- function(data, h, label_text) {
  ggplot(data, aes(x = x)) +
    geom_histogram(
      binwidth = h,
      boundary = ancla_fija, # Origen fijo
      fill = "white",    
      color = "black",   
      linewidth = 0.6
    ) +
    # Rug plot para ver los datos individuales
    geom_rug(sides = "b", alpha = 0.6, linewidth = 0.5) +
    # Fijamos el eje X para que sean comparables, eje Y libre (pues la frecuencia cambia con h)
    coord_cartesian(xlim = c(-4, 4)) +
    # Etiqueta discreta
    annotate("text", x = -3.8, y = Inf, vjust = 1.5,
             label = label_text, hjust = 0, size = 3.5, fontface = "italic") +
    labs(x = NULL, y = "Frecuencia") + 
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey92"),
      axis.text = element_text(size = 9, color = "black"),
      axis.title.y = element_text(size = 8)
    )
}

# 3. Generar 4 gráficos con anchos de ventana progresivos
# Notarás cómo la forma cambia de "ruidosa" a "bloque"
p1 <- plot_bandwidth(datos, h = 0.25, label_text = "(a) Ancho: 0.25 (Subsuavizado)")
p2 <- plot_bandwidth(datos, h = 0.50, label_text = "(b) Ancho: 0.50")
p3 <- plot_bandwidth(datos, h = 1.00, label_text = "(c) Ancho: 1.00")
p4 <- plot_bandwidth(datos, h = 2.00, label_text = "(d) Ancho: 2.00 (Sobresuavizado)")

# 4. Composición final
# Ajustamos las etiquetas de los ejes para que no sean redundantes
p1 <- p1 + labs(x = NULL)
p2 <- p2 + labs(x = NULL, y = NULL)
p3 <- p3 + labs(x = "x")
p4 <- p4 + labs(x = "x", y = NULL)

grafico_final <- (p1 + p2) / (p3 + p4)

# 5. Guardar
ggsave("histograma_ancho.png", plot = grafico_final, width = 7, height = 5, dpi = 300)



