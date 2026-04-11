#!/usr/bin/env Rscript
# =============================================================================
# Script: 01_normalization.R
# Tarea: Normalización de datos de microarreglos Illumina (diabetes tipo 2)
# Métodos: Corrección de fondo (normexp) + Normalización quantile
# Autor: Equipo de Bioinformática
# Fecha: 2025-04-11 (simulación)
# =============================================================================

# -----------------------------------------------------------------------------
# 1. Configuración inicial y carga de librerías
# -----------------------------------------------------------------------------

# Limpiar entorno
rm(list = ls())

# Cargar paquetes necesarios (se asume que están instalados)
suppressPackageStartupMessages({
  library(limma)          # para backgroundCorrect y normalizeBetweenArrays
  library(preprocessCore) # para normalize.quantiles (alternativa)
  library(dplyr)          # para manipulación de datos
  library(readr)          # para lectura rápida
})

# Establecer semilla para reproducibilidad (si se usa algún muestreo)
set.seed(123)

# Definir rutas (usando here::here o relativas)
data_dir <- "data"
raw_file <- file.path(data_dir, "raw_expression.txt")
output_rds <- file.path(data_dir, "expresion_normalizada.rds")
output_csv <- file.path(data_dir, "expresion_normalizada.csv")
log_file <- file.path(data_dir, "normalization_log.txt")

# Crear directorio data si no existe
if (!dir.exists(data_dir)) {
  dir.create(data_dir, recursive = TRUE)
}

# Iniciar archivo de log
sink(log_file, split = TRUE)
cat("===== LOG DE NORMALIZACIÓN =====\n")
cat("Fecha y hora:", as.character(Sys.time()), "\n\n")

# -----------------------------------------------------------------------------
# 2. Carga de datos crudos
# -----------------------------------------------------------------------------
# Simulación: en un caso real se leería un archivo como:
#   raw_data <- read.delim(raw_file, row.names = 1, check.names = FALSE)
# Aquí generamos datos sintéticos realistas para demostración.

cat("Cargando datos crudos...\n")
set.seed(456)
n_genes <- 15000
n_samples <- 24  # 12 controles, 12 DM2
# Simular intensidades crudas (valores positivos, con ruido)
raw_data <- matrix(
  rgamma(n_genes * n_samples, shape = 2, scale = 10),
  nrow = n_genes,
  ncol = n_samples,
  dimnames = list(
    paste0("GEN", sprintf("%05d", 1:n_genes)),
    paste0("Muestra", 1:n_samples)
  )
)
# Añadir algunos valores cercanos a cero para simular ruido de fondo
raw_data[sample(length(raw_data), 500)] <- runif(500, 0, 0.5)

cat(sprintf("Dimensiones de datos crudos: %d genes x %d muestras\n", 
            nrow(raw_data), ncol(raw_data)))
cat("Rango de intensidades crudas: [", min(raw_data), ",", 
    max(raw_data), "]\n\n")

# -----------------------------------------------------------------------------
# 3. Corrección de fondo (background correction)
# -----------------------------------------------------------------------------
# Usamos el método 'normexp' de limma, que modela la señal + ruido de fondo.
# Esto evita valores negativos posteriores.

cat("Aplicando corrección de fondo con método 'normexp'...\n")
raw_corrected <- backgroundCorrect(raw_data, method = "normexp", 
                                   offset = 16)  # offset típico para Illumina

# Verificar que no haya valores negativos
if (any(raw_corrected < 0)) {
  warning("Aún hay valores negativos después de backgroundCorrect. Se fijarán a 0.")
  raw_corrected[raw_corrected < 0] <- 0
}

cat("Corrección completada. Nuevo rango: [", min(raw_corrected), ",", 
    max(raw_corrected), "]\n\n")

# -----------------------------------------------------------------------------
# 4. Transformación logarítmica (opcional pero recomendada)
# -----------------------------------------------------------------------------
# Para estabilizar la varianza y hacer los datos más simétricos.
# Usamos log2(x + 1) para manejar ceros.

cat("Aplicando transformación log2(x+1)...\n")
expr_log <- log2(raw_corrected + 1)
cat("Rango después de log: [", min(expr_log), ",", max(expr_log), "]\n\n")

# -----------------------------------------------------------------------------
# 5. Normalización quantile (entre muestras)
# -----------------------------------------------------------------------------
# Método estándar para microarreglos: fuerza la misma distribución de 
# intensidades en todas las muestras.

cat("Ejecutando normalización quantile...\n")
expr_norm <- normalize.quantiles(as.matrix(expr_log))
# Restaurar nombres de filas y columnas
colnames(expr_norm) <- colnames(expr_log)
rownames(expr_norm) <- rownames(expr_log)

cat("Normalización completada. Rango final: [", min(expr_norm), ",", 
    max(expr_norm), "]\n\n")

# Verificación final: no debe haber NAs
if (any(is.na(expr_norm))) {
  stop("Error: Se generaron valores NA durante la normalización.")
}

# Verificar que no haya valores negativos (ahora es muy improbable)
if (any(expr_norm < 0)) {
  warning("Valores negativos detectados después de quantile. Se investigará.")
  # En un caso real aquí se podría aplicar un shift o revisar.
}

# -----------------------------------------------------------------------------
# 6. Guardar resultados
# -----------------------------------------------------------------------------
# Guardar como RDS (para scripts posteriores)
saveRDS(expr_norm, file = output_rds)
cat(sprintf("Datos normalizados guardados en: %s\n", output_rds))

# Guardar como CSV (para inspección externa)
write.csv(expr_norm, file = output_csv, row.names = TRUE)
cat(sprintf("Datos normalizados guardados en CSV: %s\n", output_csv))

# Guardar también un resumen estadístico
summary_stats <- data.frame(
  Estadística = c("Mínimo", "Primer cuartil", "Mediana", "Media", 
                  "Tercer cuartil", "Máximo"),
  Valor = round(c(min(expr_norm), quantile(expr_norm, 0.25), 
                  median(expr_norm), mean(expr_norm), 
                  quantile(expr_norm, 0.75), max(expr_norm)), 4)
)
write.csv(summary_stats, file = file.path(data_dir, "normalization_summary.csv"), 
          row.names = FALSE)

# -----------------------------------------------------------------------------
# 7. Generar un gráfico de diagnóstico (opcional, solo si se puede plotear)
# -----------------------------------------------------------------------------
# Para simulación, creamos un archivo PNG con boxplots antes/después
if (interactive() || capabilities("png")) {
  png(file.path(data_dir, "normalization_boxplots.png"), 
      width = 1200, height = 600, res = 120)
  par(mfrow = c(1, 2))
  boxplot(expr_log, main = "Antes de quantile", las = 2, cex.axis = 0.6)
  boxplot(expr_norm, main = "Después de quantile", las = 2, cex.axis = 0.6)
  dev.off()
  cat("Gráfico de diagnóstico guardado: normalization_boxplots.png\n")
}

# -----------------------------------------------------------------------------
# 8. Finalización
# -----------------------------------------------------------------------------
cat("\nNormalización completada exitosamente.\n")
cat("Resumen de la matriz final (primeras 6 filas y columnas):\n")
print(expr_norm[1:6, 1:6])
cat("\n")
sink()

# Mensaje en consola
cat("Proceso finalizado. Revisar el log en:", log_file, "\n")
