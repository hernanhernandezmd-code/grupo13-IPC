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
# Verificar la instalación del gestor de paquetes principal de Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

ruta <- "grupo13-IPC/data/expresion_GSE21232.RData"  # Ruta de carga
load(ruta)

# Descartamos sondas con NAs para asegurar la integridad biológica
datos_filt <- datos_expresion[complete.cases(datos_expresion), ] # Filtro

# Convertimos proporciones Beta a valores M para habilitar estadística
valores_m <- log2(datos_filt / (1 - datos_filt))                 # Normaliza

ruta_out <- "grupo13-IPC/data/norm_GSE21232.RData"   # Ruta exportación
save(valores_m, file = ruta_out)

# Mensaje en consola
cat("Proceso finalizado. Revisar el log en:", log_file, "\n")
