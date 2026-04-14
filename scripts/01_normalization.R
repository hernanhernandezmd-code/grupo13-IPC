#!/usr/bin/env Rscript
# =============================================================================
# Script: 01_normalization.R
# Tarea: Normalización de datos de microarreglos Illumina (diabetes tipo 2)
# Métodos:  Normalización quantile
# Autor: Equipo de Bioinformática
# Fecha: 2026-04-11
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
# Verificamos la instalación del gestor BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Instalamos y cargamos la librería estadística base
if (!require("limma")) BiocManager::install("limma")
library(limma)                                      # Dependencia

ruta <- "grupo13-IPC/data/expresion_GSE21232.RData"
load(ruta)                                          # Carga matriz

# Normalizamos directamente las proporciones Beta con cuantiles
datos_norm <- normalizeBetweenArrays(datos_expresion) # Normaliza

ruta_out <- "grupo13-IPC/data/norm_GSE21232.RData"
save(datos_norm, file = ruta_out)                   # Exportación
cat("Proceso finalizado. Revisar el log en:", log_file, "\n")
