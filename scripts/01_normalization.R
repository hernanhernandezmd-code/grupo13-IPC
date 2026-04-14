#!/usr/bin/env Rscript
# =============================================================================
# Script: 01_normalization.R
# Tarea: Normalización de datos de microarreglos Illumina (diabetes tipo 2)
# Métodos: Corrección de fondo + Normalización quantile
# Autor: Equipo de Bioinformática G13
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
# Verificar la instalación del gestor de paquetes principal de Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Instalamos el paquete principal y el manifiesto genómico para los 27k
if (!require("minfi")) BiocManager::install("minfi")
if (!require("IlluminaHumanMethylation27kmanifest")) {
  BiocManager::install("IlluminaHumanMethylation27kmanifest")
}

library(minfi)                                    # Carga algorítmica

# Definimos la ruta donde se deben alojar los archivos IDAT crudos
ruta_idat <- "grupo13-IPC/data/idat_files"        # Directorio base

# Leemos recursivamente las intensidades de fluorescencia dicromáticas
objetos_rg <- read.metharray.exp(base = ruta_idat)# Lectura primaria

# Ejecutamos la normalización robusta con el algoritmo de cuantiles
objeto_norm <- preprocessQuantile(objetos_rg)     # Ajuste estadístico

# Extraemos directamente los valores M estabilizados y homocedásticos
valores_m_minfi <- getM(objeto_norm)              # Transformación M

# Parametrizamos la exportación hacia el disco duro asegurando la ruta
ruta_out <- "grupo13-IPC/data/minfi_GSE21232.RData"
save(valores_m_minfi, file = ruta_out)            # Almacenamiento
