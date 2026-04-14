# Verificar la instalación del gestor de paquetes principal
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!require("limma")) BiocManager::install("limma")
library(limma)                                      # Herramienta base

ruta_norm <- "grupo13-IPC/data/norm_GSE21232.RData"
load(ruta_norm)                                     # Carga métricas

ruta_fen <- "grupo13-IPC/data/phenodata_diabetes.csv"
fenotipo <- read.csv(ruta_fen)                      # Carga clínica

condicion <- factor(fenotipo$group.ch1)
diseño <- model.matrix(~ 0 + condicion)             # Matriz lineal
colnames(diseño) <- levels(condicion)               # Nombres columnas

ajuste <- lmFit(datos_norm, diseño)                 # Regresión
contr <- makeContrasts(CASE - CTL, levels = diseño) # Diferencial
ajuste_c <- contrasts.fit(ajuste, contr)
ajuste_b <- eBayes(ajuste_c)                        # Moderación Bayes

ruta_fig <- "grupo13-IPC/figures"                   # Ruta de figuras
if (!dir.exists(ruta_fig)) dir.create(ruta_fig)     # Crear directorio

ruta_pdf <- file.path(ruta_fig, "volcano_limma.pdf")
pdf(ruta_pdf, width = 8, height = 6)                # Dispositivo PDF

# Invocamos la función especializada sobre el objeto estadístico
volcanoplot(ajuste_b, coef = 1, highlight = 15,     # Gráfico nativo
            names = rownames(ajuste_b),
            main = "Volcán de Metilación (limma)")

dev.off()                                           # Guardar vector
