# Análisis metilómico en diabetes tipo 2 en páncreas humano (grupo13-IPC)
Este repositorio documenta un análisis metilómico de ADN en muestras de páncreas humano asociadas a diabetes tipo 2, empleando microarreglos Illumina. El código aplica procesos de normalización y estadística bayesiana empírica para identificar genes o regiones diferencialmente metiladas, permitiendo explorar biomarcadores epigenéticos con potencial utilidad.

## Descripción del proyecto
Este repositorio contiene un análisis metilómico en muestras de páncreas humano asociadas a diabetes tipo 2, basado en datos de microarreglos Illumina.

## Propósito
Organizar y documentar un flujo de trabajo reproducible para el análisis metilómico del ADN y la identificación de patrones de metilación diferencial.

## Objetivos
- Organizar los datos y scripts del proyecto
- Aplicar un flujo de preprocesamiento y normalización
- Identificar genes o regiones diferencialmente metiladas
- Generar resultados reproducibles y debidamente documentados

## Estructura del repositorio
- `data/`: datos crudos y procesados del proyecto
- `scripts/`: scripts de preprocesamiento y análisis
- `results/`: resultados generados, incluyendo tablas y figuras
- `docs/`: documentación adicional del flujo de trabajo
- `notebooks/`: análisis exploratorio y pruebas preliminares

## Instrucciones de uso
1. Clonar el repositorio.
2. Ubicar los datos de entrada en la carpeta correspondiente dentro de `data/`.
3. Ejecutar los scripts de preprocesamiento ubicados en `scripts/`.
4. Ejecutar los scripts de análisis diferencial.
5. Revisar los resultados generados en la carpeta `results/`.

## Interpretación biológica
A partir del análisis metilómico realizado, se identifican diferencias en los patrones de metilación entre muestras control y aquellas asociadas a diabetes tipo 2. Se observa que algunos genes relacionados con la función pancreática presentan cambios en su nivel de metilación, lo cual podría influir en su regulación y en procesos metabólicos relevantes. Asimismo, se evidencian posibles alteraciones en genes asociados a procesos inflamatorios, lo cual resulta consistente con el papel de la inflamación en el desarrollo de la diabetes tipo 2. Estos hallazgos permiten identificar posibles biomarcadores epigenéticos con potencial relevancia para el estudio de la enfermedad.

## Herramientas utilizadas
- Git y GitHub
- Bash
- R/Bioconductor
- Microarreglos Illumina
- Python

## Conclusión
Este proyecto permitió estructurar y documentar un flujo de trabajo orientado al análisis metilómico en el contexto de la diabetes tipo 2. A través de la organización de datos y la aplicación de procesos de normalización y análisis diferencial, se logró simular un enfoque bioinformático utilizado en estudios reales. El repositorio facilita la reproducibilidad del análisis y contribuye a la comprensión de cómo los cambios epigenéticos pueden estar relacionados con el desarrollo de enfermedades metabólicas.

## Licencia
Este proyecto se distribuye bajo la licencia GPL-2.0.
