
library(dplyr)
library(readxl)

# Carga de datos clínicos
clin <- read_excel('C:/Users/lulim/OneDrive/Escritorio/SCAN B Data/SCANB_clinical_data_original_41523_2022_465_MOESM2_ESM.xlsx')
dim(clin) # 9206 87

# Carga de datos de conteo (FPKM)
load('C:/Users/lulim/OneDrive/Escritorio/SCAN B Data/SCANB.9206.genematrix_noNeg.Rdata')
ls()
dim(SCANB.9206.genematrix_noNeg) # 19675  9206


# Filtrado de datos: Nos quedamos sólo con muestras de TNBC.
clin_tnbc <- clin %>% filter(ClinGroup == "TNBC")

unique_patients <- unique(clin_tnbc$Patient) #715 pacientes únicos

id_tnbc <- clin_tnbc[[1]]

scanb <- as.data.frame(SCANB.9206.genematrix_noNeg)

scanb_tnbc <- scanb %>% select(all_of(id_tnbc))
 # Nos quedan 874 muestras y 19675 genes.

# Filtrado de datos: nos quedamos con solo una muestra por paciente. Eliminación de replicados técnicos.
# Contar la cantidad de ocurrencias de cada valor de la columna Patient
patient_counts <- table(clin_tnbc$Patient)

# Identificar los pacientes que aparecen más de una vez
multiple_patients <- names(patient_counts[patient_counts > 1]) #128 pacientes tienen replicados en el dataset.
repeated_patients <- clin_tnbc[clin_tnbc$Patient %in% multiple_patients, ]

# Uso de FRACTION DUPLICATION como criterio de selección de muestras
mean_duplication <- mean(clin_tnbc$FRACTION_DUPLICATION, na.rm = TRUE)

# Filtrar y seleccionar la fila más cercana a la media para cada paciente
filtered_clin_tnbc <- do.call(rbind, lapply(split(clin_tnbc, clin_tnbc$Patient), function(group) {
  if (nrow(group) > 1) {
    closest_row <- group[which.min(abs(group$FRACTION_DUPLICATION - mean_duplication)), ]
    return(closest_row)
  } else {
    return(group)
  }
}))

# Filtrar scanb_tnbc
filtered_id_tnbc <- filtered_clin_tnbc[[1]]

filtered_scanb_tnbc <-scanb_tnbc %>% select(all_of(filtered_id_tnbc))

# Filtrado de datos: eliminación de los genes con un FPKM de 0 en más del 20% de muestras.

  # zero_sum_columns <- sum(colSums(scanb_tnbc) == 0)  ### 86 columnas tienen exactamente 0 FPKM en total.

threshold <- nrow(scanb_tnbc) * 0.2

scanb_tnbc <- scanb_tnbc %>%
  select(where(~ sum(. == 0) <= threshold)) # Se me reduce en 3971 columnas. Ahora quedan 15704 genes.

# Intersección de los genes con los conjuntos de datos de validación externa.


# Variance Stabilizing Transformation (VST)


# Guardar los datos. CAMBIAR ESTO SI DECIDO FILTRAR!!!!!
saveRDS(scanb_tnbc, 'C:/Users/lulim/OneDrive/Documentos/GitHub/TFM-Lucia-Martinez-Rodriguez/data/ext_data/scanb_tnbc.rds')

saveRDS(clin_tnbc, 'C:/Users/lulim/OneDrive/Documentos/GitHub/TFM-Lucia-Martinez-Rodriguez/data/ext_data/clin_tnbc.rds')

# Prueba de normalidad para 5 genes al azar

set.seed(42)
random_cols <- sample(names(counts), 5)

par(mfrow = c(3, 2))  # Establecer un layout de gráficos
for (col in random_cols) {
  qqnorm(counts[[col]], main = paste("Q-Q plot:", col))
  qqline(counts[[col]], col = "red")
}

# Se obtienen gráficos con colas pesadas a la derecha.

par(mfrow = c(3, 2))  # Establecer un layout de gráficos
for (col in random_cols) {
  hist(counts[[col]], main = paste("Histograma:", col), xlab = "Valores", ylab = "Frecuencia", col = "lightblue", border = "black")
}

# Los datos de expresión génica suelen estar sesgados hacia la izquierda y pueden seguir una 
# distribución log-normal debido a la naturaleza de los datos biológicos y a la prevalencia de valores 
# pequeños.
