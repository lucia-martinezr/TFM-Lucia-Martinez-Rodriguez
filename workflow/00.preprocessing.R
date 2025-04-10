
library(dplyr)
library(readxl)
library(DeSeq2)

# Carga de datos clínicos
clin <- read_excel('C:/Users/lulim/OneDrive/Escritorio/SCAN B Data/SCANB_clinical_data_original_41523_2022_465_MOESM2_ESM.xlsx')


# Carga de datos de conteo (FPKM)
load('C:/Users/lulim/OneDrive/Escritorio/SCAN B Data/SCANB.9206.genematrix_noNeg.Rdata')
ls()
dim(SCANB.9206.genematrix_noNeg) # 19675  9206


# Filtrado de datos. Nos quedamos sólo con muestras de TNBC.
clin_tnbc <- clin %>% filter(ClinGroup == "TNBC")

unique_patients <- unique(clin_tnbc$Patients)

# id_tnbc <- clin_tnbc[[1]]

scanb <- as.data.frame(SCANB.9206.genematrix_noNeg)

scanb_tnbc <- scanb %>% select(all_of(id_tnbc))

scanb_tnbc <- as.data.frame(t(scanb_tnbc))

# Filtrado de datos: eliminación de los genes con un FPKM de 0 en más del 20% de muestras.

  # zero_sum_columns <- sum(colSums(scanb_tnbc) == 0)  ### 86 columnas tienen exactamente 0 FPKM en total.

threshold <- nrow(scanb_tnbc) * 0.2

filtered_scanb_tnbc <- scanb_tnbc %>%
  select(where(~ sum(. == 0) <= threshold)) # Se me reduce en 3971 columnas.

# Intersección de los genes con los conjuntos de datos de validación externa.


# Variance Stabilizing Transformation (VST)


# Guardar los datos
saveRDS(scanb_tnbc, 'C:/scanb_tnbc.rds')

saveRDS(clin_tnbc, 'C:/clin_tnbc.rds')



