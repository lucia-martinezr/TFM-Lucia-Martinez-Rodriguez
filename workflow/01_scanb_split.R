# source("requirements.R") # Carga y ejecuta el script requirements.R

# Inputpaths
metadata_inputpath <- "C:/Users/lulim/OneDrive/Documentos/GitHub/TFM-Lucia-Martinez-Rodriguez/data/ext_data/clin_tnbc.rds" #Pregunta: metadata es datos clínicos?
counts_inputpath <- "C:/Users/lulim/OneDrive/Documentos/GitHub/TFM-Lucia-Martinez-Rodriguez/data/ext_data/scanb_tnbc.rds"


# Outputpaths
outputdir <- "C:/Users/lulim/OneDrive/Documentos/GitHub/TFM-Lucia-Martinez-Rodriguez/data/data_partitions/"

# Load data
metadata <- readRDS(metadata_inputpath)
counts <- readRDS(counts_inputpath)

# Make names. Esto garantiza que los nombres de las columnas de metadata no contengan caracteres que puedan causar problemas al trabajar con ellos en R. Los hace únicos añadiendo sufijos en caso necesario.
names(metadata) <- make.names(names(metadata))

# Define partitions
smpSize <- floor(0.9 * nrow(metadata))
set.seed(55)
trainIdx <- sample(seq_len(nrow(metadata)), size = smpSize)

# Split into train and test
metadata_train <- metadata[trainIdx, ]
metadata_test <- metadata[-trainIdx, ]

counts_train <- counts[trainIdx, ] # 786 obs
counts_test <- counts[-trainIdx, ] # 88 obs

# Save partitions
saveRDS(metadata_train, file.path(outputdir, "metadata_train.rds"))
saveRDS(metadata_test, file.path(outputdir, "metadata_test.rds"))
saveRDS(counts_train, file.path(outputdir, "counts_train.rds"))
saveRDS(counts_test, file.path(outputdir, "counts_test.rds"))