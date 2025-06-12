
# Inputpaths
metadata_inputpath <- "C:/Users/lulim/OneDrive/Documentos/GitHub/TFM-Lucia-Martinez-Rodriguez/data/ext_data/clin_tnbc.rds" #Pregunta: metadata es datos clínicos?
counts_inputpath <- "C:/Users/lulim/OneDrive/Documentos/GitHub/TFM-Lucia-Martinez-Rodriguez/data/corradjust_clean_data/df_data_complete_clean.csv"


# Outputpaths
outputdir <- "C:/Users/lulim/OneDrive/Documentos/GitHub/TFM-Lucia-Martinez-Rodriguez/data/clean_data_partitions/"

# Load data
metadata <- readRDS(metadata_inputpath)
counts <- read.csv(counts_inputpath, header = TRUE, sep = ",")

# Modificar formato de counts
rownames(counts) <- counts[[1]]  
counts <- counts[ , -1]
# name_map <- setNames(dict$ensembl_gene_id, dict$external_gene_name)
# colnames(counts) <- sapply(colnames(counts), function(x) ifelse(x %in% names(name_map), name_map[x], x))

# Make names. Esto garantiza que los nombres de las columnas de metadata no contengan caracteres que puedan causar problemas al trabajar con ellos en R. Los hace únicos añadiendo sufijos en caso necesario.
names(metadata) <- make.names(names(metadata))

# Define partitions
smpSize <- floor(0.8 * nrow(metadata))
set.seed(55)
trainIdx <- sample(seq_len(nrow(metadata)), size = smpSize)

# Split into train and test
metadata_train <- metadata[trainIdx, ]
metadata_test <- metadata[-trainIdx, ]

counts_train <- counts[trainIdx, ] # 572 obs 
counts_test <- counts[-trainIdx, ] # 143 obs

# # Save partitions
# saveRDS(metadata_train, file.path(outputdir, "clean_metadata_train.rds"))
# saveRDS(metadata_test, file.path(outputdir, "clean_metadata_test.rds"))
# saveRDS(counts_train, file.path(outputdir, "clean_counts_train.rds"))
# saveRDS(counts_test, file.path(outputdir, "clean_counts_test.rds"))