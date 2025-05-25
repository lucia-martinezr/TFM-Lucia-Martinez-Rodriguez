# Feature selection with FCBF
# source("R/utils.R") Dónde está este archivo.
# source("requirements.R")

# devtools::install_github("lubianat/FCBF")
# BiocManager::install("AnnotationDbi")
# BiocManager::install("org.Hs.eg.db")
library(FCBF)
library(org.Hs.eg.db)
library(AnnotationDbi)

# Inputpaths
metadata_inputpath <- "C:/Users/lulim/OneDrive/Documentos/GitHub/TFM-Lucia-Martinez-Rodriguez/data/data_partitions/metadata_train.rds"
counts_inputpath <- "C:/Users/lulim/OneDrive/Documentos/GitHub/TFM-Lucia-Martinez-Rodriguez/data/data_partitions/counts_train.rds"

# Outputpaths
outputpath <- "C:/Users/lulim/OneDrive/Documentos/GitHub/TFM-Lucia-Martinez-Rodriguez/data/tnbc_signature/fcbf_signature_annot.rds"

# Arguments
min_su <- 0.0025

# Load data
counts_train <- readRDS(counts_inputpath)
metadata_train <- readRDS(metadata_inputpath)

# Save gene names
gene_names <- colnames(counts_train)

# Discretize data
discretized_data <- discretize_exprs(
  expression_table = counts_train,
  number_of_bins = 3,
  method = "varying_width",
  alpha = 1,
  centers = 3,
  min_max_cutoff = 0.25,
  progress_bar = TRUE
)

# Set target vector 
target_vector <- metadata_train$OS_event

# Run FCBF 
selected_features <- fcbf(
  feature_table = discretized_data,
  target_vector = target_vector,
  minimum_su = 0.0025,
  n_genes_selected_in_first_step = NULL,
  verbose = FALSE,
  samples_in_rows = TRUE, 
  balance_classes = FALSE
)

# Con minimum_su 0.25 no encuentra nada.
# Con 0.0025, que es el threshold que pone Jose, me encuentra 17 características.

# Retrieve Ensembl IDs
rn<-gene_names[selected_features$index]
rownames(selected_features)<-rn 
#Nota: ASUMO que los índices coinciden con el orden original de los IDs

# Erase Ensembl transcript information
rownames(selected_features) <- sub("\\..*", "", rownames(selected_features)) 

# Map Ensembl IDs to gene symbols using org.Hs.eg.db
signature <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = rownames(selected_features), #Nota: mis genes estaban en rows, por eso no se puede usar names()
  columns = c("ENSEMBL", "SYMBOL"),
  keytype = "ENSEMBL",
  multiVals = "first"
)

# Store signature
saveRDS(signature, file = file.path(outputpath))
