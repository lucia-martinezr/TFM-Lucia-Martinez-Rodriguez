# Feature selection

# Inputpaths
feature_scores_basal_path <- "C:/Users/lulim/OneDrive/Documentos/GitHub/TFM-Lucia-Martinez-Rodriguez/data/corradjust_clean_data/feature_scores_test_basal_tnbc_enriched_pathways.csv"

# Outputpaths
outputpath <- "C:/Users/lulim/OneDrive/Documentos/GitHub/TFM-Lucia-Martinez-Rodriguez/data/clean_tnbc_signature/"


# Load data
fs_basal <- read.csv(feature_scores_basal_path, header = TRUE, sep = ',')

# Select top ranked genes (highest enrichment + lowest padj)
fs_basal_sorted <- fs_basal[order(-fs_basal$enrichment, fs_basal$padj), ]

# Select 56 first genes
top_49_genes <- head(fs_basal_sorted$feature_id, 49)
print(top_49_genes)

clean_signature <- data.frame(SYMBOL = top_49_genes)

# # Store signature
# saveRDS(clean_signature, file.path(outputpath, "clean_signature.rds"))
