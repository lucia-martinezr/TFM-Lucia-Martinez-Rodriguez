# Creation of a gene dictionary

# Cargar/instalar paquetes necesarios usando pacman
if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman")
}
pacman::p_load(biomaRt, dplyr, readr)

# Conectar con Ensembl usando biomaRt
cat("Conectando a Ensembl...\n")

ensembl_mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                        dataset = "hsapiens_gene_ensembl",
                        host = "https://www.ensembl.org") 
cat("Conexión establecida.\n")

# 2. Obtener la información de los genes (ensembl id, símbolo, entrez id, biotype)
cat("Obteniendo datos de los genes desde biomaRt...\n")
gene_annotations <- getBM(
  attributes = c("ensembl_gene_id", 
                 "external_gene_name", 
                 "entrezgene_id", 
                 "gene_biotype"),
  mart = ensembl_mart
)
cat("Datos de genes obtenidos. Número total de entradas:", nrow(gene_annotations), "\n")

print(head(gene_annotations))
str(gene_annotations)

# Guardar la tabla de anotaciones en formato TSV y RDS
output_dir <- "C:/Users/lulim/OneDrive/Documentos/GitHub/TFM-Lucia-Martinez-Rodriguez/data/dictionary" # Directorio de salida deseado
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
output_filename_tsv <- file.path(output_dir, "ensembl_hsapiens_gene_annotations.tsv")
output_filename_rds <- file.path(output_dir, "ensembl_hsapiens_gene_annotations.rds") # También guardamos como RDS por si acaso

cat("Guardando anotaciones en:", output_filename_tsv, "\n")
readr::write_tsv(gene_annotations, output_filename_tsv)

cat("Guardando anotaciones también en formato RDS:", output_filename_rds, "\n")
saveRDS(gene_annotations, file = output_filename_rds)

cat("Proceso completado. El archivo de anotaciones se ha guardado.\n")

# Desconectar 
detach("package:biomaRt", unload = TRUE) 
