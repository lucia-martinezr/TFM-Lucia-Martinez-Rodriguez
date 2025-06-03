# Creation of a gene dictionary

# Cargar/instalar paquetes necesarios usando pacman
if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman")
}
pacman::p_load(biomaRt, dplyr, readr)

# 1. Conectar con Ensembl usando biomaRt
cat("Conectando a Ensembl...\n")
# Usar una versión específica de Ensembl puede ser bueno para la reproducibilidad
# Para encontrar la versión más reciente o una específica, puedes visitar archive.ensembl.org
# O simplemente usar la versión actual por defecto:
# ensembl_mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
# O especificar un host de archivo si la conexión directa falla o para versiones antiguas:
ensembl_mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                        dataset = "hsapiens_gene_ensembl",
                        host = "https://www.ensembl.org") # o "https://dec2023.archive.ensembl.org" para una versión específica

cat("Conexión establecida.\n")

# 2. Obtener la información de los genes
# Atributos que queremos:
# - ensembl_gene_id: ID de Ensembl (sin versión)
# - external_gene_name: Símbolo del gen (lo que queremos como 'feature_name')
# - entrezgene_id: ID de Entrez (puede ser útil)
# - gene_biotype: Tipo de gen (para filtrar por 'protein_coding' si es necesario)
cat("Obteniendo datos de los genes desde biomaRt...\n")
gene_annotations <- getBM(
  attributes = c("ensembl_gene_id", 
                 "external_gene_name", 
                 "entrezgene_id", 
                 "gene_biotype"),
  mart = ensembl_mart
)
cat("Datos de genes obtenidos. Número total de entradas:", nrow(gene_annotations), "\n")

# Inspeccionar las primeras filas y la estructura
print(head(gene_annotations))
str(gene_annotations)

# 3. Limpieza y formateo (opcional, pero recomendado)
# Renombrar columnas para mayor claridad si es necesario (los nombres actuales son buenos)
# Los nombres que da biomaRt ya son bastante descriptivos.
# Si external_gene_name está vacío para algunos, podríamos querer manejarlo (ej. reemplazar con ensembl_gene_id)
# Por ahora, los mantendremos tal cual.

# Filtrar por 'protein_coding' si solo te interesan esos (como mencionaste)
# Descomenta la siguiente línea si quieres filtrar:
# gene_annotations_protein_coding <- dplyr::filter(gene_annotations, gene_biotype == "protein_coding")
# cat("Número de entradas de 'protein_coding':", nrow(gene_annotations_protein_coding), "\n")
# print(head(gene_annotations_protein_coding))
# Si filtras, asegúrate de usar 'gene_annotations_protein_coding' en el paso de guardado.

# Por ahora, guardaremos todas las anotaciones obtenidas.
# Si quieres filtrar, cambia 'gene_annotations' por 'gene_annotations_protein_coding' en la línea readr::write_tsv

# 4. Guardar la tabla de anotaciones en un formato legible por Python (CSV o TSV)
# TSV (Tab Separated Values) suele ser mejor si los nombres de los genes pueden tener comas.
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

# Desconectar (opcional, biomaRt no mantiene una conexión persistente de esa manera)
detach("package:biomaRt", unload = TRUE) # Pacman podría manejar esto o puedes dejarlo cargado.
