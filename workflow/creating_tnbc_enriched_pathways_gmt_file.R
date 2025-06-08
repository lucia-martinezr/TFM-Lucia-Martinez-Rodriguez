
'Creación del archivo gmt con los enriched pathways de los non-basal TNBC según Aine et al. 2025.'



library(msigdbr)
library(dplyr)

# --- CONFIGURACIÓN ---
non_basal_pathways_to_include <- c(
  # KEGG
  "hsa01212", "hsa00640", "hsa04146", "hsa00280", "hsa01240", "hsa00130",
  "hsa01040", "hsa00480", "hsa00350", "hsa00071", "hsa00061", "hsa01250",
  "hsa04141", "hsa00520", "hsa00051", "hsa00983", "hsa00040", "hsa03320",
  "hsa00100", "hsa00410", "hsa00360", "hsa00062", "hsa00120",
  # GO MF (Molecular Function)
  "GO:0050660", "GO:0016627", "GO:0016616", "GO:0016614", "GO:0071949",
  "GO:0016628", "GO:0019842", "GO:0016747", "GO:0016853", "GO:0009055",
  "GO:0004032", "GO:0008603",
  # MSigDB Hallmark
  "HALLMARK_ESTROGEN_RESPONSE_EARLY", "HALLMARK_ESTROGEN_RESPONSE_LATE",
  "HALLMARK_ADIPOGENESIS", "HALLMARK_ANDROGEN_RESPONSE",
  "HALLMARK_FATTY_ACID_METABOLISM", "HALLMARK_XENOBIOTIC_METABOLISM",
  "HALLMARK_PEROXISOME", "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_BILE_ACID_METABOLISM",
  # Reactome
  "R-HSA-8978868", "R-HSA-9609507", "R-HSA-156580", "R-HSA-8963691",
  "R-HSA-9033241", "R-HSA-390918", "R-HSA-211859", "R-HSA-8957322",
  "R-HSA-156590", "R-HSA-191273", "R-HSA-196854", "R-HSA-70895"
)

# Nombre del archivo .gmt de salida 
OUTPUT_GMT_FILE_R <- "C:/Users/lulim/OneDrive/Documentos/GitHub/TFM-Lucia-Martinez-Rodriguez/data/tnbc_enriched_pathways_gmt_files/non_basal_enriched_pathways_R.gmt"

# Tipo de identificador de gen a incluir en el archivo .gmt (Opciones: "gene_symbol", "entrez_gene", "ensembl_gene")
GENE_IDENTIFIER_TYPE <- "gene_symbol"

# --- OBTENER DATA de MSigDB---

# Ver todas las colecciones y subcolecciones disponibles para Homo sapiens
collections_df <- msigdbr_collections(db_species = "Hs")
print(collections_df) 

# Definir las colecciones específicas de MSigDB
hallmark_data <- msigdbr(species = "Homo sapiens", collection = "H")
kegg_legacy_data <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG_LEGACY")
kegg_medicus_data <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG_MEDICUS")
kegg_data <- rbind(kegg_legacy_data, kegg_medicus_data)
reactome_data <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME")
go_mf_data <- msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:MF")


# --- CONSTRUIR ARCHIVO .gmt ---
pathway_gmt_list <- list()
found_pathways_count_R <- 0

for (pathway_id_to_find in non_basal_pathways_to_include) {
  genes_for_pathway <- NULL
  
  # Determinar la colección basada en el prefijo/formato del ID
  if (startsWith(pathway_id_to_find, "HALLMARK_")) {
    source_data <- hallmark_data
    filter_col <- "gs_name"
  } else if (startsWith(pathway_id_to_find, "hsa")) {
    source_data <- kegg_data
    filter_col <- "gs_exact_source"
  } else if (startsWith(pathway_id_to_find, "R-HSA-")) {
    source_data <- reactome_data
    filter_col <- "gs_exact_source"
  } else if (startsWith(pathway_id_to_find, "GO:")) {
    source_data <- go_mf_data 
    filter_col <- "gs_exact_source"
  } else {
    warning(paste("Formato de ID no reconocido para:", pathway_id_to_find, "- Se omitirá.", sep=" "))
    next # Saltar al siguiente ID
  }
  
  # Construir la condición de filtro dinámicamente
  filter_condition <- sprintf("%s == '%s'", filter_col, pathway_id_to_find)
  
  genes_for_pathway <- source_data %>%
    filter(eval(parse(text = filter_condition))) %>%
    pull({{GENE_IDENTIFIER_TYPE}}) %>%
    unique() %>%
    sort()
  
  if (length(genes_for_pathway) > 0) {
    pathway_gmt_list[[pathway_id_to_find]] <- list(
      name = pathway_id_to_find, 
      desc = "NA", 
      genes = genes_for_pathway
    )
    found_pathways_count_R <- found_pathways_count_R + 1
  } else {
    warning(paste("Pathway ID '", pathway_id_to_find, 
                  "' no encontrado en la colección MSigDB designada o sin genes.", sep=""))
  }
}

# Nota: pathways hsa01212, hsa01240, hsa00130, hsa00061, hsa01250, hsa04141, hsa00062, GO:0016616 no encontrados.


if (length(pathway_gmt_list) > 0) {
  file_conn_R <- file(OUTPUT_GMT_FILE_R, "w")
  for (pathway_entry in pathway_gmt_list) {
    line_to_write <- paste(
      pathway_entry$name, 
      pathway_entry$desc, 
      paste(pathway_entry$genes, collapse = "\t"), 
      sep = "\t"
    )
    writeLines(line_to_write, file_conn_R)
  }
  close(file_conn_R)
  message(paste("Archivo .gmt '", OUTPUT_GMT_FILE_R, "' creado con ", 
                found_pathways_count_R, "/", length(non_basal_pathways_to_include), 
                " pathways encontrados.", sep=""))
} else {
  message("No se encontraron pathways para escribir en el archivo .gmt.")
}






'Creación del archivo gmt con los enriched pathways de los basal TNBC según Aine et al. 2025.'

# --- CONFIGURACIÓN ---

basal_pathways_to_include <- c(
  # KEGG (Basal1)
  "hsa05415", "hsa04714", "hsa04723", "hsa05208",  
  "hsa04932", "hsa05022", "hsa04922", 
  "hsa05133", "hsa01200", #"hsa00190", "hsa05016", "hsa05012","hsa00630","hsa05020","hsa04310" Estos se quitan porque sólo encuentran 1 gen en la ruta y a problemas en Corradjust
  # GO MF (Basal1)
  "GO:0015453", "GO:0008137", "GO:0050136", "GO:0003954", "GO:0008955", 
  "GO:0016651", "GO:0009055", "GO:0015399", "GO:0022804", 
  # Reactome (Basal1)
  "R-HSA-611105", "R-HSA-163200", "R-HSA-1428517", "R-HSA-6799198", "R-HSA-210991",
  # MSigDB Hallmark (Basal2)
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "HALLMARK_HYPOXIA",
  "HALLMARK_UV_RESPONSE_DN", "HALLMARK_TGF_BETA_SIGNALING",
  "HALLMARK_IL2_STAT5_SIGNALING", "HALLMARK_ANGIOGENESIS",
  # MSigDB Hallmark (Basal3)
  "HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_ALLOGRAFT_REJECTION",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_COMPLEMENT",
  "HALLMARK_INFLAMMATORY_RESPONSE", # HALLMARK_IL2_STAT5_SIGNALING ya está
  "HALLMARK_IL6_JAK_STAT3_SIGNALING", "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_KRAS_SIGNALING_UP"
)
basal_pathways_to_include <- unique(basal_pathways_to_include) 

# Nombre del archivo .gmt de salida 
OUTPUT_GMT_FILE_BASAL <- "C:/Users/lulim/OneDrive/Documentos/GitHub/TFM-Lucia-Martinez-Rodriguez/data/tnbc_enriched_pathways_gmt_files/basal_enriched_pathways_R.gmt"

# Tipo de identificador de gen a incluir en el archivo .gmt
GENE_IDENTIFIER_TYPE <- "gene_symbol"

# --- OBTENER DATA de MSigDB ---

# Colecciones de MSigDB
hallmark_data <- msigdbr(species = "Homo sapiens", collection = "H")
kegg_legacy_data <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG_LEGACY")
kegg_medicus_data <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG_MEDICUS")
kegg_data <- rbind(kegg_legacy_data, kegg_medicus_data) %>% 
  distinct(gs_name, gs_exact_source, {{GENE_IDENTIFIER_TYPE}}, .keep_all = TRUE)
reactome_data <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME")
go_mf_data <- msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:MF")


# --- CONSTRUIR ARCHIVO .gmt ---
pathway_gmt_list_basal <- list()
found_pathways_count_basal <- 0

for (pathway_id_to_find in basal_pathways_to_include) {
  genes_for_pathway <- NULL
  
  if (startsWith(pathway_id_to_find, "HALLMARK_")) {
    source_data <- hallmark_data
    filter_col <- "gs_name"
  } else if (startsWith(pathway_id_to_find, "hsa")) {
    source_data <- kegg_data
    filter_col <- "gs_exact_source"
  } else if (startsWith(pathway_id_to_find, "R-HSA-")) {
    source_data <- reactome_data
    filter_col <- "gs_exact_source"
  } else if (startsWith(pathway_id_to_find, "GO:")) {
    # Aquí asumo que todos los GO son MF.
    source_data <- go_mf_data 
    filter_col <- "gs_exact_source"
  } else {
    warning(paste("Formato de ID no reconocido para:", pathway_id_to_find, "- Se omitirá.", sep=" "))
    next
  }
  
  filter_condition <- sprintf("%s == '%s'", filter_col, pathway_id_to_find)
  
  genes_for_pathway <- source_data %>%
    filter(eval(parse(text = filter_condition))) %>%
    pull({{GENE_IDENTIFIER_TYPE}}) %>%
    unique() %>%
    sort()
  
  if (length(genes_for_pathway) > 0) {
    pathway_gmt_list_basal[[pathway_id_to_find]] <- list(
      name = pathway_id_to_find, 
      desc = "NA", 
      genes = genes_for_pathway
    )
    found_pathways_count_basal <- found_pathways_count_basal + 1
  } else {
    warning(paste("Pathway ID '", pathway_id_to_find, 
                  "' no encontrado en la colección MSigDB designada o sin genes.", sep=""))
  }
}

# Nota: pathways hsa05415, hsa04714, hsa04723, hsa05208, hsa04932, hsa05022, hsa04922, hsa05133, hsa01200, GO:0050136, GO:0008955, R-HSA-163200 no encontrados.



if (length(pathway_gmt_list_basal) > 0) {
  file_conn_R_basal <- file(OUTPUT_GMT_FILE_BASAL, "w")
  for (pathway_entry in pathway_gmt_list_basal) {
    line_to_write <- paste(
      pathway_entry$name, 
      pathway_entry$desc, 
      paste(pathway_entry$genes, collapse = "\t"), 
      sep = "\t"
    )
    writeLines(line_to_write, file_conn_R_basal)
  }
  close(file_conn_R_basal)
  message(paste("Archivo .gmt '", OUTPUT_GMT_FILE_BASAL, "' creado con ", 
                found_pathways_count_basal, "/", length(basal_pathways_to_include), 
                " pathways encontrados.", sep=""))
} else {
  message("No se encontraron pathways para escribir en el archivo .gmt para Basal.")
}


