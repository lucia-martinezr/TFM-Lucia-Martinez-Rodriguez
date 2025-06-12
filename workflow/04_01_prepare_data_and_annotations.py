# 01_prepare_data_and_annotations.py

# %% --- Importo las librerías necesarias ---
import pandas as pd
import os
import numpy as np

# --- Directorio de salida para los datos procesados ---
processed_data_dir = "/home/lucia2024/TFM/processed_data"
os.makedirs(processed_data_dir, exist_ok=True)

# %% --- 1. Importar datos "crudos" (vienen del script 00_preprocesing.R)---
print("Cargando datos crudos...")
# FPKM counts data
df_counts = pd.read_csv(
    "/home/lucia2024/TFM/ext_data/scanb_tnbc.tsv",
    sep="\t", index_col=0
)

# Clinical data
df_clin = pd.read_csv(
    "/home/lucia2024/TFM/ext_data/clin_tnbc.tsv",
    sep="\t", index_col=0
)

# Setting df_counts index
df_counts.index = df_clin.index
print("Datos crudos cargados y alineados por índice.")
# print("Clinical Data (primeras filas):\n", df_clin.head())
print("Counts Data (primeras filas):\n", df_counts.head())

# Normalice with log2 transformation
print("Normalizando datos con transformación log2...")
df_counts = np.log2(df_counts + 1)
print("Datos de conteo normalizados (primeras filas):\n", df_counts.head())


# %% --- 2. Cargar anotaciones de Ensembl y crear df_features ---
print("\nCargando anotaciones de Ensembl y creando df_features...")
annotation_file_path = "/home/lucia2024/TFM/ext_data/ensembl_hsapiens_gene_annotations.tsv"
try:
    gene_annotation_df = pd.read_csv(annotation_file_path, sep='\t')
except FileNotFoundError:
    print(f"Error: Archivo de anotación no encontrado en '{annotation_file_path}'")
    print("¿Has ejecutado el script de R para generar este archivo? ¿Es la ruta correcta?")
    exit()

# print("Tabla de anotaciones cargada (primeras filas):\n", gene_annotation_df.head())

# Columnas esperadas del archivo TSV
ensembl_id_col = "ensembl_gene_id"    # Contiene los IDs SIN versión
gene_symbol_col = "external_gene_name"
gene_biotype_col = "gene_biotype"

# Creo los diccionarios de mapeo (ID Ensembl SIN versión -> Símbolo/Biotipo)
ensembl_to_symbol_dict = pd.Series(
    gene_annotation_df[gene_symbol_col].values,
    index=gene_annotation_df[ensembl_id_col]
).to_dict()

ensembl_to_biotype_dict = pd.Series(
    gene_annotation_df[gene_biotype_col].values,
    index=gene_annotation_df[ensembl_id_col]
).to_dict()

# Construyo df_features
# feature_id en df_features será el ID original CON versión de las columnas de df_counts
feature_ids_from_columns_with_version = df_counts.columns.tolist()
feature_names_list = []
feature_types_list = []

for original_id_col_with_version in feature_ids_from_columns_with_version:
    id_no_version = original_id_col_with_version.split('.')[0]
    
    symbol = ensembl_to_symbol_dict.get(id_no_version)
    feature_names_list.append(symbol if symbol and pd.notna(symbol) else id_no_version) 
    # Si no se encuentra un valor en el diccionario (ensembl_to_symbol_dict) para un identificador específico (o si el valor encontrado no es válido), se utilizará el id_no_version como valor alternativo.
    biotype = ensembl_to_biotype_dict.get(id_no_version)
    feature_types_list.append(biotype if biotype and pd.notna(biotype) else None) # Es None si no se encuentra biotipo

df_features = pd.DataFrame({
    "feature_id": feature_ids_from_columns_with_version, # ID original CON versión
    "feature_name": feature_names_list,                  # Símbolo (o id_no_version)
    "feature_type": feature_types_list                   # Biotipo (o None)
})

# Establezco feature_id como índice
df_features = df_features.set_index("feature_id")

print("df_features creado (primeras filas):\n", df_features.head())
print(f"Dimensiones de df_features: {df_features.shape}")

# %% --- 3. Dividir datos de conteo en entrenamiento y prueba ---
# Estos DataFrames (df_counts_train_orig, df_counts_test_orig)
# todavía tienen feature_id CON VERSIÓN como nombres de columna.
print("\nDividiendo datos de conteo en entrenamiento y prueba...")
df_counts_train_orig = df_counts.iloc[::2]
df_counts_test_orig = df_counts.iloc[1::2]

print("Training Set Counts Data (primeras filas):\n", df_counts_train_orig.head())
print(f"Dimensiones de df_counts_train_orig: {df_counts_train_orig.shape}")
print("Test Set Counts Data (primeras filas):\n", df_counts_test_orig.head())
print(f"Dimensiones de df_counts_test_orig: {df_counts_test_orig.shape}")

# %% --- 4. Guardar DataFrames para el siguiente script ---
print("\nGuardando DataFrames procesados...")
df_features_path = os.path.join(processed_data_dir, "df_features.pkl")
df_counts_train_orig_path = os.path.join(processed_data_dir, "df_counts_train_orig.pkl")
df_counts_test_orig_path = os.path.join(processed_data_dir, "df_counts_test_orig.pkl")

df_features.to_pickle(df_features_path)
df_counts_train_orig.to_pickle(df_counts_train_orig_path)
df_counts_test_orig.to_pickle(df_counts_test_orig_path)

print(f"DataFrame df_features guardado en: {df_features_path}")
print(f"DataFrame df_counts_train_orig guardado en: {df_counts_train_orig_path}")
print(f"DataFrame df_counts_test_orig guardado en: {df_counts_test_orig_path}")

print("\n--- Script 01 completado ---")
# %%
