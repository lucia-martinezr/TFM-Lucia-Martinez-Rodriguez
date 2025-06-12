# 02_run_corradjust.py

# %%
import pandas as pd
from corradjust import CorrAdjust # Asegúrate de que corradjust esté instalado
import os

# --- Directorio de datos procesados (entrada) y salida para CorrAdjust ---
processed_data_dir = "/home/lucia2024/TFM/processed_data"
corradjust_output_dir = "/home/lucia2024/TFM/out_data/CorrAdjust_ScanbTNBC_PCs_opt" 
os.makedirs(corradjust_output_dir, exist_ok=True)

# %% 
# --- 1. Cargar DataFrames preprocesados ---
print("Cargando DataFrames preprocesados...")
df_features_path = os.path.join(processed_data_dir, "df_features.pkl")
df_counts_train_orig_path = os.path.join(processed_data_dir, "df_counts_train_orig.pkl")
df_counts_test_orig_path = os.path.join(processed_data_dir, "df_counts_test_orig.pkl")

try:
    df_features = pd.read_pickle(df_features_path)
    df_counts_train_orig = pd.read_pickle(df_counts_train_orig_path)
    df_counts_test_orig = pd.read_pickle(df_counts_test_orig_path)
except FileNotFoundError as e:
    print(f"Error al cargar archivos .pkl: {e}")
    print("¿Has ejecutado 01_prepare_data_and_annotations.py primero?")
    exit()

print("DataFrames cargados:")
print(f"  df_features: {df_features.shape}")
print(f"  df_counts_train_orig: {df_counts_train_orig.shape}")
print(f"  df_counts_test_orig: {df_counts_test_orig.shape}")



# %%
# --- 2. Preparación de datos específicos para CorrAdjust ---
print("\nPreparando datos para CorrAdjust...")

# Cargo la anotación original (df_features.pkl)
original_annotation_df = pd.read_pickle(df_features_path)

# 2.1 Creo id_to_name_map (Ensembl ID con versión -> Símbolo de Gen)
id_source_keys = original_annotation_df["feature_id"] if "feature_id" in original_annotation_df.columns else original_annotation_df.index
id_to_name_map = pd.Series(original_annotation_df["feature_name"].values, index=id_source_keys).to_dict()
print(f"id_to_name_map creado con {len(id_to_name_map)} entradas.")

# 2.2 Renombro columnas de datos de expresión a SÍMBOLOS y filtro
def procesar_datos_expresion(df_counts, mapper):
    df_renamed = df_counts.rename(columns=mapper)
    # Filtro por símbolos válidos (no NaN, no "ENSG...", y no vacío) y únicos
    valid_symbol_cols = [
        col for col in df_renamed.columns
        if pd.notna(col) and isinstance(col, str) and not col.startswith("ENSG") and col != ""
    ]
    df_symbols = df_renamed[valid_symbol_cols]
    return df_symbols.loc[:, ~df_symbols.columns.duplicated(keep='first')]

df_data_train_symbols = procesar_datos_expresion(df_counts_train_orig, id_to_name_map)
df_data_test_symbols = procesar_datos_expresion(df_counts_test_orig, id_to_name_map)

print(f"  df_data_train (columnas son símbolos preliminares): {df_data_train_symbols.shape}")
if df_data_train_symbols.empty:
    print("ADVERTENCIA: df_data_train_symbols está vacío.")
    # exit() # Considera salir si es crítico

# 2.3 Preparar df_feature_ann (NOTA: ¡INDEXADO POR SÍMBOLOS PORQUE SINO CORRADJUST DA ERROR!)
# Debe tener columnas 'feature_name', 'feature_type' y estar alineado con los datos de expresión
df_ann_prep = original_annotation_df[["feature_name", "feature_type"]].copy()
df_ann_prep = df_ann_prep[df_ann_prep["feature_type"] == "protein_coding"] # Filtro de tipo
df_ann_prep = df_ann_prep[df_ann_prep["feature_name"].isin(df_data_train_symbols.columns)] # Símbolos presentes en datos
df_ann_prep = df_ann_prep.dropna(subset=["feature_name"]) # Quito NaN en símbolos
df_ann_prep = df_ann_prep[df_ann_prep["feature_name"] != ""] # Quito símbolos vacíos
df_ann_prep = df_ann_prep.drop_duplicates(subset=["feature_name"], keep="first").set_index("feature_name")

if df_ann_prep.empty:
    print("  Error: df_ann_prep está vacío después de filtrar. No hay anotaciones para los símbolos en los datos.")
    # exit() # Considera salir

# Creo df_feature_ann con la estructura final requerida por CorrAdjust
df_feature_ann_final = pd.DataFrame({
    "feature_name": df_ann_prep.index,
    "feature_type": df_ann_prep["feature_type"]
}, index=df_ann_prep.index)
print(f"  df_feature_ann (indexado por símbolos preliminares): {df_feature_ann_final.shape}")

# 2.4 Alineación final de df_data_train y df_feature_ann por SÍMBOLOS comunes
common_symbols = df_data_train_symbols.columns.intersection(df_feature_ann_final.index)

if common_symbols.empty:
    print("  Error: No hay símbolos comunes finales entre datos de expresión y anotaciones.")
    # exit() # Considera salir

df_data_train = df_data_train_symbols[common_symbols]
df_feature_ann = df_feature_ann_final.loc[common_symbols] # Este es el df_feature_ann final

# Alinear df_data_test con las columnas finales de df_data_train
df_data_test = df_data_test_symbols.reindex(columns=df_data_train.columns, fill_value=0)

print(f"\nDimensiones finales para CorrAdjust:")
print(f"  df_data_train: {df_data_train.shape} (Columnas: {df_data_train.columns[:3].tolist() if not df_data_train.empty else 'Vacío'}...)")
print(f"  df_data_test: {df_data_test.shape}")
print(f"  df_feature_ann: {df_feature_ann.shape} (Índice: '{df_feature_ann.index.name}', Columnas: {df_feature_ann.columns.tolist()})")

# Verificación final
if not df_data_train.empty and (df_feature_ann.shape[0] != df_data_train.shape[1] or not df_data_train.columns.equals(df_feature_ann.index)):
    print("Error CRÍTICO: La alineación final entre df_data_train y df_feature_ann falló.")
    # exit()
else:
    print("Alineación final entre df_data_train y df_feature_ann parece correcta.")

print("\nVerificación de df_data_train.head() (columnas deberían ser símbolos):")
print(df_data_train.head())

# --- Preparación de datos de expresión (train y test) ---
def prepare_expression_data_for_corradjust(df_orig_counts, id_map, target_feature_names_index):
    df_renamed = df_orig_counts.rename(columns=id_map)
    df_renamed = df_renamed.loc[:, pd.notna(df_renamed.columns)]
    df_renamed = df_renamed.loc[:, ~df_renamed.columns.duplicated(keep='first')]
    common_cols = target_feature_names_index.intersection(df_renamed.columns)
    return df_renamed[common_cols]

df_data_train = prepare_expression_data_for_corradjust(df_counts_train_orig, id_to_name_map, df_feature_ann.index)
df_feature_ann = df_feature_ann.loc[df_data_train.columns] # Alinear df_feature_ann con las columnas finales de train

df_data_test = prepare_expression_data_for_corradjust(df_counts_test_orig, id_to_name_map, df_feature_ann.index)
df_data_test = df_data_test.reindex(columns=df_data_train.columns, fill_value=0)

print(f"Dimensiones finales para CorrAdjust:")
print(f"  df_data_train: {df_data_train.shape}")
print(f"  df_data_test: {df_data_test.shape}")
print(f"  df_feature_ann: {df_feature_ann.shape}")
print(f"  Columnas en df_feature_ann: {df_feature_ann.columns.tolist()}") # Chequea columnas
print(f"  Índice de df_feature_ann: {df_feature_ann.index.name}") # Chequea índice

# %%
# --- 3. Definición de ref_feature_colls (archivos con las rutas metabólicas) ---
#ref_feature_colls = {
    #"Oncogenic Signatures": {
     #   "path": "/home/lucia2024/TFM/ext_data/GMT_files/c6.all.v2024.1.Hs.symbols.gmt",
      #  "sign": "absolute",
       # "feature_pair_types": ["protein_coding-protein_coding"],
        #"high_corr_frac": 0.01
    #},
    #"Computational Gene Sets (Cancer Modules)": {
       # "path": "/home/lucia2024/TFM/ext_data/GMT_files/c4.all.v2024.1.Hs.symbols.gmt",
        #"sign": "absolute",
        #"feature_pair_types": ["protein_coding-protein_coding"],
        #"high_corr_frac": 0.01
   # }
#}

ref_feature_colls = {
    "Non-basal TNBC enriched pathways ": {
        "path": "/home/lucia2024/TFM/ext_data/GMT_files/non_basal_enriched_pathways_R.gmt",
        "sign": "absolute",
        "feature_pair_types": ["protein_coding-protein_coding"],
        "high_corr_frac": 0.01
    },
    "Basal TNBC enriched pathways": {
        "path": "/home/lucia2024/TFM/ext_data/GMT_files/basal_enriched_pathways_R.gmt",
        "sign": "absolute",
        "feature_pair_types": ["protein_coding-protein_coding"],
        "high_corr_frac": 0.01
    }
}


# %%
# --- 4. Inicialización y ajuste del modelo CorrAdjust ---
print(f"\nInicializando Modelo CorrAdjust. Directorio de salida: {corradjust_output_dir}")

# Me aseguro de que df_feature_ann y df_data_train no están vacíos y tienen features antes de meterlos al modelo
if df_feature_ann.empty or df_data_train.empty or df_data_train.shape[1] == 0:
    print("Error: df_feature_ann o df_data_train están vacíos o sin features. No se puede continuar con CorrAdjust.")
    exit()

MIN_PAIRS_THRESHOLD = 100

model = CorrAdjust(
    df_feature_ann, # ¡Debe estar indexado por los nombres de columna de df_data_train!
    ref_feature_colls,
    corradjust_output_dir,
    min_pairs_to_score=MIN_PAIRS_THRESHOLD
)

print(f"Ajustando el modelo CorrAdjust con {df_data_train.shape[1]} features y {df_data_train.shape[0]} muestras...")
model.fit(df_data_train) 
print("Modelo ajustado.")


# %%
# --- 5. Revisión de los datos de ajuste (esto es opcional) ---
fit_data_path = os.path.join(corradjust_output_dir, "fit.tsv")
try:
    fit_data = pd.read_csv(fit_data_path, sep="\t", index_col=0, dtype={"PC": "Int64"})
    print("\nDatos de ajuste (fit.tsv):")
    print(fit_data) 
except FileNotFoundError:
    print(f"\nAdvertencia: No se encontró el archivo {fit_data_path}. Comprueba la ejecución de model.fit().")

# %%
# --- 6. Cálculo de "feature scores" y evaluación con el CONJUNTO DE PRUEBA ---
print("\nCalculando feature scores para el conjunto de entrenamiento...")
feature_scores_train = model.compute_feature_scores(df_data_train)

if not df_data_test.empty and df_data_test.shape[1] > 0:
    print("\nEvaluando el modelo con el conjunto de prueba...")
    feature_scores_test = model.compute_feature_scores(df_data_test)

    print("\nEjemplo de Feature Scores (limpios, TEST SET) para la colección 'Oncogenic Signatures':")
    if "Clean" in feature_scores_test and "Oncogenic Signatures" in feature_scores_test["Clean"]:
        print(feature_scores_test["Clean"]["Oncogenic Signatures"].head())
    else:
        print("No se pudieron encontrar los scores para 'Clean' y 'Oncogenic Signatures' en el test set.")

    # --- 6.1. Volcanoplot para el conjunto de prueba ---
    print("\nGenerando el volcanoplot para el conjunto de prueba...")
    volcano_filename = "volcano.test_samples.png" # Solo se mete el nombre del archivo, nada más
    model.make_volcano_plot(
        feature_scores_test,
        volcano_filename, # Pasar solo el nombre del archivo
        annotate_features=10,
    )
    print(f"Gráfico de volcán guardado en: {os.path.join(corradjust_output_dir, volcano_filename)}")

    # --- 6.2. Gráfico de distribución de correlaciones para el conjunto de prueba ---
    print("\nGenerando gráfico de distribución de correlaciones para el conjunto de prueba...")
    corr_distr_filename = "corr_distr.test_samples.png" # Solo el nombre del archivo
    model.make_corr_distr_plot(
        df_data_test,
        corr_distr_filename # Pasar solo el nombre del archivo
    )
    print(f"Gráfico de distribución de correlaciones guardado en: {os.path.join(corradjust_output_dir, corr_distr_filename)}")

    # --- 6.3. Obtener la matriz de datos de prueba corregida ---
    print("\nObteniendo la matriz de datos de prueba corregida (limpia)...")
    df_data_test_clean, df_rsquareds_test = model.transform(df_data_test)
    print("Matriz de datos de prueba corregida (primeras 5 filas/columnas):")
    print(df_data_test_clean.iloc[:5, :5])

else:
    print("\ndf_data_test está vacío o sin features después de la preparación, se omiten los pasos de evaluación del test set.")

# %% 
# --- 6.4. Guardar los features scores de ambas colecciones ---

def guardar_scores(clave, nombre_archivo):
    try:
        scores = feature_scores_test.get("Clean", {}).get(clave, None)
        if isinstance(scores, (pd.Series, pd.DataFrame)):
            csv_output_path = os.path.join(corradjust_output_dir, nombre_archivo)
            scores.to_csv(csv_output_path, header=True)
            return f"Scores de '{clave}' guardados en: {csv_output_path}"
        return f"Error: Los scores de '{clave}' no son un DataFrame o Series."
    except Exception as e:
        return f"Error al guardar scores de '{clave}': {e}"

if 'feature_scores_test' in locals() and feature_scores_test is not None:
    resultados = []
    # resultados.append(guardar_scores("Oncogenic Signatures", "feature_scores_test_oncogenic_signatures.csv"))
    #resultados.append(guardar_scores("Computational Gene Sets (Cancer Modules)", "feature_scores_test_comp_gene_sets_cancer_modules.csv"))
    resultados.append(guardar_scores("Non-basal TNBC enriched pathways ", "feature_scores_test_non_basal_tnbc_enriched_pathways.csv"))
    resultados.append(guardar_scores("Basal TNBC enriched pathways", "feature_scores_test_basal_tnbc_enriched_pathways.csv"))
    for resultado in resultados:
        print(resultado)
else:
    print("'feature_scores_test' no está definido o es None.")

print("\n--- Script 02 completado ---")
# %%
