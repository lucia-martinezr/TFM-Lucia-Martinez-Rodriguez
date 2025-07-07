# Inputpaths
metadata_inputpath <- "C:/Users/lulim/OneDrive/Documentos/GitHub/TFM-Lucia-Martinez-Rodriguez/data/ext_data/clin_tnbc.rds" 
counts_inputpath <- "C:/Users/lulim/OneDrive/Escritorio/Corradjust Data/05_Ouput_Corradjust_Pcs_opt_2_collection_TNBC_paper_y_log2_transformation/df_data_complete_clean.csv"

# Load data
metadata <- readRDS(metadata_inputpath) # 715 obs 87 variables
counts <- read.csv(counts_inputpath, header = TRUE, sep = ',') # 715 obs 15467 variables
rownames(counts) <- counts[,1]
counts <- counts[,-1] # 715 obs 15466 variables

# Define partitions
smpSize <- floor(0.8 * nrow(metadata))
set.seed(55)
trainIdx <- sample(seq_len(nrow(metadata)), size = smpSize)

# Split into train and test
metadata_train <- metadata[trainIdx, ]
metadata_test <- metadata[-trainIdx, ]

counts_train <- counts[trainIdx, ] # 572 obs
counts_test <- counts[-trainIdx, ] # 143 obs




# Feature selection with FCBF DATA CLEAN

# devtools::install_github("lubianat/FCBF")
# BiocManager::install("AnnotationDbi")
# BiocManager::install("org.Hs.eg.db")
library(FCBF)
library(org.Hs.eg.db)
library(AnnotationDbi)

# Arguments
min_su <- 0.0025

# # Save gene names
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

# Me encuentra 4 genes

# Retrieve gene names
rn<-gene_names[selected_features$index]
rownames(selected_features)<-rn 










# Build Cox model based on cleaned signature

library(survival)
library(dplyr)

# --- INICIO DE LA CORRECCIÓN 1: Definir los nombres de los genes ---
# Esta es la forma correcta de obtener los nombres de los genes de tu objeto FCBF
# Haz esto justo después de la sección donde creas y pones nombres de fila a "selected_features"
selected_gene_names <- rownames(selected_features)
# --- FIN DE LA CORRECCIÓN 1 ---


names(metadata_train) <- make.names(names(metadata_train))
names(metadata_test) <- make.names(names(metadata_test))
# Esta línea no hace nada, puedes eliminarla si quieres
# make.names("Age (5-year range, e.g., 35(31-35), 40(36-40), 45(41-45) etc.)")

# Arguments  
surv_vars <- c("OS_days", "OS_event")
cvrts <- c("Age..5.year.range..e.g...35.31.35...40.36.40...45.41.45..etc..") 


# Prepare train data
# La preparación de 'train' ya estaba bien en tu última versión
# colnames(counts_train) <- gsub("\\..*", "", colnames(counts_train)) # Esto ya no es necesario si tus counts ya están limpios
train <- counts_train[, selected_gene_names] # Usamos la variable definida
train <- scale(train, scale = TRUE)

train <- cbind.data.frame(
  metadata_train[, c(surv_vars, cvrts)],
  train
)
rownames(train) <- metadata_train$Patient

# Train model
cox_all <- coxph(Surv(OS_days, OS_event) ~ ., data = train, x = TRUE)
cox_cvrts <- coxph(Surv(OS_days, OS_event) ~ ., data = train[,c(surv_vars, cvrts)], x = TRUE)
cox_sign <- coxph(Surv(OS_days, OS_event) ~ ., data = train %>% dplyr::select(-all_of(cvrts)), x = TRUE)


# Prepare test data
# colnames(counts_test) <- gsub("\\..*", "", colnames(counts_test)) # Ya no es necesario
test <- counts_test[, selected_gene_names] # Usamos la misma variable
test <- scale(test, scale = TRUE)

test <- cbind.data.frame(
  metadata_test[, c(surv_vars, cvrts)],
  test
)
rownames(test) <- metadata_test$Patient

# Performances in train ----------------
pred_train_all <- survival::concordance(Surv(OS_days, OS_event) ~ predict(cox_all, train), train)
pred_train_cvrts <- survival::concordance(Surv(OS_days, OS_event) ~ predict(cox_cvrts, train[, c(surv_vars, cvrts)]), train)


pred_train_sign <- survival::concordance(
  Surv(OS_days, OS_event) ~ predict(cox_sign, train[, c(surv_vars, selected_gene_names)]), # Usamos la variable definida
  train
)

cindex_train <- data.frame(
  Study = c("signature + cvrts", "cvrts", "signature"),
  C = c(pred_train_all$concordance, pred_train_cvrts$concordance, pred_train_sign$concordance),
  ci_l = c(
    pred_train_all$concordance - 1.96 * sqrt(pred_train_all$var),
    pred_train_cvrts$concordance - 1.96 * sqrt(pred_train_cvrts$var),
    pred_train_sign$concordance - 1.96 * sqrt(pred_train_sign$var)
  ),  
  ci_u = c(
    pred_train_all$concordance + 1.96 * sqrt(pred_train_all$var),
    pred_train_cvrts$concordance + 1.96 * sqrt(pred_train_cvrts$var),
    pred_train_sign$concordance + 1.96 * sqrt(pred_train_sign$var)
  ),
  Subset = "Train"
)

# --- Hazard Ratios ---

# Creo una función para extraer HR, IC y p-valores de un modelo de Cox
extract_hr_table <- function(cox_model) {
  # Extraer el resumen del modelo
  summary_model <- summary(cox_model)
  
  # Extraer los coeficientes y sus exponenciales (HR)
  coefs <- summary_model$coefficients
  
  # Extraer los intervalos de confianza
  conf_ints <- summary_model$conf.int
  
  # Crear la tabla de resultados
  hr_table <- data.frame(
    Variable = rownames(coefs),
    Coefficient = coefs[, "coef"],
    Hazard_Ratio = coefs[, "exp(coef)"],
    CI_Lower = conf_ints[, "lower .95"],
    CI_Upper = conf_ints[, "upper .95"],
    P_Value = coefs[, "Pr(>|z|)"]
  )
  
  # Redondear para que se entienda mejor
  hr_table[, -1] <- round(hr_table[, -1], 4)
  
  # Ordenar por p-valor para ver las variables más significativas primero
  hr_table <- hr_table[order(hr_table$P_Value), ]
  
  rownames(hr_table) <- NULL
  return(hr_table)
}

# Uso la función para extraer los HR de cada uno de los modelos
hr_table_all <- extract_hr_table(cox_all)
hr_table_cvrts <- extract_hr_table(cox_cvrts)
hr_table_sign <- extract_hr_table(cox_sign)

# --- Testing (SCAN-B test) ----------------
pred_test_all <- survival::concordance(Surv(OS_days, OS_event) ~ predict(cox_all, test), test)
pred_test_cvrts <- survival::concordance(Surv(OS_days, OS_event) ~ predict(cox_cvrts, test[, c(surv_vars, cvrts)]), test)

pred_test_sign <- survival::concordance(
  Surv(OS_days, OS_event) ~ predict(cox_sign, test[, c(surv_vars, selected_gene_names)]), # Usamos la misma variable
  test
)

cindex_test <- data.frame(
  Study = c("signature + cvrts", "cvrts", "signature"),
  C = c(pred_test_all$concordance, pred_test_cvrts$concordance, pred_test_sign$concordance),
  ci_l = c(
    pred_test_all$concordance - 1.96 * sqrt(pred_test_all$var),
    pred_test_cvrts$concordance - 1.96 * sqrt(pred_test_cvrts$var),
    pred_test_sign$concordance - 1.96 * sqrt(pred_test_sign$var)
  ),  
  ci_u = c(
    pred_test_all$concordance + 1.96 * sqrt(pred_test_all$var),
    pred_test_cvrts$concordance + 1.96 * sqrt(pred_test_cvrts$var),
    pred_test_sign$concordance + 1.96 * sqrt(pred_test_sign$var)
  ), 
  Subset = "Test"
)

# Results of C-Index (train and test)
cindex_scanb_fcbf <- rbind.data.frame(cindex_train, cindex_test)





