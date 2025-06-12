# Build Cox model based on SCANB train dataset


library(survival)
library(dplyr)

# Inputpaths
clean_counts_train_inputpath <- "C:/Users/lulim/OneDrive/Documentos/GitHub/TFM-Lucia-Martinez-Rodriguez/data/clean_data_partitions/clean_counts_train.rds"
clean_counts_test_inputpath <- "C:/Users/lulim/OneDrive/Documentos/GitHub/TFM-Lucia-Martinez-Rodriguez/data/clean_data_partitions/clean_counts_test.rds"
clean_metadata_train_inputpath <- "C:/Users/lulim/OneDrive/Documentos/GitHub/TFM-Lucia-Martinez-Rodriguez/data/clean_data_partitions/clean_metadata_train.rds"
clean_metadata_test_inputpath <- "C:/Users/lulim/OneDrive/Documentos/GitHub/TFM-Lucia-Martinez-Rodriguez/data/clean_data_partitions/clean_metadata_test.rds"
signature_inputpath <- "C:/Users/lulim/OneDrive/Documentos/GitHub/TFM-Lucia-Martinez-Rodriguez/data/clean_tnbc_signature/clean_signature.rds"

# Outputpath
outputdir <- "C:/Users/lulim/OneDrive/Documentos/GitHub/TFM-Lucia-Martinez-Rodriguez/data/cox/corradjust_scanb/"

# Load data
clean_counts_train <- readRDS(clean_counts_train_inputpath)
clean_counts_test <- readRDS(clean_counts_test_inputpath)
clean_metadata_train <- readRDS(clean_metadata_train_inputpath)
clean_metadata_test <- readRDS(clean_metadata_test_inputpath)
clean_sign <- readRDS(signature_inputpath)
# clean_sign$SYMBOL <- make.names(clean_sign$SYMBOL)



# Arguments  
surv_vars <- c("OS_days", "OS_event")
cvrts <- c("Age..5.year.range..e.g...35.31.35...40.36.40...45.41.45..etc..") # Importante: quito NCN.PAM50 porque me da problemas

# Prepare train data
train <- clean_counts_train[, clean_sign$SYMBOL]
train <- scale(train, scale = TRUE)


train <- cbind.data.frame(
  clean_metadata_train[, c(surv_vars, cvrts)],
  train
)

# filas_con_na <- train[rowSums(is.na(train)) > 0, ]
# write.csv(filas_con_na, "C:/Users/lulim/OneDrive/Escritorio/filas_con_na.csv", row.names = FALSE)

# Train model
cox_all <- coxph(Surv(OS_days, OS_event) ~ ., data = train, x = TRUE)
cox_cvrts <- coxph(Surv(OS_days, OS_event) ~ ., data = train[,c(surv_vars, cvrts)], x = TRUE)
cox_sign <- coxph(Surv(OS_days, OS_event) ~ ., data = train %>% dplyr::select(-all_of(cvrts)), x = TRUE)


# Prepare test data
test <- clean_counts_test[, clean_sign$SYMBOL]
test <- scale(test, scale = FALSE)

test <- cbind.data.frame(
  clean_metadata_test[, c(surv_vars, cvrts)],
  test
)

# test$vital_status <- ifelse(test$vital_status == "Dead", TRUE, FALSE)
names(test) <- make.names(names(test))


# Performances in train ----------------
pred_train_all <- survival::concordance(Surv(OS_days, OS_event) ~ predict(cox_all, train), train)
pred_train_cvrts <- survival::concordance(Surv(OS_days, OS_event) ~ predict(cox_cvrts, train[, c(surv_vars, cvrts)]), train)
pred_train_sign<- survival::concordance(Surv(OS_days, OS_event) ~ predict(cox_sign, train[, c(surv_vars, clean_sign$SYMBOL)]), train)

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

#--- Hazard Ratios ---

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
clean_hr_table_all <- extract_hr_table(cox_all)
clean_hr_table_cvrts <- extract_hr_table(cox_cvrts)
clean_hr_table_sign <- extract_hr_table(cox_sign)

# # Guardo la tablas
# write.table(clean_hr_table_all, file = file.path(outputdir, "c_hr_table_all.tsv"),  sep = "\t", row.names = FALSE, quote = FALSE)
# write.table(clean_hr_table_cvrts, file = file.path(outputdir, "c_hr_table_cvrts.tsv"),  sep = "\t", row.names = FALSE, quote = FALSE)
# write.table(clean_hr_table_sign, file = file.path(outputdir, "c_hr_table_sign.tsv"),  sep = "\t", row.names = FALSE, quote = FALSE)

# Testing (SCAN-B test) ----------------
# Three scenarios:
#   Model with signature + cvrts
#   Model with only cvrts
#   Model with only signature

pred_test_all <- survival::concordance(Surv(OS_days, OS_event) ~ predict(cox_all, test), test)
pred_test_cvrts <- survival::concordance(Surv(OS_days, OS_event) ~ predict(cox_cvrts, test[, c(surv_vars, cvrts)]), test)
pred_test_sign <- survival::concordance(Surv(OS_days, OS_event) ~ predict(cox_sign, test[, c(surv_vars, clean_sign$SYMBOL)]), test)

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

# Results of C-Index in TCGA (train and test)
cindex_scanb_corr <- rbind.data.frame(cindex_train, cindex_test)


# # Save model
# saveRDS(train, file = file.path(outputdir, "c_train_all.rds"))
# saveRDS(test, file = file.path(outputdir, "c_test_all.rds"))
# 
# saveRDS(cox_all, file = file.path(outputdir, "c_cox_cvrts_sign.rds"))
# saveRDS(cox_cvrts, file = file.path(outputdir, "c_cox_cvrts.rds"))
# saveRDS(cox_sign, file = file.path(outputdir, "c_cox_sign.rds"))
# 
# saveRDS(cindex_scanb_corr, file = file.path(outputdir, "c_cindex_scanb.rds"))