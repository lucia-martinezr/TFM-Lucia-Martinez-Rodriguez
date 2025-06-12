# Build Cox model based on SCANB train dataset


library(survival)
library(dplyr)

# Inputpaths
counts_train_inputpath <- "C:/Users/lulim/OneDrive/Documentos/GitHub/TFM-Lucia-Martinez-Rodriguez/data/data_partitions/counts_train.rds"
counts_test_inputpath <- "C:/Users/lulim/OneDrive/Documentos/GitHub/TFM-Lucia-Martinez-Rodriguez/data/data_partitions/counts_test.rds"
metadata_train_inputpath <- "C:/Users/lulim/OneDrive/Documentos/GitHub/TFM-Lucia-Martinez-Rodriguez/data/data_partitions/metadata_train.rds"
metadata_test_inputpath <- "C:/Users/lulim/OneDrive/Documentos/GitHub/TFM-Lucia-Martinez-Rodriguez/data/data_partitions/metadata_test.rds"
signature_inputpath <- "C:/Users/lulim/OneDrive/Documentos/GitHub/TFM-Lucia-Martinez-Rodriguez/data/tnbc_signature/fcbf_signature_annot.rds"

# Outputpath
outputdir <- "C:/Users/lulim/OneDrive/Documentos/GitHub/TFM-Lucia-Martinez-Rodriguez/data/cox/fcbf_scanb/"

# Load data
counts_train <- readRDS(counts_train_inputpath)
counts_test <- readRDS(counts_test_inputpath)
metadata_train <- readRDS(metadata_train_inputpath)
metadata_test <- readRDS(metadata_test_inputpath)
fcbf_sign <- readRDS(signature_inputpath)
fcbf_sign$SYMBOL <- make.names(fcbf_sign$SYMBOL)



# Arguments  
surv_vars <- c("OS_days", "OS_event")
cvrts <- c("Age..5.year.range..e.g...35.31.35...40.36.40...45.41.45..etc..") 


# Prepare train data
colnames(counts_train) <- gsub("\\..*", "", colnames(counts_train))
train <- counts_train[, fcbf_sign$ENSEMBL]
colnames(train) <- fcbf_sign$SYMBOL
train <- scale(train, scale = TRUE)

train <- cbind.data.frame(
  metadata_train[, c(surv_vars, cvrts)],
  train
)
rownames(train) <- metadata_train$Patient # Jose usa "barcode" pero como tengo pacientes únicos, yo uso "Patient".
# train$vital_status <- ifelse(train$vital_status == "Dead", TRUE, FALSE) Yo no necesito cambiar a lógico porque tengo la variable de evento como binaria.
# Dejo las covariables como factor, porque he visto que el modelo de cox ya las gestiona solo.

# filas_con_na <- train[rowSums(is.na(train)) > 0, ]
# write.csv(filas_con_na, "C:/Users/lulim/OneDrive/Escritorio/filas_con_na.csv", row.names = FALSE)

# Train model
cox_all <- coxph(Surv(OS_days, OS_event) ~ ., data = train, x = TRUE)
cox_cvrts <- coxph(Surv(OS_days, OS_event) ~ ., data = train[,c(surv_vars, cvrts)], x = TRUE)
cox_sign <- coxph(Surv(OS_days, OS_event) ~ ., data = train %>% dplyr::select(-all_of(cvrts)), x = TRUE)


# Prepare test data
colnames(counts_test) <- gsub("\\..*", "", colnames(counts_test))
test <- counts_test[, fcbf_sign$ENSEMBL]
colnames(test) <- fcbf_sign$SYMBOL
test <- scale(test, scale = TRUE)

test <- cbind.data.frame(
  metadata_test[, c(surv_vars, cvrts)],
  test
)
rownames(test) <- metadata_test$Patient
# test$vital_status <- ifelse(test$vital_status == "Dead", TRUE, FALSE)
names(test) <- make.names(names(test))


# Performances in train ----------------
pred_train_all <- survival::concordance(Surv(OS_days, OS_event) ~ predict(cox_all, train), train)
pred_train_cvrts <- survival::concordance(Surv(OS_days, OS_event) ~ predict(cox_cvrts, train[, c(surv_vars, cvrts)]), train)
pred_train_sign<- survival::concordance(Surv(OS_days, OS_event) ~ predict(cox_sign, train[, c(surv_vars, fcbf_sign$SYMBOL)]), train)

#Nota: en este paso tuve un problema: mis objetos pred_train no tenían un componente std.err sino var, por lo que el
# código que sigue, para crear el C-index, difiere del de yaccs en que en vez de poner
# pred_train_all$concordance - 1.96 * pred_train_all$std.err
# Pone:
# pred_train_all$concordance - 1.96 * sqrt(pred_train_all$var)
# pero en esencia ambos hacen lo mismo. Es por la actualización del paquete.

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
hr_table_all <- extract_hr_table(cox_all)
hr_table_cvrts <- extract_hr_table(cox_cvrts)
hr_table_sign <- extract_hr_table(cox_sign)

# # Guardo la tablas
# write.table(hr_table_all, file = file.path(outputdir, "f_hr_table_all.tsv"),  sep = "\t", row.names = FALSE, quote = FALSE)
# write.table(hr_table_cvrts, file = file.path(outputdir, "f_hr_table_cvrts.tsv"),  sep = "\t", row.names = FALSE, quote = FALSE)
# write.table(hr_table_sign, file = file.path(outputdir, "f_hr_table_sign.tsv"),  sep = "\t", row.names = FALSE, quote = FALSE)

# Testing (SCAN-B test) ----------------
# Three scenarios:
#   Model with signature + cvrts
#   Model with only cvrts
#   Model with only signature

pred_test_all <- survival::concordance(Surv(OS_days, OS_event) ~ predict(cox_all, test), test)
pred_test_cvrts <- survival::concordance(Surv(OS_days, OS_event) ~ predict(cox_cvrts, test[, c(surv_vars, cvrts)]), test)
pred_test_sign <- survival::concordance(Surv(OS_days, OS_event) ~ predict(cox_sign, test[, c(surv_vars, fcbf_sign$SYMBOL)]), test)

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


# # Save model
# saveRDS(train, file = file.path(outputdir, "f_train_all.rds"))
# saveRDS(test, file = file.path(outputdir, "f_test_all.rds"))
# 
# saveRDS(cox_all, file = file.path(outputdir, "f_cox_cvrts_sign.rds"))
# saveRDS(cox_cvrts, file = file.path(outputdir, "f_cox_cvrts.rds"))
# saveRDS(cox_sign, file = file.path(outputdir, "f_cox_sign.rds"))
# 
# saveRDS(cindex_scanb_fcbf, file = file.path(outputdir, "f_cindex_scanb_fcbf.rds"))
