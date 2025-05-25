# Build Cox model based on SCANB train dataset
# ----------------------
# source("requirements.R")
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

# Arguments  
# Nota: no meto Ki67 porque es una columna vacía tanto en train como en test. 
# No meto ER, HER2, PR porque todos son "Negative".
# No meto ClinGroup porque todos son TNBC.

surv_vars <- c("OS_event", "OS_days")
cvrts <- c(
  "BiopsyType", "LN", "NHG", "Age..5.year.range..e.g...35.31.35...40.36.40...45.41.45..etc..", 
  "Size.mm", "pT",  "InvCa.type", "NCN.PAM50", "NCN.ROR.risk.cat"
)

# Pregunta: qué hago si en muchas de las covariables hay valores NA???




# Load data
counts_train <- readRDS(counts_train_inputpath)
counts_test <- readRDS(counts_test_inputpath)
metadata_train <- readRDS(metadata_train_inputpath)
metadata_test <- readRDS(metadata_test_inputpath)
fcbf_sign <- readRDS(signature_inputpath)
fcbf_sign$SYMBOL <- make.names(fcbf_sign$SYMBOL)

# Prepare train data
colnames(counts_train) <- gsub("\\..*", "", colnames(counts_train))
train <- counts_train[, fcbf_sign$ENSEMBL]
colnames(train) <- fcbf_sign$SYMBOL
train <- scale(train, scale = FALSE) # Centra los datos.Los ajusta para que cada columna tenga media cero al restar la media de cada columna de sus valores.
# No los normaliza.


# train <- cbind.data.frame(
#   metadata_train[,surv_vars],
#   train
# )


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
cox_all <- coxph(Surv(OS_event, OS_days) ~ ., data = train, x = TRUE)
cox_cvrts <- coxph(Surv(OS_event, OS_days) ~ ., data = train[,c(surv_vars, cvrts)], x = TRUE)
cox_sign <- coxph(Surv(OS_event, OS_days) ~ ., data = train %>% dplyr::select(-all_of(cvrts)), x = TRUE)








####################### A PARTIR DE AQUÍ SIN ADAPTAR NI USAR ##########################################


# Prepare test data
test <- counts_test[, yaccs$ENSEMBL]
colnames(test) <- yaccs$SYMBOL
test <- scale(test, scale = FALSE)

test <- cbind.data.frame(
  metadata_test[, c(surv_vars, cvrts)],
  test
)
rownames(test) <- metadata_test$barcode
test$vital_status <- ifelse(test$vital_status == "Dead", TRUE, FALSE)
names(test) <- make.names(names(test))


# Performances in train ----------------
pred_train_all <- survival::survConcordance(Surv(SurvTime, vital_status) ~ predict(cox_all, train), train)
pred_train_cvrts <- survival::survConcordance(Surv(SurvTime, vital_status) ~ predict(cox_cvrts, train[, c(surv_vars, cvrts)]), train)
pred_train_yaccs <- survival::survConcordance(Surv(SurvTime, vital_status) ~ predict(cox_yaccs, train[, c(surv_vars, yaccs$SYMBOL)]), train)

cindex_train <- data.frame(
  Study = c("yaccs + cvrts", "cvrts", "yaccs"),
  C = c(pred_train_all$concordance, pred_train_cvrts$concordance, pred_train_yaccs$concordance),
  ci_l = c(
    pred_train_all$concordance - 1.96 * pred_train_all$std.err,
    pred_train_cvrts$concordance - 1.96 * pred_train_cvrts$std.err,
    pred_train_yaccs$concordance - 1.96 * pred_train_yaccs$std.err),
  ci_u = c(
    pred_train_all$concordance + 1.96 * pred_train_all$std.err,
    pred_train_cvrts$concordance + 1.96 * pred_train_cvrts$std.err,
    pred_train_yaccs$concordance + 1.96 * pred_train_yaccs$std.err),
  Subset = "Train"
)


# Testing (TCGA-COAD test) ----------------
# Three scenarios:
#   Model with yaccs + cvrts
#   Model with only cvrts
#   Model with only yaccs

pred_test_all <- survival::survConcordance(Surv(SurvTime, vital_status) ~ predict(cox_all, test), test)
pred_test_cvrts <- survival::survConcordance(Surv(SurvTime, vital_status) ~ predict(cox_cvrts, test[, c(surv_vars, cvrts)]), test)
pred_test_yaccs <- survival::survConcordance(Surv(SurvTime, vital_status) ~ predict(cox_yaccs, test[, c(surv_vars, yaccs$SYMBOL)]), test)

cindex_test <- data.frame(
  Study = c("yaccs + cvrts", "cvrts", "yaccs"),
  C = c(pred_test_all$concordance, pred_test_cvrts$concordance, pred_test_yaccs$concordance),
  ci_l = c(
    pred_test_all$concordance - 1.96 * pred_test_all$std.err,
    pred_test_cvrts$concordance - 1.96 * pred_test_cvrts$std.err,
    pred_test_yaccs$concordance - 1.96 * pred_test_yaccs$std.err),
  ci_u = c(
    pred_test_all$concordance + 1.96 * pred_test_all$std.err,
    pred_test_cvrts$concordance + 1.96 * pred_test_cvrts$std.err,
    pred_test_yaccs$concordance + 1.96 * pred_test_yaccs$std.err),
  Subset = "Test"
)

# Results of C-Index in TCGA (train and test)
cindex_tcga <- rbind.data.frame(cindex_train, cindex_test)


# Save model
saveRDS(train, file = file.path(outputdir, "train_all.rds"))
saveRDS(test, file = file.path(outputdir, "test_all.rds"))

saveRDS(cox_all, file = file.path(outputdir, "cox_cvrts_yaccs.rds"))
saveRDS(cox_cvrts, file = file.path(outputdir, "cox_cvrts.rds"))
saveRDS(cox_yaccs, file = file.path(outputdir, "cox_yaccs.rds"))

saveRDS(cindex_tcga, file = file.path(outputdir, "cindex_tcga.rds"))