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
cvrts <- c("BiopsyType", "Age..5.year.range..e.g...35.31.35...40.36.40...45.41.45..etc..") # Importante: quito NCN.PAM50 porque me da problemas

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
cox_all <- coxph(Surv(OS_days, OS_event) ~ ., data = train, x = TRUE)
cox_cvrts <- coxph(Surv(OS_days, OS_event) ~ ., data = train[,c(surv_vars, cvrts)], x = TRUE)
cox_sign <- coxph(Surv(OS_days, OS_event) ~ ., data = train %>% dplyr::select(-all_of(cvrts)), x = TRUE)


# Prepare test data
colnames(counts_test) <- gsub("\\..*", "", colnames(counts_test))
test <- counts_test[, fcbf_sign$ENSEMBL]
colnames(test) <- fcbf_sign$SYMBOL
test <- scale(test, scale = FALSE)

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

cindex_train <- data.frame(
  Study = c("signature + cvrts", "cvrts", "signature"),
  C = c(pred_train_all$concordance, pred_train_cvrts$concordance, pred_train_sign$concordance),
  ci_l = c(
    pred_train_all$concordance - 1.96 * pred_train_all$std.err,
    pred_train_cvrts$concordance - 1.96 * pred_train_cvrts$std.err,
    pred_train_sign$concordance - 1.96 * pred_train_sign$std.err),
  ci_u = c(
    pred_train_all$concordance + 1.96 * pred_train_all$std.err,
    pred_train_cvrts$concordance + 1.96 * pred_train_cvrts$std.err,
    pred_train_sign$concordance + 1.96 * pred_train_sign$std.err),
  Subset = "Train"
)


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
    pred_test_all$concordance - 1.96 * pred_test_all$std.err,
    pred_test_cvrts$concordance - 1.96 * pred_test_cvrts$std.err,
    pred_test_sign$concordance - 1.96 * pred_test_sign$std.err),
  ci_u = c(
    pred_test_all$concordance + 1.96 * pred_test_all$std.err,
    pred_test_cvrts$concordance + 1.96 * pred_test_cvrts$std.err,
    pred_test_sign$concordance + 1.96 * pred_test_sign$std.err),
  Subset = "Test"
)

# Results of C-Index in TCGA (train and test)
cindex_scanb_fcbf <- rbind.data.frame(cindex_train, cindex_test)


# Save model
saveRDS(train, file = file.path(outputdir, "train_all.rds"))
saveRDS(test, file = file.path(outputdir, "test_all.rds"))

saveRDS(cox_all, file = file.path(outputdir, "cox_cvrts_yaccs.rds"))
saveRDS(cox_cvrts, file = file.path(outputdir, "cox_cvrts.rds"))
saveRDS(cox_yaccs, file = file.path(outputdir, "cox_yaccs.rds"))

saveRDS(cindex_tcga, file = file.path(outputdir, "cindex_tcga.rds"))