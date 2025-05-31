##################### Prueba #####################################

# Nuevas covariantes
covariates <- c("BiopsyType", "Age..5.year.range..e.g...35.31.35...40.36.40...45.41.45..etc..", 
                "VIM", "CTSB", "LTF", "TPT1", "EPHX3", "CCN2", "CXCL5", "RPS17", 
                "S100A2", "RPL10A", "CHI3L1", "ACTA2", "RPL19", "CLU", "RPL10", "COX6C", "EEF2")




# Agrupar CoreBiopsy y CoreBiopsy2nd en BiopsyType
train$BiopsyType <- factor(train$BiopsyType,
                           levels = c("CoreBiopsy", "CoreBiopsy2nd", "OP"),
                           labels = c("CoreBiopsy", "CoreBiopsy", "OP"))

# Agrupar LumA, LumB y unclassified en NCN.PAM50
# train$NCN.PAM50 <- factor(train$NCN.PAM50,
# levels = c("Basal", "Her2", "LumA", "LumB", "Normal", "unclassified"),
# labels = c("Basal", "Her2", "Lum/other", "Lum/other", "Normal", "Lum/other"))

# Verificar los cambios
table(train$BiopsyType)
# table(train$NCN.PAM50)


# Eliminar niveles no utilizados después de reagrupación
train$BiopsyType <- droplevels(train$BiopsyType)
# train$NCN.PAM50 <- droplevels(train$NCN.PAM50)

# Revisar los niveles para asegurarte de que sean los correctos
# levels(train$NCN.PAM50)
levels(train$BiopsyType)

# Volver a ajustar los modelos
univ_formulas <- sapply(covariates, function(x) as.formula(paste('Surv(OS_days, OS_event) ~', x)))

univ_models <- lapply(univ_formulas, function(x) {
  tryCatch({
    coxph(x, data = train)
  }, error = function(e) {
    message("Error en modelo para formula: ", deparse(x), "\n", e)
    return(NULL)
  })
})

# Extracción de datos
univ_results <- lapply(univ_models, function(x) {
  if (is.null(x)) return(NULL)
  x_summary <- summary(x)
  p.value <- signif(x_summary$wald["pvalue"], digits = 2)
  wald.test <- signif(x_summary$wald["test"], digits = 2)
  beta <- signif(x_summary$coef[1], digits = 2)
  HR <- signif(x_summary$coef[2], digits = 2)
  HR.confint.lower <- signif(x_summary$conf.int[,"lower .95"], 2)
  HR.confint.upper <- signif(x_summary$conf.int[,"upper .95"], 2)
  HR <- paste0(HR, " (", HR.confint.lower, "-", HR.confint.upper, ")")
  res <- c(beta, HR, wald.test, p.value)
  names(res) <- c("beta", "HR (95% CI for HR)", "wald.test", "p.value")
  return(res)
})

# Filtrar resultados consistentes
univ_results_filtered <- Filter(Negate(is.null), univ_results)

# Crear DataFrame solo con resultados consistentes
res <- t(as.data.frame(univ_results_filtered, check.names = FALSE))
as.data.frame(res)



