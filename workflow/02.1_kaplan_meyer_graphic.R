
# install.packages("ggsurvfit")
library(ggsurvfit)

# devtools::install_github("zabore/condsurv")
library(condsurv)

library(ggplot2)

# Inputpaths
metadata_inputpath <- "C:/Users/lulim/OneDrive/Documentos/GitHub/TFM-Lucia-Martinez-Rodriguez/data/data_partitions/metadata_train.rds"
counts_inputpath <- "C:/Users/lulim/OneDrive/Documentos/GitHub/TFM-Lucia-Martinez-Rodriguez/data/data_partitions/counts_train.rds"

# Outputpaths
outputpath <- "C:/Users/lulim/OneDrive/Documentos/GitHub/TFM-Lucia-Martinez-Rodriguez/figures"

# Load data
metadata_train <- readRDS(metadata_inputpath)
counts_train <- readRDS(counts_inputpath)

# Overall survival probability plot
gg<- survfit2(Surv(metadata_train$OS_days, metadata_train$OS_event) ~ 1, data = metadata_train) %>% 
  ggsurvfit() +
  labs(
    x = "Days",
    y = "Overall survival probability"
  ) + 
  add_confidence_interval() +
  add_risktable()

output_path <- file.path(outputpath, "kmplot.png")
ggsave(filename = output_path, plot = gg, width = 8, height = 6, dpi = 300)


# Creating two groups based on OS_event (1 or 0)
km_groups <- as.data.frame(metadata_train$OS_event, )


# Overall survival probability plot BY RISK GROUP
gg <- survfit2(Surv(OS_days, OS_event) ~ metadata_train$OS_event, data = metadata_train) %>% 
  ggsurvfit() +
  labs(
    x = "Days",
    y = "Overall survival probability"
  ) + 
  add_confidence_interval() +
  add_risktable() +
  ggtitle("Kaplan-Meier Plot by Risk Group")

output_path <- file.path(outputpath, "kmplot_grouped.png")
ggsave(filename = output_path, plot = gg, width = 8, height = 6, dpi = 300)

print(gg)
