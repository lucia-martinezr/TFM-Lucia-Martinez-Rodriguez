# Create forest plots 

library(forestplot)
library(dplyr)
library(viridis)


# Inputpath
fcbf_hr_path <- "C:/Users/lulim/OneDrive/Documentos/GitHub/TFM-Lucia-Martinez-Rodriguez/data/cox/fcbf_scanb/f_hr_table_all.tsv"
corradjust_hr_path <- "C:/Users/lulim/OneDrive/Documentos/GitHub/TFM-Lucia-Martinez-Rodriguez/data/cox/corradjust_scanb/c_hr_table_all.tsv"

# Outputpath
output_figures_dir <- "C:/Users/lulim/OneDrive/Documentos/GitHub/TFM-Lucia-Martinez-Rodriguez/figures"

# Load data
hr_fcbf_all <- read.table(fcbf_hr_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
hr_corradjust_all <- read.table(corradjust_hr_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
names(hr_fcbf_all)
names(hr_corradjust_all)

# Pipeline 1: fcbf

# Preparo los datos para el gráfico
fcbf_plot_data <- hr_fcbf_all %>%
  mutate(Variable = ifelse(grepl("Age..5.year.range", Variable, fixed = TRUE), "Edad", Variable)) %>% # Acorto el nombre de la variable Age
  arrange(Hazard_Ratio) %>% # Ordeno por Hazard Ratio
  mutate(
    labeltext = Variable,
    mean = Hazard_Ratio,
    lower = CI_Lower,
    upper = CI_Upper
  ) # Creación de las columnas que requiere un forest plot

# Añado el asterisco en caso de que sea significativo
fcbf_plot_data <- mutate(fcbf_plot_data, Significancia = case_when(
  P_Value <= 0.001 ~ "***",
  P_Value <= 0.01  ~ "**",
  P_Value <= 0.05  ~ "*",
  TRUE             ~ ""
))

fcbf_plot_data <- mutate(fcbf_plot_data, labeltext = paste0(Variable, " ", Significancia))

fcbf_plot_data <- arrange(fcbf_plot_data, Hazard_Ratio)

fcbf_plot_data <- mutate(fcbf_plot_data,
                         mean = Hazard_Ratio,
                         lower = CI_Lower,
                         upper = CI_Upper
                          )



# Preparo la configuración de la imagen 
png(
  filename = file.path(output_figures_dir, "forestplot_fcbf.png"),
  width = 800, # Ancho en píxeles
  height = 1200 # Alto en píxeles
)

# Genero el forest plot 
fcbf_plot_data %>%
  forestplot(
    labeltext = labeltext,
    mean = mean,
    lower = lower,
    upper = upper,
    title = "Firma FCBF: Hazard Ratios",
    xlab = "Hazard Ratio (HR)",
    zero = 1,
    boxsize = 0.35,
    lwd.ci =2,
    lwd.zero =2,
    # lineheight =  unit(1.5, "lines"), # Para configurar el espacio entre líneas
    lineheight = "auto",
    txt_gp = fpTxtGp(label = gpar(cex = 1.2),ticks = gpar(cex = 1.1), xlab  = gpar(cex = 1.3), title = gpar(cex = 1.5)),
    col = fpColors(box = "darkblue", line = "gray50", summary = "darkblue")
  ) %>%
  fp_add_header(labeltext = "Variable")

# Cerrar el archivo gráfico
dev.off()


# Pipeline 2: Corradjust

# Preparo los datos para el gráfico
corradjust_plot_data <- hr_corradjust_all %>%
  mutate(Variable = ifelse(grepl("Age..5.year.range", Variable, fixed = TRUE), "Edad", Variable)) %>%
  arrange(Hazard_Ratio) %>%
  mutate(
    labeltext = Variable,
    mean = Hazard_Ratio,
    lower = CI_Lower,
    upper = CI_Upper
  )

# Añado el asterisco en caso de que sea significativo
corradjust_plot_data <- mutate(fcbf_plot_data, Significancia = case_when(
  P_Value <= 0.001 ~ "***",
  P_Value <= 0.01  ~ "**",
  P_Value <= 0.05  ~ "*",
  TRUE             ~ ""
))

corradjust_plot_data <- mutate(corradjust_plot_data, labeltext = paste0(Variable, " ", Significancia))

corradjust_plot_data <- arrange(corradjust_plot_data, Hazard_Ratio)

corradjust_plot_data <- mutate(corradjust_plot_data,
                         mean = Hazard_Ratio,
                         lower = CI_Lower,
                         upper = CI_Upper
)

# Preparo la configuración de la imagen 
png(
  filename = file.path(output_figures_dir, "forestplot_corradjust.png"),
  width = 800,
  height = 1200
)

# Genero el forest plot
corradjust_plot_data %>%
  forestplot(
    labeltext = labeltext,
    mean = mean,
    lower = lower,
    upper = upper,
    title = "Firma CorrAdjust: Hazard Ratios",
    xlab = "Hazard Ratio (HR)",
    zero = 1,
    boxsize = 0.35,
    lwd.ci =2,
    lwd.zero =2,
    lineheight = "auto",
    txt_gp = fpTxtGp(label = gpar(cex = 1.2),ticks = gpar(cex = 1.1), xlab  = gpar(cex = 1.3), title = gpar(cex = 1.5)),
    col = fpColors(box = "darkred", line = "gray50", summary = "darkred")
  ) %>%
  fp_add_header(labeltext = "Variable")

# Cerrar el dispositivo
dev.off()