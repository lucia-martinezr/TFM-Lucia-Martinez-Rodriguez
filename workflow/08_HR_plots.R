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


# Pipeline 1: fcbf

# Preparo los datos para el gráfico
fcbf_plot_data <- hr_fcbf_all %>%
  mutate(Variable = ifelse(grepl("Age..5.year.range", Variable, fixed = TRUE), "Age", Variable)) %>% # Acorto el nombre de la variable Age
  arrange(Hazard_Ratio) %>% # Ordeno por Hazard Ratio
  mutate(
    labeltext = Variable,
    mean = Hazard_Ratio,
    lower = CI_Lower,
    upper = CI_Upper
  ) # Creación de las columnas que requiere un forest plot

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
    title = "Pipeline FCBF: Hazard Ratios",
    xlab = "Hazard Ratio (HR)",
    zero = 1,
    boxsize = 0.25,
    lineheight =  unit(1.5, "lines"), # Para configurar el espacio entre líneas
    col = fpColors(box = "darkblue", line = "gray50", summary = "darkblue")
  ) %>%
  fp_add_header(labeltext = "Variable")

# Cerrar el archivo gráfico
dev.off()


# Pipeline 2: Corradjust

# Preparo los datos para el gráfico
corradjust_plot_data <- hr_corradjust_all %>%
  mutate(Variable = ifelse(grepl("Age..5.year.range", Variable, fixed = TRUE), "Age", Variable)) %>%
  arrange(Hazard_Ratio) %>%
  mutate(
    labeltext = Variable,
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
    title = "Pipeline CorrAdjust: Hazard Ratios",
    xlab = "Hazard Ratio (HR)",
    zero = 1,
    boxsize = 0.25,
    lineheight = unit(1.5, "lines"),
    col = fpColors(box = "darkred", line = "gray50", summary = "darkred")
  ) %>%
  fp_add_header(labeltext = "Variable")

# Cerrar el dispositivo
dev.off()