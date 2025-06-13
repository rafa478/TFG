library(ggplot2)
library(dplyr)
library(readr)

leer_datos <- function(archivos) {
  datos_completos <- data.frame()
  
  for (archivo in archivos) {
    datos <- read_csv(archivo)
    datos_completos <- bind_rows(datos_completos, datos)
  }
  
  return(datos_completos)
}

args <- commandArgs(trailingOnly = TRUE)
datos <- leer_datos(args)

resumen <- datos %>%
  group_by(Tiempo) %>%
  summarise(
    Proliferacion_media = mean(Proliferacion, na.rm = TRUE),
    Proliferacion_sd = sd(Proliferacion, na.rm = TRUE),
    Metastasis_media = mean(Metastasis, na.rm = TRUE),
    Metastasis_sd = sd(Metastasis, na.rm = TRUE),
    Metastasis_Proliferacion_media = mean(MetastasisYProliferacion, na.rm = TRUE),
    Metastasis_Proliferacion_sd = sd(MetastasisYProliferacion, na.rm = TRUE),
    Ninguna_media = mean(Ninguna, na.rm = TRUE),
    Ninguna_sd = sd(Ninguna, na.rm = TRUE),
    Muertas_media = mean(Muertas, na.rm = TRUE),
    Muertas_sd = sd(Muertas, na.rm = TRUE)
  )

ggplot(resumen, aes(x = Tiempo)) +
  geom_line(aes(y = Proliferacion_media, color = 'Proliferación'), size = 1) +
  geom_line(aes(y = Metastasis_media, color = 'Metástasis'), size = 1) +
  geom_line(aes(y = Metastasis_Proliferacion_media, color = 'Metástasis y Proliferación'), size = 1) +
  geom_line(aes(y = Muertas_media, color = 'Muertas'), size = 1) +
  geom_ribbon(aes(ymin = Proliferacion_media - Proliferacion_sd, ymax = Proliferacion_media + Proliferacion_sd), 
              fill = "blue", alpha = 0.2) +
  geom_ribbon(aes(ymin = Metastasis_media - Metastasis_sd, ymax = Metastasis_media + Metastasis_sd), 
              fill = "orange", alpha = 0.2) +
  geom_ribbon(aes(ymin = Metastasis_Proliferacion_media - Metastasis_Proliferacion_sd, 
                  ymax = Metastasis_Proliferacion_media + Metastasis_Proliferacion_sd), 
              fill = "green", alpha = 0.2) +
  geom_ribbon(aes(ymin = Muertas_media - Muertas_sd, ymax = Muertas_media + Muertas_sd), 
              fill = "purple", alpha = 0.2) +
  labs(title = "Evolución de Proliferación, Metástasis, Ninguna y Muertes", 
       x = "Tiempo", 
       y = "Cantidad", 
       color = "Grupo") +
  theme_minimal() +
  theme(legend.position = "top")
