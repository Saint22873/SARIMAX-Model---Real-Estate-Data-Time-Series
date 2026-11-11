install.packages(c("readr", "dplyr", "lubridate", "janitor", "forecast", "tseries", "ggplot2", "FinTS", "lmtest", "sandwich"))

rm(list = ls())

library(readr)
library(janitor)
library(tidyr)
library(dplyr)
library(lubridate)
library(forecast)
library(ggplot2) 
library(patchwork) 
library(tseries)
library(FinTS)
library(lmtest)
library(sandwich)

ma_raw <- read_csv(
  "C:/Users/LENOVO/Desktop/Trabajos de Econometria/Real Estate Time Series Data/ma_lga_12345.csv"
)
#VERIFICACION DE LOS DATOS -----------------------------
names(ma_raw)
# Esperado: "saledate" "MA" "type" "bedrooms"

dim(ma_raw)
# Esperado: ~248  4

head(ma_raw, 3)
# Esperado: filas con 30/09/2007, 441854, house, 2, etc.


#LIMPIEZA DE LOS DATOS --------------------------------
ma_clean <- ma_raw %>%
  clean_names() %>%                                  # MA -> ma
  mutate(saledate = dmy(saledate)) %>%               # dd/mm/yyyy
  mutate(
    type     = factor(type, levels = c("house", "unit")),
    bedrooms = as.integer(bedrooms),
    ma       = as.numeric(ma)
  ) %>%
  drop_na(saledate, ma, type, bedrooms) %>%
  distinct() %>%
  distinct(saledate, type, bedrooms, .keep_all = TRUE) %>%
  arrange(type, bedrooms, saledate)



glimpse(ma_clean)
range(ma_clean$saledate)
table(ma_clean$type, ma_clean$bedrooms)




#INICIO DE LA METODOLOGIA BOX-JENKINS----------------------------

#IDENTIFICAR POR MEDIO DE UNA GRAFICA SI EXISTE LINEALIDAD O OUTLINERS
#IDENTIFICAR:

#CICLO
#TENDENCIA
#ALEATORIEDAD
#ESTACIONALIDAD


# Elige la serie que quieres analizar
filtro_tipo <- "house"
filtro_dorm <- 3

serie_df <- ma_clean %>%
  filter(type == filtro_tipo, bedrooms == filtro_dorm) %>%
  arrange(saledate)

# Convertir a ts trimestral (frecuencia = 4)
# start = c(año, trimestre) de la primera observación
y <- ts(serie_df$ma,
        start = c(year(min(serie_df$saledate)),
                  quarter(min(serie_df$saledate))),
        frequency = 4)

y


# Descomposición multiplicativa (usual cuando la varianza crece con el nivel)
desc <- decompose(y, type = "multiplicative")

# Gráfica con los 4 componentes apilados:
#   observed  -> serie original (contiene ciclo + tendencia + estacionalidad + ruido)
#   trend     -> Tendencia
#   seasonal  -> Estacionalidad
#   random    -> Aleatoriedad (residuo)
plot(desc)
# STL SEPARA TENDENCIA Y ESTACIONALIDAD SIN ASUMIR FORMA FIJA 

stl_fit <- stl(y, s.window = "periodic", robust = FALSE)
plot(stl_fit)

#--------------------------------
#GRAFICACIÓN

autoplot(y) +
  labs(
    title = paste("Serie trimestral:", filtro_tipo, "-", filtro_dorm, "dorm."),
    subtitle = "Paso 1 Box-Jenkins: identificación de media y varianza constantes",
    x = "Año", y = "MA"
  ) +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed",
              color = "red", linewidth = 0.7) +
  theme_minimal()
  
  

#Primero ubicamos si necesita alguna correccion por metodo box-cox-------------
 


  # ---- 4.1 Estimar el λ óptimo de Box-Cox ----
lambda_opt <- BoxCox.lambda(y, method = "guerrero")
  # Alternativa: method = "loglik"
lambda_opt

#INICIO DE LA DIFERENCIACIÓN --------------------------
y_t <- y 

# 3. Primera diferencia regular (d = 1)
# Esto estabiliza la media eliminando la TENDENCIA
dy <- diff(y_t, differences = 1)

# 4. Verificar la gráfica de la serie diferenciada
autoplot(dy) +
  labs(
    title = "Serie diferenciada (d = 1)",
    subtitle = "Eliminación de tendencia para estabilizar la media",
    x = "Año", 
    y = "Diferencia de MA"
  ) +
  theme_minimal()



#--------------------------------------------------
#Otra parte de la diferenciación
#En esta parte hara que las olas desaparezcan (obtenemos dy_s)
dy_s <- diff(dy, lag = 4)


autoplot(dy_s) +
  labs(
    title = "Serie con Diferencia Regular y Estacional (dy_s)",
    subtitle = "Serie estacionaria: sin tendencia y sin estacionalidad",
    x = "Año",
    y = "Diferencia de MA"
  ) +
  theme_minimal()

#---------------------------------------------------------------------------
  
#4.3 Pruebas de Raíz Unitaria sobre la serie final (dy_s)
# H0 ADF: Hay raíz unitaria (No estacionaria). Si p < 0.05, RECHAZAMOS y es estacionaria.
adf.test(dy_s, alternative = "stationary")

# H0 KPSS: Es estacionaria. Si p > 0.05, NO RECHAZAMOS y es estacionaria.
kpss.test(dy_s, null = "Level")



#GENERACION DE UN CORRELOGRAMA CON LA SERIE YA DIFERENCIADA----------------

#Correlograma de la serie final
ggAcf(dy_s, lag.max = 24) + labs(title = "ACF - dy_s")
ggPacf(dy_s, lag.max = 24) + labs(title = "PACF - dy_s")



# En la gráfica de dy_s vimos un pico gigante en 2008/2009. 
# Para que el modelo no se confunda, creamos una variable dummy (1 en crisis, 0 en el resto).
# Asumiendo datos trimestrales, time(y) da valores como 2008.00, 2008.25, etc.
crisis_2008 <- ifelse(time(y) >= 2008.50 & time(y) <= 2009.75, 1, 0)

# Agrupamos la dummy en una matriz llamada xreg (regresores externos)
xreg_dummy <- matrix(crisis_2008, ncol = 1)
colnames(xreg_dummy) <- "crisis"

modelo_auto <- auto.arima(y, d = 1, D = 1, xreg = xreg_dummy,stepwise = FALSE, approximation = FALSE)

summary(modelo_auto)                              

#---------------------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# 8. MODELACIÓN
# ------------------------------------------------------------------------------

# Candidato 3: AR(1) regular y SMA(1) estacional (lo que sugiere la teoría)
mod1 <- Arima(y, order = c(1, 1, 0), seasonal = c(0, 1, 1), xreg = xreg_dummy)

# Candidato 4: ARMA(1,1) regular y SMA(1) estacional (el clásico "Airlines" pero con AR)
mod2 <- Arima(y, order = c(1, 1, 1), seasonal = c(0, 1, 1), xreg = xreg_dummy)

# Candidato 5: Solo SMA(1) estacional (por si la parte regular no es necesaria)
mod3 <- Arima(y, order = c(0, 1, 0), seasonal = c(0, 1, 1), xreg = xreg_dummy)

# Imprimimos los AIC de todos para comparar
cat("AIC mod1 (1,1,0)(0,1,1):", mod1$aic, "\n")
cat("AIC mod2 (1,1,1)(0,1,1):", mod2$aic, "\n")
cat("AIC mod3 (0,1,0)(0,1,1):", mod3$aic, "\n")

# ------------------------------------------------------------------------------
# 9. CONCLUSIÓN Y PRONÓSTICO FINAL (Usando el modelo ganador mod1)
# ------------------------------------------------------------------------------
modelo_final <- mod1

# Ver el resumen de los coeficientes (esto va en tu trabajo escrito)
summary(modelo_final)

# Diagnóstico visual final (esto debe mostrar un ACF limpio)
checkresiduals(modelo_final)

# Pronóstico a 8 trimestres (2 años)
# Recuerda que la crisis ya pasó, así que la dummy futura es 0
futuro_dummy <- matrix(0, nrow = 8, ncol = 1)
colnames(futuro_dummy) <- "crisis"

pronostico_final <- forecast(modelo_final, h = 8, xreg = futuro_dummy)

# Gráfica final del pronóstico
autoplot(pronostico_final) +
  labs(
    title = "Pronóstico Final - Modelo Box-Jenkins",
    subtitle = "Modelo SARIMA(1,1,0)(0,1,1)[4] con intervención por crisis 2008",
    x = "Año",
    y = "Precio (MA)"
  ) +
  theme_minimal()

# ------------------------------------------------------------------------------
# 10. Evaluación del modelo
# -------------------------------------------------------------------
residuos <- residuals(modelo_final)

cat("\n--- FASE 4: DIAGNÓSTICO DE RESIDUOS ---\n")

# 4.1 Autocorrelación (Ljung-Box)
cat("\n1. Prueba de Autocorrelación (Ljung-Box):\n")
print(Box.test(residuos, lag = 8, type = "Ljung-Box", fitdf = 2))

# 4.2 Normalidad (Jarque-Bera)
cat("\n2. Prueba de Normalidad (Jarque-Bera):\n")
print(jarque.bera.test(residuos))

# 4.3 Heterocedasticidad condicional (Prueba ARCH)
cat("\n3. Prueba de Efectos ARCH (Heterocedasticidad):\n")
print(ArchTest(residuos, lags = 4))

# 4.4 Gráfico de diagnóstico integrado
checkresiduals(modelo_final)



# La prueba de Jarque-Bera no la paso por las razones de que
# el shock de la crisis inmobiliaria fue tan brutal que dejó "colas pesadas"
# (outliers) en los residuos. 
# La prueba de Jarque-Bera es extremadamente
# sensible a los valores atípicos.
# Ademas del mercado inmobiliario suele tener
# caídas abruptas y subidas más lentas,
# lo que genera una distribución asimétrica 
# (no es una campana de Gauss perfecta)


# Extraer los errores estándar robustos de tu modelo (necesitamos hacerlo manual)
# Como el modelo es ARIMA, usamos la función coeftest
coeftest(modelo_final)

# ------------------------------------------------------------------------------
# 11. PRONÓSTICO
# ------------------------------------------------------------------------------
# 5.1 Preparar matriz futura (sin crisis)
futuro_xreg <- matrix(0, nrow = 8, ncol = 1)
colnames(futuro_xreg) <- "crisis"

# 5.2 Pronosticar
pronostico <- forecast(modelo_final, h = 8, xreg = futuro_xreg)

# 5.3 Gráfica final
autoplot(pronostico) +
  labs(title = "Pronóstico a 2 Años (SARIMAX)",
       subtitle = "Incluye dummy por crisis 2008", x = "Año", y = "Precio") +
  theme_minimal()

