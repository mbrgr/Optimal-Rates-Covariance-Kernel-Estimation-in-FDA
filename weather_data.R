library(tidyverse)
library(plotly)

library(lubridate)
library(hms)

library(biLocPol)
library(future.apply)


# In produkt*.txt stehen folgende Parameter zur Verf?gung:
# STATIONS_ID Stationsidentifikationsnummer
# MESS_DATUM Zeitstempel yyyymmddhhmi 

# PP_10 Luftdruck auf Stationsh?he hPa
# TT_10 Lufttemperatur in 2m H?he ?C
# TM5_10 Lufttemperatur 5cm H?he ?C
# RF_10 relative Feuchte in 2m H?he %
# TD_10 Taupunkttemperatur in 2m H?he ?C
# eor Ende data record

# Fehlwerte sind mit -999 gekennzeichnet. Die Messungen sind vor dem Jahr 2000 einem Zeitstempel in
# MEZ zugeordnet, ab dem Jahr 2000 einem Zeitstempel in UTC. Die Taupunktstemperatur ist aus der
# Lufttemperatur in 2m H?he und der relativen Feuchtemessungen berechnet. Die Werte sind Mittelwerte ?ber
# die Minute, welche zum Zeitstempel endet. Der Zeitstempel ist vor dem Jahr 2000 in MEZ gegeben, ab dem
# Jahr 2000 in UTC

load("C:/Users/mberger.PC12599/HESSENBOX/GitHub/FDA-Covariance-Estimator/data/weather_data_nuremberg.RData")
N$DATUM = str_split_i(N$MESS_DATUM, pattern = " ", i = 1) |> ymd()

tage = c(1, 4, 8, 12, 15, 18, 22, 25, 29)

N = N |>
  as_tibble() |> 
  filter(TAG %in% tage) 

N |>  
  filter(TAG %in% c(1, 8, 15, 22, 29)) |> 
  ggplot() +
  geom_line(aes(x = UHRZEIT, y = TT_10, group = JAHR*TAG, colour = JAHR), alpha = .4) +
  facet_wrap(MONAT ~.)

N |>  
  filter(TAG %in% c(1, 15, 29), 
         MONAT == 8) |> 
  ggplot() +
  geom_line(aes(x = UHRZEIT, y = TT_10, group = JAHR*TAG, colour = JAHR), alpha = .6) + 
  labs(y = "Temp. in C°", x = "hours", title = "Temp. in August", colour = "year") + 
  theme(text = element_text(size = 18)) 

N_wide = N |> 
  filter(TAG %in% tage) |>  # Beschraenkung auf nicht aufeinanderfolgende Tage
  mutate(UHRZEIT = as.character(UHRZEIT)) |> 
  select(JAHR, MONAT, TAG, UHRZEIT, TT_10) |> 
  pivot_wider(names_from = UHRZEIT,
              values_from = TT_10)

###### August ######
Y = N_wide |> 
  filter(MONAT == 8) |> 
  select(-(1:3))

Z = Y |> 
  as.data.frame() |> # TODO: fix this in biLocPol
  observation_transformation(na.rm = T)


# h_cv = k_fold_cv(Y, H, K = 5, m = 1, h.parallel = T, h.parallel.environment = T, na.rm = T)
p.eval = 72
W = local_polynomial_weights(144, 0.2, p.eval, T, m = 1)
g_hat = eval_weights(W, Z)
var_hat = diag(g_hat)

eval_time = N$UHRZEIT[1:146][seq(2, 144, 2)]


var_est = tibble(var_hat, x = eval_time)
var_est |> 
  ggplot(aes(x = x, y = sqrt(var_hat))) + 
  geom_line(linewidth = .7) + 
  lims(y = c(0.2, 5)) + 
  labs(y = NULL, x = "hour", title = "Std. deviation of temperatur in August") + 
  theme(legend.position = "none", 
        text = element_text(size = 18)) 

cov_est_df = data.frame(x = W$x.eval, y = W$x.eval, z = g_hat)
cs2 = list(c(0, 1), c("lightblue", "darkred"))

plot_ly(cov_est_df, x = ~x, y = ~y, z = ~g_hat, size = .4) |> 
  add_surface(colorscale = cs2, alpha = .3) |> 
  layout(scene = list(xaxis = list(title = ""), 
                      yaxis = list(title = ""), 
                      zaxis = list(title = "")))

temp = matrix(diag(g_hat), p.eval, p.eval)
cor_hat = g_hat / sqrt( temp * t(temp) )


plot_ly(cov_est_df, x = ~x*24, y = ~y*24, z = ~cor_hat, size = .4) |> 
  add_surface(colorscale = cs2, alpha = .3) |> 
  layout(scene = list(xaxis = list(title = ""), 
                      yaxis = list(title = ""), 
                      zaxis = list(title = "")))

###### Januar ######
Y1 = N_wide |> 
  filter(MONAT == 1) |> 
  select(-(1:3))

Z1 = Y1 |> 
  as.data.frame() |> # TODO: fix this in biLocPol
  observation_transformation(na.rm = T)

N |>  
  filter(TAG %in% c(1, 15, 29), 
         MONAT == 1) |> 
  ggplot() +
  geom_line(aes(x = UHRZEIT, y = TT_10, group = JAHR*TAG, colour = JAHR), alpha = .7) + 
  labs(y = "Temp. in C°", x = "hours", title = "Temp. in January", colour = "year") + 
  theme(text = element_text(size = 18)) 


# h_cv = k_fold_cv(Y, H, K = 5, m = 1, h.parallel = T, h.parallel.environment = T, na.rm = T)
# p.eval = 72
# W1 = local_polynomial_weights(144, 0.2, p.eval, T, m = 2)
g_hat1 = eval_weights(W, Z1)
var_hat1 = diag(g_hat1)

var_est1 = tibble(var_hat = var_hat1, x = eval_time)
var_est1 |> 
  ggplot(aes(x = x, y = sqrt(var_hat1))) + 
  geom_line(size = .6) + 
  lims(y = c(0.2, 5)) + 
  labs(y = NULL, x = "hour", title = "Std. deviation of temperatur in January") + 
  theme(text = element_text(size = 18)) 

# both std. deviations in one picture
ggplot() + 
  geom_line(data    = var_est1, 
            mapping = aes(x = x, y = sqrt(var_hat1), lty = "January")) + 
  geom_line(data    = var_est, 
            mapping = aes(x = x, y = sqrt(var_hat), lty = "August")) + 
  geom_line(size = .6) + 
  lims(y = c(0.2, 5)) + 
  labs(y = NULL, x = "hour") + 
  theme(text = element_text(size = 18), legend.position = "top") +
  scale_linetype_discrete(name = "Std. deviation in month:  ") 


cov_est_df1 = data.frame(x = W$x.eval, y = W$x.eval, z = g_hat1)
cs2 = list(c(0, 1), c("lightblue", "darkred"))

plot_ly(cov_est_df1, x = ~x, y = ~y, z = ~g_hat1, size = .4) |> 
  add_surface(colorscale = cs2, alpha = .3) |> 
  layout(scene = list(xaxis = list(title = ""), 
                      yaxis = list(title = ""), 
                      zaxis = list(title = "")))

temp = matrix(diag(g_hat1), p.eval, p.eval)
cor_hat1 = g_hat1 / sqrt( temp * t(temp) )


plot_ly(cov_est_df1, x = ~x*24, y = ~y*24, z = ~cor_hat1, size = .4) |> 
  add_surface(colorscale = cs2, alpha = .3) |> 
  layout(scene = list(xaxis = list(title = ""), 
                      yaxis = list(title = ""), 
                      zaxis = list(title = "")))

###### Januar ######
Y1 = N_wide |> 
  filter(MONAT == 1) |> 
  select(-(1:3))

Z1 = Y1 |> 
  as.data.frame() |> # TODO: fix this in biLocPol
  observation_transformation(na.rm = T)

N |>  
  filter(TAG %in% c(1, 15, 29), 
         MONAT == 1) |> 
  ggplot() +
  geom_line(aes(x = UHRZEIT, y = TT_10, group = JAHR*TAG, colour = JAHR), alpha = .7) + 
  labs(y = "Temp. in C°", x = "hours", title = "Temp. in January", colour = "year") + 
  theme(text = element_text(size = 18)) 


# h_cv = k_fold_cv(Y, H, K = 5, m = 1, h.parallel = T, h.parallel.environment = T, na.rm = T)
# p.eval = 72
# W1 = local_polynomial_weights(144, 0.2, p.eval, T, m = 2)
g_hat1 = eval_weights(W, Z1)
var_hat1 = diag(g_hat1)

var_est1 = tibble(var_hat = var_hat1, x = eval_time)
var_est1 |> 
  ggplot(aes(x = x, y = sqrt(var_hat1))) + 
  geom_line(size = .6) + 
  lims(y = c(0.2, 5)) + 
  labs(y = NULL, x = "hour", title = "Std. deviation of temperatur in January") + 
  theme(text = element_text(size = 18)) 


cov_est_df1 = data.frame(x = W$x.eval, y = W$x.eval, z = g_hat1)
cs2 = list(c(0, 1), c("lightblue", "darkred"))

plot_ly(cov_est_df1, x = ~x, y = ~y, z = ~g_hat1, size = .4) |> 
  add_surface(colorscale = cs2, alpha = .3) |> 
  layout(scene = list(xaxis = list(title = ""), 
                      yaxis = list(title = ""), 
                      zaxis = list(title = "")))

temp = matrix(diag(g_hat1), p.eval, p.eval)
cor_hat1 = g_hat1 / sqrt( temp * t(temp) )


plot_ly(cov_est_df1, x = ~x*24, y = ~y*24, z = ~cor_hat1, size = .4) |> 
  add_surface(colorscale = cs2, alpha = .3) |> 
  layout(scene = list(xaxis = list(title = ""), 
                      yaxis = list(title = ""), 
                      zaxis = list(title = "")))

###### April ######
Y4 = N_wide |> 
  filter(MONAT == 4) |> 
  select(-(1:3))

Z4 = Y4 |> 
  as.data.frame() |> # TODO: fix this in biLocPol
  observation_transformation(na.rm = T)

N |>  
  filter(TAG %in% c(1, 15, 29), 
         MONAT == 1) |> 
  ggplot() +
  geom_line(aes(x = UHRZEIT, y = TT_10, group = JAHR*TAG, colour = JAHR), alpha = .7) + 
  labs(y = "Temp. in C°", x = "hours", title = "Temp. in January", colour = "year") + 
  theme(text = element_text(size = 18)) 

g_hat4 = eval_weights(W, Z4)
var_hat4 = diag(g_hat4)
var_est4 = tibble(var_hat = var_hat4, x = eval_time)

cov_est_df4 = data.frame(x = W$x.eval, y = W$x.eval, z = g_hat4)
cs2 = list(c(0, 1), c("lightblue", "darkred"))

temp = matrix(diag(g_hat4), p.eval, p.eval)
cor_hat4 = g_hat4 / sqrt( temp * t(temp) )


# all std deviations in one picture
ggplot() + 
  geom_line(data    = var_est1, 
            mapping = aes(x = x, y = sqrt(var_hat1), lty = "January")) + 
  geom_line(data    = var_est, 
            mapping = aes(x = x, y = sqrt(var_hat), lty = "August")) + 
  geom_line(data    = var_est4, 
            mapping = aes(x = x, y = sqrt(var_hat4), lty = "April")) + 
  geom_line(linewidth = .6) + 
  lims(y = c(0.2, 5.8)) + 
  labs(y = NULL, x = "hour") + 
  theme(text = element_text(size = 18), legend.position = "top") +
  scale_linetype_discrete(name = "Std. deviation in month:  ") 

plot_ly(cov_est_df1, x = ~x*24, y = ~y*24, z = ~cor_hat1, size = .4) |> 
  add_surface(colorscale = cs2, alpha = .3) |> 
  layout(scene = list(xaxis = list(title = ""), 
                      yaxis = list(title = ""), 
                      zaxis = list(title = "")))


##### Luftdruck #####

N_wide_pp = N |> 
  mutate(UHRZEIT = as.character(UHRZEIT)) |> 
  select(JAHR, MONAT, TAG, DATUM, UHRZEIT, PP_10) |> 
  pivot_wider(names_from = UHRZEIT,
              values_from = PP_10, ) |> 
  relocate(JAHR, MONAT, TAG, 
           .before = DATUM)

Y_pp = N_wide_pp |> 
  filter(MONAT == 8) |> 
  select(-(1:4))

# mean
Y_mean = colMeans(Y_pp, na.rm = T)
Y_mean_est = locpol::locLinSmootherC(0:143/144, Y_mean, 0:143/144, 0.08, EpaK)
Y_mean_est = data.frame(Y_mean_est, UHRZEIT = N$UHRZEIT[7:150])


ggplot() + 
  geom_line(data    = N |> filter(MONAT == 8, JAHR > 2020), 
            mapping = aes(x = UHRZEIT, y = PP_10, group = DATUM)) + 
  geom_line(data    = Y_mean_est, 
            mapping = aes(x = UHRZEIT, y = beta0), 
            col     = "red")


Z_pp = Y_pp |> 
  as.data.frame() |> # TODO: fix this in biLocPol
  observation_transformation(na.rm = T)


# h_cv = k_fold_cv(Y, H, K = 5, m = 1, h.parallel = T, h.parallel.environment = T, na.rm = T)
p.eval = 75
W_pp = local_polynomial_weights(144, 0.3, p.eval, T)
g_hat_pp = eval_weights(W_pp, Z_pp)

cov_est_pp = data.frame(x = W$x.eval, y = W$x.eval, z = g_hat_pp)
cs2 = list(c(0, 1), c("lightblue", "darkred"))

plot_ly(cov_est_df, x = ~24*x, y = ~24*y, z = ~g_hat_pp, size = .4) |> 
  add_surface(colorscale = cs2, alpha = .3) |> 
  layout(scene = list(xaxis = list(title = "x"), 
                      yaxis = list(title = "y"), 
                      zaxis = list(title = "")))

temp = matrix(diag(g_hat_pp), p.eval, p.eval)
cor_hat_pp = g_hat_pp / sqrt( temp * t(temp) )


plot_ly(cov_est_pp, x = ~x*24, y = ~y*24, z = ~cor_hat_pp, size = .4) |> 
  add_surface(colorscale = cs2, alpha = .3) |> 
  layout(scene = list(xaxis = list(title = "x"), 
                      yaxis = list(title = "y"), 
                      zaxis = list(title = "")))
