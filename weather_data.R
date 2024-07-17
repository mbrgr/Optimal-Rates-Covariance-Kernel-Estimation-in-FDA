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


load("data/weather_data_raw.RData")
N$DATUM = str_split_i(N$MESS_DATUM, pattern = " ", i = 1) |> ymd()

tage = c(1, 4, 8, 12, 15, 18, 22, 25, 29)

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
  filter(TAG %in% tage) |>  # No consecutive days
  mutate(UHRZEIT = as.character(UHRZEIT)) |> 
  select(JAHR, MONAT, TAG, UHRZEIT, TT_10) |> 
  pivot_wider(names_from = UHRZEIT,
              values_from = TT_10)

#### Plot empirical covariance for consecutive days vs. non consecutive #####

Y_all = N[-(1:6),] |> # exclude first 6 entries since the days is incomplete
  mutate(UHRZEIT = as.character(UHRZEIT)) |> 
  select(JAHR, MONAT, TAG, UHRZEIT, TT_10) |> 
  pivot_wider(names_from = UHRZEIT, values_from = TT_10) |> 
  filter(MONAT == 4) |> 
  select(-(1:3))
empirical_cov_all = Y_all |> 
  as.data.frame() |> 
  observation_transformation(na.rm = T, grid.type = "full")


empirical_cov_filtered = N_wide |> 
  filter(MONAT == 4) |> 
  select(-(1:3)) |> 
  as.data.frame() |> 
  observation_transformation(na.rm = T, grid.type = "full")


df_all = tibble(empirical_cov_all, empirical_cov_filtered, observation_grid(144, comp = "full"))
# grafic not contained in paper
plot_ly(df_all, size = .4) |>
  add_markers(x = ~Var1, y = ~Var2, z = ~empirical_cov_all) |> 
  add_markers(x = ~Var1, y = ~Var2, z = ~empirical_cov_filtered, alpha = 0.4)

#### NA Count #####
N |>  summarise(.by = c(JAHR, MONAT, TAG), 
                n = n(), mean = mean(TT_10)) |> 
  filter(n == 144) |> 
  summarise(.by = MONAT, 
            n = n()) |> 
  arrange(MONAT)

cov_estimation = function(month, weights = W, nw =  N_wide){
  y = nw |> 
    filter(MONAT == month) |> 
    select(-(1:3))
  z = y |> 
    as.data.frame() |> # TODO: fix this in biLocPol
    observation_transformation(na.rm = T)
  eval_weights(weights, z)
}


N_wide |>
  filter(MONAT == 3) |> 
  select(-(1:3)) |> is.na() |> sum()


##### calculate weights #####
p.eval = 72
W = local_polynomial_weights(144, 0.2, p.eval, T, m = 1)
Wh05 = local_polynomial_weights(144, 0.5, p.eval, T, m = 1)
Wh01 = local_polynomial_weights(144, 0.1, p.eval, T, m = 1)
g_hat = eval_weights(W, Z)
eval_time = N$UHRZEIT[1:146][seq(2, 144, 2)]

###### Januar ######

# temp curves
N |>  filter(TAG %in% c(1, 15, 29), 
             MONAT == 1) |> 
  ggplot() +
  geom_line(aes(x = UHRZEIT, y = TT_10, group = JAHR*TAG, colour = JAHR), alpha = .7) + 
  labs(y = "Temp. in C°", x = "hours", title = "Temp. in January", colour = "year") + 
  theme(text = element_text(size = 18)) 


g_hat1 = cov_estimation(1)
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



###### August ######
#est8 = 
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
sd_tibble = rbind(var_est1, 
                  var_est4, 
                  var_est, 
                  var_est10) |> 
  as_tibble() |> 
  mutate(var_hat = sqrt(var_hat), 
         month = gl(4, 72, labels = c("Jan", "Apr", "Aug", "Oct"))) |> 
  rename(sd = var_hat)
sd_tibble |> 
  ggplot(aes(x = x, y = sd, col = month, lty = month)) + 
  geom_line(linewidth = .7) + 
  lims(y = c(0.2, 5.8)) + 
  labs(y = NULL, x = "hour") + 
  theme(text = element_text(size = 18), legend.position = "top") 


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


###### std_devaition ######

sd_tibble_m1h02 = sapply(1:12, 
                         function(m){
                           est = cov_estimation(m) |> diag() |> sqrt()
                         }) |> 
  as_tibble() 
sd_tibble_m1h05 = sapply(1:12, 
                         function(m){
                           est = cov_estimation(m, weights = Wh05) |> diag() |> sqrt()
                         }) |> as_tibble() 

sd_tibble_m1h01 = sapply(1:12, 
                         function(m){
                           est = cov_estimation(m, weights = Wh01) |> diag() |> sqrt()
                         }) |> as_tibble() 

sd_tibble = sd_tibble_m1h01 |>
  rbind(sd_tibble_m1h02, sd_tibble_m1h05) |> 
  mutate(time = rep(eval_time, 3), 
         h = gl(3, 72, labels = c("144", "288", "720"))) |> 
  pivot_longer(cols = 1:12, 
               names_to = "month", 
               values_to = "sd", 
               cols_vary = "slowest") |> 
  mutate(month = gl(12, 72, 
                    labels = c("Jan", "Feb", "Mar", "Apr", 
                               "May", "Jun", "Jul", "Aug", 
                               "Sep", "Oct", "Nov", "Dec")) |> rep(each = 3))

sd_tibble$time  = sd_tibble$time |> as.POSIXct(format = "%H:%M")

sd_tibble |> 
  ggplot(aes(x = time, y = sd, lty = h, col = h)) + 
  geom_line(linewidth = .8) + 
  lims(y = c(0.2, 6)) + 
  labs(y = NULL, x = NULL) + 
  theme(text = element_text(size = 18)) + 
  facet_wrap(.~month, nrow = 3) + 
  scale_linetype_manual(values = c(2,5,4), name = "h (min)") + 
  scale_color_manual(values = 1:3, name = "h (min)") + 
  scale_x_datetime(date_breaks = "8 hours", date_labels = "%H:%M")

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

#

