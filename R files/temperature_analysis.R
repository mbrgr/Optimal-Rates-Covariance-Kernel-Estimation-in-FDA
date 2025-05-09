# Code for real data example with the weather data from Nuremberg 
# Figures 9, 10 and 11 are generated her

library(tidyverse)
library(plotly)

library(lubridate)
library(hms)

library(biLocPol) # please install this package from Github first. See "README.md" file for instructions
library(future.apply)


# Parameter description
# @STATIONS_ID identification number of the weather station 
# MESS_DATUM time of the observation given by yyyymmddhhmi 

# PP_10 air pressure at the height of the station
# TT_10 air temperature at the height of 2 metres (°C)
# TM5_10 air temperatur at the height of 5 centimetres (°C)
# RF_10 relative humidity at 2 metres
# TD_10 dew point temperature at 2 metres (°C)
# eor Ende data record

# Missing values are labelled with -999. The measurements before the year 2000 are assigned a time stamp in
# CET, from the year 2000 onwards a time stamp in UTC. The dew point temperature is calculated from the
# air temperature at an altitude of 2 metres and the relative humidity measurements. The values are mean values over
# the minute ending at the time stamp. The time stamp is given before the year 2000 in CET, from the
# year 2000 in UTC


load("data/weather_data_raw.RData")
N$DATUM = str_split_i(N$MESS_DATUM, pattern = " ", i = 1) |> ymd()

tage = c(1, 4, 8, 12, 15, 18, 22, 25, 29)
N_full = N[-(1:6),]
N = N %>% 
  filter(TAG %in% tage)


N |>  
  filter(TAG %in% c(1, 8, 15, 22, 29)) |> 
  ggplot() +
  geom_line(aes(x = UHRZEIT, y = TT_10, group = JAHR*TAG, colour = JAHR), alpha = .4) +
  facet_wrap(MONAT ~.)

# Figure 9b: Temp. curves in August
N |>  
  filter(TAG %in% c(1, 15, 29), 
         MONAT == 8) |> 
  ggplot() +
  geom_line(aes(x = UHRZEIT, y = TT_10, group = JAHR*TAG, colour = JAHR), alpha = .6) + 
  labs(y = "Temp. in C°", x = "hours", title = "Temp. in August", colour = "year") + 
  theme(text = element_text(size = 14)) 

N |>  
  filter(TAG %in% tage, 
         MONAT == 1) |> 
  ggplot() +
  geom_line(aes(x = UHRZEIT, y = TT_10, group = JAHR*TAG, colour = JAHR), alpha = .6) + 
  labs(y = "Temp. in C°", x = "hours", title = "Temp. in January", colour = "year") + 
  theme(text = element_text(size = 14)) 
ggsave("grafics/january_temp_curves.png", device = "png", width = 5, height = 4, unit = "in")

N |>  
  filter(TAG %in% tage, 
         MONAT == 7) |> 
  ggplot() +
  geom_line(aes(x = UHRZEIT, y = TT_10, group = JAHR*TAG, colour = JAHR), alpha = .6) + 
  labs(y = "Temp. in C°", x = "hours", title = "Temp. in July", colour = "year") + 
  theme(text = element_text(size = 14))
ggsave("grafics/july_temp_curves.png", device = "png", width = 5, height = 4, unit = "in")


N_wide = N |> 
  filter(TAG %in% tage) |>  # No consecutive days
  mutate(UHRZEIT = as.character(UHRZEIT)) |> 
  select(JAHR, MONAT, TAG, UHRZEIT, TT_10) |> 
  pivot_wider(names_from = UHRZEIT,
              values_from = TT_10)

#### Plot empirical covariance for consecutive days vs. non consecutive #####

Y_all = N_full |> 
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
eval_time = N$UHRZEIT[1:146][seq(2, 144, 2)]

###### January ######
# Figure 9a: temp curves in January
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


cov_est_df1 = data.frame(x = W$x.eval, y = W$x.eval, z = g_hat1)
cs2 = list(c(0, 1), c("lightblue", "darkred"))

# Plot of Covariance estimation: Not in Paper
plot_ly(cov_est_df1, x = ~x, y = ~y, z = ~g_hat1, size = .4) |> 
  add_surface(colorscale = cs2, alpha = .3) |> 
  layout(scene = list(xaxis = list(title = ""), 
                      yaxis = list(title = ""), 
                      zaxis = list(title = "")))

temp = matrix(diag(g_hat1), p.eval, p.eval)
cor_hat1 = g_hat1 / sqrt( temp * t(temp) )

# Figure 11a
plot_ly(cov_est_df1, x = ~x*24, y = ~y*24, z = ~cor_hat1, size = .4) |> 
  add_surface(colorscale = cs2, alpha = .3) |> 
  layout(scene = list(xaxis = list(title = ""), 
                      yaxis = list(title = ""), 
                      zaxis = list(title = "")))



###### August ######
g_hat = cov_estimation(8)
var_hat = diag(g_hat)
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

# Plot of Covariance Kernel: Not in Paper
plot_ly(cov_est_df, x = ~x, y = ~y, z = ~g_hat, size = .4) |> 
  add_surface(colorscale = cs2, alpha = .3) |> 
  layout(scene = list(xaxis = list(title = ""), 
                      yaxis = list(title = ""), 
                      zaxis = list(title = "")))

temp = matrix(diag(g_hat), p.eval, p.eval)
cor_hat = g_hat / sqrt( temp * t(temp) )

# Figure 11b: Plot of Correlation in August.
plot_ly(cov_est_df, x = ~x*24, y = ~y*24, z = ~cor_hat, size = .4) |> 
  add_surface(colorscale = cs2, alpha = .3) |> 
  layout(scene = list(xaxis = list(title = ""), 
                      yaxis = list(title = ""), 
                      zaxis = list(title = "")))


###### July ######
g_hat = cov_estimation(7)
var_hat = diag(g_hat)
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

# Plot of Covariance Kernel: Not in Paper
plot_ly(cov_est_df, x = ~x, y = ~y, z = ~g_hat, size = .4) |> 
  add_surface(colorscale = cs2, alpha = .3) |> 
  layout(scene = list(xaxis = list(title = ""), 
                      yaxis = list(title = ""), 
                      zaxis = list(title = "")))

temp = matrix(diag(g_hat), p.eval, p.eval)
cor_hat = g_hat / sqrt( temp * t(temp) )

# Plot of Correlation in July
plot_ly(cov_est_df, x = ~x*24, y = ~y*24, z = ~cor_hat, size = .4) |> 
  add_surface(colorscale = cs2, alpha = .3) |> 
  layout(scene = list(xaxis = list(title = ""), 
                      yaxis = list(title = ""), 
                      zaxis = list(title = "")))
###### std_deviation ######

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

# Figure 10: all standard deviations. 
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
