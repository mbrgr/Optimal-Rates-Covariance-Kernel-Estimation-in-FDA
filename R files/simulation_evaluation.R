# Code contains evaluation of 'full_vs_mirrored_estimator.R' and 'Simulations_Covariance.R' files.
# Further Figures 3, 4, 5, 6 and 8 are produced


load("data/bw_comp_OU.RData")
bw_comparison_tbl = Reduce(rbind, bw_comparison) |>
  as_tibble() |>
  rename(n = V1, p = V2, h = V3, sup.err = V4) |>
  mutate(p = as.factor(p), n = as.factor(n)) 

bw_comparison_tbl |> 
  filter(n != 50) |> 
  ggplot() + 
  geom_line(aes(x = h, y = sup.err, col = n, linetype = p), size = .6) + 
  ylim(c(0, 0.65)) +
  theme(text = element_text(size = 18)) 

# Figure 3:
bw_comparison_tbl |> 
  filter(n == 400) |> 
  ggplot() + 
  geom_line(aes(x = h, y = sup.err, col = p, linetype = p), size = .6) + 
  ylim(c(0.02, 0.52)) +
  theme(text = element_text(size = 18)) + 
  labs(title = "n = 400")

# add: dissertation
bw_comparison_tbl |> 
  filter(n == 400) |> 
  ggplot() + 
  geom_point(aes(x = h, y = sup.err, col = p, pch = p)) + 
  ylim(c(0.02, 0.52)) +
  labs(title = "n = 400")

ggsave("Grafics/cov_optimal_bw_n400.png", device = "png", width = 5, height = 3.8, units = "in")

# add: dissertation
bw_comparison_tbl |> 
  ggplot() + 
  geom_point(aes(x = h, y = sup.err, col = p, pch = p)) + 
  lims(y = c(0, .95)) +
  facet_wrap(n~., nrow = 1)
ggsave("Grafics/cov_optimal_bw_various_n.png", device = "png", width = 8, height = 4, units = "in")



bw_comparison_tbl |> 
  filter(n == 100, p == 75) |>
  ggplot() + 
  geom_line(aes(x = h, y = sup.err, linetype = p)) + 
  ylim(c(0, 0.9))

min_h_tibble = bw_comparison_tbl |>
  group_by(n, p) |> 
  slice_min(sup.err) 

bw_comparison_tbl |> 
  filter(n == 100) |> 
  ggplot( ) + 
  geom_line(aes(x = h, y = sup.err, col = p, linetype = p)) +
  ylim(c(0.03, 0.65)) + 
  labs(title = "n = 100") 

##### One -Fold CV #####
# not in paper
load("data/one_fold_cv_n400_OU.RData")
one_fold_tbl_n400 = tibble(h = one_fold_cv_n400 %>% unlist(), p = gl(3, 1000, labels = p.seq))
one_fold_tbl_n400 %>% 
  ggplot(aes(y = h, x = p, col = p)) + 
  geom_boxplot()

n400_table = one_fold_tbl_n400 %>% 
  group_by_all() %>% 
  summarise(n = n()/1000)

n400_table %>% 
  ggplot(aes(h, n)) + 
  geom_point(size = 3) + 
  ylim(c(0, 0.1))+
  facet_wrap(.~p)

one_fold_tbl_n400 %>% 
  summarise(.by = p, mean(h))



##### five fold cv #####
load("data/five_fold_cv_OU.RData")
five_fold_tbl = tibble(h = five_fold_cv %>% unlist(), p = gl(5, 1000, labels = p.seq))
five_fold_tbl %>% 
  ggplot(aes(y = h, x = p, col = p)) + 
  geom_boxplot()


five_fold_table = five_fold_tbl %>% 
  group_by_all() %>% 
  summarise(n = n()/1000)

# Figure in Appendix
five_fold_table %>% 
  ggplot(aes(h, n)) + 
  geom_point(size = 3) + 
  ylim(c(0, 0.1))+
  facet_wrap(.~p) +
  ylab(NULL) + 
  labs(title = "n = 100")

five_fold_tbl %>% 
  summarise(.by = p, mean(h))

##### five fold cv #####
# In Paper
five_fold_tbl_n400 = tibble(h = five_fold_cv_n400 %>% unlist(), p = gl(5, 1000, labels = p.seq))
# Figure 4
five_fold_tbl_n400 %>% 
  ggplot(aes(y = h, x = p, col = p)) + 
  geom_boxplot(size = .6) + 
  labs(title = "n = 400") #+
 # theme(text = element_text(size = 18)) 
ggsave("Grafics/cov_5fcv_bw_n400.png", device = "png", width = 5, height = 3.8, units = "in")



five_fold_table_n400 = five_fold_tbl_n400 %>% 
  group_by_all() %>% 
  summarise(n = n()/1000)

# Figure 5
five_fold_table_n400 %>% 
  ggplot(aes(h, n)) + 
  geom_point(size = 3) + 
  facet_wrap(.~p, nrow = 1) +
  ylab(NULL) + 
  labs(title = "n = 400") #+
#  theme(text = element_text(size = 18)) 
ggsave("Grafics/cov_5fcv_bw_n400_table.png", device = "png", width = 5, height = 3.8, units = "in")



five_fold_tbl_n400 %>% 
  summarise(.by = p, mean(h))

##### Error Decomposition #####
load("data/error_decomp.RData")
# error_decomp_arr = error_decomp_arr[, -5, ]
error_decomp_arr |> dimnames() = list(p.seq, c("eps", "dsc", "prc", "mix", "sup"), n.seq)
error_decomp_tbl = error_decomp_arr |> as_tibble() |> 
  pivot_longer(cols = everything(),
               names_to = c("term", "n"), 
               names_pattern = "(...).(.*)",
               names_transform = list(term = as.factor, n = as.double),
               values_to = "error") |> 
  mutate(p = rep(p.seq, each = 20))
         
error_decomp_tbl |> 
  filter(n == 400, term == "sup") |> 
  print(n = 100) 

# Figure 6
error_decomp_tbl |> 
  filter(p != 75) |> 
  ggplot(aes(x = n, y = error, col = term, lty = term)) + 
  geom_line(size = .6) + 
  facet_wrap(.~p, nrow = 1) +
  theme(text = element_text(size = 18)) +
  labs(y = "sup.error")

# Not in Paper
error_decomp_tbl |> 
  ggplot(aes(x = p, y = error, col = term, lty = term)) + 
  geom_line() + 
  facet_wrap(.~n)

#### Estimator Comparison ####
load("data/estimator_comparison.RData")

to_tibble = function(bw_comp, est = "mir", m = 1, Z = "OU") {
  Reduce(rbind, bw_comp) |>
    as_tibble() |>
    rename(n = V1, p = V2, h = V3, sup.err = V4) |>
    mutate(p = as.factor(p), 
           n = as.factor(n), 
           estimator = est,
           m = m,
           Z = Z)
}

# alrady contained in dataset
OU_m0 = to_tibble(bw_comparison_OU_m0, m = 0)
OU_m1 = to_tibble(bw_comparison_OU_m1)
OU_m2 = to_tibble(bw_comparison_OU_m2, m = 2)
OU_m0_full = to_tibble(bw_comparison_OU_m0_full, m = 0, est = "wd")
OU_m1_full = to_tibble(bw_comparison_OU_m1_full, est = "wd")
OU_m2_full = to_tibble(bw_comparison_OU_m2_full, m = 2, est = "wd")
rv2_m0 = to_tibble(bw_comparison_OU_m0_2rv, m = 0, Z  = "2rv")
rv2_m0_full = to_tibble(bw_comparison_OU_m0_full_2rv, m = 0, est = "wd", Z  = "2rv")
rv2_m1 = to_tibble(bw_comparison_OU_m1_2rv, Z  = "2rv")
rv2_m1_full = to_tibble(bw_comparison_OU_m1_full_2rv, est = "wd", Z  = "2rv")
rv2_m2 = to_tibble(bw_comparison_OU_m2_2rv, Z  = "2rv", m = 2)
rv2_m2_full = to_tibble(bw_comparison_OU_m2_full_2rv, est = "wd", Z  = "2rv", m = 2)

est_comp = rbind(OU_m0, OU_m1,
                 OU_m0_full, OU_m1_full,
                 OU_m2, OU_m2_full, 
                 rv2_m0, rv2_m1,
                 rv2_m0_full, rv2_m1_full,
                 rv2_m2, rv2_m2_full) |> 
  mutate(estimator = as.factor(estimator), 
         m = as.factor(m), Z = as.factor(Z))

# Figure 8a
est_comp |> 
  filter(n == 100, p == 50, Z == "OU") |> 
  ggplot(aes(x = h, y = sup.err, lty = estimator, col = m)) + 
  geom_line(size = .6) + 
  lims(y = c(0.03, 1.2)) +
  theme(text = element_text(size = 18)) + 
  labs(title = "Ornstein-Uhlenbeck", subtitle = "n = 100, p = 50")

# Figure 8b
est_comp |> 
  filter(n == 100, p == 50, Z == "2rv") |> 
  ggplot(aes(x = h, y = sup.err, lty = estimator, col = m)) + 
  geom_line(size = .6) + 
  lims(y = c(0.03, 1.1))+
  theme(text = element_text(size = 18)) + 
  labs(title = "Process 2", subtitle = "n = 100, p = 50")

# Figure 8a
est_comp |> 
  filter(n == 100, p == 50, Z == "OU") |> 
  ggplot(aes(x = h, y = sup.err, pch = estimator, col = m)) + 
  geom_point() + 
  lims(y = c(0.03, 1.2)) +
  labs(title = "Ornstein-Uhlenbeck", subtitle = "n = 100, p = 50")
ggsave("Grafics/est_comp_OU_points.png", device = "png", width = 5, height = 3.8, units = "in")


# Figure 8b
est_comp |> 
  filter(n == 100, p == 50, Z == "2rv") |> 
  ggplot(aes(x = h, y = sup.err, pch = estimator, col = m)) + 
  geom_point() + 
  lims(y = c(0.03, 1.1))+
  labs(title = "Process 2", subtitle = "n = 100, p = 50")
ggsave("Grafics/est_comp_2rv_points.png", device = "png", width = 5, height = 3.8, units = "in")



