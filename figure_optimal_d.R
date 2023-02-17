if (!require("pacman")) {
  install.packages("pacman")
}
pacman::p_load(
  tidyverse
)

# Dense data set whole
load("./data/dens_factorization.Rda")

data_dens_whole <- eval_grid %>% select(r2, rmse) %>%
  mutate(dims = sapply(eval_grid$out, function(x)
    x$P %>% dim() %>% .[2])) %>%
  mutate(ratio = 0.05,
         dataset = "Dense")

# Dense data set sub
load("./data/sparse_exp.Rda")

data_dens_sub <- eval_grid %>% select(ratio, r2, rmse) %>%
  mutate(dims = sapply(eval_grid$out, function(x)
    x$P %>% dim() %>% .[2])) %>%
  mutate(dataset = "Dense")

# Sparse data set
load("./data/sparse_factorization.Rda")

data_sparse <- eval_grid %>% select(r2, rmse) %>%
  mutate(dims = sapply(eval_grid$out, function(x)
    x$P %>% dim() %>% .[2])) %>%
  mutate(ratio = 0.05) %>%
  mutate(dataset = "Sparse")

# Intense data set
load("./data/caravan_mf_opt_para.Rda")
load("./data/caravan_mf_eval_gof.Rda")

data_intense <- eval_grid2 %>%
  mutate(dims = sapply(opt_para, function(x) x$dim)) %>%
  mutate(ratio = 0.05) %>%
  mutate(dataset = "Intense")

# combine -----------------------------------------------------------------

data_plot <- data_dens_whole %>%
  bind_rows(data_dens_sub) %>%
  bind_rows(data_sparse) %>%
  bind_rows(data_intense) %>% 
  select(dataset, everything())

save(data_plot, "./data/optimal_d.Rda")
