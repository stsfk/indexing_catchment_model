if (!require("pacman")) {
  install.packages("pacman")
}

pacman::p_load(
  tidyverse
)

# Sparse experiments on Dense data set
load("data/sparse_exp_summary.Rda")

data_dens_sub <- eval_grid_summary %>% transmute(
  density = paste0(ratio * 100, "%"),
  r2 = paste0(
    sprintf("%.3f", round(r2_mean, 3)),
    "(",
    formatC(r2_sd, format = "e", digits = 2),
    ")"
  ),
  rmse = paste0(
    sprintf("%.3f", round(rmse_mean, 3)),
    "(",
    formatC(rmse_sd, format = "e", digits = 2),
    ")"
  )
) %>%
  mutate(dataset="Dense")

# Dense data set whole
load("./data/dens_factorization.Rda")

eval_grid <- eval_grid %>%
  filter(!is.na(r2))

eval_grid$r2 %>% mean(na.rm = T)
eval_grid$r2 %>% var(na.rm = T)

eval_grid$rmse %>% mean(na.rm = T)
eval_grid$rmse %>% var(na.rm = T)

data_dens_whole <- tibble(
  r2_mean = eval_grid$r2 %>% mean(),
  r2_sd = eval_grid$r2 %>% sd(),
  rmse_mean = eval_grid$rmse %>% mean(),
  rmse_sd = eval_grid$rmse %>% sd()
) %>% transmute(
  r2 = paste0(
    sprintf("%.3f", round(r2_mean, 3)),
    "(",
    formatC(r2_sd, format = "e", digits = 2),
    ")"
  ),
  rmse = paste0(
    sprintf("%.3f", round(rmse_mean, 3)),
    "(",
    formatC(rmse_sd, format = "e", digits = 2),
    ")"
  )
) %>%
  mutate(density = "100%",
         dataset="Dense")

# Caravan data set
load("./data/caravan_mf_eval_gof.Rda")

data_intense <- tibble(
  r2_mean = eval_grid2$r2 %>% mean(),
  r2_sd = eval_grid2$r2 %>% sd(),
  rmse_mean = eval_grid2$rmse %>% mean(),
  rmse_sd = eval_grid2$rmse %>% sd()
) %>% transmute(
  r2 = paste0(
    sprintf("%.3f", round(r2_mean, 3)),
    "(",
    formatC(r2_sd, format = "e", digits = 2),
    ")"
  ),
  rmse = paste0(
    sprintf("%.3f", round(rmse_mean, 3)),
    "(",
    formatC(rmse_sd, format = "e", digits = 2),
    ")"
  )
) %>%
  mutate(density = "5%",
         dataset="Intense")

# Sparse data set

load("./data/sparse_factorization.Rda")

eval_grid <- eval_grid %>%
  filter(!is.na(r2))

data_sparse <- tibble(
  r2_mean = eval_grid$r2 %>% mean(),
  r2_sd = eval_grid$r2 %>% sd(),
  rmse_mean = eval_grid$rmse %>% mean(),
  rmse_sd = eval_grid$rmse %>% sd()
) %>% transmute(
  r2 = paste0(
    sprintf("%.3f", round(r2_mean, 3)),
    "(",
    formatC(r2_sd, format = "e", digits = 2),
    ")"
  ),
  rmse = paste0(
    sprintf("%.3f", round(rmse_mean, 3)),
    "(",
    formatC(rmse_sd, format = "e", digits = 2),
    ")"
  )
) %>%
  mutate(density = "5%",
         dataset="Sparse")

# Combine -----------------------------------------------------------------

data_dens_sub %>%
  bind_rows(data_dens_whole) %>%
  bind_rows(data_sparse) %>%
  bind_rows(data_intense) %>% 
  select(dataset, everything()) %>%
  view()



