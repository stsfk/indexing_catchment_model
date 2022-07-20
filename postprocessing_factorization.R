if (!require("pacman")) {
  install.packages("pacman")
}

pacman::p_load(tidyverse,
               lubridate,
               zeallot,
               recosystem,
               hydroGOF,
               caret,
               tidymodels,
               mlrMBO,
               ModelMetrics)


# data --------------------------------------------------------------------

load("./data/dens_factorization.Rda")

dense_gof <- eval_grid %>%
  select(r2, rmse, repeats) %>%
  mutate(case = "dense")


load("./data/sparse_factorization.Rda")

sparse_gof <- eval_grid %>%
  select(r2, rmse, repeats) %>%
  mutate(case = "dense")

sparse_gof$r2 %>% mean()
sparse_gof$r2 %>% var()

sparse_gof$rmse %>% mean()
sparse_gof$rmse %>% var()
