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
               matlib)

# seed --------------------------------------------------------------------

set.seed(1234)
nthread <- 10 # number of CPU thread

# data --------------------------------------------------------------------

similarity_m <- read_csv("./data/SM.csv",
                         col_names = c("KGE", "NSE", "RMSE"))

model_list  <- c(
  'm_01_collie1_1p_1s',
  'm_05_ihacres_7p_1s',
  'm_07_gr4j_4p_2s',
  'm_13_hillslope_7p_2s',
  'm_18_simhyd_7p_3s',
  'm_22_vic_10p_3s',
  'm_27_tank_12p_4s',
  'm_28_xinanjiang_12p_4s',
  'm_34_flexis_12p_5s',
  'm_37_hbv_15p_5s'
)

n_catchment <- 533
n_paras <- nrow(similarity_m)/n_catchment/length(model_list) # 533 catchments, n_paras sets of paras for each model
catchment_id <- rep(1:n_catchment, each = n_paras)

data_raw <- similarity_m %>%
  bind_cols(expand_grid(model_list, catchment_id)) %>%
  rename(model_name = model_list) %>%
  mutate(
    para_set_id = rep(1:n_paras, n() / n_paras),
    model_id = paste(model_name, para_set_id, sep = "_"),
    model_id = factor(model_id, levels = unique(model_id)),
    model_id = as.integer(model_id) # assign a unique id to each model
  ) %>%
  select(catchment_id, model_name, model_id, KGE, NSE, RMSE)

# Recommender system ------------------------------------------------------

data_process <- data_raw %>%
  mutate(NNSE = 1 / (2 - NSE)*10) %>% 
  dplyr::select(
    catchment_id,
    model_id,
    NNSE) %>%
  rename(rating = NNSE)

data_process <- data_process %>%
  sample_frac(0.05)

# splitting data

data_split <- initial_split(data_process, prop = 4/5, strata = NULL)

dtrain <- analysis(data_split)
dtest <- assessment(data_split)

# get training and test fold
training_set <- data_memory(
  user_index = dtrain$catchment_id,
  item_index = dtrain$model_id,
  rating = dtrain$rating,
  index1 = T
)

test_set <- data_memory(
  user_index = dtest$catchment_id,
  item_index = dtest$model_id,
  rating = dtest$rating,
  index1 = T
)

# construct recommender model
r = Reco()

opts = r$tune(training_set,
              opts = list(
                dim = c(10,25,50), #c(1:25) * 2,
                lrate = c(0.001, 0.005, 0.1, 0.2),
                nthread = 10,
                niter = 20
              ))

r$train(training_set, opts = c(opts$min, nthread = nthread, niter = 20))

c(P, Q) %<-% r$output(out_memory(), out_memory())

# Postprocessing ----------------------------------------------------------

# retrieve data and model)

# predicting test set
pred_rvec <- r$predict(test_set)

# goodness-of-fit
rmse <- ModelMetrics::rmse(actual = dtest$rating, predicted = pred_rvec)
r2 <- cor(dtest$rating, pred_rvec)^2
mae <- hydroGOF::mae(sim = pred_rvec, obs = dtest$rating)



