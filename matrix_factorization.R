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

# spliting data
v <- 5
cv_folds <- vfold_cv(data_process, v = v)

# check if all catchments included in the training and test sets
flag <- 0 # flag = 0 means the number of classes in training and test folds are the same

for (i in 1:nrow(cv_folds)){
  training_set <- analysis(cv_folds$splits[[i]])
  test_set <- assessment(cv_folds$splits[[i]])
  
  temp1 <- training_set %>% count(catchment_id) %>% nrow()
  temp2 <- test_set %>% count(catchment_id) %>% nrow()
  
  temp3 <- training_set %>% count(model_id) %>% nrow()
  temp4 <- test_set %>% count(model_id) %>% nrow()
  
  flag <- flag + temp1 - temp2 + temp3 - temp4
}

flag # pass the test if flag == 0

# storage for experiment results in CV iterations
Ps <- vector("list", nrow(cv_folds))
Qs <- vector("list", nrow(cv_folds))
Optim_paras <- vector("list", nrow(cv_folds))
Optim_res <- vector("list", nrow(cv_folds))
rs <- vector("list", nrow(cv_folds))

# iterate over CV folds
for (i in 1:nrow(cv_folds)) {
  # get training and test fold
  dtrain <- analysis(cv_folds$splits[[i]])
  training_set <- data_memory(
    user_index = dtrain$catchment_id,
    item_index = dtrain$model_id,
    rating = dtrain$rating,
    index1 = T
  )
  
  dtest <- assessment(cv_folds$splits[[i]])
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
                  dim = c(1:25) * 2,
                  lrate = c(0.001, 0.005, 0.1, 0.2),
                  nthread = 10,
                  niter = 20
                ))
  
  r$train(training_set, opts = c(opts$min, nthread = nthread, niter = 20))
  
  c(P, Q) %<-% r$output(out_memory(), out_memory())
  
  # save result
  Ps[[i]] <- P
  Qs[[i]] <- Q
  Optim_paras[[i]] <- opts$min
  Optim_res[[i]] <- opts$res
  rs[[i]] <- r
}

save(cv_folds, Ps, Qs, Optim_paras, Optim_res, rs, file = "./data/mf_results.Rda")

# Postprocessing ----------------------------------------------------------


dim(P)
dim(Q)

pred_rvec <- r$predict(test_set)

gof(pred_rvec, dtest$rating)


as_tibble()%>%
  group_by(dim) %>%
  summarise(loss = min(loss_fun)) %>%
  ggplot(aes(dim, loss)) +
  geom_point()+
  geom_line()+
  labs(x = "dimensions",
       y = "loss")
# code length vs. prediction errors
opts$res %>%
  as_tibble()%>%
  group_by(dim) %>%
  summarise(loss = min(loss_fun)) %>%
  ggplot(aes(dim, loss)) +
  geom_point()+
  geom_line()+
  labs(x = "dimensions",
       y = "loss")
