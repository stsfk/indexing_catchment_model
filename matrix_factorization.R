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
  mutate(NNSE = 1 / (2 - NSE),# normalized NSE
         catchment_id = factor(catchment_id, level = 1:length(unique(catchment_id))),
         model_id = factor(model_id, level = 1:length(unique(model_id))),
  ) %>% 
  dplyr::select(
    catchment_id,
    model_id,
    NNSE) %>%
  rename(rating = NNSE)

# spliting data
v <- 5
cv_folds <- vfold_cv(data_process, v = 5)

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

flag

# start training
i <- 1

dtrain <- analysis(cv_folds$splits[[i]])
dtest <- assessment(cv_folds$splits[[i]])

# construct recommender model
r = Reco()

training_set <- data_memory(
  user_index = dtrain$catchment_id,
  item_index = dtrain$model_id,
  rating = dtrain$rating,
  index1 = T
)

opts = r$tune(training_set, opts = list(
  dim = c(10),
  lrate = c(0.001, 0.005, 0.1, 0.2),
  nthread = 10,
  niter = 20
)) # 1:10

r$train(training_set, opts = c(opts$min, nthread = 10, niter = 20))

c(P, Q) %<-% r$output(out_memory(), out_memory())

dim(P)
dim(Q)

test_set <- data_memory(
  user_index = dtest$user_index,
  item_index = dtest$item_index,
  rating = dtest$rating
)
pred_rvec <- r$predict(test_set)

gof(pred_rvec, dtest$rating)

# code length vs. prediction errors
opts$res %>%
  as_tibble()%>%
  group_by(dim) %>%
  summarise(loss = min(loss_fun)) %>%
  ggplot(aes(dim, loss)) +
  geom_point()+
  geom_line()+
  scale_x_continuous(breaks = c(1:10)) +
  labs(x = "dimensions",
       y = "loss")
