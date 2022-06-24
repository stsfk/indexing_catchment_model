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
               mlrMBO
               )

# seed --------------------------------------------------------------------

set.seed(1234)
nthread <- 10 # number of CPU thread

# data --------------------------------------------------------------------

similarity_m <- read_csv("./data/SM.csv",
                         col_names = c("KGE", "NSE", "RMSE"))

model_class  <- c(
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

n_catchments <- 533
n_model_classes <- length(model_class)
n_model_instances <-  nrow(similarity_m)/n_catchments/n_model_classes # number of instance of each model class

catchment_id <- rep(1:n_catchments, each = n_model_instances)
data_raw <- similarity_m %>%
  bind_cols(expand_grid(model_class, catchment_id)) %>%
  mutate(
    instance_id = rep(1:n_model_instances, n()/n_model_instances),
    model_id = paste(model_class, instance_id, sep = "_"),
    model_id = factor(model_id, levels = unique(model_id)),
    model_id = as.integer(model_id) # assign a unique id to each model
  ) %>%
  select(model_class, catchment_id, model_id, KGE, NSE, RMSE)

data_process <- data_raw %>%
  mutate(NNSE = 1 / (2 - NSE) * 10) %>%
  dplyr::select(catchment_id,
                model_id,
                NNSE) %>%
  rename(rating = NNSE) %>%
  mutate(record_id = 1:n())


# Experiments -------------------------------------------------------------

frac <- 0.1
train_portion <- 0.6
val_portion <- 0.2
test_portion <- 0.2

data_sample <- data_process %>%
  group_by(model_id) %>% 
  sample_frac(frac) %>%
  ungroup()

data_train_val <- data_sample %>%
  group_by(model_id) %>% 
  sample_frac(train_portion + val_portion) %>%
  ungroup()

data_train <- data_train_val %>%
  group_by(model_id) %>% 
  sample_frac(train_portion/(val_portion+train_portion)) %>%
  ungroup()

data_val <- data_train_val %>%
  filter(record_id %in% setdiff(data_train_val$record_id, data_train$record_id))

data_test <- data_sample %>%
  filter(record_id %in% setdiff(data_sample$record_id, data_train_val$record_id))

list(
  data_train = data_train %>% select(-record_id),
  data_val = data_val %>% select(-record_id),
  data_test = data_test %>% select(-record_id),
  record_id_train = data_train$record_id,
  record_id_val = data_val$record_id,
  record_id_test = data_test$record_id
)











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

# experiment
scoringFunction <- function(x){
  
  dim <- x["dim"] %>% unlist()
  lrate <- x["lrate"] %>% unlist()
  niter <- x["niter"] %>% unlist()
  
  costp_l1 <- x["costp_l1"] %>% unlist()
  costp_l2 <- x["costp_l2"] %>% unlist()
  costq_l1 <- x["costq_l1"] %>% unlist()
  costq_l2 <- x["costq_l2"] %>% unlist()
  
  dim <- 5 * dim
  lrate <- 10^lrate
  niter <- 10 * niter
  
  r = Reco()
  
  r$train(training_set, opts = c(opts$min, nthread = nthread, niter = niter, verbose = F))
  
  pred_rvec <- r$predict(test_set)
  r2 <- cor(dtest$rating, pred_rvec)^2
  
  score <- r2
  return(score)
} 




obj_fun <- makeSingleObjectiveFunction(
  fn = scoringFunction,
  par.set = makeParamSet(
    makeIntegerParam("dim", lower= 1, upper = 20),
    makeNumericParam("lrate",  lower= -4,   upper = -0.5),
    makeIntegerParam("niter", lower = 1,  upper = 20),
    makeNumericParam("costp_l1",  lower= 0,   upper = 1),
    makeNumericParam("costp_l2",  lower= 0,   upper = 1),
    makeNumericParam("costq_l1",  lower= 0,   upper = 1),
    makeNumericParam("costq_l2",  lower= 0,   upper = 1)
  ),
  has.simple.signature = FALSE,
  minimize = FALSE
)

des = generateDesign(
  n = 4 * getNumberOfParameters(obj_fun),
  par.set = getParamSet(obj_fun),
  fun = lhs::randomLHS
)

des$y = apply(des, 1, obj_fun)

control <- makeMBOControl() %>%
  setMBOControlTermination(., iters = 100 - 4 * getNumberOfParameters(obj_fun))

run <- mbo(
  fun = obj_fun,
  design = des,
  control = control,
  show.info = TRUE
)

run$best.ind
run$y

plot(run)






# recycle


# construct recommender model
r = Reco()

opts = r$tune(training_set,
              opts = list(
                dim = c(5, 10,20,50), #c(1:25) * 2,
                lrate = c(0.001, 0.005, 0.1, 0.2),
                nthread = 10,
                niter = 30,
                verbose = T
              ))

r$train(training_set, opts = c(opts$min, nthread = nthread, niter = 1000))
c(P, Q) %<-% r$output(out_memory(), out_memory())

# Postprocessing ----------------------------------------------------------
pred_rvec <- r$predict(test_set)
cor(dtest$rating, pred_rvec)^2

# retrieve data and model)

# predicting test set
pred_rvec <- r$predict(test_set)

# goodness-of-fit
rmse <- ModelMetrics::rmse(actual = dtest$rating, predicted = pred_rvec)
r2 <- cor(dtest$rating, pred_rvec)^2
mae <- hydroGOF::mae(sim = pred_rvec, obs = dtest$rating)

