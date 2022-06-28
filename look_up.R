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
               DiceKriging,
               rgenoud,
               GA,
               ModelMetrics,
               doParallel)

# seed --------------------------------------------------------------------

#set.seed(1234)
nthread <- 14 # number of CPU thread

# data --------------------------------------------------------------------

weights <- read_csv("./data/SM.csv",
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
n_model_instances <-  nrow(weights)/n_catchments/n_model_classes # number of instance of each model class

catchment_id <- rep(1:n_catchments, each = n_model_instances)
data_raw <- weights %>%
  bind_cols(expand_grid(model_class, catchment_id)) %>%
  mutate(
    instance_id = rep(1:n_model_instances, n() / n_model_instances),
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

# keep 500 catchments for model building and 30 catchments for look-up
catchment_id_building <- sample(1:n_catchments, 500) %>% sort()
catchment_id_lookup <- setdiff(1:n_catchments, catchment_id_building)

data_process_building <- data_process %>%
  filter(catchment_id %in% catchment_id_building)

data_process_lookup <- data_process %>%
  filter(catchment_id %in% catchment_id_lookup)


# Building model ----------------------------------------------------------

train_portion = 0.8
val_portion = 0.2

data_train <- data_process_building %>%
  group_by(model_id) %>%
  sample_frac(train_portion) %>%
  ungroup()

data_val <- data_process_building %>%
  filter(record_id %in% setdiff(data_process_building$record_id, data_train$record_id))

record_id_train <- data_train$record_id
record_id_val <- data_val$record_id

train_set <- data_memory(
  user_index = data_train$catchment_id,
  item_index = data_train$model_id,
  rating = data_train$rating,
  index1 = T
)

val_set <- data_memory(
  user_index = data_val$catchment_id,
  item_index = data_val$model_id,
  rating = data_val$rating,
  index1 = T
)

train_val_set <- data_memory(
  user_index = data_process_building$catchment_id,
  item_index = data_process_building$model_id,
  rating = data_process_building$rating,
  index1 = T
)

scoringFunction <- function(x) {
  dim <- x["dim"] %>% unlist()
  lrate <- x["lrate"] %>% unlist()
  niter <- x["niter"] %>% unlist()
  
  costp_l1 <- x["costp_l1"] %>% unlist()
  costp_l2 <- x["costp_l2"] %>% unlist()
  costq_l1 <- x["costq_l1"] %>% unlist()
  costq_l2 <- x["costq_l2"] %>% unlist()
  
  r = Reco()
  
  r$train(
    train_set,
    opts = list(
      dim = dim,
      costp_l1 = costp_l1,
      costp_l2 = costp_l2,
      costq_l1 = costq_l1,
      costq_l2 = costq_l2,
      lrate = lrate,
      niter = niter,
      verbose = F,
      nthread = nthread
    )
  )
  
  rmse <-
    ModelMetrics::rmse(actual = data_val$rating,
                       predicted = r$predict(val_set))
  
  return(rmse)
}

obj_fun <- makeSingleObjectiveFunction(
  fn = scoringFunction,
  par.set = makeParamSet(
    makeIntegerParam("dim", lower = 1, upper = 100),
    makeNumericParam("lrate",  lower = 1e-04,   upper = 0.2),
    makeIntegerParam("niter", lower = 5,  upper = 1000),
    makeNumericParam("costp_l1",  lower = 0,   upper = 0.1),
    makeNumericParam("costp_l2",  lower = 0,   upper = 0.1),
    makeNumericParam("costq_l1",  lower = 0,   upper = 0.1),
    makeNumericParam("costq_l2",  lower = 0,   upper = 0.1)
  ),
  has.simple.signature = FALSE,
  minimize = TRUE
)

des <- generateDesign(
  n = 3 * getNumberOfParameters(obj_fun),
  par.set = getParamSet(obj_fun),
  fun = lhs::randomLHS
)

des$y = apply(des, 1, obj_fun)

control <- makeMBOControl() %>%
  setMBOControlTermination(., iters = 30 - 4 * getNumberOfParameters(obj_fun))

run <- mbo(
  fun = obj_fun,
  design = des,
  control = control,
  show.info = TRUE
)

r = Reco()

opts <- run$x
opts$nthread <- nthread
opts$verbose <- F
r$train(train_val_set, opts = opts)

c(P,Q) %<-% r$output(out_memory(), out_memory())


# Estimating latent variable of new users ---------------------------------

n_evaluation <- 50

i <- 2
model_ind_evaluation <- sample(1:(n_model_instances*n_model_classes), n_evaluation) %>%
  sort()
model_ind_test <- setdiff(1:(n_model_instances*n_model_classes), model_ind_evaluation)

Q_evaluation <- Q[model_ind_evaluation,]
rating_evaluation <- data_process_lookup %>%
  filter(catchment_id == catchment_id_lookup[i],
         model_id %in% model_ind_evaluation) %>%
  arrange(model_id) %>%
  pull(rating)

fn <- function(x){
  
  m1 <- matrix(x, nrow = 1)
  m2 <- t(Q_evaluation)
  
  pred <- m1%*%m2 %>%
    as.vector()
  
  -ModelMetrics::rmse(actual = rating_evaluation, predicted = pred)
}

LB <- rep(range(P,na.rm = T)[1]*1.25, dim(P)[2])
UB <- rep(range(P,na.rm = T)[2]*1.25, dim(P)[2])

GA <- ga(type = "real-valued", fitness = fn, lower = LB, upper = UB, maxiter = 5000)

P_new_catchment <- GA@solution[1,]
pred <- P_new_catchment %*%
  t(Q_evaluation) %>%
  as.vector()
ModelMetrics::rmse(actual = rating_evaluation, pred)
cor(rating_evaluation, pred)^2
plot(rating_evaluation, pred)



# evaluation
model_ind_test <- setdiff(1:(n_model_instances*n_model_classes), model_ind_evaluation) %>%
  sort()
Q_test<- Q[model_ind_test,]
rating_test <- data_process_lookup %>%
  filter(catchment_id == catchment_id_lookup[i],
         model_id %in% model_ind_test) %>%
  arrange(model_id) %>%
  pull(rating)

pred_test <-P_new_catchment %*%
  t(Q_test) %>%
  as.vector()
ModelMetrics::rmse(actual = rating_test, pred_test)
cor(rating_test, pred_test)^2
plot(rating_test, pred_test)




# LM experiments ----------------------------------------------------------

data_lm <- Q_evaluation %>%
  as.data.frame() %>%
  mutate(rating = rating_evaluation)

model <- lm(rating ~ . - 1, data = data_lm)
preds <- predict(model, data_lm)
plot(preds, data_lm$rating)
cor(data_lm$rating, preds)^2


data_lm_test <- Q[model_ind_test,] %>%
  as.data.frame() %>%
  mutate(rating = rating_test)
preds <- predict(model, data_lm_test)

plot(preds, data_lm_test$rating)

cor(data_lm_test$rating, preds)^2
plot(data_lm_test$rating, preds)


