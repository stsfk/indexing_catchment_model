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
               ModelMetrics)

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
  n = 5 * getNumberOfParameters(obj_fun),
  par.set = getParamSet(obj_fun),
  fun = lhs::randomLHS
)

des$y = apply(des, 1, obj_fun)

control <- makeMBOControl() %>%
  setMBOControlTermination(., iters = 40 - 5 * getNumberOfParameters(obj_fun))

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



i <- 1
j <- 4023

retrieve_from_model <- function(i, j){
  P[i,] %*% Q[j,] %>% as.vector()
}
n_evaluated <- 1000
js <- sample(1:(n_model_instances*n_model_classes), n_evaluated)

build_lm_data_set <- function(i, js){
  
  ratings <- data_process %>%
    filter(catchment_id == i,
           model_id %in% js) %>%
    pull(rating)
  
  Q[js,] %>%
    as.data.frame() %>%
    mutate(rating = ratings) %>%
    as_tibble()
}

data_lm <- build_lm_data_set(i, js)

model <- lm(rating ~ . - 1, data = data_lm)
preds <- predict(model, data_lm)

plot(preds, data_lm$rating)



js <- sample(1:(n_model_instances*n_model_classes), 100) %>% sort()
ratings <- data_process %>%
  filter(catchment_id == i,
         model_id %in% js) %>%
  pull(rating)

js <- sort(js)

i <- catchment_id_lookup[[2]]
ratings <- data_process %>%
  filter(catchment_id == i,
         model_id %in% js) %>%
  pull(rating)

fn <- function(x){
  
  m1 <- matrix(x, nrow = 1)
  m2 <- t(Q[js,])
  
  pred <- m1%*%m2 %>%
    as.vector()
  
  -ModelMetrics::rmse(actual = ratings, predicted = pred)
}

GA <- ga(type = "real-valued", fitness = fn, lower = rep(-1, 43), upper = rep(1, 43), maxiter = 2000)


pred <- GA@solution %*%
  t(Q[js,]) %>%
  as.vector()
actual = ratings

ModelMetrics::rmse(actual, pred)
plot(actual, pred)


js <- sample(1:(n_model_instances*n_model_classes), 1000) %>% sort()
ratings <- data_process %>%
  filter(catchment_id == i,
         model_id %in% js) %>%
  pull(rating)
pred <- GA@solution %*%
  t(Q[js,]) %>%
  as.vector()
actual = ratings

ModelMetrics::rmse(actual, pred)
cor(actual, pred)^2

plot(actual, pred)











  
  
data_lm_all <- build_lm_data_set(i, 1:100)
preds <- predict(model, data_lm_all)

plot(preds, data_lm_all$rating)














sparse_gof_wrapper <- function(frac = 0.1) {
  c(data_train,
    data_val,
    data_test,
    record_id_train,
    record_id_val,
    record_id_test) %<-% prepare_modeling_data(frac = frac)
  
  # training
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
  
  data_train_val <- data_train %>%
    bind_rows(data_val)
  train_val_set <- data_memory(
    user_index = data_train_val$catchment_id,
    item_index = data_train_val$model_id,
    rating = data_train_val$rating,
    index1 = T
  )
  
  test_set <- data_memory(
    user_index = data_test$catchment_id,
    item_index = data_test$model_id,
    rating = data_test$rating,
    index1 = T
  )
  
  # experiment
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
    n = 5 * getNumberOfParameters(obj_fun),
    par.set = getParamSet(obj_fun),
    fun = lhs::randomLHS
  )
  
  des$y = apply(des, 1, obj_fun)
  
  control <- makeMBOControl() %>%
    setMBOControlTermination(., iters = 100 - 5 * getNumberOfParameters(obj_fun))
  
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
  pred_rvec <- r$predict(test_set)
  r2 <- cor(data_test$rating, pred_rvec) ^ 2
  rmse <- ModelMetrics::rmse(actual = data_test$rating,
                             predicted = pred_rvec)
  
  c(P,Q) %<-% r$output(out_memory(), out_memory())
  
  out <- list(
    r2 = r2,
    rmse = rmse,
    run = run,
    des = des,
    P = P,
    Q = Q,
    record_id_train = record_id_train,
    record_id_val = record_id_val,
    record_id_test = record_id_test
  )
}

eval_grid <- expand_grid(
  ratio = c(0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.2, 0.5, 0.8, 1),
  r2 = 0,
  rmse = 0,
  out = vector("list",1),
  repeats = c(1)
)

sta_time <- Sys.time()

for (i in 1:nrow(eval_grid)){
  frac <- eval_grid$ratio[i]
  eval_grid$out[[i]] <- sparse_gof_wrapper(frac)
  eval_grid$r2[[i]] <- eval_grid$out[[i]]$r2
  eval_grid$rmse[[i]] <- eval_grid$out[[i]]$rmse
  
  gc()
}

end_time <- Sys.time()

end_time - sta_time

save(eval_grid, file = "sparse_exp.Rda")

# recycle -----------------------------------------------------------------
