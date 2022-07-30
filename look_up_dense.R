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
               doParallel,
               Rfast,
               Rfast2)

# seed --------------------------------------------------------------------

set.seed(1234)

nthread <- detectCores()
# cl <- makeCluster(nthread-2)
# registerDoParallel(cl)
# 
# clusterEvalQ(cl, library(Rfast))
# clusterEvalQ(cl, library(Rfast2))

# data --------------------------------------------------------------------

weights <- read_csv("./data/SM.csv",
                    col_names = c("KGE", "NSE", "RMSE"))

# 10 model classes in the dense data set
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
n_instances_per_class <-  nrow(weights)/n_catchments/n_model_classes # number of instance of each model class
n_instances <- n_model_classes *n_instances_per_class

# assign catchment id and model instance id to each row
catchment_id <- rep(1:n_catchments, each = n_instances_per_class)
data_raw <- weights %>%
  bind_cols(expand_grid(model_class, catchment_id)) %>%
  mutate(
    instance_id = rep(1:n_instances_per_class, n() / n_instances_per_class),
    model_id = paste(model_class, instance_id, sep = "_"),
    model_id = factor(model_id, levels = unique(model_id)),
    model_id = as.integer(model_id) # assign a unique id to each model
  ) %>%
  select(model_class, catchment_id, model_id, KGE, NSE, RMSE)

# normalizing NSE on a scale from 0 to 10
data_process <- data_raw %>%
  mutate(NNSE = 1 / (2 - NSE) * 10) %>%
  dplyr::select(catchment_id,
                model_id,
                NNSE) %>%
  rename(rating = NNSE) %>%
  mutate(record_id = 1:n())

# Function ----------------------------------------------------------------

#  Splitting with important groups, i.e., catchment_id
train_folds <- groupKFold(data_process$catchment_id, k = 10)

base_look_up_split <- function(train_fold){
  # This function splits data_process into a data set for building recommender system,
  # and a data set for look up experiments
  
  data_base <- data_process[train_fold,]
  data_look_up <- data_process[-train_fold,]
  
  # get the catchment_ids for look up experiments
  catchment_id_look_ups <- data_look_up %>% pull(catchment_id) %>% unique() %>% sort()
  catchment_id_bases <- setdiff(1:n_catchments, catchment_id_look_ups)
  
  list(
    data_base = data_base,
    data_look_up = data_look_up,
    catchment_id_bases = catchment_id_bases,
    catchment_id_look_ups = catchment_id_look_ups
  )
}

base_further_split <- function(data_base, train_portion = 0.8) {
  # This function further splits a subset of data_base in to train_portion and val_portion,
  # to optimize the hyperparameters of matrix factorization
  # train_portion and 1-train_portion of the subset of the data are used for different roles.
  
  # split the subset into train and validation sets
  data_train <- data_base %>%
    group_by(model_id) %>%
    slice_sample(prop = train_portion) %>%
    ungroup() %>%
    arrange(catchment_id, model_id)
  
  data_val <- data_base %>%
    filter(record_id %in% setdiff(data_base$record_id, data_train$record_id)) %>%
    arrange(catchment_id, model_id)
  
  # return the data set and the row id
  list(
    data_train = data_train %>% select(-record_id),
    data_val = data_val %>% select(-record_id),
    record_id_train = data_train$record_id,
    record_id_val = data_val$record_id
  )
}

derive_PQ <- function(data_train, data_val) {
  # This function derives the P and Q value from data_base
  
  # prepare the data sets in the format required by recosystem
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
    user_index = data_base$catchment_id,
    item_index = data_base$model_id,
    rating = data_base$rating,
    index1 = T
  )
  
  scoringFunction <- function(x) {
    # 7 hyperparameters are used
    
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
        nthread = nthread # nbin = 4*nthread^2+1
      )
    )
    
    rmse <-
      ModelMetrics::rmse(actual = data_val$rating,
                         predicted = r$predict(val_set))
    
    return(rmse)
  }
  
  # create objective function in mlrMBO required format
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
  
  c(P, Q) %<-% r$output(out_memory(), out_memory())
  
  list(P = P, Q = Q)
}

derive_PQ_wrapper <- function(data_base){
  
  # further split data_base into data_train and data_val for derive optimal parameters for matrix factorization
  c(data_train, data_val, record_id_train, record_id_val) %<-% base_further_split(data_base)
  
  # derive P and Q
  c(P, Q) %<-% derive_PQ(data_train, data_val)
  
  # output
  list(
    P = P,
    Q = Q,
    record_id_train = record_id_train, 
    record_id_val = record_id_val
  )
}

fn_factory <- function(Q, rating) {
  function(x) {
    # This function computes the predicted
    pred <- Rfast::eachrow(Q, x, "*")
    pred <- Rfast::rowsums(pred)
    
    - ModelMetrics::rmse(actual = rating, predicted = pred)
  }
}

p_rmse <- function(p, Q, rating) {
  # This function computes the RMSE of predicted ratings of models specified by Q
  # the target rating is `rating`
  ModelMetrics::rmse(predicted = p %*% t(Q) %>% as.vector(),
                     actual = rating)
}

p_r2 <- function(p, Q, rating) {
  # This function computes the R-squared of predicted ratings of models specified by Q
  # the target rating is `rating`
  cor(p %*% t(Q) %>% as.vector(), rating) ^ 2
}

p_at_k <- function(pred, actual, k = 1) {
  # compute precision at K
  sum(pred$model_id[1:k] %in% actual$model_id[1:k]) / k
}

mean_diff_at_k <- function(pred, actual, k = 1) {
  # difference in mean ranking
  mean(pred$rating[1:k]) - mean(actual$rating[1:k])
}

prepare_Qs_for_look_up  <- function(Q, n_probed, train_portion, catchment_id_look_up){
  # This function split Q into Q_probed, Q_train, Q_val, and Q_test
  # weights of the links between the look-up catchment and the models associated with Q_probed is known
  # Q_probed is further split into Q_train and Q_val for deriving the optimal number of iterations in GA
  # Q_test is the remaining models with unknown weight links with the look-up catchment.
  
  # n_evalution: the number of links between Q and the look-up catchment with known weights
  
  model_id_probed <- sample(1:n_instances, n_probed) %>%
    sort()
  
  Q_probed <- Q[model_id_probed,]
  
  rating_probed <- data_process %>%    # data_process is from global
    filter(catchment_id == catchment_id_look_up,
           model_id %in% model_id_probed) %>%
    arrange(model_id) %>%
    pull(rating)
  
  # Q_test and rating_test
  model_id_test <- setdiff(1:n_instances, model_id_probed)
  Q_test <- Q[model_id_test,]
  rating_test <- data_look_up %>%
    filter(catchment_id == catchment_id_look_up,
           model_id %in% model_id_test) %>%
    arrange(model_id) %>%
    pull(rating)
  
  # train and validation split
  n_train <- round(length(model_id_probed)*train_portion)
  
  sample_id <- sample(seq_along(model_id_probed), size = n_train) %>% sort()
  model_id_train <- model_id_probed[sample_id]
  model_id_val <- model_id_probed[-sample_id]
  Q_train <- Q[model_id_train,]
  Q_val <- Q[model_id_val,]
  rating_train <- data_look_up %>%
    filter(catchment_id == catchment_id_look_up,
           model_id %in% model_id_train) %>%
    arrange(model_id) %>%
    pull(rating)
  rating_val <- data_look_up %>%
    filter(catchment_id == catchment_id_look_up,
           model_id %in% model_id_val) %>%
    arrange(model_id) %>%
    pull(rating)
  
  # output
  list(
    Q_probed = Q_probed,
    rating_probed = rating_probed,
    Q_test = Q_test,
    rating_test = rating_test,
    model_id_test = model_id_test,
    Q_train = Q_train,
    rating_train = rating_train,
    Q_val = Q_val,
    rating_val = rating_val
  )
}

derive_p <- function(Q, n_probed, train_portion, catchment_id_look_up){
  # train_portion: the portion used in training during GA hyperparameter optimization
  
  # prepare look up data sets
  c(
    Q_probed, rating_probed,
    Q_test, rating_test, model_id_test,
    Q_train, rating_train, 
    Q_val, rating_val
  ) %<-%
    prepare_Qs_for_look_up(Q, n_probed, train_portion, catchment_id_look_up)
  
  # get model_id of the models used in the evaluation and the retrieval experiments
  # run ga with Q_train and rating_train
  fn <- fn_factory(Q_train, rating_train)
  
  LB <- apply(P[catchment_id_bases,], 2, min, na.rm = T)
  LB <- LB - abs(LB) * 0.25
  
  UB <- apply(P[catchment_id_bases,], 2, max, na.rm = T)
  UB <- UB + abs(UB) * 0.25
  
  GA <-
    ga(
      type = "real-valued",
      popSize = length(LB)*2,
      fitness = fn,
      lower = LB,
      upper = UB,
      maxiter = 5000,
      parallel = F,
      keepBest = T,
      monitor = F
    )
  
  # validate and find the optimal iteration using Q_val and rating_val
  intermediate_solution <- GA@bestSol %>%
    lapply(function(x) x[1,]) %>%
    unlist() %>%
    matrix(byrow = T, ncol = length(LB))
  
  out <-
    apply(intermediate_solution, 1, function(x)
      p_rmse(x, Q_val, rating_val))
  
  optimal_iter <- which.min(out)
  
  # run ga with Q_probed and rating_probed
  fn <- fn_factory(Q_probed, rating_probed)
  
  GA <-
    ga(
      type = "real-valued",
      popSize = length(LB)*2,
      fitness = fn,
      lower = LB,
      upper = UB,
      maxiter = optimal_iter,
      parallel = F,
      keepBest = F,
      monitor = F
    )
  
  # estimated p
  p <- GA@solution[1,] %>% unname()
  
  # evaluate the quality of p
  pred <- p %*%
    t(Q_test) %>%
    as.vector()
  
  # identify top 100 models specified by Q_test and p
  top100_pred <- tibble(
    model_id = model_id_test,
    rating = pred,
    case = "pred"
  ) %>%
    arrange(desc(rating)) %>%
    slice(1:100)
  
  top100_actual <- tibble(
    model_id = model_id_test,
    rating = rating_test,
    case = "actual"
  ) %>%
    arrange(desc(rating)) %>%
    slice(1:100)
  
  # output
  list(
    p = p,
    rmse = p_rmse(p, Q_test, rating_test),
    r2 = p_r2(p, Q_test, rating_test),
    top100_actual = top100_actual,
    top100_pred = top100_pred
  )
}

# Derive P and Q ----------------------------------------------------------

# iterate over the train_folds, i.e., repeat the look-up experiment on different splits of the data
eval_grids <- vector("list", length(train_folds))
Ps  <- vector("list", length(train_folds))
Qs  <- vector("list", length(train_folds))

for (i in seq_along(train_folds)){
  
  # split data set into data_base for estimating P and Q, and data_look_up for look up experiments
  c(data_base, data_look_up, catchment_id_bases, catchment_id_look_ups) %<-% 
    base_look_up_split(train_fold = train_folds[[i]])
  
  # derive P and Q
  c(P, Q, record_id_train, record_id_val) %<-% derive_PQ_wrapper(data_base)

  # n_probed: the numbers of edges evaluated for estimating p
  
  eval_grid <- expand_grid(
    n_probed = round(dim(P)[2]*0.5*1:4), # 0.5, 1, 1.5, and 2 times of the latent dimension
    catchment_id_look_up = catchment_id_look_ups,
    out = vector("list", 1)
  )
  
  # iterate over eval_grid
  for (j in 1:nrow(eval_grid)){
    
    n_probed <- eval_grid$n_probed[[j]]
    catchment_id_look_up <- eval_grid$catchment_id_look_up[[j]]
    train_portion <- 0.8
    
    # evaluation
    eval_grid$out[[j]] <- derive_p(Q, n_probed, train_portion, catchment_id_look_up)
  }
  
  eval_grid <- eval_grid %>%
    mutate(train_fold_id = i)
  
  # output
  eval_grids[[i]] <- eval_grid
  Ps[[i]] <- P
  Qs[[i]] <- Q
}

eval_grid <- eval_grids %>%
  bind_rows()

save(Ps, Qs, train_folds, eval_grid, file = "./data/dense_look_up.Rda")


# Recycle -----------------------------------------------------------------



# LM experiments ----------------------------------------------------------

data_lm <- Q_probed %>%
  as.data.frame() %>%
  mutate(rating = rating_probed)

model <- lm(rating ~ . - 1, data = data_lm)
preds <- predict(model, data_lm)
plot(preds, data_lm$rating)
cor(data_lm$rating, preds)^2


data_lm_test <- Q[model_id_test,] %>%
  as.data.frame() %>%
  mutate(rating = rating_test)
preds <- predict(model, data_lm_test)

plot(preds, data_lm_test$rating)

cor(data_lm_test$rating, preds)^2
plot(data_lm_test$rating, preds)


