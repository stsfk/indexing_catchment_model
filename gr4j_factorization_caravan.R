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
               parallel)

# seed --------------------------------------------------------------------

set.seed(1234)
nthread <- parallel::detectCores() - 1 # number of CPU thread

# data --------------------------------------------------------------------

load("./data/gr4j_carvan.Rda")

catchment_ids <- results$catchment_id %>% unique()

# normalizing NSE on a scale from 0 to 10
data_process <- results %>%
  mutate(catchment_id = factor(catchment_id, levels = catchment_ids) %>% as.numeric()) 

# Functions ---------------------------------------------------------------


prepare_modeling_data <-
  function(frac = 1,
           train_portion = 0.6,
           val_portion = 0.2,
           test_portion = 0.2) {
    # This function splits a subset of the data set in to train_portion, val_portion, and test_portion.
    # frac: the fraction of samples in the data set used in modeling
    # train_portion, val_portion, and test_portion of the subset of the data are used for different roles.
    
    # Creating a subset of "data_process" for modeling
    if (frac == 1){
      data_sample <- data_process
    } else {
      data_sample <- data_process %>%
        group_by(model_id) %>%
        sample_frac(frac) %>%
        ungroup() %>%
        arrange(model_id)
    }
    
    # split the subset into train, validation, and test sets
    data_train_val <- data_sample %>%
      group_by(model_id) %>%
      sample_frac(train_portion + val_portion) %>%
      ungroup()
    
    data_train <- data_train_val %>%
      group_by(model_id) %>%
      sample_frac(train_portion / (val_portion + train_portion)) %>%
      ungroup()
    
    data_val <- data_train_val %>%
      filter(record_id %in% setdiff(data_train_val$record_id, data_train$record_id))
    
    data_test <- data_sample %>%
      filter(record_id %in% setdiff(data_sample$record_id, data_train_val$record_id))
    
    # return the data set and the row id
    list(
      data_train = data_train %>% select(-record_id),
      data_val = data_val %>% select(-record_id),
      data_test = data_test %>% select(-record_id),
      record_id_train = data_train$record_id,
      record_id_val = data_val$record_id,
      record_id_test = data_test$record_id
    )
  }

factorization_wrapper <- function(i, frac = 1) {
  # splitting
  c(data_train,
    data_val,
    data_test,
    record_id_train,
    record_id_val,
    record_id_test) %<-% prepare_modeling_data(frac = frac)
  
  # training set
  train_set <- data_memory(
    user_index = data_train$catchment_id,
    item_index = data_train$model_id,
    rating = data_train$rating,
    index1 = T
  )
  
  # validation set
  val_set <- data_memory(
    user_index = data_val$catchment_id,
    item_index = data_val$model_id,
    rating = data_val$rating,
    index1 = T
  )
  
  # training and validation set combined
  data_train_val <- data_train %>%
    bind_rows(data_val)
  train_val_set <- data_memory(
    user_index = data_train_val$catchment_id,
    item_index = data_train_val$model_id,
    rating = data_train_val$rating,
    index1 = T
  )
  
  # test set
  test_set <- data_memory(
    user_index = data_test$catchment_id,
    item_index = data_test$model_id,
    rating = data_test$rating,
    index1 = T
  )
  
  # experiment
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
        nthread = nthread,
        nbin = 2*nthread# nbin = 4*nthread^2+1
      )
    )
    
    # rmse is used as the optimization objective
    rmse <- ModelMetrics::rmse(actual = data_val$rating,
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
  
  # statistical design 
  des <- generateDesign(
    n = 5 * getNumberOfParameters(obj_fun),
    par.set = getParamSet(obj_fun),
    fun = lhs::randomLHS
  )
  
  des$y = apply(des, 1, obj_fun)
  
  # set number of generation
  control <- makeMBOControl() %>%
    setMBOControlTermination(., iters = 100 - 5 * getNumberOfParameters(obj_fun))
  
  # run Bayesian optimization 
  run <- mbo(
    fun = obj_fun,
    design = des,
    control = control,
    show.info = TRUE
  )
  
  # get the optimal hyperparameters and run on the data set combining training and validation
  r = Reco()
  
  opts <- run$x
  opts$nthread <- nthread
  opts$verbose <- F
  r$train(train_val_set, opts = opts)
  
  # evaluating on test set
  pred_rvec <- r$predict(test_set)
  r2 <- cor(data_test$rating, pred_rvec) ^ 2
  rmse <- ModelMetrics::rmse(actual = data_test$rating,
                             predicted = pred_rvec)
  
  # outputting
  # c(P,Q) %<-% r$output(out_memory(), out_memory())
  P_file = out_file(paste0("./data/mat_P_",i,".txt"))
  Q_file = out_file(paste0("./data/mat_Q_",i,".txt"))
  
  r$output(out_P = P_file, out_Q = Q_file)
  
  out <- list(
    r2 = r2,
    rmse = rmse,
    run = run,
    des = des,
    record_id_train = record_id_train,
    record_id_val = record_id_val,
    record_id_test = record_id_test
  )
}

# Run ---------------------------------------------------------------------

eval_grid <- expand_grid(
  r2 = 0,
  rmse = 0,
  out = vector("list",1),
  repeats = c(1:1)
)

sta_time <- Sys.time()

for (i in 1:nrow(eval_grid)){
  
  # Matrix factorization can stuck in local minimal, in which case r2 value is NA
  # check if r2 is NA; if so, run matrix factorization again
  okay_result <- FALSE
  while (!okay_result){
    temp <- factorization_wrapper(i, frac = 1)
    
    if (!is.na(temp$r2)){
      okay_result <- TRUE
    }
    
    gc()
  }
  
  eval_grid$out[[i]] <- temp
  eval_grid$r2[[i]] <- temp$r2
  eval_grid$rmse[[i]] <- temp$rmse
  
}

end_time <- Sys.time()

end_time - sta_time

# save(eval_grid, file = "./data/dens_factorization.Rda")


# Post processing ---------------------------------------------------------
# load("./data/dens_factorization.Rda")

eval_grid <- eval_grid %>%
  filter(!is.na(r2))

eval_grid$r2 %>% mean(na.rm = T)
eval_grid$r2 %>% var(na.rm = T)

eval_grid$rmse %>% mean(na.rm = T)
eval_grid$rmse %>% var(na.rm = T)

# plot --------------------------------------------------------------------



# Speed test --------------------------------------------------------------


c(data_train,
  data_val,
  data_test,
  record_id_train,
  record_id_val,
  record_id_test) %<-% prepare_modeling_data(frac = frac)

# training set
train_set <- data_memory(
  user_index = data_train$catchment_id,
  item_index = data_train$model_id,
  rating = data_train$rating,
  index1 = T
)

# validation set
val_set <- data_memory(
  user_index = data_val$catchment_id,
  item_index = data_val$model_id,
  rating = data_val$rating,
  index1 = T
)

# training and validation set combined
data_train_val <- data_train %>%
  bind_rows(data_val)
train_val_set <- data_memory(
  user_index = data_train_val$catchment_id,
  item_index = data_train_val$model_id,
  rating = data_train_val$rating,
  index1 = T
)

# test set
test_set <- data_memory(
  user_index = data_test$catchment_id,
  item_index = data_test$model_id,
  rating = data_test$rating,
  index1 = T
)



sta_time <- Sys.time()
r = Reco()
r$train(
  train_set,
  opts = list(
    dim = 100,
    costp_l1 = 0,
    costp_l2 = 0,
    costq_l1 = 0,
    costq_l2 = 0,
    lrate = 0.005,
    niter = 10,
    verbose = T,
    nthread = nthread, 
    nbin = 20 # nbin = 4*nthread^2+1
  )
)
end_time <- Sys.time()
dif1 <- end_time - sta_time




sta_time <- Sys.time()
r = Reco()
r$train(
  train_set,
  opts = list(
    dim = 100,
    costp_l1 = 0,
    costp_l2 = 0,
    costq_l1 = 0,
    costq_l2 = 0,
    lrate = 0.01,
    niter = 150,
    verbose = T,
    nthread = nthread, 
    nbin = nthread * 2 # nbin = 4*nthread^2+1
  )
)
end_time <- Sys.time()
dif2 <- end_time - sta_time

# evaluating on test set
pred_rvec <- r$predict(test_set)
r2 <- cor(data_test$rating, pred_rvec) ^ 2
rmse <- ModelMetrics::rmse(actual = data_test$rating,
                           predicted = pred_rvec)

# outputting
r$output(out_P = out_file("mat_P.txt"), out_Q = out_file("mat_Q.txt"))

read.table("mat_P.txt") %>%
  as.matrix()













sta_time <- Sys.time()
r = Reco()
r$train(
  train_set,
  opts = list(
    dim = 100,
    costp_l1 = 0,
    costp_l2 = 0,
    costq_l1 = 0,
    costq_l2 = 0,
    lrate = 0.005,
    niter = 10,
    verbose = T,
    nthread = nthread, 
    nbin = 4*nthread^2+1
  )
)
end_time <- Sys.time()
dif3 <- end_time - sta_time
