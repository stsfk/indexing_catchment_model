if (!require("pacman")) {
  install.packages("pacman")
}

pacman::p_load(tidyverse,
               lubridate,
               zeallot,
               recosystem,
               hydroGOF,
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

# splitting
c(data_train,
  data_val,
  data_test,
  record_id_train,
  record_id_val,
  record_id_test) %<-% prepare_modeling_data(frac = 1)

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
    lrate = 0.01,
    niter = 10,
    verbose = T,
    nthread = nthread, 
    nbin = nthread * 2 # nbin = 4*nthread^2+1
  )
)
end_time <- Sys.time()
dif <- end_time - sta_time






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
    niter = 10,
    verbose = T,
    nthread = nthread, 
    nbin = 4*nthread^2+1 # nbin = 4*nthread^2+1
  )
)
end_time <- Sys.time()
dif2 <- end_time - sta_time