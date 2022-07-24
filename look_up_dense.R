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

cl <- makeCluster(nthread)
registerDoParallel(cl)

clusterEvalQ(cl, library(Rfast))
clusterEvalQ(cl, library(Rfast2))

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
  list(
    data_base = data_process[train_fold,],
    data_look_up = data_process[-train_fold,]
  )
}

prepare_modeling_data <-
  function(data_base,
           train_portion = 0.8,
           val_portion = 0.2) {
    # This function splits a subset of data_base in to train_portion and val_portion.
    # train_portion and val_portion of the subset of the data are used for different roles.
    
    # split the subset into train and validation sets
    data_train <- data_base %>%
      group_by(model_id) %>%
      slice_sample(prop = train_portion) %>%
      ungroup() %>%
      arrange(model_id)
    
    data_val <- data_base %>%
      filter(record_id %in% setdiff(data_base$record_id, data_train$record_id)) %>%
      arrange(model_id)
    
    # return the data set and the row id
    list(
      data_train = data_train %>% select(-record_id),
      data_val = data_val %>% select(-record_id),
      record_id_train = data_train$record_id,
      record_id_val = data_val$record_id
    )
  }


derive_PQ <- function(data_train, data_val){
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
  
  c(P,Q) %<-% r$output(out_memory(), out_memory())
  
  list(P = P, Q = Q)
}



# Derive P and Q ----------------------------------------------------------

i <- 1

c(data_base, data_look_up) %<-% base_look_up_split(train_fold = train_folds[[i]])
c(data_train, data_val, record_id_train = data_train$record_id, record_id_val) %<-% prepare_modeling_data(data_base)
# derive P and Q
c(P, Q) %<-% derive_PQ(data_train, data_val)

# get the catchment_ids for look up experiments
catchment_id_look_ups <- data_look_up %>% pull(catchment_id) %>% unique() %>% sort()
catchment_id_bases <- setdiff(1:n_catchments, catchment_id_look_ups)

save(P,Q, train_folds, data_base, data_look_up, file = "./data/temp_look_up.Rda")


# Look up -----------------------------------------------------------------


load("./data/temp_look_up.Rda")

# set the number of edges evaluated
n_evaluation <- dim(P)[2]
# set the portion used in training due GA hyperparameter optimization
train_portion <- 0.8

# Functions
fn_factory <- function(Q, rating){
  function(x){
    # This function computes the predicted 
    pred <- Rfast::eachrow(Q, x, "*")
    pred <- Rfast::rowsums(pred)
    
    -ModelMetrics::rmse(actual = rating, predicted = pred)
  }
}

p_rmse <- function(p, Q, rating){
  # This function computes the RMSE of predicted ratings of models specified by Q
  # the target rating is `rating`
  ModelMetrics::rmse(
    predicted = p %*% t(Q) %>% as.vector(),
    actual = rating
  )
}

p_r2 <- function(p, Q, rating){
  # This function computes the R-squared of predicted ratings of models specified by Q
  # the target rating is `rating`
  cor(p %*% t(Q) %>% as.vector(), rating)^2
}

split_Q <- function(Q, n_evaluation, train_portion, catchment_id_look_up){
  # This function split Q into Q_evaluation, Q_train, Q_val, and Q_test
  # weights of the links between the look-up catchment and the models associated with Q_evaluation is known
  # Q_evaluation is further split into Q_train and Q_val for deriving the optimal number of iterations in GA
  # Q_test is the remaining models with unknown weight links with the look-up catchment.
  
  # n_evalution: the number of links between Q and the look-up catchment with known weights
  # data_look_up is from global env
  
  model_ind_evaluation <- sample(1:n_instances, n_evaluation) %>%
    sort()
  
  Q_evaluation <- Q[model_ind_evaluation,]
  
  rating_evalution <- data_look_up %>%
    filter(catchment_id == catchment_id_look_up,
           model_id %in% model_ind_evaluation) %>%
    arrange(model_id) %>%
    pull(rating)
  
  # Q_test and rating_test
  model_ind_test <- setdiff(1:n_instances, model_ind_evaluation)
  Q_test <- Q[model_ind_test,]
  rating_test <- data_look_up %>%
    filter(catchment_id == catchment_id_look_up,
           model_id %in% model_ind_test) %>%
    arrange(model_id) %>%
    pull(rating)
  
  # train and validation split
  n_train <- round(length(model_ind_evaluation)*train_portion)
  
  sample_ind <- sample(seq_along(model_ind_evaluation), size = n_train) %>% sort()
  model_ind_train <- model_ind_evaluation[sample_ind]
  model_ind_val <- model_ind_evaluation[-sample_ind]
  Q_train <- Q_evaluation[sample_ind,]
  Q_val <- Q_evaluation[-sample_ind,]
  rating_train <- data_look_up %>%
    filter(catchment_id == catchment_id_look_up,
           model_id %in% model_ind_train) %>%
    arrange(model_id) %>%
    pull(rating)
  rating_val <- data_look_up %>%
    filter(catchment_id == catchment_id_look_up,
           model_id %in% model_ind_val) %>%
    arrange(model_id) %>%
    pull(rating)
  

  # output
  list(
    Q_evaluation = Q_evaluation,
    rating_evalution = rating_evalution,
    Q_test = Q_test,
    rating_test = rating_test,
    Q_train = Q_train,
    rating_train = rating_train,
    Q_val = Q_val,
    rating_val = rating_val
  )
}


# iterate over catchment_id_look_ups
# catchment id used in current experiments
j <- 1
catchment_id_look_up <- catchment_id_look_ups[[j]]

c(
  Q_evaluation, rating_evalution,
  Q_test, rating_test,
  Q_train, rating_train, 
  Q_val, rating_val
) %<-%
  split_Q(Q, n_evaluation, train_portion, catchment_id_look_up)


# get model_ind of the models used in the evaluation and the retrieval experiments
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

# run ga with Q_evaluation and rating_evaluation
fn <- fn_factory(Q_evaluation, rating_evaluation)

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

# testing
p_look_up <- GA@solution[1,]
pred <- p_look_up %*%
  t(Q_test) %>%
  as.vector()

p_rmse(p_look_up, Q_test, rating_test)

cor(rating_test, pred)^2
plot(rating_test, pred)




# Plotting ----------------------------------------------------------------

n_evaluation <- 200
max_rating <- rep(0,33)
max_rating_lookup <- rep(0,33)
r2s <- rep(0,33)
rmses <- rep(0,33)
rankings <- rep(0,33)

for (i in 1:33){
  model_ind_evaluation <- sample(1:(n_instances_per_class*n_model_classes), n_evaluation) %>%
    sort()
  model_ind_test <- setdiff(1:(n_instances_per_class*n_model_classes), model_ind_evaluation)
  
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
  
  GA <- ga(type = "real-valued", fitness = fn, lower = LB, upper = UB, maxiter = 200)
  
  P_new_catchment <- GA@solution[1,]
  
  # evaluation
  model_ind_test <- setdiff(1:(n_instances_per_class*n_model_classes), model_ind_evaluation) %>%
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
  rmses[[i]] <- ModelMetrics::rmse(actual = rating_test, pred_test)
  r2s[[i]] <- cor(rating_test, pred_test)^2
  #plot(rating_test, pred_test)
  
  max_rating[[i]] <- max(rating_test)
  max_rating_lookup[[i]] <- rating_test[which.max(pred_test)]
  rankings[[i]] <- rank(rating_test)[which.max(pred_test)]/length(rating_test)
}

plot(max_rating_lookup, max_rating)
hist(max_rating_lookup - max_rating)
sum(max_rating_lookup - max_rating)

data_plot <- tibble(
  max_rating,
  max_rating_lookup,
  rankings
)

ggplot(data_plot, aes(max_rating_lookup, max_rating)) +
  geom_point()+
  theme_bw()+
  labs(x = "Rating of the model retrived using catchment embeddings",
       y = "Rating of the best model in the data set")

ggplot(data_plot, aes(rankings * 100))+
  geom_histogram()+
  theme_bw()




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


