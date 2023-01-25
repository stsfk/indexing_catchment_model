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
               rgenoud,
               GA,
               ModelMetrics,
               doParallel,
               Rfast,
               Rfast2,
               airGR)

# seed --------------------------------------------------------------------

set.seed(1234)

nthread <- detectCores()
cl <- makeCluster(nthread-2)
registerDoParallel(cl)

clusterEvalQ(cl, library(lubridate))
clusterEvalQ(cl, library(airGR))
clusterEvalQ(cl, library(tidyverse))

# data --------------------------------------------------------------------

# load gr4j parameters
load("./data/gr4j_carvan.Rda")
rm(results)

n_instances <- nrow(paras)

# load PQ data
P <- read.table("./data/mat_P_1.txt", header = FALSE, sep = " ")
Q <- read.table("./data/mat_Q_1.txt", header = FALSE, sep = " ")

# load CAMELS data
data_raw <- read_csv("./data/CAMELS_US.csv")

data_process <- data_raw %>%
  transmute(
    catchment_id = catchment_id,
    DatesR = as.POSIXct(date, format = "%Y-%m-%d", tz = "UTC"),
    P = P,
    `T` = `T`,
    E = PET,
    Qmm = Q
  )

catchments <- data_process$catchment_id %>% unique()


# Function ----------------------------------------------------------------

get_catchment_gof <- function(selected_catchment_id, selected_para_ids, selected_paras){
  
  # subsetting by catchment
  BasinObs <- data_process %>%
    dplyr::filter(catchment_id == selected_catchment_id)
  
  InputsModel <- CreateInputsModel(FUN_MOD = RunModel_GR4J, DatesR = BasinObs$DatesR,
                                   Precip = BasinObs$P, PotEvap = BasinObs$E)
  
  Ind_Run <- seq(which(format(BasinObs$DatesR, format = "%Y-%m-%d") == "1988-10-01"), 
                 which(format(BasinObs$DatesR, format = "%Y-%m-%d") == "1998-09-30"))
  
  RunOptions <- CreateRunOptions(FUN_MOD = RunModel_GR4J,
                                 InputsModel = InputsModel, IndPeriod_Run = Ind_Run,
                                 IniStates = NULL, IniResLevels = NULL, IndPeriod_WarmUp = NULL)
  
  InputsCrit <- CreateInputsCrit(FUN_CRIT = ErrorCrit_NSE, InputsModel = InputsModel, 
                                 RunOptions = RunOptions, VarObs = "Q", Obs = BasinObs$Qmm[Ind_Run])
  
  
  OutputsCrit_wrapper <- function(i){
    Param = selected_paras[i,]
    
    OutputsModel <- RunModel_GR4J(InputsModel = InputsModel, RunOptions = RunOptions, Param = Param)
    OutputsCrit <- ErrorCrit(InputsCrit = InputsCrit, OutputsModel = OutputsModel, verbose = FALSE)
    OutputsCrit$CritValue
  }
  
  out <- foreach(x=1:nrow(selected_paras)) %dopar%
    OutputsCrit_wrapper(x) %>%
    unlist()
  
  tibble(
    para_id = selected_para_ids,
    nse = out
  )
}

prepare_Qs_for_look_up  <- function(n_probed, train_portion, selected_catchment_id){
  # This function split Q into Q_probed, Q_train, Q_val, and Q_test
  # weights of the links between the look-up catchment and the models associated with Q_probed is known
  # Q_probed is further split into Q_train and Q_val for deriving the optimal number of iterations in GA
  # Q_test is the remaining models with unknown weight links with the look-up catchment.
  
  # n_probed: the number of links between Q and the look-up catchment with known weights
  
  selected_para_ids <- sample(1:n_instances, n_probed) %>% sort()
  selected_paras <- paras[selected_para_ids,]
  
  rating_probed <-
    get_catchment_gof(selected_catchment_id, selected_para_ids, selected_paras) %>% 
    mutate(rating = 1 / (2 - nse) * 10) %>% pull(rating)
  
  Q_probed <- Q[selected_para_ids,]
  
  
  # Q_test and rating_test
  model_id_test <- setdiff(1:n_instances, selected_para_ids)
  Q_test <- Q[model_id_test,]
  
  # train and validation split
  n_train <- round(length(selected_para_ids)*train_portion)
  
  sample_id <- sample(seq_along(selected_para_ids), size = n_train) %>% sort()
  model_id_train <- selected_para_ids[sample_id]
  model_id_val <- selected_para_ids[-sample_id]
  Q_train <- Q[model_id_train,]
  Q_val <- Q[model_id_val,]
  rating_train <- rating_probed[sample_id]
  rating_val <- rating_probed[-sample_id]
  
  # output
  list(
    Q_probed = Q_probed,
    rating_probed = rating_probed,
    Q_test = Q_test,
    model_id_test = model_id_test,
    Q_train = Q_train,
    rating_train = rating_train,
    Q_val = Q_val,
    rating_val = rating_val
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

derive_p <- function(Q, n_probed, train_portion, catchment_id_look_up){
  # train_portion: the portion used in training during GA hyperparameter optimization
  
  # prepare look up data sets
  c(
    Q_probed, rating_probed,
    Q_test, model_id_test,
    Q_train, rating_train, 
    Q_val, rating_val
  ) %<-%
    prepare_Qs_for_look_up(n_probed, train_portion, selected_catchment_id)
  
  # get model_id of the models used in the evaluation and the retrieval experiments
  # run ga with Q_train and rating_train
  fn <- fn_factory(Q_train, rating_train)
  
  LB <- apply(P, 2, min, na.rm = T)
  LB <- LB - abs(LB) * 0.25
  
  UB <- apply(P, 2, max, na.rm = T)
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
    probed_ratio = 0.5*1:4, # 0.5, 1, 1.5, and 2 times of the latent dimension
    catchment_id_look_up = catchment_id_look_ups,
    out = vector("list", 1)
  ) %>%
    mutate(n_probed = round(dim(P)[2]*probed_ratio))
  
  out_wrapper <- function(x){
    n_probed <- eval_grid$n_probed[[x]]
    catchment_id_look_up <- eval_grid$catchment_id_look_up[[x]]
    train_portion <- 0.8
    
    # output
    derive_p(Q, n_probed, train_portion, catchment_id_look_up)
  }
  
  # iterate over eval_grid
  outs <- foreach(x=1:nrow(eval_grid)) %dopar% 
    out_wrapper(x)
  
  # save result to eval_grid
  eval_grid <- eval_grid %>%
    mutate(out = outs) %>%
    mutate(train_fold_id = i) # register fold id
  
  # output
  eval_grids[[i]] <- eval_grid
  Ps[[i]] <- P
  Qs[[i]] <- Q
}

eval_grid <- eval_grids %>%
  bind_rows()

save(Ps, Qs, train_folds, eval_grid, file = "./data/dense_look_up.Rda")

# stop
stopCluster(cl)



# Analysis ----------------------------------------------------------------

load("./data/dense_look_up.Rda")

catchment_lookup_groups <- lapply(train_folds,
                                  function(x)
                                    data_process[-x, ] %>% pull(catchment_id) %>% unique())

eval_summary <- eval_grid %>%
  mutate(r2 = map_dbl(out, function(x) x$r2),
         rmse = map_dbl(out, function(x) x$rmse)) %>%
  group_by(probed_ratio) %>%
  summarise(r2_mean = mean(r2),
            r2_sd = sd(r2),
            rmse_mean = mean(rmse),
            rmse_sd = sd(rmse))


eval_summary %>% transmute(
  r2 = paste0(
    round(r2_mean, 3),
    "(",
    formatC(r2_sd, format = "e", digits = 2),
    ")"
  ),
  rmse = paste0(
    round(rmse_mean, 3),
    "(",
    formatC(rmse_sd, format = "e", digits = 2),
    ")"
  )
)


# compare top retrieved model with actual top models

get_predicted_model_rating <- function(catchment_id_look_up, out){
  predicted_model_rating <- data_process %>%
    filter(catchment_id == catchment_id_look_up,
           model_id %in% out$top100_pred$model_id) %>%
    select(model_id, actual_rating=rating)
  
  out$top100_pred %>% 
    left_join(predicted_model_rating, by = "model_id") %>% # to keep the predicted model ranking
    pull(actual_rating)
}

eval_grid_expand <- eval_grid %>%
  mutate(
    predicted_model_rating = map2(catchment_id_look_up, out, get_predicted_model_rating)
  ) 

eval_grid_expand %>%
  mutate(catchment_group = map_dbl(catchment_id_look_up, function(y) which(
    sapply(catchment_lookup_groups, function(x)
      y %in% x)
  )))

eval_grid_expand %>%
  mutate(
    actual_best_model_rating = map_dbl(out, function(x) x$top100_actual$rating %>% max),
    hit1 = map_dbl(out, function(x) x$top100_actual$model_id[[1]] %in% x$top100_pred$model_id[1]),
    hit10 = map_dbl(out, function(x) x$top100_actual$model_id[[1]] %in% x$top100_pred$model_id[1:10]),
    hit25 = map_dbl(out, function(x) x$top100_actual$model_id[[1]] %in% x$top100_pred$model_id[1:25]),
    hit100 = map_dbl(out, function(x) x$top100_actual$model_id[[1]] %in% x$top100_pred$model_id),
    diff1 = map2_dbl(predicted_model_rating, actual_best_model_rating, function(x,y) y-x[1]),
    diff10 = map2_dbl(predicted_model_rating, actual_best_model_rating, function(x,y) y-max(x[1:10])),
    diff25 = map2_dbl(predicted_model_rating, actual_best_model_rating, function(x,y) y-max(x[1:25])),
    diff100 = map2_dbl(predicted_model_rating, actual_best_model_rating, function(x,y) y-max(x[1:100])),
    diff1p = diff1/actual_best_model_rating*100,
    diff10p = diff10/actual_best_model_rating*100,
    diff25p = diff25/actual_best_model_rating*100,
    diff100p = diff100/actual_best_model_rating*100
  ) %>%
  group_by(train_fold_id,probed_ratio) %>%
  summarise(
    hit1 = mean(hit1),
    hit10 = mean(hit10),
    hit100 = mean(hit100),
    diff1 = mean(diff1),
    diff10 = mean(diff10),
    diff100 = mean(diff100),
    diff1p = mean(diff1p),
    diff10p = mean(diff10p),
    diff100p = mean(diff100p),
  ) %>%
  group_by(probed_ratio) %>%
  summarise(
    hit1_sd = sd(hit1),
    hit10_sd = sd(hit10),
    hit100_sd = sd(hit100),
    hit1 = mean(hit1),
    hit10 = mean(hit10),
    hit100 = mean(hit100),
    diff1_sd = sd(diff1),
    diff10_sd = sd(diff10),
    diff100_sd = sd(diff100),
    diff1 = mean(diff1),
    diff10 = mean(diff10),
    diff100 = mean(diff100),
    diff1p_sd = sd(diff1p),
    diff10p_sd = sd(diff10p),
    diff100p_sd = sd(diff100p),
    diff1p = mean(diff1p),
    diff10p = mean(diff10p),
    diff100p = mean(diff100p),
  ) %>% 
  transmute(
    probed_ratio = probed_ratio,
    hit1 = paste0(
      formatC(hit1, format = "f", digits = 3),
      "(",
      formatC(hit1_sd, format = "f", digits = 3),
      ")"
    ),
    hit10 = paste0(
      formatC(hit10, format = "f", digits = 3),
      "(",
      formatC(hit10_sd, format = "f", digits = 3),
      ")"
    ),
    hit100 = paste0(
      formatC(hit100, format = "f", digits = 3),
      "(",
      formatC(hit100_sd, format = "f", digits = 3),
      ")"
    ),
    diff1 = paste0(
      formatC(diff1, format = "f", digits = 3),
      "(",
      formatC(diff1_sd, format = "f", digits = 3),
      ")"
    ),
    diff10 = paste0(
      formatC(diff10, format = "f", digits = 3),
      "(",
      formatC(diff10_sd, format = "f", digits = 3),
      ")"
    ),
    diff100 = paste0(
      formatC(diff100, format = "f", digits = 3),
      "(",
      formatC(diff100_sd, format = "f", digits = 3),
      ")"
    ),
    diff1p = paste0(
      formatC(diff1p, format = "f", digits = 3),
      "(",
      formatC(diff1p_sd, format = "f", digits = 3),
      ")"
    ),
    diff10p = paste0(
      formatC(diff10p, format = "f", digits = 3),
      "(",
      formatC(diff10p_sd, format = "f", digits = 3),
      ")"
    ),
    diff100p = paste0(
      formatC(diff100p, format = "f", digits = 3),
      "(",
      formatC(diff100p_sd, format = "f", digits = 3),
      ")"
    )
  ) %>%
  view()











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


