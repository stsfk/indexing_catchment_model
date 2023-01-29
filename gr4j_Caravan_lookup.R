if (!require("pacman")) {
  install.packages("pacman")
}

pacman::p_load(tidyverse,
               lubridate,
               zeallot,
               recosystem,
               GA,
               ModelMetrics,
               doParallel,
               Rfast,
               airGR)

# seed --------------------------------------------------------------------

set.seed(1234)

nthread <- detectCores()
cl <- makeCluster(nthread-2)
registerDoParallel(cl)

clusterEvalQ(cl, {
  pacman::p_load(tidyverse,
                 lubridate,
                 zeallot,
                 recosystem,
                 GA,
                 ModelMetrics,
                 doParallel,
                 Rfast,
                 airGR)
})


# data --------------------------------------------------------------------

# load gr4j parameters
paras <- read.table("./data/gr4jparas.txt", header = FALSE, sep = " ") %>% data.matrix()

n_instances <- nrow(paras)

# load PQ data
P <- read.table("./data/mat_P_1.txt", header = FALSE, sep = " ") %>% data.matrix()
Q <- read.table("./data/mat_Q_1.txt", header = FALSE, sep = " ") %>% data.matrix()

latent_dim <- dim(Q)[2]

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

# calibrated NSEs
load("./data/OutputsCalibs.Rda")

calibrated_NSE <- tibble(
  catchment_id = catchments,
  nse = sapply(OutputsCalibs, function(x) x$CritFinal)
)

# Function ----------------------------------------------------------------

p_rmse <- function(p, Q, rating) {
  # This function computes the RMSE of predicted ratings of models specified by Q
  # the target rating is `rating`
  ModelMetrics::rmse(predicted = p %*% t(Q) %>% as.vector(),
                     actual = rating)
}

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
  
  out <- foreach(x=1:nrow(selected_paras)) %do%
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
    pred <- Rfast::rowsums(pred, parallel = F)
    
    - ModelMetrics::rmse(actual = rating, predicted = pred)
  }
}

derive_p <- function(n_probed, train_portion, selected_catchment_id){
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
  GA@solution[1,] %>% unname()
}

top_n_nse <- function(p, n, selected_catchment_id){
  
  pred <- p %*% t(Q) %>% as.vector()
  
  top_n <- tibble(
    model_id = 1:n_instances,
    rating = pred,
    case = "pred"
  ) %>%
    arrange(desc(rating)) %>%
    slice(1:n) %>%
    pull(model_id)
  
  top_n_nse <- get_catchment_gof(selected_catchment_id,
                                  selected_para_ids = top_n,
                                  selected_paras = paras[top_n, ])%>% 
    pull(nse)
  
  # output
  tibble(
    para_id = top_n,
    nse = top_n_nse,
    rank = 1:n
  )
}

eval_grid <- expand_grid(
  selected_catchment_id = catchments,
  n_probed = round(latent_dim *c(0.5,1,2,4)),
  repeats = 1:5,
  results = vector("list",1)
)

catchment_top_n_nse_wrapper <- function(i, n_retrieved = 200) {
  
  n_probed <- eval_grid$n_probed[[i]]
  selected_catchment_id <- eval_grid$selected_catchment_id[[i]]
  
  p <- derive_p(n_probed, train_portion=0.8, selected_catchment_id)
  
  top_n_nse(p, n_retrieved, selected_catchment_id)
}

# Modeling ----------------------------------------------------------------

eval_grid <- expand_grid(
  selected_catchment_id = catchments,
  n_probed = round(latent_dim *c(0.5,1,2,4)),
  repeats = 1:10,
  results = vector("list",1)
)

n_eval <- nrow(eval_grid)

sta_time <- Sys.time()
out <- foreach(x=1:n_eval) %dopar%
  catchment_top_n_nse_wrapper(x)
end_time <- Sys.time()
dif1 <- end_time - sta_time

eval_grid$results[1:n_eval] <- out

# Stop --------------------------------------------------------------------
stopCluster(cl)

# Result analysis ---------------------------------------------------------
eval_grid <- eval_grid %>%
  rename(catchment_id = selected_catchment_id)

data_plot <- eval_grid %>%
  dplyr::filter(!is.null(results)) %>%
  unnest(cols = c(results)) %>%
  group_by(catchment_id, repeats, n_probed) %>%
  summarise(retr_nse = max(nse)) 

data_plot <- data_plot %>%
  group_by(catchment_id, n_probed) %>%
  summarise(
    max_nse = max(retr_nse),
    mean_nse = mean(retr_nse),
    min_nse = min(retr_nse)
  ) %>%
  left_join(calibrated_NSE, by="catchment_id")

ggplot(data_plot, aes(mean_nse, nse))+
  geom_point(color = "steelblue")+
  geom_linerange(aes(xmin = min_nse, xmax = max_nse))+
  geom_abline(slope = 1)+
  coord_cartesian(xlim = c(0,1), ylim = c(0,1))+
  labs(x = "NSE of retrieved models",
       y = "NSE of calibrated models")+
  facet_wrap(~n_probed, nrow = 1)+
  theme_bw()
