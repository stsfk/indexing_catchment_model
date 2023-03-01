# setting experiment environment ------------------------------------------

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
               airGR,
               cowplot)

set.seed(12345)

# data --------------------------------------------------------------------

# load gr4j parameters
paras <- read.table("./data/gr4jparas.txt", header = FALSE, sep = " ") %>% data.matrix()

n_instances <- nrow(paras)

# load PQ data
P <- read.table("./data/mat_P_1.txt", header = FALSE, sep = " ") %>% data.matrix()
Q <- read.table("./data/mat_Q_1.txt", header = FALSE, sep = " ") %>% data.matrix()

latent_dim <- dim(Q)[2]

# load L0123001 data from airGR package
data(L0123001)

# Model calibration -------------------------------------------------------
# The code in this section is from the "Get Started with airGR" user guide included in the
# airGR package.

InputsModel <- CreateInputsModel(FUN_MOD = RunModel_GR4J, DatesR = BasinObs$DatesR,
                                 Precip = BasinObs$P, PotEvap = BasinObs$E)

Ind_Run <- seq(which(format(BasinObs$DatesR, format = "%Y-%m-%d") == "1990-01-01"), 
               which(format(BasinObs$DatesR, format = "%Y-%m-%d") == "1999-12-31"))

RunOptions <- CreateRunOptions(FUN_MOD = RunModel_GR4J,
                               InputsModel = InputsModel, IndPeriod_Run = Ind_Run,
                               IniStates = NULL, IniResLevels = NULL, IndPeriod_WarmUp = NULL)

InputsCrit <- CreateInputsCrit(FUN_CRIT = ErrorCrit_NSE, InputsModel = InputsModel, 
                               RunOptions = RunOptions, VarObs = "Q", Obs = BasinObs$Qmm[Ind_Run])

CalibOptions <-
  CreateCalibOptions(
    FUN_MOD = RunModel_GR4J,
    FUN_CALIB = Calibration_Michel,
    SearchRanges = matrix(c(1, -10, 1, 1, 2000, 15, 300, 15), nrow = 2, byrow = T)
  )

OutputsCalib <- Calibration_Michel(InputsModel = InputsModel, RunOptions = RunOptions,
                                   InputsCrit = InputsCrit, CalibOptions = CalibOptions,
                                   FUN_MOD = RunModel_GR4J)

calibrated_para <- OutputsCalib$ParamFinalR # calibrated parameter
calibrated_nse <- OutputsCalib$CritFinal

# Function ----------------------------------------------------------------

p_rmse <- function(p, Q, rating) {
  # This function computes the RMSE of predicted ratings of models specified by Q
  ModelMetrics::rmse(predicted = p %*% t(Q) %>% as.vector(),
                     actual = rating)
}

para_gof <- function(para){
  # compute the gof metric for the parameter set "para"
  # InputsModel and RunOptions are defined in the global environment
  
  OutputsModel <- RunModel_GR4J(InputsModel = InputsModel, RunOptions = RunOptions, Param = para)
  OutputsCrit <- ErrorCrit(InputsCrit = InputsCrit, OutputsModel = OutputsModel, verbose = FALSE)
  OutputsCrit$CritValue
}

para_ids_gof <- function(selected_para_ids){
  # get the ratings for the model instances specified by selected_para_ids
  selected_paras <- paras[selected_para_ids,]
  
  nses <- foreach(i=1:nrow(selected_paras)) %do%
    para_gof(selected_paras[i,]) %>%
    unlist()
  
  tibble(
    para_id = selected_para_ids,
    nse = nses,
    rating = 1 / (2 - nses) * 10
  )
}

prepare_Qs_for_retrieval  <- function(n_probed, train_portion){
  # This function split Q into Q_probed, Q_train, Q_val
  # weights of the links between the look-up catchment and the models associated with Q_probed is known
  # Q_probed is further split into Q_train and Q_val for deriving the optimal number of iterations in GA

  # n_probed: the number of links between Q and the look-up catchment with known weights
  
  selected_para_ids <- sample(1:n_instances, n_probed) %>% sort()
  rating_probed <- para_ids_gof(selected_para_ids) %>% pull(rating)
  Q_probed <- Q[selected_para_ids,]
  
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
    pred <- Rfast::rowsums(pred, parallel = T)
    
    - ModelMetrics::rmse(actual = rating, predicted = pred)
  }
}

derive_p <- function(n_probed, train_portion){
  # train_portion: the portion used in training during GA hyperparameter optimization
  
  # prepare data sets for retrieval
  c(Q_probed,
    rating_probed,
    Q_train,
    rating_train,
    Q_val,
    rating_val) %<-%
    prepare_Qs_for_retrieval(n_probed, train_portion)
  
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

top_n_nse <- function(p, n_retrieved){
  
  pred <- p %*% t(Q) %>% as.vector()
  
  top_n <- tibble(
    model_id = 1:n_instances,
    rating = pred
  ) %>%
    arrange(desc(rating)) %>%
    slice(1:n_retrieved) %>%
    pull(model_id)
  
  gof <- para_ids_gof(top_n)
  
  # output
  tibble(
    para_id = top_n,
    rating = gof$rating,
    nse = gof$nse,
    rank = 1:n_retrieved
  )
}

# Modeling ----------------------------------------------------------------

n_probed <- latent_dim * 4
n_retrieved <- 200
train_portion <- 0.8

p <- derive_p(n_probed, train_portion)
retrieval_result <- top_n_nse(p, n_retrieved)

retrieved_nse <- retrieval_result$nse %>% max()
retrieved_para <- paras[retrieval_result$para_id[which.max(retrieval_result$nse)],] %>% unname()

# Result analysis ---------------------------------------------------------

Q_calibrated <- RunModel_GR4J(InputsModel = InputsModel, RunOptions = RunOptions, Param = calibrated_para)$Qsim
Q_retrieved <- RunModel_GR4J(InputsModel = InputsModel, RunOptions = RunOptions, Param = retrieved_para)$Qsim
Q_actual <- BasinObs[Ind_Run,] %>% pull(Qmm)

data_plot <- tibble(
  date = seq(from = ymd("1990-01-01"), to = ymd("1999-12-31"), by = 1),
  calibrated = Q_calibrated,
  retrieved = Q_retrieved,
  actual = Q_actual
)

p1 <- data_plot %>% 
  gather(item, value, -date) %>%
  ggplot(aes(date, value, color = item, linetype = item))+
  geom_line(size = 0.25)+
  labs(color = "",
       linetype = "",
       y = "Flow [mm/day]") +
  theme_bw(base_size = 10)+
  theme(legend.position = "top",
        axis.title.x = element_blank())
  
p1 +
  cowplot::draw_text(
    x = ymd("1990-01-01"),
    y = 19,
    size = 10,
    hjust = 0,
    text = paste0("NSE of calibrated model = ", round(calibrated_nse,3))
  )+
  cowplot::draw_text(
    x = ymd("1990-01-01"),
    y = 17.5,
    size = 10,
    hjust = 0,
    text = paste0("NSE of retrieved model = ", round(retrieved_nse,3))
  )

ggsave(filename = "./data/plot/fig_tutorial.pdf", width = 15, height = 8, unit = "cm")

