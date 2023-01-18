if (!require("pacman")) {
  install.packages("pacman")
}

pacman::p_load(tidyverse,
               lubridate,
               doParallel,
               airGR,
               recosystem)


# Parallel setting --------------------------------------------------------

nthread <- detectCores()
cl <- makeCluster(nthread-2)
registerDoParallel(cl)

clusterEvalQ(cl, library(lubridate))
clusterEvalQ(cl, library(airGR))
clusterEvalQ(cl, library(tidyverse))


# Data --------------------------------------------------------------------

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


# Functions ---------------------------------------------------------------

random_para_gen <- function(n_para = 4, LB = c(1, -10, 1, 1), UB = c(2000, 15, 300, 15)){
  runif(n_para) * (UB - LB) + LB
}

get_catchment_gof <- function(selected_catchment_id, paraset){
  
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
    Param = paraset[i,]
    
    OutputsModel <- RunModel_GR4J(InputsModel = InputsModel, RunOptions = RunOptions, Param = Param)
    OutputsCrit <- ErrorCrit(InputsCrit = InputsCrit, OutputsModel = OutputsModel, verbose = FALSE)
    OutputsCrit$CritValue
  }
  
  out <- foreach(x=1:nrow(paraset)) %dopar%
    OutputsCrit_wrapper(x) %>%
    unlist()
  
  tibble(
    para_id = 1:nrow(paraset),
    nse = out
  )
}

# Processing --------------------------------------------------------------

# generate data

paraset <- lapply(1:1000, function(x) random_para_gen()) %>%
  unlist() %>%
  matrix(ncol = 4, byrow = T)

catchment_ids <- data_process$catchment_id %>% unique()

results <- tibble(
  catchment_id = catchment_ids,
  out = vector("list", 1)
)

for (i in 1:nrow(results)){
  selected_catchment_id = catchment_ids[[i]]
  results$out[[i]] <- get_catchment_gof(selected_catchment_id, paraset)
}

results <-results %>% unnest(out) %>%
  mutate(rating = 1 / (2 - nse) * 10) %>%
  dplyr::select(catchment_id,
                para_id,
                rating) %>%
  mutate(record_id = 1:n()) # give each row a unique ID

stopCluster(cl)

# matrix factorization



data_train <- tibble(
  user_index = factor(results$catchment_id) %>% as.numeric(),
  item_index = results$para_id,
  rating = results$rating
)

train_set <- data_memory(
  user_index = data_train$user_index,
  item_index = data_train$item_index,
  rating = data_train$rating,
  index1 = T
)

r = Reco()

opts = r$tune(train_set, opts = list(dim = c(50), lrate = c(0.01),
                                     costp_l1 = 0.01, costq_l1 = 0.01,
                                     costp_l2 = 0, costq_l2 = 0,
                                     nthread = 10, niter = 200))
opts


r$train(
  train_set,
  opts = list(
    dim = 100,
    costp_l1 = 0.01,
    costp_l2 = 0,
    costq_l1 = 0.01,
    costq_l2 = 0,
    lrate = 0.005,
    niter = 100,
    verbose = F,
    nthread = 10
  )
)

cor(r$predict(train_set), data_train$rating)^2


# Stop --------------------------------------------------------------------

stopCluster(cl)




# Recycle -----------------------------------------------------------------


selected_catchment_id = catchment_ids[[1]]

start_time <- Sys.time()
out <- get_catchment_gof(selected_catchment_id, paraset)
Sys.time() - start_time

out %>% range()


catchment_ids <- data_process$catchment_id %>% unique()

selected_catchment_id = catchment_ids[[1]]

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


CalibOptions <- CreateCalibOptions(FUN_MOD = RunModel_GR4J, FUN_CALIB = Calibration_Michel)

OutputsCalib <- Calibration_Michel(InputsModel = InputsModel, RunOptions = RunOptions,
                                   InputsCrit = InputsCrit, CalibOptions = CalibOptions,
                                   FUN_MOD = RunModel_GR4J)
Param <- OutputsCalib$ParamFinalR

OutputsModel <- RunModel_GR4J(InputsModel = InputsModel, RunOptions = RunOptions, Param = Param)
str(OutputsModel)

plot(OutputsModel, Qobs = BasinObs$Qmm[Ind_Run])










## loading catchment data
data(L0123001)

## loading generalist parameter sets
data(Param_Sets_GR4J)


Param_Sets_GR4J$X4 <- Param_Sets_GR4J$X4u / 5.995 * BasinInfo$BasinArea^0.3
Param_Sets_GR4J$X4u <- NULL
Param_Sets_GR4J <- as.matrix(Param_Sets_GR4J)

summary(Param_Sets_GR4J)

## preparation of the InputsModel object
InputsModel <- CreateInputsModel(FUN_MOD = RunModel_GR4J, DatesR = BasinObs$DatesR, 
                                 Precip = BasinObs$P, PotEvap = BasinObs$E)


Ind_Cal <- seq(which(format(BasinObs$DatesR, format = "%d/%m/%Y %H:%M")=="01/01/1990 00:00"), 
               which(format(BasinObs$DatesR, format = "%d/%m/%Y %H:%M")=="28/02/1990 00:00"))

## preparation of the RunOptions object for the calibration period
RunOptions_Cal <- CreateRunOptions(FUN_MOD = RunModel_GR4J,
                                   InputsModel = InputsModel, IndPeriod_Run = Ind_Cal)

InputsCrit_Cal  <- CreateInputsCrit(FUN_CRIT = ErrorCrit_NSE, InputsModel = InputsModel, 
                                    RunOptions = RunOptions_Cal, Obs = BasinObs$Qmm[Ind_Cal])


Ind_Val <- seq(which(format(BasinObs$DatesR, format = "%d/%m/%Y %H:%M")=="01/03/1990 00:00"), 
               which(format(BasinObs$DatesR, format = "%d/%m/%Y %H:%M")=="31/12/1999 00:00"))

## preparation of the RunOptions object for the validation period
RunOptions_Val <- CreateRunOptions(FUN_MOD = RunModel_GR4J,
                                   InputsModel = InputsModel, IndPeriod_Run = Ind_Val)

## efficiency criterion (Nash-Sutcliffe Efficiency) on the validation period
InputsCrit_Val  <- CreateInputsCrit(FUN_CRIT = ErrorCrit_NSE, InputsModel = InputsModel, 
                                    RunOptions = RunOptions_Val, Obs = BasinObs$Qmm[Ind_Val])

OutputsCrit_Loop <- apply(Param_Sets_GR4J, 1, function(iParam) {
  OutputsModel_Cal <- RunModel_GR4J(InputsModel = InputsModel, RunOptions = RunOptions_Cal,
                                    Param = iParam)
  OutputsCrit <- ErrorCrit_NSE(InputsCrit = InputsCrit_Cal, OutputsModel = OutputsModel_Cal)
  return(OutputsCrit$CritValue)
})


CalibOptions <- CreateCalibOptions(FUN_MOD = RunModel_GR4J, FUN_CALIB = Calibration_Michel)

Ind_Run <- seq(which(format(BasinObs$DatesR, format = "%Y-%m-%d") == "1990-01-01"), 
               which(format(BasinObs$DatesR, format = "%Y-%m-%d") == "1999-12-31"))
str(Ind_Run)


RunOptions <- CreateRunOptions(FUN_MOD = RunModel_GR4J,
                               InputsModel = InputsModel, IndPeriod_Run = Ind_Run,
                               IniStates = NULL, IniResLevels = NULL, IndPeriod_WarmUp = NULL)

InputsCrit <- CreateInputsCrit(FUN_CRIT = ErrorCrit_NSE, InputsModel = InputsModel, 
                               RunOptions = RunOptions, VarObs = "Q", Obs = BasinObs$Qmm[Ind_Run])
str(InputsCrit)


Param <- OutputsCalib$ParamFinalR


OutputsCalib <- Calibration_Michel(InputsModel = InputsModel, RunOptions = RunOptions,
                                   InputsCrit = InputsCrit, CalibOptions = CalibOptions,
                                   FUN_MOD = RunModel_GR4J)

OutputsModel <- RunModel_GR4J(InputsModel = InputsModel, RunOptions = RunOptions, Param = Param)
str(OutputsModel)

plot(OutputsModel, Qobs = BasinObs$Qmm[Ind_Run])

