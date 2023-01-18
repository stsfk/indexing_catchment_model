if (!require("pacman")) {
  install.packages("pacman")
}

pacman::p_load(tidyverse,
               lubridate,
               doParallel,
               airGR)


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
  
  
  out <- vector("double", nrow(paraset)) * NA
  for (i in 1:nrow(paraset)){
    Param = paraset[i,]
    
    OutputsModel <- RunModel_GR4J(InputsModel = InputsModel, RunOptions = RunOptions, Param = Param)
    OutputsCrit <- ErrorCrit(InputsCrit = InputsCrit, OutputsModel = OutputsModel, verbose = FALSE)
    out[[i]] = OutputsCrit$CritValue
  }
  
  out
}

catchment_ids <- data_process$catchment_id %>% unique()

paraset <- lapply(1:1000, function(x) random_para_gen()) %>%
  unlist() %>%
  matrix(ncol = 4, byrow = T)

get_catchment_gof_wrapper <- function(i){
  get_catchment_gof(catchment_ids[[i]], paraset)
}

start_time <- Sys.time()
outs <- foreach(x=1:14) %dopar% 
  get_catchment_gof_wrapper(x)
Sys.time() - start_time


start_time <- Sys.time()
outs <- foreach(x=1:14) %do% 
  get_catchment_gof_wrapper(x)
Sys.time() - start_time




out <- get_catchment_gof("01022500", paraset)
range(out)

# Stop --------------------------------------------------------------------

stopCluster(cl)




# Recycle -----------------------------------------------------------------


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


out = vector("list", 10000)
for (i in 1:10000){
  Param = random_para_gen()
  
  OutputsModel <- RunModel_GR4J(InputsModel = InputsModel, RunOptions = RunOptions, Param = Param)
  OutputsCrit <- ErrorCrit(InputsCrit = InputsCrit, OutputsModel = OutputsModel, verbose = FALSE)
  out[[i]] = OutputsCrit$CritValue
}

out %>% unlist %>% range()












Param = random_para_gen()

OutputsModel <- RunModel_GR4J(InputsModel = InputsModel, RunOptions = RunOptions, Param = Param)
OutputsCrit <- ErrorCrit(InputsCrit = InputsCrit, OutputsModel = OutputsModel)
OutputsCrit$CritValue




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

