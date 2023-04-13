if (!require("pacman")) {
  install.packages("pacman")
}

pacman::p_load(tidyverse,
               lubridate,
               doParallel,
               airGR)


# seed --------------------------------------------------------------------

set.seed(1234)


# data --------------------------------------------------------------------

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


# function ----------------------------------------------------------------

catchment_ids <- data_process$catchment_id %>% unique()
OutputsCalibs <- vector("list", length(catchment_ids))

selected_catchment_id <- catchment_ids[[3]]

BasinObs <- data_process %>%
  dplyr::filter(catchment_id == selected_catchment_id)

N <- 500
sta_inds <- sample(x=c(365:3635), N, T) %>% sort()
end_inds <- sta_inds + 400
InputsModels <- vector("list", N) 
Ind_Runs <- vector("list", N)
RunOptionss <- vector("list", N)
InputsCrits <- vector("list", N)
OutputsCalibss <- vector("list", N)
paras <- vector("list", N)
sims <- vector("list", N)
  
for (i in 1:N){
  InputsModels[[i]] <- CreateInputsModel(FUN_MOD = RunModel_GR4J, DatesR = BasinObs$DatesR,
                                         Precip = BasinObs$P, PotEvap = BasinObs$E)
  
  Ind_Runs[[i]] <- sta_inds[[i]]:end_inds[[i]]
  
  RunOptionss[[i]] <- CreateRunOptions(FUN_MOD = RunModel_GR4J,
                                       InputsModel = InputsModels[[i]], IndPeriod_Run = Ind_Runs[[i]],
                                       IniStates = NULL, IniResLevels = NULL, IndPeriod_WarmUp = NULL)
  
  InputsCrits[[i]] <- CreateInputsCrit(FUN_CRIT = ErrorCrit_NSE, InputsModel = InputsModels[[i]], 
                                       RunOptions = RunOptionss[[i]], VarObs = "Q", Obs = BasinObs$Qmm[Ind_Runs[[i]]])
  
  
  CalibOptions <-
    CreateCalibOptions(
      FUN_MOD = RunModel_GR4J,
      FUN_CALIB = Calibration_Michel,
      SearchRanges = matrix(c(1, -10, 1, 1, 2000, 15, 300, 15), nrow = 2, byrow = T)
    )
  
  OutputsCalibss[[i]] <- Calibration_Michel(InputsModel = InputsModels[[i]], RunOptions = RunOptionss[[i]],
                                            InputsCrit = InputsCrits[[i]], CalibOptions = CalibOptions,
                                            FUN_MOD = RunModel_GR4J)
  
  paras[[i]] <- OutputsCalibss[[i]]$ParamFinalR
}

for (i in 1:N){
  Param = paras[[i]]
  
  RunOptions <- CreateRunOptions(FUN_MOD = RunModel_GR4J,
                                 InputsModel = InputsModels[[i]], IndPeriod_Run = 4000:8000,
                                 IniStates = NULL, IniResLevels = NULL, IndPeriod_WarmUp = NULL)
  
  OutputsModel <- RunModel_GR4J(InputsModel = InputsModel, RunOptions = RunOptions, Param = Param)
  
  sims[[i]] <- OutputsModel$Qsim
}

pred <- sims %>%
  unlist() %>%
  matrix(nrow = N, byrow = T) %>%
  colMeans()

hydroGOF::NSE(sim = pred, obs = BasinObs$Qmm[4000:8000])


-0.01239293
0.063689
0.08169528
0.08517038
0.08201663

0.6146744

# conventional

BasinObs <- data_process %>%
  dplyr::filter(catchment_id == selected_catchment_id)

InputsModel <- CreateInputsModel(FUN_MOD = RunModel_GR4J, DatesR = BasinObs$DatesR,
                                 Precip = BasinObs$P, PotEvap = BasinObs$E)

Ind_Run <- 365:3635

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

OutputsCalibs <- Calibration_Michel(InputsModel = InputsModel, RunOptions = RunOptions,
                                         InputsCrit = InputsCrit, CalibOptions = CalibOptions,
                                         FUN_MOD = RunModel_GR4J)

Param <-OutputsCalibs$ParamFinalR

RunOptions <- CreateRunOptions(FUN_MOD = RunModel_GR4J,
                               InputsModel = InputsModels[[i]], IndPeriod_Run = 4000:8000,
                               IniStates = NULL, IniResLevels = NULL, IndPeriod_WarmUp = NULL)

OutputsModel <- RunModel_GR4J(InputsModel = InputsModel, RunOptions = RunOptions, Param = Param)

hydroGOF::NSE(sim = OutputsModel$Qsim, obs = BasinObs$Qmm[4000:8000])










for (i in seq_along(catchment_ids)){
  selected_catchment_id <- catchment_ids[[i]]
  
  BasinObs <- data_process %>%
    dplyr::filter(catchment_id == selected_catchment_id)
  
  InputsModel <- CreateInputsModel(FUN_MOD = RunModel_GR4J, DatesR = BasinObs$DatesR,
                                   Precip = BasinObs$P, PotEvap = BasinObs$E)
  
  Ind_Run <- 365:3635
  
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
  
  OutputsCalibs[[i]] <- Calibration_Michel(InputsModel = InputsModel, RunOptions = RunOptions,
                                     InputsCrit = InputsCrit, CalibOptions = CalibOptions,
                                     FUN_MOD = RunModel_GR4J)
}

save(OutputsCalibs, file = "data/OutputsCalibs.Rda")


# check results -----------------------------------------------------------

sapply(OutputsCalibs, function(x) x$CritFinal) %>% hist()

sapply(OutputsCalibs, function(x) x$NRuns) %>% hist()
