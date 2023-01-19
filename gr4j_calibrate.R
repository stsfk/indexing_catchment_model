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

for (i in seq_along(catchment_ids)){
  selected_catchment_id <- catchment_ids[[i]]
  
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

save(OutputsCalibs, file = "OutputsCalibs.Rda")


sapply(OutputsCalibs, function(x) x$CritFinal) %>% hist()

sapply(OutputsCalibs, function(x) x$NRuns) %>% hist()


