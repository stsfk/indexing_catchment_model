if (!require("pacman")) {
  install.packages("pacman")
}

pacman::p_load(tidyverse,
               lubridate,
               doParallel,
               airGR,
               recosystem)


# seed --------------------------------------------------------------------

set.seed(1234)

# Parallel setting --------------------------------------------------------

nthread <- detectCores()
cl <- makeCluster(nthread-2)
registerDoParallel(cl)

clusterEvalQ(cl, library(lubridate))
clusterEvalQ(cl, library(airGR))
clusterEvalQ(cl, library(tidyverse))


# Data --------------------------------------------------------------------
# This file contains processed P, T, PET, Q data for 2116 catchments outside the US
data_raw <- read_csv("./data/processed_caravan_data.csv") 

# Some PET is negative, causing problems in GR4J modeling
data_raw$PET[data_raw$PET < 0] = 0

# Rename variables, following GR4J convention
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
  # This function generates one set of random parameters within the ranges specified by LB and UB
  # n_para: number of parameters
  # LB and UB: lower and upper bound of parameters
  runif(n_para) * (UB - LB) + LB
}


get_catchment_gof <- function(selected_catchment_id, selected_para_ids, selected_paras){
  # This function computes the goodness-of-fit metrics for selected catchment and  selected paras
  # "selected_para_ids" give the id of the selected parameter set 
  
  # Extract observation of the selected catchment to create "InputsModel" object
  BasinObs <- data_process %>%
    dplyr::filter(catchment_id == selected_catchment_id)
  
  InputsModel <- CreateInputsModel(FUN_MOD = RunModel_GR4J, DatesR = BasinObs$DatesR,
                                   Precip = BasinObs$P, PotEvap = BasinObs$E)
  
  # Specify run period
  Ind_Run <- seq(which(format(BasinObs$DatesR, format = "%Y-%m-%d") == "1990-01-01"), 
                 which(format(BasinObs$DatesR, format = "%Y-%m-%d") == "2010-12-31"))
  
  # Using default run options
  RunOptions <- CreateRunOptions(FUN_MOD = RunModel_GR4J,
                                 InputsModel = InputsModel, IndPeriod_Run = Ind_Run,
                                 IniStates = NULL, IniResLevels = NULL, IndPeriod_WarmUp = NULL)
  
  # Using NSE as goodness-of-fit metrics
  InputsCrit <- CreateInputsCrit(FUN_CRIT = ErrorCrit_NSE, InputsModel = InputsModel, 
                                 RunOptions = RunOptions, VarObs = "Q", Obs = BasinObs$Qmm[Ind_Run])
  
  # This function computes the NSE for the i-th parameter set
  OutputsCrit_wrapper <- function(i){
    Param = selected_paras[i,]
    
    OutputsModel <- RunModel_GR4J(InputsModel = InputsModel, RunOptions = RunOptions, Param = Param)
    OutputsCrit <- ErrorCrit(InputsCrit = InputsCrit, OutputsModel = OutputsModel, verbose = FALSE)
    OutputsCrit$CritValue
  }
  
  # Run OutputsCrit_wrapper in parallel for computing NSE of selected_paras
  out <- foreach(x=1:nrow(selected_paras)) %dopar%
    OutputsCrit_wrapper(x) %>%
    unlist()
  
  # Outputting selected para ids and the corresponding NSE
  tibble(
    para_id = selected_para_ids,
    nse = out
  )
}

# Processing --------------------------------------------------------------

# 1,000,000 model instances, and 5% of the catchment-model instance associations are explicitly assessed
n_paras <- 1000000
sample_frac <- 0.05

# generate random parameters of 1,000,000 model instances
paras <- lapply(1:n_paras, function(x) random_para_gen()) %>%
  unlist() %>%
  matrix(ncol = 4, byrow = T)

# list all unique catchment-model instance pairs, for each instance, select 5% of the instances
results <- expand_grid(catchment_id = data_process$catchment_id %>% unique(),
                       para_id = 1:n_paras) %>%
  group_by(para_id) %>%
  sample_frac(sample_frac) %>%
  ungroup() %>% 
  group_by(catchment_id) %>% # group and summarize by "catchment_id" to make the tibble more concise
  summarise(para_ids = list(para_id)) %>%
  mutate(out = vector("list", 1))

# run model simulations using "get_catchment_gof" function
for (i in 1:nrow(results)){
  selected_catchment_id <- results$catchment_id[[i]]
  selected_para_ids <- results$para_ids[[i]]
  selected_paras <- paras[selected_para_ids,]
  
  results$out[[i]] <- get_catchment_gof(selected_catchment_id, selected_para_ids, selected_paras)
}

# convert NSE to NNSE
results <-results %>% unnest(out) %>%
  mutate(rating = 1 / (2 - nse) * 10) %>%
  dplyr::select(catchment_id,
                model_id = para_id,
                rating) %>%
  mutate(record_id = 1:n()) # give each row a unique ID

# saving results
save(paras, results, file = "./data/gr4j_caravan.Rda")

write.table(paras, quote = F, sep = " ", row.names = F, col.names = F, file = "data/gr4jparas.txt")

# Stop --------------------------------------------------------------------

stopCluster(cl)
