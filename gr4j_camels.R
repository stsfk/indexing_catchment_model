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

get_catchment_gof <- function(selected_catchment_id, selected_para_ids, selected_paras, ver){
  
  # subsetting by catchment
  BasinObs <- data_process %>%
    dplyr::filter(catchment_id == selected_catchment_id)
  
  InputsModel <- CreateInputsModel(FUN_MOD = RunModel_GR4J, DatesR = BasinObs$DatesR,
                                   Precip = BasinObs$P, PotEvap = BasinObs$E)
  
  if (ver == 1){
    Ind_Run <- seq(which(format(BasinObs$DatesR, format = "%Y-%m-%d") == "1988-10-01"), 
                   which(format(BasinObs$DatesR, format = "%Y-%m-%d") == "1998-09-30"))
  } else{
    Ind_Run <- seq(which(format(BasinObs$DatesR, format = "%Y-%m-%d") == "1998-10-01"), 
                   which(format(BasinObs$DatesR, format = "%Y-%m-%d") == "2009-09-30"))
  }

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

# Processing --------------------------------------------------------------

n_paras <- 10000
sample_frac <- 0.05

# generate data

paras <- lapply(1:n_paras, function(x) random_para_gen()) %>%
  unlist() %>%
  matrix(ncol = 4, byrow = T)

results <- expand_grid(catchment_id = data_process$catchment_id %>% unique(),
                       para_id = 1:n_paras, ver = c(1:2)) %>%
  group_by(para_id)%>%
  sample_frac(sample_frac) %>%
  ungroup() %>% 
  group_by(catchment_id, ver) %>%
  summarise(para_ids = list(para_id)) %>%
  ungroup() %>%
  mutate(out = vector("list", 1)) 


for (i in 1:nrow(results)){
  
  selected_catchment_id <- results$catchment_id[[i]]
  selected_para_ids <- results$para_ids[[i]]
  selected_paras <- paras[selected_para_ids,]
  ver <- results$ver[[i]]
  
  results$out[[i]] <- get_catchment_gof(selected_catchment_id, selected_para_ids, selected_paras, ver)
}

results <-results %>% unnest(out) %>%
  mutate(rating = 1 / (2 - nse) * 10) %>%
  dplyr::select(catchment_id,
                ver,
                model_id = para_id,
                rating) %>%
  mutate(record_id = 1:n()) # give each row a unique ID

# saving results
save(paras, results, file = "./data/gr4j_camels.Rda")

write.table(paras, quote = F, sep = " ", row.names = F, col.names = F, file = "./data/gr4jparas_camels.txt")

# Stop --------------------------------------------------------------------

stopCluster(cl)
