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

# data --------------------------------------------------------------------

# load gr4j parameters
paras <- read.table("./data/gr4jparas.txt", header = FALSE, sep = " ") %>% data.matrix()

n_instances <- nrow(paras)

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

# NSE of random sampling method -------------------------------------------

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
  
  # adjustment for only 1 samples
  if (!is.matrix(selected_paras)){
    selected_paras <- matrix(selected_paras, nrow = 1)
  }
  
  out <- foreach(x=1:nrow(selected_paras)) %do%
    OutputsCrit_wrapper(x) %>%
    unlist()
  
  tibble(
    para_id = selected_para_ids,
    nse = out
  )
}

eval_grid_random_sample <- expand_grid(
  n_samples = c(1, 10, 100),
  catchment_id = catchments,
  selected_paras = vector("list", 1),
  results = vector("list", 1)
)

n_repeats <- 10
for (i in 1:nrow(eval_grid_random_sample)){
  
  selected_catchment_id <- eval_grid_random_sample$catchment_id[[i]]
  n_samples <- eval_grid_random_sample$n_samples[[i]]
  
  selected_para_ids_list <- lapply(1:n_repeats, function(x) sample(1:1e6, n_samples))
  
  temp_results <- vector("list", n_repeats) # temporary store results
  for (j in 1:n_repeats){
    selected_para_ids <- selected_para_ids_list[[j]]
    selected_paras <- paras[selected_para_ids,]
    
    temp_results[[j]] <- get_catchment_gof(
      selected_catchment_id = selected_catchment_id,
      selected_paras = selected_paras,
      selected_para_ids = selected_para_ids
    )
  }
  
  eval_grid_random_sample$selected_paras[[i]] <-  selected_para_ids_list %>% unlist()
  eval_grid_random_sample$results[[i]] <- sapply(temp_results, function(x) max(x$nse))
}

save(eval_grid_random_sample, file = "data/caravan_random_sample_retrieval.Rda")


# Plot --------------------------------------------------------------------

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

# load calibration, retrial, and random results
load("./data/OutputsCalibs.Rda")
load("./data/carvan_look_up.Rda")
load("data/caravan_random_sample_retrieval.Rda")

calibrated_NSE <- tibble(
  catchment_id = catchments,
  nse = sapply(OutputsCalibs, function(x) x$CritFinal)
) %>%
  mutate(rating = 1 / (2 - nse) * 10)


prepare_data_plot <- function(n_retrieved){
  data_plot <- eval_grid %>%
    dplyr::filter(!is.null(results)) %>%
    unnest(cols = c(results)) %>%
    filter(rank %in% c(1:n_retrieved)) %>%
    mutate(rating = 1 / (2 - nse) * 10) %>%
    group_by(catchment_id, repeats, n_probed) %>%
    summarise(retr_rating = max(rating)) 
  
  data_plot %>%
    group_by(catchment_id, n_probed) %>%
    summarise(
      max_rating = max(retr_rating),
      mean_rating = mean(retr_rating),
      min_rating = min(retr_rating)
    ) %>%
    left_join(calibrated_NSE, by="catchment_id") %>%
    mutate(n_retrieved = n_retrieved)
}


eval_grid_random_sample_processed <- eval_grid_random_sample %>%
  mutate(
    max_rating_random = map_dbl(results, function(x)
      max(x)),
    mean_rating_random = map_dbl(results, function(x)
      mean(x))
  ) %>%
  transmute(
    catchment_id = catchment_id,
    n_retrieved = n_samples,
    max_rating_random = 1 / (2 - max_rating_random) * 10,
    mean_rating_random = 1 / (2 - mean_rating_random) * 10
  )
  
data_plot<-lapply(c(1,10,100), prepare_data_plot) %>%
  bind_rows()

data_plot <- data_plot %>%
  left_join(eval_grid_random_sample_processed, by = c("catchment_id", "n_retrieved"))

data_plot <- data_plot %>%
  mutate(n_probed = factor(
    n_probed,
    labels = paste0("n/d=", c(0.5, 1, 2, 4), "\n(", c(18,37,74,148), " instances sampled)")
  )) %>%
  mutate(n_retrieved = factor(
    n_retrieved
  )) %>%
  ungroup()

data_plot_text <- data_plot %>% 
  slice(1) %>%
  mutate(text = "reference line y=x",
         mean_rating = 3,
         rating = 2.25)

ggplot(data_plot, aes(mean_rating_random, mean_rating, color = n_retrieved, shape = n_retrieved))+
  geom_point(alpha = 0.8, size = 0.75, stroke = 0.35)+
  geom_abline(slope = 1)+
  geom_text(data = data_plot_text, aes(mean_rating,rating, label= text), size = 2.5, hjust = 0, color = "black")+
  facet_wrap(~n_probed)+
  scale_color_manual(values = c("#377eb8", "#e41a1c", "#4daf4a"))+
  scale_shape(solid = FALSE)+
  labs(color = "Number of model(s) retrieved/randomly sampled",
       shape = "Number of model(s) retrieved/randomly sampled") +
  xlab((bquote(Assocation~r["i,j"]~value~between~randomly~sampled~model~instance~and~a~catchment))) +
  ylab((bquote(Assocation~r["i,j"]~value~between~retrieved~model~instance~and~a~catchment))) + 
  theme_bw(9)+
  theme(legend.position = "top")

ggsave(filename = "data/plot/fig_retrieve_baseline.png", width = 17, height = 12, units = "cm")
ggsave(filename = "data/plot/fig_retrieve_baseline.pdf", width = 17, height = 12, units = "cm")



# baseline vs. calibrated
ggplot(data_plot, aes(mean_rating_random, rating, color = n_retrieved, shape = n_retrieved))+
  geom_point(alpha = 0.8, size = 0.75, stroke = 0.35)+
  geom_abline(slope = 1)+
  geom_text(data = data_plot_text, aes(mean_rating,rating, label= text), size = 2.5, hjust = 0, color = "black")+
  facet_wrap(~n_probed)+
  scale_color_manual(values = c("#377eb8", "#e41a1c", "#4daf4a"))+
  scale_shape(solid = FALSE)+
  labs(color = "Number of model(s) randomly sampled",
       shape = "Number of model(s) randomly sampled") +
  xlab((bquote(Assocation~r["i,j"]~value~between~randomly~sampled~model~instance~and~a~catchment))) +
  ylab((bquote(Assocation~r["i,j"]~value~between~calibrated~model~instance~and~a~catchment))) + 
  theme_bw(9)+
  theme(legend.position = "top")

ggsave(filename = "data/plot/fig_retrieve_baseline_cali.png", width = 17, height = 12, units = "cm")
ggsave(filename = "data/plot/fig_retrieve_baseline_cali.pdf", width = 17, height = 12, units = "cm")



# baseline max vs. retrieved max

ggplot(data_plot, aes(max_rating_random, max_rating, color = n_retrieved, shape = n_retrieved))+
  geom_point(alpha = 0.8, size = 0.75, stroke = 0.35)+
  geom_abline(slope = 1)+
  geom_text(data = data_plot_text, aes(mean_rating,rating, label= text), size = 2.5, hjust = 0, color = "black")+
  facet_wrap(~n_probed)+
  scale_color_manual(values = c("#377eb8", "#e41a1c", "#4daf4a"))+
  scale_shape(solid = FALSE)+
  labs(color = "Number of model(s) randomly sampled",
       shape = "Number of model(s) randomly sampled") +
  xlab((bquote(Maximum~assocation~r["i,j"]~value~between~randomly~sampled~model~instance~and~a~catchment))) +
  ylab((bquote(Maximum~assocation~r["i,j"]~value~between~retrieved~model~instance~and~a~catchment))) + 
  theme_bw(9)+
  theme(legend.position = "top")

ggsave(filename = "data/plot/fig_retrieve_max.png", width = 17, height = 13, units = "cm")
ggsave(filename = "data/plot/fig_retrieve_max.pdf", width = 17, height = 13, units = "cm")



