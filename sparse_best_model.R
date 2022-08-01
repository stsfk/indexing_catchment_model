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
               mlrMBO,
               ModelMetrics)

# seed --------------------------------------------------------------------

set.seed(1234)

model_classes  <- c(
  'm_01_collie1_1p_1s',
  'm_02_wetland_4p_1s',
  'm_03_collie2_4p_1s',
  'm_04_newzealand1_6p_1s',
  'm_05_ihacres_7p_1s',
  'm_06_alpine1_4p_2s',
  'm_07_gr4j_4p_2s',
  'm_08_us1_5p_2s',
  'm_09_susannah1_6p_2s',
  'm_10_susannah2_6p_2s',
  'm_11_collie3_6p_2s',
  'm_12_alpine2_6p_2s',
  'm_13_hillslope_7p_2s',
  'm_14_topmodel_7p_2s',
  'm_15_plateau_8p_2s',
  'm_16_newzealand2_8p_2s',
  'm_17_penman_4p_3s',
  'm_18_simhyd_7p_3s',
  'm_19_australia_8p_3s',
  'm_20_gsfb_8p_3s',
  'm_21_flexb_9p_3s',
  'm_22_vic_10p_3s',
  'm_24_mopex1_5p_4s',
  'm_25_tcm_6p_4s',
  'm_26_flexi_10p_4s',
  'm_27_tank_12p_4s',
  'm_28_xinanjiang_12p_4s',
  'm_29_hymod_5p_5s',
  'm_30_mopex2_7p_5s',
  'm_31_mopex3_8p_5s',
  'm_32_mopex4_10p_5s',
  'm_34_flexis_12p_5s',
  'm_35_mopex5_12p_5s',
  'm_36_modhydrolog_15p_5s',
  'm_37_hbv_15p_5s',
  'm_40_smar_8p_6s'
)


# data --------------------------------------------------------------------

model_index <- read_csv("./data/sparse_index.csv",
                        col_names = c("model_class", "catchment_id"))

weights <- read_csv("./data/sparse_result.csv",
                    col_names = c("KGE", "NSE", "RMSE"))


n_catchments <- 533
n_model_classes <- (model_index$model_class %>% unique() %>% length()) 

n_model_instances <- 1000 # n_model_instances
n_catchment_per_instance <-  nrow(weights)/n_model_classes/n_model_instances # n_catchment_per_instance 

# assign instance_id to model_index data base
model_index <- model_index %>%
  group_by(model_class) %>%
  mutate(instance_id = rep(1:n_model_instances, each = n_catchment_per_instance)) %>%
  ungroup()

data_raw <- model_index %>%
  bind_cols(weights) %>%
  filter(model_class != "m_40_smar_8p_6s") %>% # evaluation errors for "m_40_smar_8p_6s"
  mutate(
    model_id = paste(model_class, instance_id, sep = "_"),
    model_id = factor(model_id, levels = unique(model_id)),
    model_id = as.integer(model_id) # assign a unique id to each model
  ) %>%
  select(model_class, catchment_id, model_id, KGE, NSE, RMSE)

data_process <- data_raw %>%
  mutate(NNSE = 1 / (2 - NSE) * 10) %>%
  dplyr::select(catchment_id,
                model_id,
                NNSE) %>%
  rename(rating = NNSE) %>%
  mutate(record_id = 1:n())


load("./data/sparse_factorization.Rda")


# Experiments -------------------------------------------------------------

load_PQ <- function(i){
  list(P = eval_grid$out[[i]]$P,
       Q = eval_grid$out[[i]]$Q)
}

c(P,Q) %<-% load_PQ(1)

temp <- P %*% t(Q)

temp2 <- temp %>%
  as.data.frame() %>%
  as_tibble() %>%
  mutate(catchment_id = 1:n()) %>%
  gather(model_id, rating, -catchment_id) %>%
  mutate(model_id = str_sub(model_id, start = 2) %>% as.numeric())


temp2 <- temp2 %>%
  group_by(catchment_id) %>%
  arrange(desc(rating)) %>%
  slice(1:100) %>%
  ungroup()

write_csv(temp2, "./data/catchment_model_pair.csv")


