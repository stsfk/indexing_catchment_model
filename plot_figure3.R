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
               ModelMetrics,
               ggsci)

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


# sparse recomemded model
catchment_model_pair <- read_csv("./data/catchment_model_pair.csv", col_types =  list(col_double()))
sparse_best <- read_csv("./data/best_fit.csv", col_types =  list(col_double()), col_names = F) %>%
  set_names("KGE", "NSE", "RMSE")

sparse_best <- catchment_model_pair %>%
  bind_cols(sparse_best) %>%
  group_by(catchment_id) %>%
  summarise(best = max(NSE)) %>%
  mutate(best = 1 / (2 - best) * 10,
         case = "Sparse data set (recommended)")


# load sparse data set
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

data_sparse <- data_raw %>%
  mutate(NNSE = 1 / (2 - NSE) * 10) %>%
  dplyr::select(catchment_id,
                model_id,
                NNSE) %>%
  rename(rating = NNSE) %>%
  mutate(record_id = 1:n())


# load sparse data set
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

data_sparse <- data_raw %>%
  mutate(NNSE = 1 / (2 - NSE) * 10) %>%
  dplyr::select(catchment_id,
                model_id,
                NNSE) %>%
  rename(rating = NNSE) %>%
  mutate(record_id = 1:n())


# load dense dataset

weights <- read_csv("./data/SM.csv",
                    col_names = c("KGE", "NSE", "RMSE"))

# 10 model classes in the dense data set
model_class  <- c(
  'm_01_collie1_1p_1s',
  'm_05_ihacres_7p_1s',
  'm_07_gr4j_4p_2s',
  'm_13_hillslope_7p_2s',
  'm_18_simhyd_7p_3s',
  'm_22_vic_10p_3s',
  'm_27_tank_12p_4s',
  'm_28_xinanjiang_12p_4s',
  'm_34_flexis_12p_5s',
  'm_37_hbv_15p_5s'
)

n_catchments <- 533
n_model_classes <- length(model_class)
n_model_instances <-  nrow(weights)/n_catchments/n_model_classes # number of instance of each model class

# assign catchment id and model instance id to each row
catchment_id <- rep(1:n_catchments, each = n_model_instances)
data_raw <- weights %>%
  bind_cols(expand_grid(model_class, catchment_id)) %>%
  mutate(
    instance_id = rep(1:n_model_instances, n() / n_model_instances),
    model_id = paste(model_class, instance_id, sep = "_"),
    model_id = factor(model_id, levels = unique(model_id)),
    model_id = as.integer(model_id) # assign a unique id to each model
  ) %>%
  select(model_class, catchment_id, model_id, KGE, NSE, RMSE)

# normalizing NSE on a scale from 0 to 10
data_dense <- data_raw %>%
  mutate(NNSE = 1 / (2 - NSE) * 10) %>%
  dplyr::select(catchment_id,
                model_id,
                NNSE) %>%
  rename(rating = NNSE) %>%
  mutate(record_id = 1:n()) # give each row a unique ID


# combine
data_dense <- data_dense %>%
  group_by(catchment_id) %>%
  summarise(best = max(rating)) %>%
  mutate(case = "Dense data set")


data_sparse <- data_sparse %>%
  group_by(catchment_id) %>%
  summarise(best = max(rating)) %>%
  mutate(case = "Sparse data set (actual)")


data_plot <- data_sparse %>%
  rbind(data_dense) %>%
  rbind(sparse_best)


data_plot <- data_plot %>%
  mutate(case = factor(
    case,
    levels = c(
      "Dense data set",
      "Sparse data set (actual)",
      "Sparse data set (recommended)"
    )
  ))

data_plot %>%
  ggplot(aes(best, color = case, linetype = case)) +
  stat_ecdf(geom = "step")+
  scale_color_aaas()+
  scale_linetype_manual(values = c("solid", "dashed", "dotdash")) +
  labs(y = "Cumulative distribution", 
       x="Model rating",
       color = "",
       linetype = "")+
  theme_bw(base_size = 8)+
  theme(legend.position = c(0.01,0.99),
        legend.justification = c(0,1),
        legend.background = element_rect(fill = "transparent", color = NA))


ggsave(filename = "./data/plot/gof.pdf",   width = 3.5,
       height = 3.25,
       units = "in")
