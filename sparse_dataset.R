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
               mlrMBO)

# seed --------------------------------------------------------------------

#set.seed(1234)
nthread <- 10 # number of CPU thread

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


# Experiments -------------------------------------------------------------

prepare_modeling_data <-
  function(frac = 1,
           train_portion = 0.6,
           val_portion = 0.2,
           test_portion = 0.2) {
    data_sample <- data_process %>%
      group_by(model_id) %>%
      sample_frac(frac) %>%
      ungroup()
    
    data_train_val <- data_sample %>%
      group_by(model_id) %>%
      sample_frac(train_portion + val_portion) %>%
      ungroup()
    
    data_train <- data_train_val %>%
      group_by(model_id) %>%
      sample_frac(train_portion / (val_portion + train_portion)) %>%
      ungroup()
    
    data_val <- data_train_val %>%
      filter(record_id %in% setdiff(data_train_val$record_id, data_train$record_id))
    
    data_test <- data_sample %>%
      filter(record_id %in% setdiff(data_sample$record_id, data_train_val$record_id))
    
    list(
      data_train = data_train %>% select(-record_id),
      data_val = data_val %>% select(-record_id),
      data_test = data_test %>% select(-record_id),
      record_id_train = data_train$record_id,
      record_id_val = data_val$record_id,
      record_id_test = data_test$record_id
    )
  }


sparse_gof_wrapper <- function(frac = 1) {
  c(data_train,
    data_val,
    data_test,
    record_id_train,
    record_id_val,
    record_id_test) %<-% prepare_modeling_data(frac = frac)
  
  # training
  train_set <- data_memory(
    user_index = data_train$catchment_id,
    item_index = data_train$model_id,
    rating = data_train$rating,
    index1 = T
  )
  
  val_set <- data_memory(
    user_index = data_val$catchment_id,
    item_index = data_val$model_id,
    rating = data_val$rating,
    index1 = T
  )
  
  data_train_val <- data_train %>%
    bind_rows(data_val)
  train_val_set <- data_memory(
    user_index = data_train_val$catchment_id,
    item_index = data_train_val$model_id,
    rating = data_train_val$rating,
    index1 = T
  )
  
  test_set <- data_memory(
    user_index = data_test$catchment_id,
    item_index = data_test$model_id,
    rating = data_test$rating,
    index1 = T
  )
  
  # experiment
  scoringFunction <- function(x) {
    dim <- x["dim"] %>% unlist()
    lrate <- x["lrate"] %>% unlist()
    niter <- x["niter"] %>% unlist()
    
    costp_l1 <- x["costp_l1"] %>% unlist()
    costp_l2 <- x["costp_l2"] %>% unlist()
    costq_l1 <- x["costq_l1"] %>% unlist()
    costq_l2 <- x["costq_l2"] %>% unlist()
    
    r = Reco()
    
    r$train(
      train_set,
      opts = list(
        dim = dim,
        costp_l1 = costp_l1,
        costp_l2 = costp_l2,
        costq_l1 = costq_l1,
        costq_l2 = costq_l2,
        lrate = lrate,
        niter = niter,
        verbose = F,
        nthread = nthread
      )
    )
    
    rmse <-
      ModelMetrics::rmse(actual = data_val$rating,
                         predicted = r$predict(val_set))
    
    return(rmse)
  }
  
  obj_fun <- makeSingleObjectiveFunction(
    fn = scoringFunction,
    par.set = makeParamSet(
      makeIntegerParam("dim", lower = 1, upper = 100),
      makeNumericParam("lrate",  lower = 1e-04,   upper = 0.2),
      makeIntegerParam("niter", lower = 5,  upper = 1000),
      makeNumericParam("costp_l1",  lower = 0,   upper = 0.1),
      makeNumericParam("costp_l2",  lower = 0,   upper = 0.1),
      makeNumericParam("costq_l1",  lower = 0,   upper = 0.1),
      makeNumericParam("costq_l2",  lower = 0,   upper = 0.1)
    ),
    has.simple.signature = FALSE,
    minimize = TRUE
  )
  
  des <- generateDesign(
    n = 5 * getNumberOfParameters(obj_fun),
    par.set = getParamSet(obj_fun),
    fun = lhs::randomLHS
  )
  
  des$y = apply(des, 1, obj_fun)
  
  control <- makeMBOControl() %>%
    setMBOControlTermination(., iters = 100 - 5 * getNumberOfParameters(obj_fun))
  
  run <- mbo(
    fun = obj_fun,
    design = des,
    control = control,
    show.info = TRUE
  )
  
  r = Reco()
  
  opts <- run$x
  opts$nthread <- nthread
  opts$verbose <- F
  r$train(train_val_set, opts = opts)
  pred_rvec <- r$predict(test_set)
  r2 <- cor(data_test$rating, pred_rvec) ^ 2
  rmse <- ModelMetrics::rmse(actual = data_test$rating,
                             predicted = pred_rvec)
  
  c(P,Q) %<-% r$output(out_memory(), out_memory())
  
  out <- list(
    r2 = r2,
    rmse = rmse,
    run = run,
    des = des,
    P = P,
    Q = Q,
    record_id_train = record_id_train,
    record_id_val = record_id_val,
    record_id_test = record_id_test
  )
}

eval_grid <- expand_grid(
  ratio = c(0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.2, 0.5, 0.8, 1),
  r2 = 0,
  rmse = 0,
  out = vector("list",1),
  repeats = c(1)
)


for (i in 1:nrow(eval_grid)){
  frac <- eval_grid$ratio[i]
  eval_grid$out[[i]] <- sparse_gof_wrapper(frac)
  eval_grid$r2[[i]] <- eval_grid$out[[i]]$r2
  eval_grid$rmse[[i]] <- eval_grid$out[[i]]$rmse
  
  gc()
}


save(eval_grid, file = "sparse_exp.Rda")




predict_score <- function(i){
  # i: catchment id
  
  tibble(
    catchment_id = i,
    model_id = 1:nrow(Q),
    pred = P[i,] %*% t(Q) %>% as.vector())
}

predict_scores <- vector("list", n_catchments)
for (i in 1:n_catchments){
  predict_scores[[i]] <- predict_score(i)
}

predict_scores <- bind_rows(predict_scores)

data_plot <- predict_scores %>%
  group_by(catchment_id) %>%
  summarise(best = max(pred))

data_plot <- data_plot %>%
  mutate(nse = (2*best-10)/best)

data_plot$nse %>% hist()

data_plot$nse %>% mean()

# plot --------------------------------------------------------------------

data_plot <- eval_grid %>%
  select(ratio, r2, rmse) %>%
  transmute(ratio = ratio * 100,
            `R-squared` = r2,
            RMSE = rmse) %>%
  gather(item, value, -ratio) 


ggplot(data_plot, aes(ratio, value))+
  geom_point(color = "steelblue")+
  geom_line()+
  facet_wrap(~item, ncol = 1, scales = "free") +
  scale_x_continuous(breaks = c(0:5)*20) +
  labs(x = "Density of data set (percentage of available weights)") +
  theme_bw()+
  theme(axis.title.y = element_blank())

# recycle -----------------------------------------------------------------
model_class  <- c(
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

n_model_classes <- n_model_classes - 1 # errors in the last model class 
