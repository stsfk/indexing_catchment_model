if (!require("pacman")) {
  install.packages("pacman")
}
pacman::p_load(
  tidyverse
)


# data --------------------------------------------------------------------

# Dense data set whole
load("./data/dens_factorization.Rda")

data_dens_whole <- eval_grid %>% select(r2, rmse) %>%
  mutate(dims = sapply(eval_grid$out, function(x)
    x$P %>% dim() %>% .[2])) %>%
  mutate(ratio = 1,
         dataset = "Dense")

# Dense data set sub
load("./data/sparse_exp.Rda")

data_dens_sub <- eval_grid %>% select(ratio, r2, rmse) %>%
  mutate(dims = sapply(eval_grid$out, function(x)
    x$P %>% dim() %>% .[2])) %>%
  mutate(dataset = "Dense")

# Sparse data set
load("./data/sparse_factorization.Rda")

data_sparse <- eval_grid %>% select(r2, rmse) %>%
  mutate(dims = sapply(eval_grid$out, function(x)
    x$P %>% dim() %>% .[2])) %>%
  mutate(ratio = 0.05) %>%
  mutate(dataset = "Sparse")

# Intense data set
load("./data/caravan_mf_opt_para.Rda")
load("./data/caravan_mf_eval_gof.Rda")

data_intense <- eval_grid2 %>%
  mutate(dims = sapply(opt_para, function(x) x$dim)) %>%
  mutate(ratio = 0.05) %>%
  mutate(dataset = "Intense")

# combine -----------------------------------------------------------------

data_plot <- data_dens_whole %>%
  bind_rows(data_dens_sub) %>%
  bind_rows(data_sparse) %>%
  bind_rows(data_intense) %>% 
  select(dataset, everything())

save(data_plot, file = "./data/optimal_d.Rda")


# Plot --------------------------------------------------------------------

load("./data/optimal_d.Rda")

data_plot <- data_plot %>%
  mutate(dataset = factor(dataset, levels = c("Dense", "Sparse", "Intense")),
         density = factor(ratio, levels = c(0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.2, 0.4, 0.8, 1), labels =
                            c("1%","2%","3%","4%","5%","10%","20%","40%","80%","100%")))

data_plot2 <- data_plot %>%
  mutate(case = paste(dataset, density, sep="\n")) %>%
  group_by(case) %>%
  mutate(r2_max = max(r2),
         r2_norm = r2/r2_max) %>%
  ungroup() %>%
  mutate(case = replace(case, case == "Sparse\n5%", "Sparse"),
         case = replace(case, case == "Intense\n5%", "Intense"),
         case = factor(case, levels=c(
           "Dense\n1%",
           "Dense\n2%" ,
           "Dense\n3%",
           "Dense\n4%",
           "Dense\n5%",
           "Dense\n10%",
           "Dense\n20%" ,
           "Dense\n40%"  ,
           "Dense\n80%",
           "Dense\n100%",
           "Sparse",
           "Intense"
         )))

ggplot(data_plot2, aes(dims, r2_norm, color = dataset)) +
  geom_point(shape = 1) +
  facet_wrap(~ case, nrow = 2) +
  scale_x_continuous(breaks = c(0,20,40,60,80,100))+
  scale_color_manual(values = c("#377eb8", "#e41a1c", "#4daf4a"))+
  labs(y = "Relative performance among 10 random experiments",
       x = "Optimal number of dimensions of latent factor vectors",
       color = "Data set") +
  theme_bw(base_size = 8) +
  theme(legend.position = "top")

ggsave(filename = "./data/plot/fig_optimal_d.pdf", width = 19, height = 10, units = "cm")


# other plots -------------------------------------------------------------

ggplot(data_plot, aes(dataset, dims, fill = density)) +
  geom_boxplot(position = position_dodge2(preserve = "single"), alpha=0.8) +
  scale_fill_manual(
    values = c(
      "#a50026",
               "#d73027",
               "#f46d43",
               "#fdae61",
               "#fee08b",
               "#d9ef8b",
               "#a6d96a",
               "#66bd63",
               "#1a9850",
               "#006837"
    )
  )+
  labs(x="Data set",
       y = "Optimal dimension of latent factor vectors",
       fill = "Density")+
  theme_bw()
