if (!require("pacman")) {
  install.packages("pacman")
}

pacman::p_load(tidyverse,
               lubridate,
               vegan,
               recosystem,
               hydroGOF,
               caret,
               tidymodels,
               data.table,
               doParallel,
               sf,
               cowplot,
               stringi)

cl <- makeCluster(detectCores()-2)
registerDoParallel(cl)
clusterEvalQ(cl, library(vegan))

# data --------------------------------------------------------------------
load("./data/mf_results.Rda")

P <- Ps[[1]]
Q <- Qs[[1]]

# geophysio
geophysio <- st_read("./data/physio_shp/physio.shp") %>%
  select(DIVISION) # PROVINCE or DIVISION

# catchment locations
data_camels <- read_csv("./data/CAMELS_US.csv")
catchment_ids <- data_camels$catchment_id %>% unique()

camels_topo <-
  read_delim(
    "./data/camels_topo.txt",
    delim = ";"
  ) %>%
  filter(gauge_id %in% catchment_ids)%>%
  select(gauge_id, gauge_lat, gauge_lon)

catchment_points <- tibble(id = catchment_ids) %>%
  left_join(camels_topo, by = c("id" = "gauge_id")) %>%
  st_as_sf(coords = c("gauge_lon","gauge_lat"), remove = T)

st_crs(catchment_points) <- st_crs(geophysio)

sf::sf_use_s2(FALSE)
catchment_points <- catchment_points %>%
  st_join(geophysio, join = st_intersects)

# functions ---------------------------------------------------------------

compute_distance <- function(x, y, LB, UB, n_samples = 1000000){
  input_dim <- length(x)
  random_numbers <- runif(n = n_samples*input_dim)
  random_matrix <- LB + (UB - LB) * matrix(random_numbers, nrow = input_dim)
  
  mean(abs(colSums((y-x)*random_matrix)))
}

# compute distance metrics ------------------------------------------------

LB <- apply(Q, 2, min)*1.1
UB <- apply(Q, 2, max)*1.1

eval_grid <- expand_grid(x=1:nrow(P), y=1:nrow(P)) # nrow(P)

n_catchments <- 533 
eval_grid <- expand_grid(x=1:n_catchments, y=1:n_catchments) %>%
  filter(x < y)

compute_distance_wrapper <- function(i){
  compute_distance(
    x = P[eval_grid$x[[i]],],
    y = P[eval_grid$y[[i]],],
    LB,
    UB
  )
}

dist <- foreach(i=1:nrow(eval_grid)) %dopar% 
  compute_distance_wrapper(i) %>%
  unlist()

eval_grid <- eval_grid %>%
  mutate(dist = dist)

eval_grid2 <- tibble(
  x = eval_grid$y,
  y = eval_grid$x,
  dist = eval_grid$dist
)

eval_grid <- eval_grid %>%
  rbind(eval_grid2)

dist_m <- as.matrix(dcast.data.table(
  data = as.data.table(eval_grid),
  x ~ y,
  value.var = "dist",
  fill = 0
)[, -1, with = FALSE])


# NMDS --------------------------------------------------------------------

NMDS_result <-
  metaMDS(
    dist_m,
    k = 2,
    try = 20,
    trymax = 5000,
    autotransform = F,
    noshare = F,
    wascores = F,
    parallel = detectCores()-2
  )

stressplot(NMDS_result)
plot(NMDS_result)


# Plotting ----------------------------------------------------------------

data_plot <- tibble(
  x = NMDS_result$points[,1],
  y = NMDS_result$points[,2],
  DIVISION = stri_trans_totitle(catchment_points$DIVISION)
)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

p1 <- ggplot(data_plot, aes(x,y, color = DIVISION))+
  geom_point(size = 0.9)+
  scale_color_manual(values=cbPalette)+
  theme_bw()+
  labs(x = "NMDS Axis 1",
       y = "NMDS Axis 2",
       color = "Physiographic regions")+
  theme_bw(9)+
  theme(
    legend.position = "right",
    legend.key.height = unit(0.16, "inches"),
    plot.margin = unit(c(0.1, 0.1, 0.25, 0.1), "inches")
  )

p2 <- ggplot(data_plot, aes(x,y, color = DIVISION))+
  geom_point(size = 0.9)+
  scale_color_manual(values=cbPalette)+
  facet_wrap(~DIVISION, ncol = 4)+
  labs(x = "NMDS Axis 1",
       y = "NMDS Axis 2")+
  theme_bw(base_size = 9)+
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "grey80",color = NA),
    panel.border = element_blank(),
    panel.background = element_rect(fill = "grey95")
  )

plot_grid(
  p1,
  p2,
  labels = c('(a)', '(b)'),
  label_fontface = "plain",
  label_size = 10,
  ncol = 1,
  rel_heights = c(1, 1.4)
)

ggsave(filename = "./data/plot/NMDS.pdf",   width = 7,
       height = 5,
       units = "in")

ggsave(filename = "./data/plot/NMDS.png",   width = 7,
       height = 5,
       units = "in",
       dpi = 400)

# Stop --------------------------------------------------------------------

stopCluster(cl)




# Recycle -----------------------------------------------------------------

nmdsresults$conf[[1]] %>% plot()
nmdsresults$conf[[2]] %>% plot()

example_NMDS=monoMDS(dist_m, k=2)
example_NMDS=metaMDS(dist_m, k = 2)

stressplot(example_NMDS)
plot(example_NMDS)

ordiplot(example_NMDS,type="n")




metaMDSiter(island.spp_distmat,
            distance = "bray",
            k = 3,
            maxit = 999, 
            trymax = 500,
            wascores = TRUE)

fit <- cmdscale(dist_m,eig=TRUE, k=2) # k is the number of dim
fit # view results

# plot solution
x <- fit$points[,1]
y <- fit$points[,2]
plot(x, y)
text(x, y, labels = row.names(mydata), cex=.7)




dist_m <- matrix(nrow = 533, ncol = 533)

x <- c(0,0,0)
y <- c(1,2,3)
LB <- c(-1,-10,-100)
UB <- c(1,10,100)

compute_distance(x,y,LB,UB)


x <- c(0,0)
y <- c(1,1)
LB <- c(-1,-1)
UB <- c(1,1)
compute_distance(x,y,LB,UB)




input_dim <- length(x)
random_numbers <- runif(n = n_samples*input_dim)
random_matrix <- LB + (UB - LB) * matrix(random_numbers, nrow = input_dim)

((y-x)*random_matrix) %>%
  colSums() %>%
  abs() %>%
  mean()

# 

a1 = 1
a2 = 1
n_samples = 10000

mean(abs(a1*runif(n_samples,-1,1) - a2*runif(n_samples,-1,1)))


mean(abs(2*runif(n_samples,-1,1) - 2*runif(n_samples,-1,1)))
mean(abs(1*runif(n_samples,-1,1) - 3*runif(n_samples,-1,1)))


mean(abs(a1*runif(n_samples,-1,1) - a2*runif(n_samples,-1,1)))



# point a: (0,0)
# point b: (1,0)
# point c: (1,1)
# distance ab:
mean(abs(1*runif(n_samples,-1,1) + 0*runif(n_samples,-1,1)))

# distance bc:
mean(abs(0*runif(n_samples,-1,1) + 1*runif(n_samples,-1,1)))

# distance ac:
mean(abs(1*runif(n_samples,-1,1) + 1*runif(n_samples,-1,1)))
# ab + bc > ac






# point a: (0,0)
# point b: (1,1)
# point c: (10,5)
# distance ab:
mean(abs(1*runif(n_samples,-1,1) + 1*runif(n_samples,-1,1)))

# distance bc:
mean(abs(9*runif(n_samples,-1,1) + 4*runif(n_samples,-1,1)))

# distance ac:
mean(abs(10*runif(n_samples,-1,1) + 5*runif(n_samples,-1,1)))
# ab + bc > ac


for (i in 1:nrow(dist_m)){
  for (j in 1:ncol(dist_m)){
    
    if (i!=j){
      dist_m[i,j] <- eval_grid %>% 
        filter(x == min(i,j), y == max(i,j)) %>%
        pull(dist)
    }
  }
}



