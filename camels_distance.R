if (!require("pacman")) {
  install.packages("pacman")
}

pacman::p_load(tidyverse,
               zeallot,
               data.table,
               doParallel,
               Rfast,
               Rfast2,
               sf)

cl <- makeCluster(detectCores()-2)
registerDoParallel(cl)

clusterEvalQ(cl, library(Rfast))
clusterEvalQ(cl, library(Rfast2))

# data --------------------------------------------------------------------

P <- read.table("./data/camels_mat_P_1.txt") %>% data.matrix()
Q <- read.table("./data/camels_mat_Q_1.txt") %>% data.matrix()

# functions ---------------------------------------------------------------

compute_distance <- function(x, y, LB, UB, n_samples){
  # This function computes the expected distance between x and y, where each dimension is scaled by a random vector
  # the Manhattan distance is computed
  
  # dimension of vector x
  input_dim <- length(x)
  range <- UB - LB
  variation <- x - y
  
  # generate uniformly distributed random numbers of size n_samples*input_dim
  M <- matrix(Rfast2::Runif(n = n_samples*input_dim), ncol = input_dim)
  # scale the random matrix
  M <- Rfast::eachrow(M, range*variation, "*")
  # adding a base
  M <- Rfast::eachrow(M, variation*LB, "+")
  
  # compute the sum of each row to get sample distance, mean and abs to get the average distance
  mean(abs(Rfast::rowsums(M)))
}

# compute distance metrics ------------------------------------------------

LB <- apply(Q, 2, min)
UB <- apply(Q, 2, max)

n_catchments <- 533*2 
x_y_pairs <- expand_grid(x=1:n_catchments, y=1:n_catchments) %>%
  filter(x < y)

compute_distance_wrapper <- function(i, n_samples=500000){
  compute_distance(
    x = P[x_y_pairs$x[[i]],],
    y = P[x_y_pairs$y[[i]],],
    LB,
    UB,
    n_samples = n_samples
  )
}

dist <- foreach(i=1:nrow(x_y_pairs)) %dopar% 
  compute_distance_wrapper(i) %>%
  unlist()

x_y_pairs <- x_y_pairs %>%
  mutate(dist = dist)

x_y_pairs2 <- tibble(
  x = x_y_pairs$y,
  y = x_y_pairs$x,
  dist = x_y_pairs$dist
)

x_y_pairs <- x_y_pairs %>%
  rbind(x_y_pairs2)

dist_m <- as.matrix(dcast.data.table(
  data = as.data.table(x_y_pairs),
  x ~ y,
  value.var = "dist",
  fill = 0
)[, -1, with = FALSE])

save(dist_m, file = "./data/camels_random_start_dist_m.Rda")

# Stop --------------------------------------------------------------------

stopCluster(cl)


# Result analysis ---------------------------------------------------------

# data
load(file = "./data/gr4j_camels.Rda")

results <- results %>%
  mutate(old_catchment_id = catchment_id,
         catchment_id = paste(catchment_id, ver, sep = "_"))

catchment_ids <- results$catchment_id %>% unique()

load("./data/camels_random_start_dist_m.Rda")


catchment_points <- st_read("./data/catchment_points/points.shp")
geophysio <- st_read("./data/physio_shp/physio.shp") %>%
  st_make_valid(geophysio) %>%
  select(DIVISION)

# catchment_points <- catchment_points %>% st_transform(st_crs(geophysio))
# 
# polygon_ids <- catchment_points %>% st_within(geophysio) %>% unlist()
# polygon_DIVISION <- geophysio %>% st_drop_geometry() %>% pull(DIVISION)
# catchment_points <- catchment_points %>% 
#   mutate(DIVISION = polygon_DIVISION[polygon_ids])


# process
data_process <- tibble(
  old_catchment_id = str_sub(catchment_ids, 1,-3),
  catchment_id = 1:1066,
  ranks = vector("list", 1)
)

for (i in 1:nrow(dist_m)){
  data_process$ranks[[i]] <- tibble(
    catchment_id = 1:1066,
    ranks= (dist_m[i,] %>% rank()) - 1)
}

data_process <- data_process %>%
  mutate(
    pair_catchment = floor((catchment_id-1)/2)*2 + c(2,1),
    pair_rank = map2_dbl(ranks, pair_catchment, function(x,y) x$ranks[[y]])
  ) 

data_plot <- data_process %>%
  group_by(old_catchment_id) %>%
  summarise(rank = mean(pair_rank)) %>%
  rename(gauge_id = old_catchment_id) %>%
  mutate(rank_interval = cut(rank, breaks = c(0,1,5,10,100,500,1066), include.lowest=F) %>% as.character(),
         rank_interval = replace(rank_interval, rank_interval == "(500,1.07e+03]", "(500,1065]"))

catchment_points_with_rank <- catchment_points %>%
  left_join(data_plot, by = "gauge_id")

data_plot %>% count(rank_interval)

st_write(catchment_points_with_rank, dsn="./data/catchment_points/points_with_rank.shp",delete_dsn = T)

