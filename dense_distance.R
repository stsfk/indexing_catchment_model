if (!require("pacman")) {
  install.packages("pacman")
}

pacman::p_load(tidyverse,
               zeallot,
               data.table,
               doParallel)

cl <- makeCluster(detectCores()-2)
registerDoParallel(cl)

# data --------------------------------------------------------------------
load("./data/dens_factorization.Rda")

# functions ---------------------------------------------------------------

load_P_Q <- function(i){
  # This function is to load P and Q from factorization matrix
  P <- eval_grid$out[[i]]$P
  Q <- eval_grid$out[[i]]$Q
  
  list(P = P, Q = Q)
}

compute_distance <- function(x, y, LB, UB, n_samples = 500000){
  input_dim <- length(x)
  random_numbers <- runif(n = n_samples*input_dim)
  random_matrix <- LB + (UB - LB) * matrix(random_numbers, nrow = input_dim)
  
  mean(abs(colSums((y-x)*random_matrix)))
}

# compute distance metrics ------------------------------------------------

dists <- vector("list", nrow(eval_grid))

for (i in 1:10){
  c(P, Q) %<-% load_P_Q(i)
  
  LB <- apply(Q, 2, min)
  UB <- apply(Q, 2, max)
  
  n_catchments <- 533 
  x_y_pairs <- expand_grid(x=1:n_catchments, y=1:n_catchments) %>%
    filter(x < y)
  
  compute_distance_wrapper <- function(i){
    compute_distance(
      x = P[x_y_pairs$x[[i]],],
      y = P[x_y_pairs$y[[i]],],
      LB,
      UB
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
  
  dists[[i]] <- dist_m
  
  save(dists, file = "./data/dense_dist.Rda")
}

# Stop --------------------------------------------------------------------

stopCluster(cl)

