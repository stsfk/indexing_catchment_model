# The same random model instances were used in the computation of all catchment pairs.

if (!require("pacman")) {
  install.packages("pacman")
}

pacman::p_load(tidyverse,
               zeallot,
               data.table,
               doParallel,
               Rfast,
               Rfast2)

cl <- makeCluster(8)
# registerDoParallel(cl)
# 
# clusterEvalQ(cl, library(Rfast))
# clusterEvalQ(cl, library(Rfast2))

# data --------------------------------------------------------------------
load("./data/dens_factorization.Rda")

# functions ---------------------------------------------------------------

load_P_Q <- function(i){
  # This function is to load P and Q from factorization matrix
  P <- eval_grid$out[[i]]$P
  Q <- eval_grid$out[[i]]$Q
  
  list(P = P, Q = Q)
}

# compute distance metrics ------------------------------------------------

dists <- vector("list", nrow(eval_grid))
random_Ms <- vector("list", nrow(eval_grid))

for (i in 1:10){
  c(P, Q) %<-% load_P_Q(i)
  
  n_samples <- 500000
  input_dim <- dim(P)[2]
  
  LB <- apply(Q, 2, min)
  UB <- apply(Q, 2, max)
  
  n_catchments <- 533 
  
  random_M <- matrix(Rfast2::Runif(n = n_samples*input_dim), ncol = input_dim)
  random_Ms[[i]] <- random_M
  
  R <- P %*% t(random_M)
  
  x_y_pairs <- expand_grid(x=1:n_catchments, y=1:n_catchments) %>%
    filter(x < y)
  
  compute_distance <- function(i){
    
    x = x_y_pairs$x[[i]]
    y = x_y_pairs$y[[i]]
    
    mean(abs(R[x,]-R[y,]))
  }
  
  clusterExport(cl, "compute_distance")
  clusterExport(cl, "x_y_pairs")
  clusterExport(cl, "R")
  
  dist <-
    parLapply(cl, 1:nrow(x_y_pairs), function(x)
      compute_distance(x)) %>% unlist()
  
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
  gc()
}

save(random_Ms, dists, file = "./data/dense_dist_fixed.Rda")

# Stop --------------------------------------------------------------------

stopCluster(cl)


# Result analysis ---------------------------------------------------------

load("./data/dense_dist_fixed.Rda")

dists <- lapply(dists, function(x) x[lower.tri(x)])

out <- matrix(nrow = 10, ncol = 10)

x_y_pairs <- expand_grid(x=1:10, y=1:10) %>%
  filter(x < y)

compute_distance_wrapper <- function(i){
  a <- dists[[x_y_pairs$x[i]]]
  b <- dists[[x_y_pairs$y[i]]]
  
  cor(a,b)
}

dist <- foreach(i=1:nrow(x_y_pairs)) %do% 
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


dist_m %>%
  view()


# recycle -----------------------------------------------------------------


# the various distance functions below are slower or much slower.

compute_distance <- function(x, y, LB, UB, n_samples = 500000){
  # This function computes the expected distance between x and y, where each dimension is scaled by a random vector
  # the Manhattan distance is computed
  
  # dimension of vector x
  input_dim <- length(x)
  
  # generate uniformly distributed random numbers of size n_samples*input_dim
  random_numbers <- runif(n = n_samples*input_dim)
  
  # arrange the random numbers into a matrix of shape = input_dim * n_samples
  random_matrix <- LB + (UB - LB) * matrix(random_numbers, nrow = input_dim)
  
  mean(abs(colSums((y-x)*random_matrix)))
}

compute_distance2 <- function(x, y, LB, UB, n_samples = 500000){
  input_dim <- length(x)
  
  M <- matrix(Rfast2::Runif(n = n_samples*input_dim), ncol = input_dim)
  M <- Rfast::eachrow(M, (UB - LB), "*")
  M <- Rfast::eachrow(M, LB, "+")
  M <- Rfast::eachrow(M, y - x, "*")
  
  mean(abs(Rfast::rowsums(M)))
}

compute_distance3 <- function(x, y, LB, UB, n_samples = 500000){
  input_dim <- length(x)
  
  M <- matrix(Rfast2::Runif(n = n_samples*input_dim), ncol = input_dim)
  M <- Rfast::eachrow(M, (UB - LB)*(x - y), "*")
  M <- Rfast::eachrow(M, (x - y)*LB, "+")
  
  mean(abs(Rfast::rowsums(M)))
}

compute_distance4 <- function(x, y, LB, UB, n_samples = 500000){
  input_dim <- length(x)
  
  v <- Rfast2::Runif(n = n_samples*input_dim)
  v <- v*(UB - LB)*(x - y) + (x - y)*LB
  
  mean(abs(Rfast::rowsums(matrix(v, ncol = input_dim, byrow = T))))
}

xs <- list(x1 = x, x2 = x, x3 = x, x4 = x, x5 = x)
ys <- list(y1 = y, y2 = y, y3 = y, y4 = y, y5 = y)

compute_distance_list <- function(xs, ys, LB, UB, n_samples = 500000){
  input_dim <- length(xs[[1]])
  n_pairs <- length(xs)
  
  x_long <- unlist(xs)
  y_long <- unlist(ys)
  
  LB_long <- rep(LB, n_pairs)
  UB_long <- rep(UB, n_pairs)
  
  M <- matrix(Rfast2::Runif(n = n_samples*input_dim*n_pairs), ncol = input_dim*n_pairs)
  M <- Rfast::eachrow(M, (UB_long - LB_long)*(x_long - y_long), "*")
  M <- Rfast::eachrow(M, (x_long - y_long)*LB_long, "+")
  
  out <- rep(0, n_pairs)
  for (i in 1:n_pairs){
    out[[i]] <- mean(abs(Rfast::rowsums(M[,(i-1)*input_dim+1:input_dim])))
  }
  out
}

Rfast2::benchmark(compute_distance(x,y,LB,UB), times = 100)

abs(Rfast::rowsums(M))

v <- matrix(unlist(M), ncol = input_dim)

M2 <- matrix(c(1:4),byrow = T, nrow = 2)
matrix(unlist(M2), nrow = 4)

c(P,Q) %<-% load_P_Q(1)
x <- P[1,]
y <- P[2,]

LB <- apply(Q, 2, min)
UB <- apply(Q, 2, max)

Rfast2::benchmark(compute_distance(x,y,LB,UB), times = 100)

