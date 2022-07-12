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
               ecodist)

cl <- makeCluster(detectCores()-2)
registerDoParallel(cl)

# data --------------------------------------------------------------------
load("./data/mf_results.Rda")

P <- Ps[[1]]
Q <- Qs[[1]]

# functions ---------------------------------------------------------------

compute_distance <- function(x, y, LB, UB, n_samples = 1000000){
  input_dim <- length(x)
  random_numbers <- runif(n = n_samples*input_dim)
  random_matrix <- LB + (UB - LB) * matrix(random_numbers, nrow = input_dim)
  
  mean(abs(colSums((y-x)*random_matrix)))
}

# computation -------------------------------------------------------------

LB <- apply(Q, 2, min)*1.1
UB <- apply(Q, 2, max)*1.1

eval_grid <- expand_grid(x=1:nrow(P), y=1:nrow(P)) # nrow(P)

n_catchments <- 50 
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

nmdsresults <- nmds(as.dist(dist_m),mindim=2,maxdim=2)

plot(nmdsresults)
nmdsresults$conf[[1]] %>% plot()
nmdsresults$conf[[2]] %>% plot()




stopCluster(cl)




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



