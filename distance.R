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
               GA,
               ModelMetrics,
               doParallel)


# data --------------------------------------------------------------------

library(tidyverse)

compute_distance <- function(x, y, LB, UB, n_samples = 100000){
  input_dim <- length(x)
  random_numbers <- runif(n = n_samples*input_dim)
  random_matrix <- LB + (UB - LB) * matrix(random_numbers, nrow = input_dim)
  
  ((y-x)*random_matrix) %>%
    colSums() %>%
    abs() %>%
    mean()
}

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






