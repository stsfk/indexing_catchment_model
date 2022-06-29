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


data_process <- tibble(
  user_index = c(rep(1, 5), rep(2, 3), rep(3, 4), rep(4, 3)),
  item_index = c(1:5, 1, 2, 4, 1, 2, 3, 5, 1:3),
  rating = c(6, 3, 8, 2, 7, 7, 3, 7, 6, 2, 7, 8, 4, 2, 1)
)


train_set <- data_memory(
  user_index = data_process$user_index,
  item_index = data_process$item_index,
  rating = data_process$rating,
  index1 = T
)

r = Reco()
r$train(train_set, opts = list(dim = 3, niter = 1000))

c(P,Q) %<-% r$output(out_memory(), out_memory())

P2 <- round(P,1)
Q2 <- round(Q, 1) %>% t()

save(P,Q, file = "illustration_factorization.Rda")
