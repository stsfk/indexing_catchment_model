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
               fpc,
               pvclust,
               mclust)


# data --------------------------------------------------------------------

load("./data/mf_results.Rda")
i <- 1

P <- Ps[[i]]

# Clustering --------------------------------------------------------------

# k-means
mydata <- P
wss <- (nrow(mydata) - 1) * sum(apply(mydata, 2, var))
for (i in 2:15)
  wss[i] <- sum(kmeans(mydata,
                       centers = i)$withinss)

plot(1:15, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")


# K-Means Cluster Analysis
fit <- kmeans(mydata, 5) # 5 cluster solution
# get cluster means 
aggregate(mydata,by=list(fit$cluster),FUN=mean)
# append cluster assignment
mydata <- data.frame(mydata, fit$cluster)

pamk(mydata)


# Ward Hierarchical Clustering
d <- dist(mydata, method = "euclidean") # distance matrix
fit <- hclust(d, method="ward") 
plot(fit) # display dendogram
groups <- cutree(fit, k=5) # cut tree into 5 clusters
# draw dendogram with red borders around the 5 clusters 
rect.hclust(fit, k=5, border="red")


fit <- pvclust(t(mydata), method.hclust="ward",
               method.dist="euclidean")
plot(fit) # dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fit, alpha=.95)


# model based
fit <- Mclust(mydata)
plot(fit) # plot results 
summary(fit) # display the best model


# K-Means Clustering with 5 clusters
fit <- kmeans(mydata, 5)

# Cluster Plot against 1st 2 principal components

# vary parameters for most readable graph
library(cluster) 
clusplot(mydata, fit$cluster, color=TRUE, shade=TRUE, 
         labels=2, lines=0)

# Centroid Plot against 1st 2 discriminant functions
library(fpc)
plotcluster(mydata, fit$cluster)








df <- P[-1,] %>%
  as.data.frame()
df <- scale(df)

km <- kmeans(df, centers = 5, nstart = 25)

# Plot --------------------------------------------------------------------

# catchment locations
camels_topo <-
  read_delim(
    "./data/CAMELS_US/camels_attributes_v2.0/camels_attributes_v2.0/camels_topo.txt",
    delim = ";"
  ) %>%
  filter(gauge_id %in% catchment_ids)%>%
  select(gauge_id, gauge_lat, gauge_lon)

# catchment_id
data_camels <- read_csv("./data/CAMELS_US/CAMELS_US.csv")
catchment_ids <- data_camels$catchment_id %>% unique()

# geophysio

sf_use_s2(FALSE)

geophysio <- st_read("./data/physio_shp/physio.shp") %>%
  select(PROVINCE) %>% # PROVINCE , DIVISION
  rename(DIVISION = PROVINCE)

catchment_points <- tibble(id = catchment_ids) %>%
  left_join(camels_topo, by = c("id" = "gauge_id")) %>%
  mutate(group = factor(km$cluster)) %>%
  st_as_sf(coords = c("gauge_lon","gauge_lat"), remove = T) %>%
  select(-id)

st_crs(catchment_points) <- st_crs(geophysio)

st_write(catchment_points,dsn = "./data/catchment_points/points.shp")

plot(geophysio, border = 'grey', axes = T)
plot(catchment_points["group"], add = TRUE)


# analysis ----------------------------------------------------------------

results <- catchment_points %>%
  st_join(geophysio)%>%
  group_by(DIVISION) %>%
  count(group)%>%
  st_drop_geometry() %>%
  group_by(DIVISION) %>%
  mutate(perc = round(n/sum(n)*100, 2)) %>%
  slice_max(perc , n = 2) %>%
  ungroup(DIVISION) %>%
  arrange(DIVISION, desc(perc)) 

results

catchment_points %>%
  st_join(geophysio)%>%
  group_by(DIVISION) %>%
  count(group)%>%
  st_drop_geometry() %>%
  group_by(DIVISION) %>%
  mutate(perc = round(n/sum(n)*100, 2)) %>%
  slice_max(perc , n = 1) %>%
  ungroup(DIVISION) %>%
  arrange(DIVISION, desc(perc)) 


# Applications ------------------------------------------------------------

# model database

data_process2 <- data_process %>%
  filter(item_index <= 3900)

data_process3 <- data_process %>%
  filter(item_index > 3900) # new models


train_ind <- sample(1:nrow(data_process2), round(0.5*nrow(data_process2)))
test_ind <- setdiff(1:nrow(data_process2), train_ind)

dtrain <- data_process2[train_ind,]
dtest <- data_process2[test_ind,]

r = Reco()

train_set <- data_memory(user_index = dtrain$user_index,
                         item_index = dtrain$item_index,
                         rating = dtrain$rating)

opts = r$tune(train_set, opts = list(dim = c(10), lrate = c(0.001, 0.005, 0.1, 0.2),
                                     nthread = 10, niter = 20)) # 1:10

r$train(train_set, opts = c(opts$min, nthread = 10, niter = 20))

c(P,Q) %<-% r$output(out_memory(), out_memory())

dim(P)
dim(Q)

test_set <- data_memory(user_index = dtest$user_index,
                        item_index = dtest$item_index,
                        rating = dtest$rating)
pred_rvec <- r$predict(test_set)

gof(pred_rvec, dtest$rating)


P_df <- P[-1,] %>%
  as.data.frame() %>%
  as_tibble()%>%
  mutate(catchment_id = 1:n())

item_index_of_interest <- 3948
eval_grid <- data_process3 %>%
  filter(item_index == item_index_of_interest) %>%
  sample_n(10) %>%
  select(-item_index) %>%
  arrange(user_index)

A <- P_df %>%
  filter(catchment_id %in% eval_grid$user_index) %>%
  select(-catchment_id) %>%
  as.matrix()

b <- eval_grid$rating
showEqn(A, b)

Q_of_interest <- solve(A, b)


ob <- data_process3 %>%
  filter(item_index == item_index_of_interest) %>%
  select(-item_index)

pred <- P_df %>%
  select(-catchment_id) %>%
  as.matrix() %*% Q_of_interest



## predict parameter values

model_parameters <- read_csv("./data/parameters.csv", col_names = F) %>% unlist() %>% unname()


Q[2:501,] %>%
  as.data.frame() %>%
  as_tibble() %>%
  set_names(paste0("q", 1:10)) %>%
  mutate(model_parameters) %>%
  gather(item, value, -model_parameters)%>%
  ggplot(aes(model_parameters, value)) +
  geom_point()+
  facet_wrap(~item, nrow = 2)+
  labs(x = "theta value",
       y = "q value")






# plot map

MainStates <- map_data("state")

ggplot() + 
  geom_polygon(data=MainStates, aes(x=long, y=lat, group=group),
               color="black", fill="lightblue" )+
  geom_point(data = data_plot, aes(x = gauge_lon, y = gauge_lat, color = group))





# Draw decision boundaries

MainStates$long %>% range()
MainStates$lat %>% range()

cgrid <- expand.grid(x = c(-125:-67), y = c(25:50))

dtrain <- tibble(
  x = data_plot$gauge_lon,
  y = data_plot$gauge_lat,
  class = as.factor(km$cluster)
)

knnModel <- train(class ~., 
                  data = dtrain, 
                  method = 'knn',
                  tuneGrid = expand.grid(k = c(1:11)*2-1))

cgrid$class <- predict(knnModel, newdata = cgrid)
cgrid <- as_tibble(cgrid)

ggplot() +   
  geom_point(data = cgrid, aes(x=x, y=y, colour=as.factor(class)), size = 0.8) +
  geom_contour_filled(data = cgrid, aes(x=x, y=y, z=as.numeric(class)), color = "black", bins = 3, breaks = c(0:6), alpha = 0.3) +
  geom_polygon(data=MainStates, aes(x=long, y=lat, group=group),
               color="black", fill="lightblue", size = 1.5, alpha = 0.4)+
  scale_color_discrete(rainbow(5))+
  scale_fill_discrete(rainbow(5))+
  theme_bw()


ggplot() + 
  geom_polygon(data=MainStates, aes(x=long, y=lat, group=group),
               color="black", fill="lightblue" )+
  geom_point(data = data_plot, aes(x = gauge_lon, y = gauge_lat, color = as.factor(km$cluster)))+
  scale_size(range = c(1,10))+
  labs(color = "group")


# Recycle -----------------------------------------------------------------

compute_cor <- function(catchment_id_selected1, catchment_id_selected2){
  temp1 <- data_raw %>%
    filter(catchment_id %in% catchment_id_selected1) %>%
    pull(NNSE)
  
  temp2 <- data_raw %>%
    filter(catchment_id %in% catchment_id_selected2) %>%
    pull(NNSE)
  
  cor(temp1, temp2)
}


compute_dist <- function(catchment_id_selected1, catchment_id_selected2){
  temp1 <- data_raw %>%
    filter(catchment_id %in% catchment_id_selected1) %>%
    pull(NNSE)
  
  temp2 <- data_raw %>%
    filter(catchment_id %in% catchment_id_selected2) %>%
    pull(NNSE)
  
  abs(temp1 - temp2) %>% 
    mean()
}

out <- combn(1:533, 2) %>%
  t() %>%
  setNames(c("V1","V2")) %>%
  as_tibble() %>%
  mutate(cors = map2_dbl(V1, V2, compute_cor))



ggplot(out, aes(V1,V2, fill = cors))+
  geom_tile()+
  scale_fill_continuous(type = "viridis", values = c(0,1))


temp <- data_raw

data_raw %>%
  dplyr::select(catchment_id, NNSE)

dat <- matrix(data_raw$NNSE, ncol = 100, byrow = T)

dij <- dist(scale(dat, center = TRUE, scale = TRUE))
clust <- hclust(dij, method = "average")


ord <- order(cutree(clust, k = 3))
coph <- cophenetic(clust)


layout(matrix(1:4, ncol = 2))
image(as.matrix(dij)[ord, ord], main = "Original distances")
image(as.matrix(coph)[ord, ord], main = "Cophenetic distances")
image((as.matrix(coph) - as.matrix(dij))[ord, ord], 
      main = "Cophenetic - Original")
plot(coph ~ dij, ylab = "Cophenetic distances", xlab = "Original distances",
     main = "Shepard Plot")
abline(0,1, col = "red")
box()
layout(1)

data_raw %>%
  group_by(model_list, catchment_id) %>%
  summarise(mean_nnse = mean(NNSE)) %>%
  ggplot(aes(mean_nnse, fill = model_list)) +
  geom_histogram(position="dodge")
