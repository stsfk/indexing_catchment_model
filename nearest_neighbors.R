if (!require("pacman")) {
  install.packages("pacman")
}

pacman::p_load(tidyverse,
               lubridate,
               zeallot,
               recosystem,
               vegan,
               hydroGOF,
               caret,
               tidymodels,
               RANN,
               sf,
               sfheaders,
               data.table)

# data --------------------------------------------------------------------

# load latent factor vector of catchments
load("./data/mf_results.Rda")
P <- Ps[[1]]

# catchment locations
data_camels <- read_csv("./data/CAMELS_US.csv")
catchment_ids <- data_camels$catchment_id %>% unique()

camels_topo <-
  read_delim(
    "./data/camels_topo.txt",
    delim = ";"
  ) %>%
  filter(gauge_id %in% catchment_ids)%>%
  select(gauge_id, gauge_lat, gauge_lon)%>%
  mutate(catchment_id = 1:n())


# Uniform neighbor --------------------------------------------------------

load("./data/dist_m.Rda")

# find the nearest neighbors for each catchment
n_catchment <- nrow(P) # number of catchments
data_process<- tibble(
  catchment = 1:n_catchment,
  n_neighbor = apply(dist_m, 1, function(x) which(rank(x)==2)) # nearest neighbor except itself
)

# keep only distinct nearest catchment pairs
data_process <- data_process %>%
  transmute(catchment_a = map2_int(n_neighbor, catchment, min),  # smaller value of the catchment IDs in a pair
         catchment_b = map2_int(n_neighbor, catchment, max)) %>% # larger value of the catchment IDs in a pair
  distinct() # keep only distinct pairs

# find the coordinate of the catchments
df <- tibble(
  id = 1:nrow(data_process),
  lon1 = 0,
  lat1= 0,
  lon2 = 0,
  lat2 = 0
)

# https://stackoverflow.com/a/51922422/3361298; Create linestring from two points in same row in dataframe
for (i in 1:nrow(data_process)){
  catchment_a <- data_process$catchment_a[i]
  c(df$lat1[i], df$lon1[i]) %<-% (camels_topo %>% filter(catchment_id == catchment_a) %>% select(gauge_lat, gauge_lon) %>% unlist())
  
  catchment_b <- data_process$catchment_b[i]
  c(df$lat2[i], df$lon2[i]) %<-% (camels_topo %>% filter(catchment_id == catchment_b) %>% select(gauge_lat, gauge_lon) %>% unlist())
}

dt <- as.data.table(df)

## To use `sfheaders` the data needs to be in long form
dt1 <- dt[, .(id, lon = lon1, lat = lat1)]
dt2 <- dt[, .(id, lon = lon2, lat = lat2)]

## Add on a 'sequence' variable so we know which one comes first
dt1[, seq := 1L ]
dt2[, seq := 2L ]

## put back together
dt <- rbindlist(list(dt1, dt2), use.names = TRUE)
setorder(dt, id, seq)

sf <- sfheaders::sf_linestring(
  obj = dt
  , x = "lon"
  , y = "lat"
  , linestring_id = "id"
) 

geophysio <- st_read("./data/physio_shp/physio.shp") %>%
  select(PROVINCE) %>% # PROVINCE , DIVISION
  rename(DIVISION = PROVINCE)
st_crs(sf) <- st_crs(geophysio)
st_write(sf,dsn = "./data/catchment_links/links_uniform.shp")







# euclidean distance ------------------------------------------------------

# searching the nearest neighbors; k=2 is to exclude itself
c(nn.idx, nn.dists) %<-% nn2(data = P, query = P, k = 2)

# creating a tibble storing the "nearest" catchment pairs and their distances
data_process <- tibble(
  catchment_a = 1:nrow(P),
  catchment_b = nn.idx[,2],
  distance = nn.dists[,2]
) %>%
  group_by(distance) %>%
  mutate(ind = 1:n())%>%
  ungroup() %>%
  filter(ind == 1)
  
# creating lines connecting neighboring catchments
camels_topo <- camels_topo %>%
  mutate(catchment_id = 1:n())

df <- tibble(
  id = 1:nrow(data_process),
  lon1 = 0,
  lat1= 0,
  lon2 = 0,
  lat2 = 0
)

# https://stackoverflow.com/a/51922422/3361298
for (i in 1:nrow(data_process)){
  catchment_a <- data_process$catchment_a[i]
  c(df$lat1[i], df$lon1[i]) %<-% (camels_topo %>% filter(catchment_id == catchment_a) %>% select(gauge_lat, gauge_lon) %>% unlist())
  
  catchment_b <- data_process$catchment_b[i]
  c(df$lat2[i], df$lon2[i]) %<-% (camels_topo %>% filter(catchment_id == catchment_b) %>% select(gauge_lat, gauge_lon) %>% unlist())
}


dt <- as.data.table(df)

## To use `sfheaders` the data needs to be in long form

dt1 <- dt[, .(id, lon = lon1, lat = lat1)]
dt2 <- dt[, .(id, lon = lon2, lat = lat2)]

## Add on a 'sequence' variable so we know which one comes first
dt1[, seq := 1L ]
dt2[, seq := 2L ]

## put back together
dt <- rbindlist(list(dt1, dt2), use.names = TRUE)
setorder(dt, id, seq)

sf <- sfheaders::sf_linestring(
  obj = dt
  , x = "lon"
  , y = "lat"
  , linestring_id = "id"
) 

geophysio <- st_read("./data/physio_shp/physio.shp") %>%
  select(PROVINCE) %>% # PROVINCE , DIVISION
  rename(DIVISION = PROVINCE)
st_crs(sf) <- st_crs(geophysio)
st_write(sf,dsn = "./data/catchment_links/links.shp")

# plot --------------------------------------------------------------------







