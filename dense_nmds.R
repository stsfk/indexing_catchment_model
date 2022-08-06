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
load("./data/dense_dist.Rda")

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

# NMDS --------------------------------------------------------------------

dist_m <- dists[[1]]

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

save(NMDS_result, file = "NMDS_result.Rda")
# Plotting ----------------------------------------------------------------

data_plot <- tibble(
  x = NMDS_result$points[,1],
  y = NMDS_result$points[,2],
  DIVISION = stri_trans_totitle(catchment_points$DIVISION)
)
#F0E442
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#A6761D", "#0072B2", "#D55E00", "#CC79A7")

p1 <- ggplot(data_plot, aes(x,y, color = DIVISION))+
  geom_point(size = 0.75, shape = 21)+
  scale_color_manual(values=cbPalette)+
  theme_bw()+
  labs(x = "NMDS Axis 1",
       y = "NMDS Axis 2",
       color = "Physiographic regions")+
  theme_bw(10)+
  theme(
    legend.position = "right",
    legend.key.height = unit(0.16, "inches"),
    plot.margin = unit(c(0.1, 0.1, 0.25, 0.1), "inches")
  )

p2 <- ggplot(data_plot, aes(x,y, color = DIVISION))+
  geom_point(size = 0.75, shape = 21)+
  scale_color_manual(values=cbPalette)+
  facet_wrap(~DIVISION, ncol = 4)+
  labs(x = "NMDS Axis 1",
       y = "NMDS Axis 2")+
  theme_bw(base_size = 10)+
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "grey80",color = NA),
    panel.border = element_blank(),
    panel.background = element_rect(fill = "grey96")
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
