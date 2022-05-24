require(ncdf4)
require(tidyverse)
require(ncdump)
require(lubridate)
library(qgraph)
library(tripack)

# function ----------------------------------------------------------------


NKGE <- function(KGE){
  1/(2 - KGE)
}


# Data --------------------------------------------------------------------


metadata <- NetCDF("./data/m_01_collie1_1p_1s_modelResults.nc")
data_raw <- nc_open("./data/m_01_collie1_1p_1s_modelResults.nc")


# Explore -----------------------------------------------------------------


ncvar_get(data_raw, "Gauge_ID")
ncvar_get(data_raw, "Calibrated_parameters")
ncvar_get(data_raw, "Objective_function_cal")
ncvar_get(data_raw, "Objective_function_eval")

ncvar_get(data_raw, "Store_warmup")
ncvar_get(data_raw, "Sim_q")
ncvar_get(data_raw, "Sim_ea")
ncvar_get(data_raw, "Sim_flux_noInternalFluxes")

ncvar_get(data_raw, "Sim_store_S1")


# Processing --------------------------------------------------------------


catchment_ids <- ncvar_get(data_raw, "Gauge_ID")
date_index <-
  seq(ymd(from = "1989-01-01"),
      to = ymd("2009-12-31"),
      by = 1)


ncvar_get(data_raw, "Sim_q")



ncvar_get(data_raw, "Calibrated_parameters")

load("./data/similarity_m.RDA")
qgraph(similarity_m, layout='spring', vsize=3, threshold = 0.5)





plot(voronoi.mosaic(similarity_m))



temp <- data_raw

data_raw %>%
  dplyr::select(catchment_id, NNSE)

dat <- similarity_m

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



