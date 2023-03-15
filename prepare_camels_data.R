# This code uses the method presented in Gnann, S. (2022, April). Sebastiangnann/camels_matlab: First release. Zenodo. 
# Rerieved from https://doi.org/10.5281/zenodo.6462821 doi: 10.5281/ zenodo.6462821

if (!require("pacman")) {
  install.packages("pacman")
}

pacman::p_load(
  tidyverse,
  lubridate,
  zeallot
)

# Processing basin attributes ---------------------------------------------

catchment_id <- read_delim(
  "./data/CAMELS_US/camels_attributes_v2.0/camels_attributes_v2.0/camels_clim.txt",
  delim = ";"
) %>%
  pull(gauge_id)

camels_topo <-
  read_delim(
    "./data/CAMELS_US/camels_attributes_v2.0/camels_attributes_v2.0/camels_topo.txt",
    delim = ";"
  ) 


# Read Q ------------------------------------------------------------------

read_streamflow <- function(fpath) {
  read_table(
    fpath,
    col_names = c(
      "catchment_id",
      "year",
      "month",
      "day",
      "discharge",
      "quality_flag"
    ),
    col_types = cols(
      catchment_id = col_character(),
      year = col_character(),
      month = col_character(),
      day = col_character(),
      discharge = col_double(),
      quality_flag = col_character()
    )
  ) %>%
    mutate(date = ymd(paste(year, month, day))) %>%
    dplyr::select(catchment_id, date, discharge, quality_flag) %>%
    arrange(date)
}

# read catchment discharge, there are 674 catchments with Q series, however only 671 included in camels 
catchment_discharge <-
  tibble(
    subfolder = dir(
      './data/CAMELS_US/basin_timeseries_v1p2_metForcing_obsFlow/basin_dataset_public_v1p2/usgs_streamflow/'
    ),
    mainfoler = "./data/CAMELS_US/basin_timeseries_v1p2_metForcing_obsFlow/basin_dataset_public_v1p2/usgs_streamflow/"
  ) %>%
  transmute(folder = paste0(mainfoler, subfolder),
            fname = map(folder, dir)) %>%
  unnest(fname) %>%
  transmute(
    catchment_id = map_chr(fname, str_extract, "[0-9]+"),
    fpath = paste0(folder, "/", fname),
    discharge = map(fpath, read_streamflow)
  ) %>%
  dplyr::select(discharge) %>%
  unnest(discharge)

# normalize by catchment area
Q <- catchment_discharge %>%
  left_join(camels_topo %>% select(catchment_id = gauge_id, area_gages2), by = c("catchment_id")) %>%
  mutate(discharge = discharge * 24 * 3600 * 0.0283168466/(area_gages2 * 1000000) * 1000) %>%# CFS to to ft3 to m3 divided by m2 to mm
  select(-area_gages2, Q = discharge)

# P,  T,  PET -------------------------------------------------------------

read_forcing <- function(fpath) {
  read_table(
    fpath,
    col_types = cols(
      YR = col_character(),
      MNTH = col_character(),
      DY = col_character(),
      HR = col_character(),
      SWE = col_double(),
      PRCP = col_double(),
      RAIM = col_double(),
      TAIR = col_double(),
      PET = col_double(),
      ET = col_double(),
      MOD_RUN = col_double(),
      OBS_RUN = col_double()
    )
  ) %>%
    mutate(date = ymd(paste(YR, MNTH, DY))) %>%
    dplyr::select(date, P = PRCP, T = TAIR, PET) %>%
    arrange(date)
}


forcing <- tibble(
  subfolder = dir(
    './data/CAMELS_US/basin_timeseries_v1p2_modelOutput_daymet/model_output/flow_timeseries/daymet/'
  ),
  mainfoler = './data/CAMELS_US/basin_timeseries_v1p2_modelOutput_daymet/model_output/flow_timeseries/daymet/'
) %>%
  transmute(folder = paste0(mainfoler, subfolder),
            fname = map(folder, dir)) %>%
  unnest(fname) %>%
  filter(str_detect(fname, "_05_model_output.txt")) %>%
  transmute(
    catchment_id = map_chr(fname, str_extract, "[0-9]+"),
    fpath = paste0(folder, "/", fname)
  ) %>%
  group_by(catchment_id) %>%
  summarise(fpath = fpath[1]) %>% # some _model_output.txt are repeated
  mutate(forcing = map(fpath, read_forcing)) %>%
  unnest(forcing) %>%
  select(catchment_id, date, P, T, PET)


# Adjust PET --------------------------------------------------------------

read_pet_coef <- function(fname) {
  fname %>%
    read_tsv(
      col_names = c("item", "value"),
      col_types = cols(item = col_character(), value = col_double())
    ) %>%
    filter(item == "pet_coef") %>%
    pull(value)
}


pet_coef_df <- tibble(
  subfolder = dir(
    './data/CAMELS_US/basin_timeseries_v1p2_modelOutput_daymet/model_output/flow_timeseries/daymet/'
  ),
  mainfoler = './data/CAMELS_US/basin_timeseries_v1p2_modelOutput_daymet/model_output/flow_timeseries/daymet/'
) %>%
  transmute(folder = paste0(mainfoler, subfolder),
            fname = map(folder, dir)) %>%
  unnest(fname) %>%
  filter(str_detect(fname, "_05_model_parameters.txt")) %>%
  transmute(
    catchment_id = map_chr(fname, str_extract, "[0-9]+"),
    fpath = paste0(folder, "/", fname),
    pet_coef = map_dbl(fpath, read_pet_coef)
  ) %>%
  select(-fpath)

forcing <- forcing %>%
  left_join(pet_coef_df, by = "catchment_id") %>%
  mutate(PET = 1.26/pet_coef * PET ) %>%  # adjust PET to standard value 1.26
  select(-pet_coef)

# Join data ---------------------------------------------------------------

all_data <- forcing %>%
  left_join(Q, by = c("catchment_id", "date"))


# Save data ---------------------------------------------------------------
save(all_data, file = "./data/CAMELS_US/all_data.Rda")


# Selected catchments and period ------------------------------------------

selected_catchments <-
  read_csv("./data/CAMELS_US/selected_catchments.csv",
           col_types = cols(ID = col_character())) %>%
  pull(ID) %>%
  str_c("00", .) %>%
  str_sub(start = -8)

interval_start <- ymd("1987-10-01")
interval_end <- ymd("2009-09-30")
selected_interval <- interval(interval_start, interval_end)
interval_length <- (interval_end - interval_start) %>% as.double()

data_process <- all_data %>%
  filter(catchment_id %in% selected_catchments,
         date %within% selected_interval)

# keep catchments contains 8035 days of record
selected_catchments <- data_process %>%
  group_by(catchment_id) %>%
  summarise(first_day = min(date),
            last_day = max(date)) %>%
  mutate(length = as.numeric(last_day - first_day)) %>%
  filter(length == interval_length) %>%
  pull(catchment_id)

data_process <- data_process %>%
  filter(catchment_id %in% selected_catchments)

# filter catchment with missing Q record
not_selected_catchments <- data_process %>%
  filter(Q < 0) %>%
  pull(catchment_id) %>%
  unique()

data_process <- data_process %>%
  filter(!(catchment_id %in% not_selected_catchments))

# check regular interval
data_process$date %>%
  diff() %>%
  table()

# check missing data
data_process %>%
  sapply(function(x) sum(is.na(x))/length(x))

# check negative
data_process %>%
  select(P, PET) %>%
  lapply(function(x) sum(x < 0)/length(x))

# Save processed data -----------------------------------------------------
data_process <- data_process %>% 
  select(-quality_flag)
save(data_process, file = "./data/CAMELS_US/data_processed.Rda")


# save to csv -------------------------------------------------------------

write_csv(data_process, file = "./data/CAMELS_US/CAMELS_US.csv")

# Basic info --------------------------------------------------------------

data_process %>% count(catchment_id)
