if (!require("pacman")) {
  install.packages("pacman")
}

pacman::p_load(
  tidyverse,
  lubridate,
  zeallot,
  sf,
  remotes
)

# USAboundaries and USAboundariesData were removed from CRAN as of 31 May 2023
if (!require("USAboundaries")){
  remotes::install_github("ropensci/USAboundaries")
  library(USAboundaries)
}

if (!require("USAboundariesData")){
  remotes::install_github("ropensci/USAboundariesData")
  library(USAboundariesData)
}

# Data --------------------------------------------------------------------
# Read contiguous US map
contiguous_us <- us_states() %>% filter(!(name %in% c("Alaska", "Hawaii", "Puerto Rico"))) %>% st_union()
# Read non-US catchement ID
not_us_catchments <- st_read("./data/Caravan/caravan_not_us_catchments.shp") %>% pull(gauge_id)

# filter out catchments contained in CAMELS
collection_names <- dir("./data/Caravan/timeseries/csv/")
ts_filename <- lapply(collection_names, function(x) paste0("./data/Caravan/timeseries/csv/", x)) %>%
  lapply(dir)

# give name to each catchment, i.e., "data collection"_"catchment name"
all_catchments <- lapply(collection_names, function(x) paste0("./data/Caravan/timeseries/csv/", x)) %>%
  lapply(dir) %>%
  unlist() %>%
  str_sub(1, -5)

(not_us_catchments %in% all_catchments) %>% sum() - length(not_us_catchments) # names of the shape file matches


data_ts <- tibble(
  file_name = paste0(not_us_catchments, ".csv"), # read only catchment not in the us
  collection_name = str_extract(file_name, ".*(?=[_])"),
  catchment_id = str_extract(file_name, "(?<=[_]).*(?=[\\.])"),
  path = paste0(
    "./data/Caravan/timeseries/csv/",
    collection_name,
    '/',
    file_name),
  data = vector("list", 1)
) 

data_ts %>% count(collection_name) # num of catchments per collection

data_ts %>% count(catchment_id) %>% arrange(desc(n)) # catchment 208009 and 27001 have the same name

data_ts %>%
  filter(catchment_id %in% c("208009", "27001")) # not the same catchments, just the same name

# read time series data
data_ts <- data_ts %>%
  mutate(data = map(path, read_csv, show_col_types = FALSE))

data_ts <- data_ts %>%
  select(-file_name)


# Filtering meteorological forcing  ---------------------------------------

# select and rename the variables
data_raw <- data_ts %>%
  mutate(data = purrr::map(data, function(x)
    x %>% dplyr::select(
      date,
      P = total_precipitation_sum,
      T = temperature_2m_mean,
      PET = potential_evaporation_sum,
      Q = streamflow
    )))

data_process <- data_raw  %>%
  mutate(
    catchment_id = paste(collection_name, catchment_id, sep = "-")
  ) %>%
  select(catchment_id, data)

# many missing data, some catchment only have 365 records
n_complete_record <- data_process %>%
  mutate(
    missingness = purrr::map_dbl(
      data, function(x)
        x %>% complete.cases() %>% sum()
    )
  ) %>%
  pull(missingness)

# all the catchment's records are of the same length = 14609
data_process %>%
  mutate(
    record_length = purrr::map_dbl(
      data, function(x)
        nrow(x)
    )
  ) %>%
  pull(record_length) %>%
  unique()

# the start date is different, some from 1981-01-01, some from 1981-01-02
data_process %>%
  mutate(
    start_date = purrr::map(
      data, function(x)
        x$date[[1]]
    )
  ) %>%
  unnest(start_date) %>%
  pull(start_date) %>%
  table()

# the end date is different, some 2020-12-30, some from 2020-12-31
data_process %>%
  mutate(
    start_date = purrr::map(
      data, function(x)
        last(x$date)
    )
  ) %>%
  unnest(start_date) %>%
  pull(start_date) %>%
  table()


# Splitting data ----------------------------------------------------------

rm(data_ts)
gc()

data_process <- data_process %>%
  unnest(data)

# use record from 1981-01-02 to 2020-12-30, as the data has different starting and ending dates
data_process <- data_process %>%
  filter(date > ymd("1981-01-01"),
         date < ymd("2020-12-31"))

# data from 1981-01-01 to 2010-12-31 is interested
# all the forcing data is available, some of the flow data is missing
# catchments with missing Q records is stored in `incomplete_catchments`

minimal_required_Q_length = 365*10 # at least 10 years of data is needed

incomplete_catchments <- data_process %>%
  filter(date >= ymd("1989-01-01"),
         date <= ymd("2010-12-31")) %>%
  group_by(catchment_id) %>%
  summarise(data = list(tibble(Q))) %>%
  mutate(
    n_complete_record = map_dbl(
      data, function(x) complete.cases(x) %>% sum()
    )
  ) %>%
  filter(n_complete_record < minimal_required_Q_length) %>%
  pull(catchment_id)


data_process %>%
  filter(!(catchment_id %in% incomplete_catchments)) %>% pull(catchment_id) %>% unique() %>% length() # 2116 catchments

data_process <- data_process %>%
  filter(!(catchment_id %in% incomplete_catchments)) %>%
  filter(date >= ymd("1989-01-01"),
         date <= ymd("2010-12-31"))

# 2116 catchments left

# save data ---------------------------------------------------------------

data_process %>%
  arrange(catchment_id, date) %>%
  write_csv(file = "./data/processed_caravan_data.csv")
