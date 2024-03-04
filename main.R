


# 1. Setup -------------------------------------------------------------------

# renv::activate()
# renv::snapshot()
renv::restore()
# renv::install()

box::purge_cache()
box::use(
  R/pull,
  R/preparation,
)
box::reload(pull)
box::reload(preparation)


# 2. Pull Datasets -----------------------------------------------------------


pull$pull_conditions(pull_from_raw = TRUE)
pull$pull_evidences(pull_from_raw = TRUE)
train <- pull$pull_training(pull_from_raw = TRUE) # adds patient ID to dataset


# 3. Preparation -------------------------------------------------------------

train_arrow <- arrow::open_dataset("data/processed/train.arrow", 
                                   format = "arrow")

box::help(preparation$unnest_evidences)


train_sample <- train_arrow %>% 
  slice_head(n = 10) %>% 
  collect() %>% 
  preparation$add_patient_ID()

train <- train_arrow %>% 
  collect() %>% 
  preparation$add_patient_ID()


# 3.1. Unnest evidences ----------------------------------------------------

start <- Sys.time()
preparation$unnest_evidences_batch(
  train, 
  chunk_size = 100000,
  result_path = "data/processed/prep1_unnestevidences.arrow"
)
end <- Sys.time()
time <- difftime(end,start,units = "mins")
time
# Note: 
# Currently, it takes 7-9 mins for all train rows
# First attempt, I was able to reduce it to 5 mins
# With chunking the data, I was able to reduce it to 3 mins



# 3.2. Unnest differential diganosis --------------------------------------

box::use(arrow[read_feather,write_feather],
         purrr[iwalk])

prep1 <- 
  read_feather(file =
                 "data/processed/prep1_unnestevidences.arrow")
prep1_sample <- prep1 %>% 
  slice_sample(n = 100000)


start <- Sys.time()
test <- 
  preparation$unnest_differential(prep1_sample, count_diagnosis_keep = 100)
end <- Sys.time()
time <- difftime(end,start,units = "mins")
time


batches <- 
  preparation$divide_data_to_batches(prep1, rows_per_batch = 100000)
rm(prep1)
gc()
length(batches)
start <- Sys.time()
iwalk(batches, 
      .f = function(x, idx){
        output <- preparation$unnest_differential(x, 
                                                  count_diagnosis_keep = 15)
        arrow::write_feather(output, 
                             sink = 
                               paste0("data/processed/unnested_differential_batch",
                                      idx,".arrow"))
        rm(output)
        gc()
        
      }, 
      .progress = TRUE)
end <- Sys.time()
time <- difftime(end,start,units = "mins")
time

