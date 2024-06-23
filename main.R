


# 1. Setup -------------------------------------------------------------------

renv::restore()

# box::purge_cache()
box::use(
  arrow[read_feather,write_feather],
  purrr[iwalk],
  magrittr[...],
  dplyr[...],
  R/pull,
  R/preparation,
  R/coda,
  R/evaluation
)
# box::reload(pull)
# box::reload(preparation)

modeling_objects <-
  pins::board_folder("modeling_objects")


# 2. Pull Datasets -----------------------------------------------------------


pull$pull_conditions(pull_from_raw = TRUE)
pull$pull_evidences(pull_from_raw = TRUE)
train <- pull$pull_training(pull_from_raw = TRUE) # adds patient ID to dataset


# 3. Preparation -------------------------------------------------------------

train_arrow <- arrow::open_dataset("data/processed/train.arrow", 
                                   format = "arrow")

box::help(preparation$unnest_evidences)


train_sample <- train_arrow %>% 
  slice_head(n = 100000) %>% 
  collect() %>% 
  preparation$add_patient_ID()

train <- train_arrow %>% 
  collect() %>% 
  preparation$add_patient_ID()


# 3.1. Unnest evidences ----------------------------------------------------

#> Description:
#> In this section, we unnest the evidences.
#> Since evidences are contained in one column,
#> we unnest them wider such that each evidence has
#> its own column.This is so each evidence can 
#> become its own predictor.

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

#> Description:
#> In this section, we unnest the differential diagnosis.
#> This means separating all the diagnosis per patient
#> and giving each one its own row. This will make it 
#> easier later on to manipulate the data.



prep1 <- 
  read_feather(file =
                 "data/processed/prep1_unnestevidences.arrow")


batches <- 
  preparation$divide_data_to_batches(prep1, rows_per_batch = 100000)
rm(prep1)
gc()
length(batches)
start <- Sys.time()
iwalk(batches, 
      .f = function(x, idx){
        output <- preparation$unnest_differential(x)
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


preparation$dd_pivot_wider(
  path_unnested_differential =
    "data/processed/unnested_differential/",
  save_path =
    "data/processed/unnested_differential_wide/"
)



# 4. Training -------------------------------------------------------------

coda$feature_engineering(
  training_data_path = "data/processed/unnested_differential_wide",
  sink = "data/processed/coda_steps/training_fe.rds"
  )

training_fe <- 
  readr::read_rds("data/processed/coda_steps/training_fe.rds")

coda$model_fitting(
  training_fe = training_fe,
  # Save pin to modeling_objects board
  board = modeling_objects, 
  name = "coda_fitted_model",
  versioned = TRUE,
  metadata =
    list(
      training_data = "release_train_patients.csv",
      size = "1,025,602 patients"
    )
)

fitted_model <-
  modeling_objects %>% 
  pins::pin_read(
    name = "coda_fitted_model"
  )

coda$prediction(
  fitted_model = fitted_model,
  data_list = training_fe,
  sink = "data/processed/coda_steps/predictions.rds"
)



# 5. Evaluation --------------------------------------------------------------

pull$pull_testing(pull_from_raw = FALSE)
test_arrow <- arrow::open_dataset("data/processed/test.arrow", 
                                   format = "arrow")

test <- test_arrow %>% 
  dplyr::collect() %>% 
  preparation$add_patient_ID()

test %>% 
  filter(patientId == 298) %>% 
  pull(DIFFERENTIAL_DIAGNOSIS)
  


# Unnest evidences
start <- Sys.time()
preparation$unnest_evidences_batch(
  test, 
  chunk_size = 100000,
  result_path = "data/processed/validation/prep1_unnestevidences.arrow"
)
end <- Sys.time()
time <- difftime(end,start,units = "mins")
time

# Unnest differential diagnosis
test_prep1 <- 
  read_feather(file =
                 "data/processed/validation/prep1_unnestevidences.arrow")


# glimpse(test_prep1)
# test_prep1$E_204 %>% unique()
# evidences %>% 
#   filter(name == "E_204") %>% 
#   pull(possible_values)

batches <- 
  preparation$divide_data_to_batches(test_prep1, rows_per_batch = 100000)
rm(test_prep1)
gc()

start <- Sys.time()
iwalk(batches, 
      .f = function(x, idx){
        output <- preparation$unnest_differential(x)
        arrow::write_feather(
          output, 
          sink = 
            paste0("data/processed/validation/unnested_differential/unnested_differential_batch",
                   idx,".arrow"))
        rm(output)
        gc()
        
      }, 
      .progress = TRUE)
end <- Sys.time()
time <- difftime(end,start,units = "mins")
time

# Feature engineering
coda$feature_engineering(
  training_data_path = "data/processed/validation/unnested_differential",
  sink = "data/processed/validation/test_fe.rds"
)

test_fe <-
  readr::read_rds("data/processed/validation/test_fe.rds")
fitted_model <-
  modeling_objects %>% 
  pins::pin_read(name = "coda_fitted_model")


coda$prediction(
  fitted_model = fitted_model,
  data_list = test_fe,
  sink = "data/processed/validation/predictions.rds"
)

predictions_tidy <-
  readr::read_rds("data/processed/validation/predictions.rds")

# Metrics


prediction_evaluation_list <-
  evaluation$pivot_longer_dd(
    predictions_tidy,
    differential_prob_threshold_actual = 0.01,
    differential_prob_threshold_pred = 0.01
  )


metrics <-
  evaluation$compute_differential_metrics(
    prediction_evaluation_list,
    workers = 4,
    save_path = "data/processed/validation/"
    )

# Save metrics
readr::write_rds(metrics,
                 "data/processed/validation/metrics.rds")

metrics <- 
  readr::read_rds("data/processed/validation/metrics.rds")
metrics$raw_scores %>% 
  arrange(ddr)


evaluation$manual_evaluate(
  predictions_tidy,
  differential_prob_threshold_actual = 0.01,
  differential_prob_threshold_pred = 0.01,
  patientId = 91166
)


# 6. Using adaptive questionnaire --------------------------------------------

test_adaptive_raw <- 
  readr::read_csv("data/input/questionnaire_df.csv")

evidence

# Reference
sample <- 
  test_arrow %>% 
  head() %>% 
  collect()

test_adaptive <-
  test_arrow %>% 
  dplyr::collect() %>% 
  mutate(EVIDENCES = test_adaptive_raw$binary_evidences) %>% 
  preparation$add_patient_ID()

# Unnest evidences
start <- Sys.time()
preparation$unnest_evidences_batch(
  test_adaptive, 
  # chunk_size = 100000,
  result_path = "data/processed/validation_adaptive/prep1_unnestevidences.arrow"
)
end <- Sys.time()
time <- difftime(end,start,units = "mins")
time



# Unnest differential diagnosis
test_adaptive_prep1 <- 
  read_feather(file =
                 "data/processed/validation_adaptive/prep1_unnestevidences.arrow")

batches <- 
  preparation$divide_data_to_batches(test_adaptive_prep1, rows_per_batch = 100000)



start <- Sys.time()
iwalk(batches, 
      .f = function(x, idx){
        output <- preparation$unnest_differential(x)
        arrow::write_feather(
          output, 
          sink = 
            paste0("data/processed/validation_adaptive/unnested_differential/unnested_differential_batch",
                   idx,".arrow"))
        rm(output)
        gc()
        
      }, 
      .progress = TRUE)
end <- Sys.time()
time <- difftime(end,start,units = "mins")
time




# Feature engineering
coda$feature_engineering(
  training_data_path = "data/processed/validation_adaptive/unnested_differential",
  sink = "data/processed/validation_adaptive/test_fe.rds"
)

test_fe <-
  readr::read_rds("data/processed/validation_adaptive/test_fe.rds")

fitted_model <-
  modeling_objects %>% 
  pins::pin_read(name = "coda_fitted_model")


coda$prediction(
  fitted_model = fitted_model,
  data_list = test_fe,
  sink = "data/processed/validation_adaptive/predictions.rds"
)

predictions_tidy <-
  readr::read_rds("data/processed/validation_adaptive/predictions.rds")



# Metrics
prediction_evaluation_list <-
  evaluation$pivot_longer_dd(
    predictions_tidy,
    differential_prob_threshold_actual = 0.01,
    differential_prob_threshold_pred = 0.01
  )


metrics_adaptive <-
  evaluation$compute_differential_metrics(
    prediction_evaluation_list,
    workers = 4,
    save_path = "data/processed/validation_adaptive/"
  )

# Save metrics
readr::write_rds(metrics_adaptive,
                 "data/processed/validation_adaptive/metrics.rds")

metrics <- 
  readr::read_rds("data/processed/validation/metrics.rds")
metrics_adaptive <- 
  readr::read_rds("data/processed/validation_adaptive/metrics.rds")


metrics$metrics_summary
metrics_adaptive$metrics_summary

evaluation$manual_evaluate(
  predictions_tidy,
  differential_prob_threshold_actual = 0.01,
  differential_prob_threshold_pred = 0.01,
  # patientId = 91166
)



# Notes -------------------------------------------------------------------

# FIXME
# fix coda data prep
# * remove pathology as predictor
# * dummify process includes all levels. should only get n - 1
# to this for all categorical evidences. not anymore for multiclass
