
# 1. Setup -------------------------------------------------------------------

renv::restore()

# box::purge_cache()
box::use(
  arrow[read_feather,write_feather],
  purrr[iwalk, map],
  magrittr[...],
  dplyr[...],
  R/pull,
  R/preparation,
  R/coda,
  R/evaluation
)
# box::reload(pull)
# box::reload(preparation)
# box::reload(evaluation)

modeling_objects <-
  pins::board_folder("modeling_objects")


# 2. Pull Datasets -----------------------------------------------------------

pull$pull_conditions(pull_from_raw = TRUE)
pull$pull_evidences(pull_from_raw = TRUE)
train <- pull$pull_training(pull_from_raw = TRUE) # adds patient ID to dataset


# 3. Preparation -------------------------------------------------------------

train_arrow <- arrow::open_dataset("data/processed/train.arrow", 
                                   format = "arrow")

# box::help(preparation$unnest_evidences)

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

#' In this section, we train the CoDa model on the prepared training data.

#' Do feature engineering
coda$feature_engineering(
  training_data_path = "data/processed/unnested_differential_wide",
  sink = "data/processed/coda_steps/training_fe.rds"
  )

training_fe <- 
  readr::read_rds("data/processed/coda_steps/training_fe.rds")

# Fit the model
coda$model_fitting(
  training_fe = training_fe,
  # Save pin to modeling_objects board
  board = modeling_objects, 
  name = "coda_fitted_model",
  versioned = TRUE,
  metadata =
    list(
      training_data = "release_train_patients.csv",
      size = "1,025,602 patients",
      scaled_predictors = FALSE
    )
)

# Read the fitted model
fitted_model <-
  modeling_objects %>% 
  pins::pin_read(
    name = "coda_fitted_model"
  )

# Predict using the fitted model
coda$prediction(
  fitted_model = fitted_model,
  data_list = training_fe,
  sink = "data/processed/coda_steps/predictions.rds"
)



# 5. Evaluation --------------------------------------------------------------

#' In this section, we test the performance of the CoDa model by predicting
#' on unseen data.

# Get the DDXPlus test dataset
pull$pull_testing(pull_from_raw = FALSE)
test_arrow <- arrow::open_dataset("data/processed/test.arrow", 
                                   format = "arrow")

# Add patient IDs
test <- test_arrow %>% 
  dplyr::collect() %>% 
  preparation$add_patient_ID()

# Sample patients
set.seed(31)
sample_patients <- 
  test$patientId %>% 
  sample(size = 132448)
test <- 
  test %>% 
  filter(patientId %in% sample_patients)

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

# Divide to batches for easier computation
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
rm(batches, test, test_arrow)
gc()

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

# Compute Metrics
prediction_evaluation_list <-
  evaluation$pivot_longer_dd(
    predictions_tidy,
    differential_prob_threshold_actual = 0.01,
    differential_prob_threshold_pred = 0.01,
    high_severity_ddprob_threshold = 0.01
  )

rm(predictions_tidy, fitted_model, test_fe)
gc()

readr::write_rds(
  prediction_evaluation_list, 
  "data/processed/validation/prediction_evaluation_list.rds"
)

prediction_evaluation_list <-
  readr::read_rds(
    "data/processed/validation/prediction_evaluation_list.rds"
  )

# Filter PATHOLOGY
sample_pathologies <- 
  c("Tuberculosis","GERD","SLE","HIV (initial infection)",
    "Pulmonary neoplasm")

prediction_evaluation_list_filtered <- 
  prediction_evaluation_list
prediction_evaluation_list_filtered$actual_pathology <- 
  prediction_evaluation_list$actual_pathology %>% 
  filter(PATHOLOGY == sample_pathologies[5])

filtered_patients <- 
  prediction_evaluation_list_filtered$actual_pathology %>% 
  pull(patientId)
prediction_evaluation_list_filtered$assessment_df <- 
  prediction_evaluation_list$assessment_df %>% 
  filter(patientId %in% filtered_patients)
prediction_evaluation_list_filtered$prediction_df <- 
  prediction_evaluation_list$prediction_df %>% 
  filter(patientId %in% filtered_patients)
# *******************

metrics <-
  evaluation$compute_differential_metrics(
    prediction_evaluation_list_filtered,
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

#' In the following section, we repeat step 5 but using the test dataset with
#' adaptive questionnaire used. 

test_adaptive_raw <- 
  readr::read_csv("data/input/questionnaire_df.csv")

# Get differential diagnosis
test_arrow <- arrow::open_dataset("data/processed/test.arrow", 
                                  format = "arrow")
test <- 
  test_arrow %>% 
  select(DIFFERENTIAL_DIAGNOSIS, AGE, SEX, PATHOLOGY) %>% 
  dplyr::collect() %>% 
  preparation$add_patient_ID()

evidence

# # Reference
# sample <- 
#   test_adaptive_raw %>% 
#   head() %>% 
#   collect()

# Sample patients
set.seed(31)
sample_patients <- 
  test$patientId %>% 
  sample(size = 132448)

test_adaptive <-
  test_adaptive_raw %>% 
  # dplyr::collect() %>% 
  mutate(EVIDENCES = test_adaptive_raw$binary_evidences) %>% 
  preparation$add_patient_ID() %>% 
  select(-`...1`) %>% 
  # get DD
  left_join(test) %>% 
  # remove unecessary fields
  select(-binary_evidences_all, -binary_evidences, -missed_evidence)

# get sample of patients (needed if more than 800k patients)
# test_adaptive <- 
#   test_adaptive %>% 
#   filter(patientId %in% sample_patients)

# Unnest evidences
start <- Sys.time()
preparation$unnest_evidences_batch(
  test_adaptive, 
  # chunk_size = 100000,
  result_path = "data/processed/validation_train_adaptive/prep1_unnestevidences.arrow"
)
end <- Sys.time()
time <- difftime(end,start,units = "mins")
time



# Unnest differential diagnosis
test_adaptive_prep1 <- 
  read_feather(file =
                 "data/processed/validation_train_adaptive/prep1_unnestevidences.arrow")

batches <- 
  preparation$divide_data_to_batches(test_adaptive_prep1, rows_per_batch = 10000)



start <- Sys.time()
iwalk(batches, 
      .f = function(x, idx){
        output <- preparation$unnest_differential(x)
        arrow::write_feather(
          output, 
          sink = 
            paste0("data/processed/validation_train_adaptive/unnested_differential/unnested_differential_batch",
                   idx,".arrow"))
        rm(output)
        gc()
        
      }, 
      .progress = TRUE)
end <- Sys.time()
time <- difftime(end,start,units = "mins")
time

rm(batches)
rm(test_adaptive_raw, test_adaptive)
rm(test, test_adaptive_prep1, test_arrow, test_prep1)
gc()

# Feature engineering
coda$feature_engineering(
  training_data_path = "data/processed/validation_train_adaptive/unnested_differential",
  sink = "data/processed/validation_train_adaptive/test_fe.rds"
)

test_fe <-
  readr::read_rds("data/processed/validation_train_adaptive/test_fe.rds")

fitted_model <-
  modeling_objects %>% 
  pins::pin_read(name = "coda_fitted_model")


coda$prediction(
  fitted_model = fitted_model,
  data_list = test_fe,
  sink = "data/processed/validation_train_adaptive/predictions.rds"
)

predictions_tidy <-
  readr::read_rds("data/processed/validation_train_adaptive/predictions.rds")



# Metrics
prediction_evaluation_list <-
  evaluation$pivot_longer_dd(
    predictions_tidy,
    differential_prob_threshold_actual = 0.01,
    differential_prob_threshold_pred = 0.01,
    high_severity_ddprob_threshold = 0.01
  )

rm(predictions_tidy, fitted_model, test_fe)
gc()

readr::write_rds(
  prediction_evaluation_list, 
  "data/processed/validation_train_adaptive/prediction_evaluation_list.rds"
)

# Filter PATHOLOGY
sample_pathologies <- 
  c("Tuberculosis","GERD","SLE","HIV (initial infection)",
    "Pulmonary neoplasm")
prediction_evaluation_list_filtered <- 
  prediction_evaluation_list
prediction_evaluation_list_filtered$actual_pathology <- 
  prediction_evaluation_list$actual_pathology %>% 
  filter(PATHOLOGY == sample_pathologies[5])

filtered_patients <- 
  prediction_evaluation_list_filtered$actual_pathology %>% 
  pull(patientId)
prediction_evaluation_list_filtered$assessment_df <- 
  prediction_evaluation_list$assessment_df %>% 
  filter(patientId %in% filtered_patients)
prediction_evaluation_list_filtered$prediction_df <- 
  prediction_evaluation_list$prediction_df %>% 
  filter(patientId %in% filtered_patients)
# *******************

metrics_adaptive <-
  evaluation$compute_differential_metrics(
    prediction_evaluation_list_filtered,
    workers = 4,
    save_path = "data/processed/validation_train_adaptive/"
  )

# # Save metrics
# readr::write_rds(metrics_adaptive,
#                  "data/processed/validation_train_adaptive/metrics.rds")

metrics <- 
  readr::read_rds("data/processed/validation_train_adaptive/metrics.rds")
metrics_adaptive <- 
  readr::read_rds("data/processed/validation_train_adaptive/metrics.rds")


metrics$metrics_summary
metrics_adaptive$metrics_summary

evaluation$manual_evaluate(
  predictions_tidy,
  differential_prob_threshold_actual = 0.01,
  differential_prob_threshold_pred = 0.01,
  # patientId = 91166
)



# Prototype  --------------------------------------------------------------

# To create a prototype results to put in the capstone paper

fitted_model <-
  modeling_objects %>% 
  pins::pin_read(
    name = "coda_fitted_model"
  )

evidences <- 
  pull$pull_evidences()


test_arrow <- arrow::open_dataset("data/processed/test.arrow", 
                                  format = "arrow")

test <- 
  test_arrow %>% 
  slice_head(n = 99) %>% 
  dplyr::collect() %>% 
  preparation$add_patient_ID()

# Create a csv file and edit the evidences and demographic profile
# write.csv(test, "data/input/prototype.csv")

prototype_sample <- 
  readr::read_csv("data/input/prototype.csv",
           col_types = "idccccc") %>% 
  bind_rows(test) %>% 
  as_tibble()

# Unnest evidences
prot_evidences <-
  preparation$unnest_evidences(data = prototype_sample)

# Unnest differential diagnosis
prot_dd <- 
  preparation$unnest_differential(prot_evidences)

prep <- 
  preparation$coda_prediction_preparation(prot_dd)

prot_pred <- 
  coda$prediction_raw(
  fitted_model = fitted_model,
  predictors_df = prep$predictors[100,]
)

prot_pred_sample <- 
  prot_pred[1,]

res <- 
  tibble(
  diagnosis = names(prot_pred_sample),
  prob = prot_pred_sample
) %>% 
  arrange(desc(prob)) %>% 
  filter(prob >= 0.01)

clipr::write_clip(res)

map(res$diagnosis[1:22],
    .f = function(x){
      
      # x <- res$diagnosis[14]
      
      coef_model <- 
        coef(fitted_model) %>% 
        as.data.frame() %>% 
        tibble::rownames_to_column(var = "predictor") %>% 
        as_tibble() %>% 
        select(predictor, contains(x))

      
      samplepat_evidences <- 
        prep$predictors %>% 
        slice(100) %>% 
        pivot_longer(everything(),
                     names_to = "predictor",
                     values_to = "value") %>% 
        mutate(
          predictor = str_replace(predictor,
                                  pattern = "@",
                                  replacement = "\\.")
        )
      samplepat_evidences
      
      coef_model2 <- 
        coef_model %>% 
        # filter(predictor != "(Intercept)") %>% 
        left_join(samplepat_evidences,
                  by = join_by(predictor)) %>% 
        # Get the yes diagnoses
        mutate(
          value = 
            if_else(predictor == "(Intercept)",
                    1, value),
          evidence_yes = 
            if_else(predictor %in% c("E_173",
                                     "E_215",
                                     "E_217",
                                     "E_53",
                                     "E_148"),
                    "yes_diag",
                    "no_diag"),
          evidence_yes = 
            if_else(predictor %in% c("AGE","SEX.M"),
                    "yes_diag",
                    evidence_yes)
        )
      
      # coef_model2 %>% 
      #   mutate(coef_value = diag_GERD * value) %>% 
      #   filter(abs(coef_value) > 0) %>% 
      #   # filter(is.na(coef_value)) %>% 
      #   summarize(pred = sum(coef_value, na.rm = TRUE))
      
      xvar <- sym(x)
      
      coef_model2 %>% 
        filter(evidence_yes == "yes_diag") %>% 
        mutate(coef_value = !!xvar * value,
               abs_coef = abs(coef_value)) %>% 
        arrange(desc(coef_value)) %>% 
        left_join(
          evidences %>% 
            select(question_en, name),
          by = join_by(predictor == name)
        ) %>% 
        mutate(
          question_en = ifelse(is.na(question_en),
                 predictor,
                 question_en),
          question_en = ifelse(question_en == "AGE",
                               paste0("AGE.",value),
                               question_en)) %>% 
        filter(coef_value > 0) %>% 
        select(question_en, coef_value)
      
      
    }) %>% system.time()




# Notes -------------------------------------------------------------------

# FIXME
# fix coda data prep
# * remove pathology as predictor
# * dummify process includes all levels. should only get n - 1
# to this for all categorical evidences. not anymore for multiclass
