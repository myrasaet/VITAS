

# renv::activate()
# renv::snapshot()
renv::restore()
# renv::install()

box::purge_cache()
box::use(
  R/pull,
  R/preparation
)
box::reload(pull)
box::reload(preparation)

pull$pull_conditions(pull_from_raw = TRUE)
pull$pull_evidences(pull_from_raw = TRUE)
train <- pull$pull_training(pull_from_raw = TRUE) # adds patient ID to 
train <- arrow::open_dataset("data/processed/train.arrow", format = "arrow")

box::help(preparation$unnest_evidences)

train_sample <- train %>% 
  slice_head(n = 10000) %>% 
  collect() %>% 
  preparation$add_patient_ID()

(test <- unnest_evidences(train_sample)) %>% system.time()
