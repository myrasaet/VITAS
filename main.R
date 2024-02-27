


# Setup -------------------------------------------------------------------

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


# Pull Datasets -----------------------------------------------------------


pull$pull_conditions(pull_from_raw = TRUE)
pull$pull_evidences(pull_from_raw = TRUE)
train <- pull$pull_training(pull_from_raw = TRUE) # adds patient ID to 


# Preparation -------------------------------------------------------------

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

(test <- preparation$unnest_evidences(train)) %>% system.time()

test


## TODO
# - speed up code, use dtplyr where possible