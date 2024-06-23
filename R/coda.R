

# 1. Setup -------------------------------------------------------------------

box::use(
  stats[...],
  magrittr[...],
  dplyr[...],
  tidyr[...],
  caret[...],
  arrow[...],
  dtplyr[...],
  compositions[...],
  pins[...],
  stringr[...],
  purrr[...],
  R/pull,
  R/utils,
  R/preparation,
  R/evaluation
)



# 2. Fetch Data --------------------------------------------------------------

evidences <- 
  pull$pull_evidences(pull_from_raw = FALSE)


# 2.5. Sampling (Only for initial run) -----------------------------------------


# 3. Data Preparation --------------------------------------------------------



# 4. Feature Engineering -----------------------------------------------------

#> Prepare the data for compositional regression

#' Feature engineer the variables for compositional regression fitting
#' 
#' @param training_data_path A string. The path to the training data. 
#' The dataset should have unnested evidences and unnested differential.
#' The differential diagnosis should be spread to 49 columns.
#' 
#' @export
feature_engineering <-
  function(training_data_path,
           sink){
    
    train_arrow <- 
      open_dataset(training_data_path,
                   format = "feather")
    
    train_prep1 <- 
      as_tibble(train_arrow)
    
    rlang::inform("Feature engineering...")
    train_prep2 <-
      preparation$coda_prediction_preparation(data = train_prep1)
    
    
    readr::write_rds(train_prep2,
                     file = sink)
    rlang::inform(
      paste0("Saved feature engineered dataset to ",
             "data/processed/coda_steps/training_fe.rds")
    )
    rm(train_prep1)
    gc()
    
  }





# 5. Modeling -------------------------------------------------------------

#' Fit the model
#' 
#' @param training_fe A list. Contains the elements needed for training or
#' prediction namely:
#' predictors, id, and outcome. Each of them is a tibble. The predictors and 
#' outcome datasets are already feature engineered.
#' 
#' @export
model_fitting <-
  function(training_fe,
           ...){
    
    # training_fe$outcome <- 
    #   training_fe$outcome %>% 
    #   slice_head(n = 2000)
    
    
    Y <- acomp(training_fe$outcome)
    X <- as.matrix(training_fe$predictors)
    train_data <- 
      data.frame(Y,X)
    
    # xTest <- X[2,]
    
    diffdiag <- 
      names(train_data) %>% 
      str_subset("diag_")
    
    # Construct the formula dynamically
    comp_formula <- paste("Y ~ . -", paste(diffdiag, collapse = " - "))
    
    rlang::inform("Training model...")
    mod <- 
      lm(as.formula(comp_formula), data = train_data)
    
    # # TEST predict
    # test_pred <-
    #   predict(mod,
    #           newdata = train_data[4,])
    
    rlang::inform("Saving model...")
    pins::pin_write(x = mod,
                    ...
    )
    
    
  }

#' Do predictions
#' 
#' @param fitted_model A CoDa fitted model.
#' @param data_list A list. Contains the elements needed for training or
#' prediction namely:
#' predictors, id, and outcome. Each of them is a tibble. The predictors and 
#' outcome datasets are already feature engineered.
#' @param sink A string. The path where to save the predictions
#' 
#' @export
prediction <-
  function(fitted_model,
           data_list,
           sink){
    
    # fitted_model <- mod

    
    Y <- acomp(data_list$outcome)
    X <- as.matrix(data_list$predictors)
    test_data <-
      data.frame(Y,X)
    
    rlang::inform("Predicting...")

    predictions <- 
      predict(fitted_model,
              newdata = test_data)
    
    
    predictions_tidy <-
      c(data_list,
        list(predicted_outcome = as_tibble(predictions)))
    
    rlang::inform("Saving predictions...")
    readr::write_rds(
      predictions_tidy,
      file = sink
    )
    
    predictions_tidy
    
  }


# test[1,]$DIFFERENTIAL_DIAGNOSIS
# test_data[1,]
# test_pred <- 
#   predict(fitted_model,
#           newdata = train_data[1,])
# 
# pred <- 
#   test_pred %>% 
#   as_tibble() %>% 
#   pivot_longer(
#     everything()
#   ) %>% 
#   arrange(desc(value))
# 
# actual <- 
#   train_data[1,] %>% 
#   as_tibble() %>% 
#   select(starts_with("diag_")) %>% 
#   pivot_longer(
#     everything()
#   ) 
# 
# round(pred$value,15) == round(actual$value, 15)
# 
# 
# pred %>% 
#   left_join(actual %>% 
#               rename(pred = value),
#             by = join_by(name)) %>% 
#   summarize(
#     across(.cols = c(value, pred),
#            .fns = ~sum(., na.rm = TRUE))
#   )
# 
# test_data[1,] %>% 
#   as_tibble() %>% 
#   select(starts_with("diag_")) %>% 
#   pivot_longer(
#     everything()
#   ) 

# readr::write_rds(
#   predictions_tidy,
#   file = "data/processed/coda_steps/predictions_tidy.rds"
# )


# 6. Evaluation -----------------------------------------------------------


# evaluation$manual_evaluate(
#   predictions_tidy,
#   differential_prob_threshold_actual = 0.01,
#   differential_prob_threshold_pred = 0.01,
#   # patientId = 678267
# )
# 
# 
# prediction_evaluation_list <-
#   evaluation$pivot_longer_dd(
#     predictions_tidy,
#     differential_prob_threshold_actual = 0.01,
#     differential_prob_threshold_pred = 0.01
#   )
# 
# 
# metrics <-
#   evaluation$compute_differential_metrics(prediction_evaluation_list)
# 
# modeling_objects %>% 
#   pins::pin_write(x = metrics, 
#                   name = "model_metrics_coda", 
#                   versioned = TRUE, 
#                   description = paste0("100,000 patients were used to",
#                                        " train this model."))
# 
# modeling_objects %>% 
#   pins::pin_read(name = "model_metrics_coda")
