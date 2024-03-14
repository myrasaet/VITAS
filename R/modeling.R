

# Setup -------------------------------------------------------------------

library(tidymodels); tidymodels_prefer()

box::use(
  R/pull,
  R/preparation
)

box::use(
  here[here],
  pins[...],
  stringr[...],
  furrr[...],
  arrow[...],
  vetiver[...],
  readr[...]
)

modeling_objects_board <-
  board_folder(path = here("modeling_objects"), versioned = TRUE)

set.seed(147)


# Import Data -------------------------------------------------------------



train_processed <- arrow::open_dataset(
  here("data",
       "processed",
       "unnested_differential"),
  format = "arrow"
)

train <- arrow::open_dataset(
  here("data",
       "processed",
       "train.arrow"),
  format = "arrow"
)


evidences <-
  pull$pull_evidences(pull_from_raw = FALSE)


# Sampling ----------------------------------------------------------------

# We should make sure that the sample contain all pathologies
train_collect <- collect(train)

# Sample 200 patients for each pathology
set.seed(130)
sample_patients <- train_collect |> 
  preparation$add_patient_ID() |> 
  group_by(PATHOLOGY) |> 
  slice_sample(n = 200) |> 
  ungroup() |> 
  pull(patientId)

train_sample <- 
  train_processed |> 
  filter(patientId %in% sample_patients) |> 
  collect()

rm(sample_patients,train_collect,train_processed,train)
gc()


# Preparation -------------------------------------------------------------

# Add importance weighting
# remove unecessary variables

ddxplus_sample <- train_sample |> 
  select(-differential_score,
         -differential_rank,
         -PATHOLOGY,
         -DIFFERENTIAL_DIAGNOSIS) |> 
  mutate(differential_case_weight = 
           importance_weights(differential_case_weight)) |> 
  arrange(patientId)
ddxplus_sample

# Splitting ---------------------------------------------------------------

set.seed(502)
ddxplus_split <- group_initial_split(ddxplus_sample,
                                     group = patientId,
                                     prop = 0.70)
ddxplus_split

rm(ddxplus_sample, train_sample)
gc()


# modeling_objects_board |> 
#   pin_write(x = ddxplus_split,
#             name = "training_split",
#             title = "Train Split Object")

ddxplus_split <-
  modeling_objects_board |> 
  pin_read("training_split")


## 2.2. Specify model (Parsnip) --------------------------------------------------------

# parsnip_addin()

# decision_tree_rpart_spec <-
#   decision_tree(tree_depth = tune(), 
#                 min_n = tune(), 
#                 cost_complexity = tune()) %>%
#   set_engine('rpart') %>%
#   set_mode('classification')
# 
# decision_tree_rpart_spec <-
#   decision_tree() %>%
#   set_engine('rpart') %>%
#   set_mode('classification')

library(ranger)

# rand_forest_ranger_spec <-
#   rand_forest(mtry = 15, min_n = 5, trees = 500) %>%
#   set_engine('ranger') %>%
#   set_mode('classification')

rand_forest_ranger_spec <-
  rand_forest(mtry = tune(), min_n = tune()) %>%
  set_engine('ranger') %>%
  set_mode('classification')

# rand_forest_ranger_spec_2 <-
#   rand_forest(mtry = tune(), min_n = tune(), trees = tune()) %>%
#   set_engine('ranger') %>%
#   set_mode('classification')



## 2.3. Feature engineering (recipe) ------------------------------------------

# The features of the model are:
# 1. patientId (ID)
# 2. AGE
# 3. differential_case_weight (case weight)
# 4. SEX
# 5. INITIAL_EVIDENCE
# 6. All evidences

ddxplus_split <-
  modeling_objects_board |> 
  pin_read(name = "training_split")
ddxplus_train <- training(ddxplus_split)

recipe <-
  recipe(diagnosis ~ ., 
         data = ddxplus_train) |> 
  update_role(patientId, new_role = "Id") |> 
  step_other(INITIAL_EVIDENCE, threshold = 0.01, skip = FALSE) |> 
  step_novel(all_string_predictors()) |> 
  step_mutate(across(where(is.character),
                     .fns = as.factor)) |> 
  step_zv(all_predictors()) |> 
  step_YeoJohnson(all_numeric_predictors()) |> 
  step_normalize(all_numeric_predictors()) |> 
  step_pca(all_numeric_predictors(), threshold = 0.40, skip = FALSE)

# check_recipe <-
#   recipe |> 
#   prep() |> 
#   bake(new_data = ddxplus_train)

## 2.4. Create the workflow --------------------------------------------------

wf <-
  workflow() %>% 
  add_case_weights(differential_case_weight) %>% 
  add_recipe(recipe) %>% 
  add_model(rand_forest_ranger_spec)

modeling_objects_board |> 
  pin_write(x = wf,
            name = "workflow",
            title = "Tidymodels Workflow Object",
            description = "Tidymodels workflow object with 
            case weights, recipe, and random forest model.")

rm(recipe)

## 2.5. Initial Model fitting -------------------------------------------------

training_start <- Sys.time()

fit <- 
  wf %>% 
  fit(ddxplus_train)

training_end <- Sys.time()
training_duration <- training_end - training_start
training_duration
# 41 seconds for E_201 rpart
# 28.5 seconds for random forest default hyperparameters and 0.40 threshold pca

# readr::write_rds(fit,
#                  file = here("data",
#                              "output",
                             # "workflow_sample_rpart_fit_E201.Rds"))

training_start <- Sys.time()

lastfit <- 
  wf %>% 
  last_fit(split = ddxplus_split)

training_end <- Sys.time()
training_duration <- training_end - training_start
training_duration

rm(wf, recipe, ddxplus_train)

collect_metrics(lastfit)


fit

class(fit)

all_predictions <- bind_rows(
  augment(fit, new_data = ddxplus_test) %>% 
    mutate(type = "test")
)

collect_metrics(all_predictions)


all_predictions %>%
  group_by(type) %>%
  metrics(nearness, .pred)

accuracy_vec(truth = as.factor(ddxplus_test$diagnosis),
             estimate = all_predictions$.pred_class)
# accuracy was 0.1168 E_201 rpart

table(all_predictions$.pred_class)


# Grid search -------------------------------------------------------------

ranger_param <- extract_parameter_set_dials(rand_forest_ranger_spec)
ranger_param %>% 
  extract_parameter_dials("min_n")
ranger_param <- 
  rand_forest_ranger_spec %>% 
  extract_parameter_set_dials() |> 
  update(
    mtry = mtry(c(1,20)),
    min_n = min_n(c(2,40))
  )

grid_regular(ranger_param, levels = 2)
dials::finalize(object = ranger_param)

# Random grid

set.seed(1301)
ranger_param %>% 
  grid_random(size = 1000) %>% # 'size' is the number of combinations
  summary()

# Latin hypercube 

set.seed(1301)
ranger_param %>% 
  grid_latin_hypercube(size = 20, original = FALSE) 

# Cross validation

set.seed(1304)
ddxplus_train_folds <- vfold_cv(ddxplus_train, v = 3)

rpart_param <- 
  wf %>% 
  extract_parameter_set_dials()


roc_res <- metric_set(roc_auc)
set.seed(1305)
rpart_reg_tune <-
  wf %>%
  tune_grid(
    ddxplus_train_folds,
    grid = ranger_param %>% 
      grid_latin_hypercube(size = 10, original = FALSE),
    metrics = roc_res
  )
rpart_reg_tune
collect_metrics(rpart_reg_tune)
show_best(rpart_reg_tune)
autoplot(rpart_reg_tune)

readr::write_rds(rpart_reg_tune,
                 here::here("data",
                            "output",
                            "rpart_reg_tune"))

# TODO 
# * Check if there are no missing rows (no missing rows   )
# * Try oversampling or undersampling (Done - we did undersampling)
# * Experiment with other data processing steps like normalization



# Re-fit with best fitting metrics ----------------------------------------


rpart_reg_tune <-
  readr::read_rds(here::here("data",
                             "output",
                             "rpart_reg_tune"))

keep_pred <- control_resamples(save_pred = TRUE)

collect_metrics(rpart_reg_tune) |> 
  arrange(desc(mean))
show_best(rpart_reg_tune)

# Get top 5 best tune parameters.
# We will save predictions then evaluate them later with
# differential diagnosis metrics
best_tune <- 
  show_best(rpart_reg_tune) |> 
  select(mtry, min_n)


# Re-fit with improved metrics (testing) ----------------------------------


wf <- modeling_objects_board |> 
  pin_read(name = "workflow")

set.seed(1304)
ddxplus_train_folds <- vfold_cv(ddxplus_train, v = 2)

# modeling_objects_board |> 
#   pin_write(x = ddxplus_train_folds,
#             name = "train_cv_folds",
#             title = "Tidymodels Cross-validation Folds",
#             description = "Tidymodels Cross-validation Folds. 
#             2 folds")

ddxplus_train_folds <- 
  modeling_objects_board |> 
  pin_read(name = "train_cv_folds")

keep_pred <- control_resamples(save_pred = TRUE)

roc_res <- metric_set(roc_auc)


ranger_param_5 <- 
  rand_forest_ranger_spec %>% 
  extract_parameter_set_dials() |> 
  update(
    mtry = mtry(c(2,3)),
    min_n = min_n(c(400,1000))
  )

set.seed(1305)
grid5 <- 
  ranger_param_5 |> 
  grid_latin_hypercube(size = 10)


set.seed(1305)
rpart_reg_tune5 <-
  wf %>%
  tune_grid(
    ddxplus_train_folds,
    grid = grid5,
    metrics = roc_res,
    control = keep_pred
  )
collect_metrics(rpart_reg_tune5) |> 
  arrange(desc(mean))


readr::write_rds(rpart_reg_tune5,
                 here::here("data",
                            "output",
                            "rpart_reg_tune5"))


preds <- collect_predictions(rpart_reg_tune5)
names(preds)



# Metrics evaluation ------------------------------------------------------

rpart_reg_tune <-
  readr::read_rds(here::here("data",
                             "output",
                             "rpart_reg_tune5"))


modeling_objects_board |> 
  pin_write(rpart_reg_tune5, 
            name = "rpart_reg_tune",
            title = "Random Forest Tuning Results", 
            description = "Tuning results with grid_latin_hypercube function
            used for tuning grid search. Correct split count.",
            versioned = T)

rpart_reg_tune <-
  modeling_objects_board |> 
  pin_read("rpart_reg_tune")

# Do last fit
collect_metrics(rpart_reg_tune)
best_fit_hparam <-
  tibble(
    mtry = 3,
    min_n = 531
  )

# rand_forest_last_fit <-
#   wf |> 
#   finalize_workflow(best_fit_hparam) |> 
#   last_fit(split = ddxplus_split)

# modeling_objects_board |>
#   pin_write(rand_forest_last_fit,
#             name = "model_last_fit",
#             title = "Model Last Fit")

rand_forest_final_workflow <-
  extract_workflow(x = rand_forest_last_fit)

modeling_objects_board |>
  pin_write(rand_forest_final_workflow,
            name = "model_final_workflow",
            title = "Model Final Workflow",
            metadata = 
              list(metrics = collect_metrics(rand_forest_last_fit)))

# collect_metrics(rand_forest_last_fit)

# library(vetiver)
# 
# v <- vetiver_model(model = rand_forest_final_workflow,
#               model_name = "rand_forest_final_workflow")
# modeling_objects_board |>
#   vetiver::vetiver_pin_write(v)



pred_prob <- collect_predictions(rand_forest_last_fit)

pred_prob1 <-
  pred_prob |>
  # filter(.config == "Preprocessor1_Model07",
  #        id == "Fold1") |>
  select(matches(".pred")) |> 
  select(-.pred_class)
  # relocate(diagnosis)


ddxplus_train_folds

pred_prob1
# analysis_data <- 
#   ddxplus_train_folds |> 
#   filter(id == "Fold1") |> 
#   pull(splits) |>
#   _[[1]] |> 
#   analysis()

assessment_data <-
  ddxplus_train_folds |> 
  filter(id == "Fold1") |> 
  pull(splits) |>
  _[[1]] |> 
  assessment()

assessment_data <-
  ddxplus_split |> 
  testing()

test_data_prob_predict <- pred_prob1
test_data <- assessment_data
diff_prob_minimum <- 0.05
diff_rank_maximum <- 5
top_n_from_test_data <- 5

gather_differential_metrics <-
  function(test_data, test_data_prob_predict,
           diff_prob_minimum = 0.05, diff_rank_maximum = 5,
           top_n_from_test_data = 15,
           filter_rule = c("or","and")){
    # collect the metrics of predicted probabilities
    # Args:
    # 1. test_data: The test data set with each differential diagnosis 
    # having its own row. Also, evidences are pivoted wider.
    # 2. test_data_prob_predict: test data with probability predictions
    # All conditions in the differential diagnosis that meet any of the two conditions
    # below are kept
    # 3. diff_prob_minimum: differential probability of the condition greater than this value is kept
    # 4. diff_rank_maximum: differential rank is below this value
    # Output:
    # a df with the following columns
    # metric1: Proportions of diagnoses in the predicted differential diagnoses 
    # that are in the actual differential diagnoses
    # metric2: Proportions of diagnoses in the predicted differential diagnoses
    # that are not in the actual differential diagnoses
    # correlation: Pearson correlation coefficient of order/rank of probabilities
    # of the predicted differential diagnoses and the actual differential diagnoses
    
    metric1_computation_df <- 
      test_data %>% 
      select(patientId, 
             diagnosis,
             differential_case_weight) %>% 
      mutate(prob_rank = 
               rank(desc(differential_case_weight)),
             .by = patientId)
    
    # Limit diagnosis in test_data
    metric1_computation_df <-
      metric1_computation_df |> 
      filter(prob_rank <= top_n_from_test_data)
    
    if (nrow(test_data_prob_predict) == nrow(test_data)){
      pred_df <- test_data_prob_predict %>% 
        mutate(patientId = test_data$patientId) %>% 
        relocate(patientId) %>% 
        slice_head(n = 1, by = patientId) # only get 1 sample per patient
    } else {
      pred_df <- test_data_prob_predict
    }
    

    pred_df_longer <- pred_df %>%
      pivot_longer(-patientId,
                   names_to = "diagnosis",
                   values_to = "prob",
                   values_ptypes = numeric()) %>%
      filter(prob > 0) %>% 
      group_by(patientId) %>% 
      arrange(desc(prob), .by_group = TRUE) %>% 
      ungroup()

    
    # only keep diagnosis that has probability more than 5%
    # or if they are within top 5 diagnoses
    pred_df_longer <- pred_df_longer %>% 
      mutate(prob_rank = rank(-prob),
             .by = patientId)
    
    filter_rule <- filter_rule[1]
    if (filter_rule == "and"){
      pred_df_longer <- pred_df_longer |> 
        filter(prob > diff_prob_minimum | prob_rank <= diff_rank_maximum) 
    } else if (filter_rule == "or"){
      pred_df_longer <- pred_df_longer |> 
        filter(prob > diff_prob_minimum & prob_rank <= diff_rank_maximum) 
    } else {
      stop("choose filter_rule 'and' or 'or'")
    }
    
    
    # remove .pred_ from diagnosis
    pred_df_longer <- pred_df_longer %>% 
      mutate(diagnosis = str_remove(diagnosis,
                                    pattern = ".pred_"))
    
    
    # Checking metric #1
    
    matching_diagnosis_prop <-
      function(truth, estimate){
        
        matches <- intersect(truth, estimate)
        truth_diag_counnt <- length(truth)
        
        metric <- length(matches)/truth_diag_counnt
        
        return(metric)
        
      }
    
    patients_testSet <- unique(test_data$patientId)
    
    # require(furrr)
    # plan(multisession, workers = 4)
    matching_diag_prop_metric <- map_dbl(patients_testSet,
                                                .f = function(patientId){
                                                  
                                                  truth <- metric1_computation_df$diagnosis[
                                                    metric1_computation_df$patientId == patientId
                                                  ]
                                                  estimate <-
                                                    pred_df_longer$diagnosis[
                                                      pred_df_longer$patientId == patientId
                                                    ]
                                                  
                                                  matching_diagnosis_prop(truth, estimate)
                                                  
                                                }, .progress = TRUE)
    
    metric1_df <- tibble(patientId = patients_testSet,
                         metric1 = matching_diag_prop_metric)
    
    message(paste0("Message: Differential Diagnosis Recall is ",
                   mean(metric1_df$metric1) |> round(3),
                   "...(1/3)"))
    
    # Checking metric #2 
    
    # plan(multisession, workers = 8)
    setdiff_diag_prop_metric <-
      map_dbl(patients_testSet,
                     .f = function(patientId){
                       
                       truth <- metric1_computation_df$diagnosis[
                         metric1_computation_df$patientId == patientId
                       ]
                       estimate <-
                         pred_df_longer$diagnosis[
                           pred_df_longer$patientId == patientId
                         ]
                       
                       different_n <- length(setdiff(estimate, truth))
                       estimate_n <- length(estimate)
                       
                       output <- 1 - (different_n/estimate_n)
                       output
                       
                     }, .progress = TRUE)
    
    metric2_df <- tibble(patientId = patients_testSet,
                         metric2 = setdiff_diag_prop_metric)
    
    message(paste0("Message: Differential Diagnosis Precision is ",
                   mean(metric2_df$metric2) |> round(3),
                   "...(2/3)"))
    
    # Checking metric #3 
    
    truth_join_pred <- metric1_computation_df %>% 
      inner_join(pred_df_longer,
                 by = c("patientId", "diagnosis")) %>% 
      mutate(prob_rank.y = rank(desc(prob)),
             .by = patientId)
    rankCorrelation_diffDiag_metric_df <- truth_join_pred %>% 
      summarize(correlation =
                  if (length(prob_rank.x)  == 1){
                    warning("Warning: Rank correlation with only one diagnosis converted from NA to 1")
                    1
                  } else {
                    cor(prob_rank.x,
                        prob_rank.y)
                  }
                ,
                .by = patientId) |> 
      suppressWarnings()
    # filter(!is.na(correlation))

    
    metric3_df <- rankCorrelation_diffDiag_metric_df |> 
      mutate(correlation = ifelse(is.na(correlation), 0.5, correlation))
    
    # message("Message: Done with metric 3...(3/3)")
    message(paste0("Message: Differential Diagnosis Spearman Correlation is ",
                   mean(metric3_df$correlation) |> round(3),
                   "...(3/3)"))
    message("Message: Tidying metrics...")
    
    # list(metric_df = metric1_df %>% 
    #        left_join(metric2_df) %>% 
    #        left_join(metric3_df),
    #      predictions = pred_df_longer)
    
    result <- metric1_df %>% 
      left_join(metric2_df) %>% 
      left_join(metric3_df)
    result
    
  }

result <-
  gather_differential_metrics(test_data, test_data_prob_predict,
                              diff_prob_minimum = diff_prob_minimum, 
                              diff_rank_maximum = diff_rank_maximum,
                              top_n_from_test_data = top_n_from_test_data,
                              filter_rule = "and")

result |> 
  summarize(across(-patientId,
                   .fns = ~mean(., na.rm = TRUE)))



# Test evaluation metrics for different values
# Which hyperparameters performed best?

fit_metrics$.config
pred_df

hparam_results <- 
  map(fit_metrics$.config,
    .f = function(config){
      
      pred_prob1 <- 
        pred_prob |> 
        filter(.config == config,
               id == "Fold1") |> 
        select(matches(".pred"))  |> 
        select(-.pred_class)
      
      assessment_data <-
        ddxplus_train_folds |> 
        filter(id == "Fold1") |> 
        pull(splits) |>
        _[[1]] |> 
        assessment()
      
      test_data_prob_predict <- pred_prob1
      test_data <- assessment_data
      diff_prob_minimum <- 0.05
      diff_rank_maximum <- 15
      top_n_from_test_data <- 15
      
      result <-
        gather_differential_metrics(test_data, test_data_prob_predict,
                                    diff_prob_minimum = diff_prob_minimum, 
                                    diff_rank_maximum = diff_rank_maximum,
                                    top_n_from_test_data = top_n_from_test_data)
      
      result |> 
        summarize(across(-patientId,
                         .fns = ~mean(., na.rm = TRUE))) |> 
        mutate(config = config, .before = metric1)
      
    }, .progress = TRUE)


bind_rows(hparam_results)

# modeling_objects_board |> 
#   pin_write(bind_rows(hparam_results), 
#             name = "hparam_project_metrics",
#             title = "Random Forest Tuning Results for Project Metrics", 
#             description = "Tuning results with grid_latin_hypercube function
#             used for tuning grid search for DDR, DDP, and Spearman Correlation.",
#             versioned = T)

dd_rules_set <-
  expand_grid(
    # config = fit_metrics$.config[6:8],
    filter_rule = c("and","or"),
    diff_prob_minimum = seq(0.01, 0.07, 0.02),
    diff_rank_maximum = c(5,7,10,12,15),
    top_n_from_test_data = c(5,7,10,12,15),
  )

# limit to top 3 models with good project-defined metric scores
dd_rule_set_results <- 
  map(1:nrow(dd_rules_set),
      .f = function(idx){
        
        pred_prob1 <- 
          pred_prob |> 
          # filter(.config == dd_rules_set$config[idx],
          #        id == "Fold1") |> 
          select(matches(".pred")) |> 
          select(-.pred_class)
        
        # assessment_data <-
        #   ddxplus_train_folds |> 
        #   filter(id == "Fold1") |> 
        #   pull(splits) |>
        #   _[[1]] |> 
        #   assessment()
        
        assessment_data <-
          ddxplus_split |> 
          testing()
        
        test_data_prob_predict <- pred_prob1
        test_data <- assessment_data
        diff_prob_minimum <- 
          dd_rules_set$diff_prob_minimum[idx]
        diff_rank_maximum <- 
          dd_rules_set$diff_rank_maximum[idx]
        top_n_from_test_data <- 
          dd_rules_set$top_n_from_test_data[idx]
        filter_rule <-
          dd_rules_set$filter_rule[idx]
        
        result <-
          gather_differential_metrics(test_data, test_data_prob_predict,
                                      diff_prob_minimum = diff_prob_minimum, 
                                      diff_rank_maximum = diff_rank_maximum,
                                      top_n_from_test_data = top_n_from_test_data,
                                      filter_rule = filter_rule)
        
        output <- result |> 
          summarize(across(-patientId,
                           .fns = ~mean(., na.rm = TRUE)))
        
        message(paste0("Done with iteration ", idx,"\nTime: ", Sys.time()))
        
        output
        
      }, .progress = TRUE)

dd_rule_set_results
dd_rules_set_results <-
  dd_rules_set |> 
  bind_cols(bind_rows(dd_rule_set_results))

  

# readr::write_rds(bind_rows(dd_rule_set_results), 
#                  file = here("modeling_objects"))

modeling_objects_board |> 
  pin_write(
    dd_rules_set_results, 
    name = "dd_rules_testing_res",
    title = "Results of Testing Different Filter Rules for Differential Diagnosis", 
    description = "Results of Testing Different Filter Rules for Differential Diagnosis
    while checking for DDR, DDP, and Spearman Correlation. 
    We check both 'and' and 'or' filter_rule",
    versioned = T)



dd_rules_set_results <-
  modeling_objects_board |> 
  pin_read("dd_rules_testing_res")

dd_rules_set_results |> 
  rowwise() |> 
  mutate(mean_score = mean(c(metric1, metric2, correlation))) |>
  ungroup() |> 
  relocate(mean_score) |> 
  arrange(desc(mean_score))
  

# TODO 
# Let's experiment with different PCA tunings (specifically higher PCs)
# We can also try increasing the sample size from 200 patients per pathology



