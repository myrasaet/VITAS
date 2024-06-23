

# Setup -------------------------------------------------------------------

box::use(
  magrittr[...],
  dplyr[...],
  tidyr[...],
  stringr[...],
  purrr[...],
  future[...],
  furrr[...],
  R/pull
)


# Functions ---------------------------------------------------------------

#' Manually check the actual and predicted differential diagnosis
#' 
#' @param idx A numeric. An index from 1 to number of unique patients
#' in the assessment dataset
#' @param predictions_tidy A list. This is a list that contains the following
#' elements:
#' id - a dataframe containing the patientId column
#' predictors - a dataframe containing all the predictors. AGE, evidences, SEX
#' outcome - a dataframe with 49 columns. Each is a diagnosis. The values
#' are the actual probabilities or differential scores
#' predicted_outcome -  a dataframe with 49 columns. Each is a diagnosis. The
#' values are the predicted probabilities.
#' @param differential_prob_threshold_actual A numeric. Is the filtering rule for 
#' the thresholds. ex. `filter(differential_prob > differential_prob_threshold)`
#' @param differential_prob_threshold_pred A numeric. Is the filtering rule for 
#' the thresholds. ex. `filter(differential_prob > differential_prob_threshold)`
#' @param patientId a character. provider the patientID if want to specify 
#' specific patientId.
#' @param severity_threshold a numeric. if a predicted diagnosis has severity
#' higher than this threhold, then override differential_prob_threshold_pred 
#' and keep the diagnosis. Note: The lower the number, the more severe the 
#' disease.
#' @param high_severity_ddprob_threshold a numeric. The dd prob threshold for
#' highly sever diseases
#' 
#' @export
manual_evaluate <-
  function(predictions_tidy,
           idx = sample(size = 1,
                        x = seq_along(unique(predictions_tidy$id$patientId))),
           patientId = NULL,
           differential_prob_threshold_actual = 0.01,
           differential_prob_threshold_pred = 0.01,
           severity_threshold = 2,
           high_severity_ddprob_threshold = 0.001){
    
    
    # idx <- 4
    if (is.null(patientId)){
      filter_idx <- 
        predictions_tidy$id %>% 
        slice(idx) %>% 
        pull(patientId)
    } else if (!is.null(patientId)){
      filter_idx <- patientId
    }
    
    
    conditions <- 
      pull$pull_conditions() %>% 
      select(condition_name, severity)
    
    pivot_longer_dd <- 
      function(x = c("outcome","predicted_outcome")){
        
        out <- 
          predictions_tidy$id %>% 
          bind_cols(predictions_tidy[[x]]) %>% 
          filter(patientId == filter_idx) %>% 
          pivot_longer(c(-patientId, -PATHOLOGY),
                       names_to = "diagnosis",
                       values_to = "differential_prob") %>% 
          group_by(patientId) %>% 
          arrange(desc(differential_prob), .by_group = TRUE) %>% 
          ungroup() %>% 
          # Add the severity
          mutate(diagnosis = str_remove(diagnosis,
                                         pattern = "diag_")) %>% 
          left_join(conditions, 
                    by = join_by(diagnosis == condition_name))
        
        if (x == "predicted_outcome"){
          out <- 
            out %>% 
            filter(differential_prob >= differential_prob_threshold_pred |
                     severity <= severity_threshold) %>% 
            filter(differential_prob >= high_severity_ddprob_threshold)
        } else if (x == "outcome"){
          
          out <- 
            out %>% 
            filter(differential_prob >= differential_prob_threshold_actual)
          
        }
       
        out
        
      }
    
    
    rlang::inform("Actual data")
    .assessment_data <-
      pivot_longer_dd("outcome")
    
    .assessment_data %>% 
      print()
    
    rlang::inform("\n\nPredicted Data")
    .predictions <-
      pivot_longer_dd("predicted_outcome")
    print(.predictions)
    
  }


#' Pivot longer the assessment data and predictions
#'
#' @param predictions_tidy A list. This is a list that contains the following
#' elements:
#' id - a dataframe containing the patientId column
#' predictors - a dataframe containing all the predictors. AGE, evidences, SEX
#' outcome - a dataframe with 49 columns. Each is a diagnosis. The values
#' are the actual probabilities or differential scores
#' predicted_outcome -  a dataframe with 49 columns. Each is a diagnosis. The
#' values are the predicted probabilities.
#' @param differential_prob_threshold_actual A numeric. Is the filtering rule for 
#' the thresholds. ex. `filter(differential_prob > differential_prob_threshold)`
#' @param differential_prob_threshold_pred A numeric. Is the filtering rule for 
#' the thresholds. ex. `filter(differential_prob > differential_prob_threshold)`
#' @param severity_threshold a numeric. if a predicted diagnosis has severity
#' higher than this threhold, then override differential_prob_threshold_pred 
#' and keep the diagnosis. Note: The lower the number, the more severe the 
#' disease.
#' @param high_severity_ddprob_threshold a numeric. The dd prob threshold for
#' highly sever diseases
#' 
#' @export
pivot_longer_dd <- 
  function(predictions_tidy,
           differential_prob_threshold_actual = 0.01,
           differential_prob_threshold_pred = 0.01,
           severity_threshold = 2,
           high_severity_ddprob_threshold = 0.001){
    
    
    
    conditions <- 
      pull$pull_conditions() %>% 
      select(condition_name, severity)
    
    pivot_longer_dd <- 
      function(x = c("outcome","predicted_outcome")){
        
        out <- 
          predictions_tidy$id %>% 
          select(patientId) %>% 
          bind_cols(predictions_tidy[[x]]) %>% 
          pivot_longer(-patientId,
                       names_to = "diagnosis",
                       values_to = "differential_prob") %>% 
          group_by(patientId) %>% 
          arrange(desc(differential_prob), .by_group = TRUE) %>% 
          ungroup() %>% 
          # Add the severity
          mutate(diagnosis = str_remove(diagnosis,
                                        pattern = "diag_")) %>% 
          left_join(conditions, 
                    by = join_by(diagnosis == condition_name))
        
        if (x == "predicted_outcome"){
          out <- 
            out %>% 
            filter(differential_prob >= differential_prob_threshold_pred |
                     severity <= severity_threshold) %>% 
            filter(differential_prob >= high_severity_ddprob_threshold)
        } else if (x == "outcome"){
          
          out <- 
            out %>% 
            filter(differential_prob >= differential_prob_threshold_actual)
          
        }
        
        out
        
      }
    
    
    .assessment_data <-
      pivot_longer_dd("outcome")
    .predictions <-
      pivot_longer_dd("predicted_outcome")
    
    
    list(
      actual_pathology = predictions_tidy$id,
      assessment_df = .assessment_data,
      prediction_df = .predictions
    )
    
  }



#' Compute the Metrics
#' 
#' @param prediction_evaluation_list A list. The output from 
#' pivot_longer_dd.
#' 
#' @export
compute_differential_metrics <-
  function(prediction_evaluation_list,
           workers = 1,
           save_path
           ){
    
    # workers <- 4
    
    if (workers > 1){
      plan(multisession, workers = workers)
      mapp_func <- future_map_dbl
    } else if (workers == 1){
      
      mapp_func <- map_dbl
      
    }
    
    # Setup
    
    .actual <- 
      prediction_evaluation_list$assessment_df
    .prediction <- 
      prediction_evaluation_list$prediction_df
    
    pathology_per_patient <- 
      prediction_evaluation_list$actual_pathology
    
    patient_set <- 
      unique(pathology_per_patient$patientId)
    
    ## Computation of Metrics ##
    
    # 1. Compute GTPA@1
    rlang::inform("Computing Metric 1: GTPA@1")
    gtpa1 <- 
      mapp_func(patient_set,
              .f = function(x){
                
                # x <- patient_set[1]
                
                .pathology <- 
                  pathology_per_patient %>% 
                  filter(patientId == x) %>% 
                  dplyr::pull(PATHOLOGY)
                
                .pred_pathology <- 
                  .prediction %>% 
                  filter(patientId == x) %>% 
                  slice_max(differential_prob) %>% 
                  dplyr::pull(diagnosis)
                
                as.numeric(.pathology == .pred_pathology)
                
              },
              .progress = TRUE
      )
    readr::write_rds(gtpa1, paste0(save_path,"gtpa1.rds"))
    
    # 2. Compute GTPA
    rlang::inform("Computing Metric 2: GTPA")
    gtpa <- 
      mapp_func(patient_set,
              .f = function(x){
                
                # x <- patient_set[1]
                
                .pathology <- 
                  pathology_per_patient %>% 
                  filter(patientId == x) %>% 
                  dplyr::pull(PATHOLOGY)
                
                .pred_dd <- 
                  .prediction %>% 
                  filter(patientId == x) %>% 
                  dplyr::pull(diagnosis)
                
                as.numeric(.pathology %in% .pred_dd)
                
                
              },
              .progress = TRUE
      )
    readr::write_rds(gtpa, paste0(save_path,"gtpa.rds"))
    
    
    # 3. Compute the Differential Diagnosis Recall (DDR)
    rlang::inform("Computing Metric 3: DDR")
    compute_ddr <-
      function(actual_dd, pred_dd){
        
        matches <- intersect(actual_dd, pred_dd)
        denominator <- length(actual_dd)
        
        metric <- length(matches)/denominator
        
        return(metric)
        
      }
    
    ddr <- 
      mapp_func(patient_set,
              .f = function(x){
                
                # x <- patient_set[1]
                
                .actual_dd <- 
                  .actual %>% 
                  filter(patientId == x) %>% 
                  dplyr::pull(diagnosis)
                
                .pred_dd <- 
                  .prediction %>% 
                  filter(patientId == x) %>% 
                  dplyr::pull(diagnosis)
                
                compute_ddr(actual_dd = .actual_dd, 
                            pred_dd = .pred_dd)
                
                
              },
              .progress = TRUE
      )
    readr::write_rds(ddr, paste0(save_path,"ddr.rds"))
    
    # 4. Compute the Differential Diagnosis Precision (DDP)
    rlang::inform("Computing Metric 4: DDP")
    compute_ddp <-
      function(actual_dd, pred_dd){
        
        matches <- intersect(actual_dd, pred_dd)
        denominator <- length(pred_dd)
        
        metric <- length(matches)/denominator
        
        return(metric)
        
      }
    
    ddp <- 
      mapp_func(patient_set,
              .f = function(x){
                
                # x <- patient_set[1]
                
                .actual_dd <- 
                  .actual %>% 
                  filter(patientId == x) %>% 
                  dplyr::pull(diagnosis)
                
                .pred_dd <- 
                  .prediction %>% 
                  filter(patientId == x) %>% 
                  dplyr::pull(diagnosis)
                
                compute_ddp(actual_dd = .actual_dd, 
                            pred_dd = .pred_dd)
                
                
              },
              .progress = TRUE
      )
    readr::write_rds(ddp, paste0(save_path,"ddp.rds"))
    
    rlang::inform("Computing Metric 5: Spearman Correlation")
    # 5. Compute Spearman Correlation
    spear_corr <- 
      mapp_func(patient_set,
              .f = function(x){
                
                # x <- patient_set[300]
                
                .actual_dd <- 
                  .actual %>% 
                  filter(patientId == x)
                
                
                .pred_dd <- 
                  .prediction %>% 
                  filter(patientId == x)
                
                res <- 
                  .actual_dd %>% 
                  inner_join(.pred_dd,
                             by = join_by(patientId, diagnosis)) %>% 
                  mutate(actual_dd_rank = rank(desc(differential_prob.x)),
                         pred_dd_rank = rank(desc(differential_prob.y)))
                
                output <- 
                  if (nrow(res) > 1){
                  stats::cor(res$actual_dd_rank, 
                             res$pred_dd_rank, 
                             method = "pearson")
                } else if (nrow(res) == 1){
                  1
                } else if (nrow(res) == 0){
                  
                  NA
                }
                
                output
                
                
                
              },
              .progress = TRUE
      )
    readr::write_rds(spear_corr, paste0(save_path,"spear_corr.rds"))
    
    
    
    raw_scores <- 
      bind_cols(patientId = patient_set,
                gtpa1 = gtpa1, 
                gtpa = gtpa,
                ddr = ddr,
                ddp = ddp,
                spear_corr = spear_corr)
    
    metrics_summary <- 
      raw_scores %>% 
      select(-patientId) %>% 
      summarize(
        across(everything(),
               .fns = ~mean(., na.rm = TRUE))
      )
    
    metrics <- 
      list(
      raw_scores = raw_scores,
      metrics_summary = metrics_summary
    )
    
      
    readr::write_rds(metrics, paste0(save_path, "metrics.rds"))
    rlang::inform("saved metrics to ",save_path)
    
    metrics
    
    
  }
