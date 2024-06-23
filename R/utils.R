
box::use(
  stats[...],
  dplyr[...],
  tibble[...],
  caret[...],
  R/pull
)





# Functions ---------------------------------------------------------------

#' Load the dataset to memory
#' 
#' @param sample.size A numeric. Indicate a sample size if you want to 
#' only load a sample dataset to memory. Put NULL if load the whole dataset.
#' @param seed A numeric. Set seed.
#' @param train_arrow_object An ArrowObject. an output from
#' preparation$unnest_differential
#'
#' @export
load_df_to_memory <- 
  function(sample.size = NULL, train_arrow_object,
           seed = 123){
    
    # Get only a sample of patients for initial run
    if (!is.null(sample.size)){
      
      sample_patients <-
        train_arrow_object %>% 
        distinct(patientId) %>% 
        pull(patientId, as_vector = TRUE)
      set.seed(seed)
      sample_patients <- 
        sample_patients %>% 
        sample(size = sample.size)
      
      train_arrow_object <- 
        train_arrow_object %>% 
        filter(patientId %in% sample_patients)
      
    }
    
    
    output <- 
      train_arrow_object %>% 
      as_tibble() 
    
    output
    
  }



#' @export
get_dummy_categories_from_pred <-
  function(predictors){
    
    
    evidences <- 
      pull$pull_evidences(pull_from_raw = FALSE)
    
    get_category_evidence_levels <-
      function(evidence_code){
        
        # evidence_code <- "E_130"
        
        output <- 
          evidences %>% 
          dplyr::filter(name == evidence_code) 
        
        # rlang::inform(output$question_en)
        
        list(value_meaning =
               output %>% 
               pull(value_meaning) %>% 
               unlist(),
             possible_values = 
               output %>% 
               pull(possible_values) %>% 
               unlist()
        )
      }
    
    
    # We create a dummy table just to make sure we capture all 
    # levels from the categorical evidences
    
    e130 <- 
      get_category_evidence_levels("E_130")$possible_values %>% 
      c(replicate(6, "V_11")) %>% 
      factor()
    # levels(e130)[c(1,2)] <- levels(e130)[c(2,1)]
    
    e135 <- 
      get_category_evidence_levels("E_135")$possible_values %>% 
      c(replicate(10, "V_10")) %>% 
      factor()
    
    e131 <- 
      get_category_evidence_levels("E_131")$possible_values %>% 
      c(replicate(10, "V_10")) %>% 
      factor()
    
    e204 <- 
      get_category_evidence_levels("E_204")$possible_values %>% 
      factor()
    
    # Create a dataframe including all levels
    
    dummy_factor_df <- 
      tibble(
        E_130 =e130,
        E_135 = e135,
        E_131 = e131,
        E_204 = e204
      )
    
    
    # Prepare factor predictors
    pred_fct <- 
      predictors %>% 
      select(-where(is.numeric))
    
    pred_fct2 <- 
      bind_rows(
        dummy_factor_df,
        pred_fct
      )
    
    pred_fct2 <- 
      pred_fct2 %>% 
      mutate(
        across(everything(),
               .fns = as.character)
      ) %>% 
      mutate(
        E_130 = ifelse(E_130 == "V_11",
                       "NA",
                       E_130),
        E_135 = ifelse(E_135 == "V_10",
                       "NA",
                       E_135),
        E_131 = ifelse(E_131 == "V_10",
                       "NA",
                       E_131),
        E_204 = ifelse(E_204 == "V_10",
                       "NA",
                       E_204),
      ) %>% 
      mutate(
        across(everything(),
               .fns = as.factor)
      )
    
    
    
    # dummify the data
    dmy <- dummyVars(" ~ .", data = pred_fct2, fullRank = TRUE)
    pred_fct_dummified <- 
      data.frame(predict(dmy, newdata = pred_fct2)) %>% 
      as_tibble()
    
    # Remove the dummy rowss
    pred_fct_dummified <- 
      pred_fct_dummified %>% 
      slice(-1:-12)
    
    pred_fct_dummified
    
    
  }


