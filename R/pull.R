


# Overview ----------------------------------------------------------------

# This script contains functions for pulling required datasets. 
# Also includes minor pre-processing to convert files to table.


# Setup -----------------------------------------------------------------

box::use(
  here[here],
  readr[read_csv],
  arrow[...],
  dplyr[...],
  magrittr[...],
  purrr[...],
  tidyr[...],
  jsonlite[...],
  stringr[...]
)

processed_path <- here("data","processed")
input_path <- here("data","input")


# Functions ---------------------------------------------------------------




#' @export
pull_training <- 
  function(pull_from_raw = FALSE){
    
    if (pull_from_raw == TRUE){
      release_train_patients_raw <- 
        read_csv(file  = here(input_path,
                              "release_train_patients"))
      
      arrow::write_feather(x = release_train_patients_raw, 
                           sink = here(processed_path, 
                                       "train.arrow"))
      release_train_patients <- 
        arrow::open_dataset("data/processed/train.arrow", format = "arrow")
      
      message(
        paste0("Saved to ", 
               here(processed_path, 
                    "train.arrow"))
      )
    } else {
      
      release_train_patients <- 
        arrow::open_dataset("data/processed/train.arrow", format = "arrow")
      
    }
    
    return(release_train_patients)
    
  }



#' @export
pull_evidences <-
  function(pull_from_raw = FALSE){
    
    if (pull_from_raw == TRUE){
      
      evidences_json <- jsonlite::read_json(path = here(input_path, 
                                                        "release_evidences.json"))

      evidence_classes <- arrow::read_feather(
        file = here("configs",
                    "evidence_classes.arrow")
      )
      
      evidences <- tibble(evidences_json) %>% 
        hoist(evidences_json, 
              "question_en",
              "name", 
              "data_type",
              "is_antecedent",
              "possible-values",
              "value_meaning") %>% 
        select(-evidences_json) %>% 
        left_join(evidence_classes)
      
      pos_values_unnested <- 
        map(evidences$`possible-values`,
            .f = function(x){
              output <- unlist(x, recursive = TRUE)
                output %>% 
                  as.character()
            })

      value_meaning_unnested <- map(evidences$value_meaning,
          .f = function(x){
            output <- unlist(x, recursive = FALSE)
            if (!is.null(output)){
              output <- output[str_detect(names(output), 
                                ".en")]
              unname(output)
            } else character(0)
          })
      
      evidences_output <- evidences %>% 
        mutate(value_meaning = value_meaning_unnested,
               possible_values = pos_values_unnested) %>%
        rowwise() %>% 
        mutate(value_meaning = list(unlist(value_meaning)),
               possible_values = list(unlist(possible_values))) %>% 
        ungroup() %>% 
        select(-`possible-values`) %>% 
        relocate(possible_values, .before = value_meaning)


      arrow::write_feather(x = evidences_output, 
                           sink = here(processed_path, 
                                       "evidences.arrow"))
      message(
        paste0("Saved to ", 
               here(processed_path, 
                    "evidences.arrow"))
      )
    } else {
      
      evidences_output <-
        arrow::read_feather(file = here(
          processed_path,
          "evidences.arrow"
        ))
      
    }
    
    return(evidences_output)
    
  }



#' @export
pull_conditions <-
  function(pull_from_raw = FALSE){
    
    if (pull_from_raw == TRUE){
      conditions_json <- 
        jsonlite::read_json(
          path = here(input_path, 
                      "release_conditions.json"))
      
      # unnest columns
      conditions <- tibble(conditions_json) %>% 
        hoist(conditions_json, 
              "condition_name",
              "icd10-id",
              "symptoms",
              "antecedents",
              "severity") %>% 
        select(-conditions_json)
      
      # get values from symptoms and antecedents
      conditions <- conditions %>% 
        group_by(condition_name) %>% 
        mutate(
          symptoms = list(names(symptoms[[1]])),
          antecedents = list(names(antecedents[[1]]))
        ) %>% 
        ungroup()
      
      arrow::write_feather(x = conditions, 
                           sink = here(processed_path, 
                                       "conditions.arrow"))
      message(
        paste0("Saved to ", 
               here(processed_path, 
                    "conditions.arrow"))
      )
    } else {
      
      conditions <- 
        arrow::read_feather(file = here(
          "data",
          "processed",
          "conditions.arrow"
        ))
      
    }
    
    
    return(conditions)
    
  }

