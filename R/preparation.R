

# Description -------------------------------------------------------------

# The goal of this script is for feature engineering. We prepare the features
# to be used for model training. 


# Setup -------------------------------------------------------------------

box::use(
  magrittr[...],
  arrow[...],
  stringr[...],
  dplyr[...],
  here[here],
  tidytext[unnest_tokens],
  tidyr[...],
  purrr[...],
  tidyselect[...],
  furrr[...],
  future[...],
  R/pull
)

# box::reload(pull)

# Functions ---------------------------------------------------------------

#' Helper function for unnest_evidences function
#' 
#' @description
#' helper function for unnest_evidences
#' create an empty tibble to be binded with the df with unnested evidences
#' this is to ensure that in the new df, all evidences are accounted for
#' 
#' @param data the output from unnest_evidences
#' @param json the json vector
#' @param code_names the evidences
#' @examples
#' extract_vector(x)
#' @return the modified number
#' @export
cook_empty_tibble <-
  function(data){
    # helper function for unnest_evidences
    # create an empty tibble to be binded with the df with unnested evidences
    # this is to ensure that in the new df, all evidences are accounted for
    # Args:
    # 1. data: the output from unnest_evidences
    
    # Note: 
    # evidence_classes_df: a df containing evidence with their
    # supposed classes
    
    evidence_classes_df <-
      arrow::read_feather(
        file = here("configs","evidence_classes.arrow")
      )
    
    missing_cols <- setdiff(evidence_classes_df$name,
                            colnames(data))
    
    empty_matrix <- matrix(0,
                           ncol = length(missing_cols),
                           nrow = nrow(data))
    colnames(empty_matrix) <- missing_cols
    empty_tibble <- as_tibble(empty_matrix)
    
    empty_tibble
    
  }




#' Mutate patient ID column
#' 
#' @description
#' Although not a necessary step for training. This is 
#' helpful for model validation
#' 
#' @param data the training dataset
#' @examples
#' add_patient_ID(training)
#' @return the training dataset with a new column "patientId"
#' @export
add_patient_ID <-
  function(data){
    
    # Add patient Id
    output <- data %>% 
      mutate(patientId = 1:nrow(.)) %>% 
      relocate(patientId)
    
    output
  }


#' Unnest the initial evidences column from the DDXPlus dataset
#' 
#' @description
#' Unnest the INITIAL_EVIDENCE column from the DDXPlus dataset.
#' This will "pivot wider" each evidence from the single column
#' to multiple columns. Each binary evidence will have its own column.
#' Presence of binary evidence is represented by numeric "1". 
#' For categorical evidence, each will also have its own column,
#' but this time, the values will depend on the possible values
#' of the evidence which may either be numeric or character. 
#' Lastly, for multiclass evidences, all possible values for each 
#' evidence will have its own column. For example, 'E_55' will have
#' columns E_55_@_V_123, E_55_@_V_14, and so on.
#' 
#' @param data the training dataset
#' @examples
#' unnest_evidences(training)
#' @return the training dataset with unnested evidences
#' @export
unnest_evidences <-
  function(data){
    
    # data <- train_sample
    
    evidence_classes_df <-
      arrow::read_feather(
        file = here("configs","evidence_classes.arrow")
      )
    evidences <- pull$pull_evidences(pull_from_raw = FALSE)
    
    prep1 <- data %>% 
      unnest_tokens(token = "regex", 
                    output = "evidence_raw",
                    input = EVIDENCES,
                    pattern = ', ',
                    to_lower = FALSE) %>% 
      mutate(
        # Clean the evidences. Remove "[", "'", and "]"
        evidence_raw = 
          str_remove_all(evidence_raw,
                         pattern = "\\[|\\'|\\]"),
        # Extract the value after the evidence e.g. E_12_@_[value]
        evidence_value = str_extract(evidence_raw,
                                     pattern = "(?<=_@_).*"),
        # re-convert all NA to 1. Since binary variables
        # from previous step were converted to NA
        evidence_value = ifelse(is.na(evidence_value),
                                1,
                                evidence_value)) 
    
    # Get data type of all evidences. 
    # For categorical variables, only keep 
    # evidence name. So during pivot wider,
    # only evidence name is kept as the column.
    # For multiclass variables, keep the value
    # within the evidence name so that each value 
    # will be given its own column
    prep1.5 <- prep1 %>% 
      mutate(evidence_name = str_extract(evidence_raw,
                                         pattern = "E_\\d+")) %>% 
      left_join(evidences %>% 
                  select(name, data_type),
                by = join_by(evidence_name == name)) %>% 
      mutate(
        evidence_raw = ifelse(data_type == "C",
                              evidence_name,
                              evidence_raw),
        evidence_value = ifelse(data_type == "M",
                                1,
                                evidence_value)
      ) %>% 
      select(-evidence_name, -data_type)

    
    prep_test <- prep1.5 %>% 
      pivot_wider(names_from = evidence_raw, 
                  values_from = evidence_value)

    profile_cols <- str_subset(colnames(prep_test),
                               pattern = "E_",
                               negate = TRUE)
    prep_profile <- prep_test[profile_cols]
    
    prep_numeric_col <- prep_test %>% 
      select(any_of(
        with(evidence_classes_df, 
             name[data_type %in% c("M","B") | class == "numeric"]))) %>% 
      mutate(across(everything(), as.numeric))
    prep_numeric_col[is.na(prep_numeric_col)] <- 0
    
    prep_char_col <-
      prep_test %>% 
      select(any_of(
        with(evidence_classes_df, 
             name[data_type == "C" & class == "character"])))
    prep_char_col[is.na(prep_char_col)] <- "NA"
    
    
    prep2 <- prep_profile %>% 
      bind_cols(prep_numeric_col) %>% 
      bind_cols(prep_char_col)
    
    # Re-add missing evidence columns
    emp <- cook_empty_tibble(data = prep2)
    
    out <- bind_cols(prep2, emp)
      
    
    # suppose column classes
    walk(evidence_classes_df$name,
         .f = function(colname){
           class(out[[colname]]) <<- 
             evidence_classes_df$class[evidence_classes_df$name == colname]
         }, .progress = TRUE
    )
    
    out
    
  }

