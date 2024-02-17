

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
  R/pull
)

box::reload(pull)

# Functions ---------------------------------------------------------------

evidences <- pull$pull_evidences()
# pull$pull_training()
train <- arrow::open_dataset("data/processed/train.arrow", format = "arrow")

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


#' @export
unnest_evidences <-
  function(data){
    # Unnest the evidences from the EVIDENCES column
    # Each evidence will have its own column
    
    # Note: 
    # evidence_classes_df: a df containing evidence with their
    # supposed classes
    evidence_classes_df <-
      arrow::read_feather(
        file = here("configs","evidence_classes.arrow")
      )
    
    # data <- train_sample # TEST
    test <- data %>% 
      unnest_tokens(token = "regex", 
                    output = "evidence_raw",
                    input = EVIDENCES,
                    pattern = ', ',
                    to_lower = FALSE) %>% 
      mutate(evidence_raw = str_remove_all(evidence_raw,
                                           pattern = "\\[|\\'|\\]")) %>% 
      mutate(evidence = str_extract(evidence_raw,
                                    pattern = "[_[:alnum:]]*(?=_@_)"),
             evidence_value = str_extract(evidence_raw,
                                          pattern = "(?<=_@_).*"),
             evidence_value = ifelse(is.na(evidence_value),
                                     1,
                                     evidence_value),
             evidence = if_else(is.na(evidence),
                                evidence_raw,
                                evidence)) %>% 
      relocate(evidence_raw, evidence, evidence_value)
    
    
    test2 <- test %>% 
      select(-evidence_raw) %>% 
      distinct(evidence, patientId, .keep_all = TRUE) %>% #only get first example from multiple evidence variable
      pivot_wider(names_from = evidence,
                  values_from = evidence_value) 
    
    # Re-add missing evidence columns
    emp <- cook_empty_tibble(data = test2)
    
    out <- bind_cols(test2, emp)
    
    # suppose column classes
    walk(evidence_classes_df$evidence,
         .f = function(colname){
           class(out[[colname]]) <<- 
             evidence_classes_df$class[evidence_classes_df$evidence == colname]
         }
    )
    
    out %>% 
      relocate(all_of(evidence_classes_df$evidence),
               .after = INITIAL_EVIDENCE)
    
  }


#### TESTING ####


train_sample <- train %>% 
  slice_sample(n = 100) %>% 
  collect()


test <- unnest_evidences(train_sample)
names(test)

evidence_classes_df
evidences


# There are 676 possible values for all multiclass evidences

multiclass_evidences <- evidences %>% 
  filter(data_type == "M") %>% 
  unnest_longer(col = c(possible_values, value_meaning)) %>% 
  mutate(name = paste0(name,"_@_",possible_values)) %>% 
  select(name, data_type, class)

singleclass_evidences <- evidences %>% 
  filter(data_type %in% c("B","C")) %>% 
  select(name, data_type, class)

all_evidence_classes <- 
  bind_rows(multiclass_evidences,singleclass_evidences)



prep1 <- train_sample %>% 
  slice(1:2) %>% 
  unnest_tokens(token = "regex", 
                output = "evidence_raw",
                input = EVIDENCES,
                pattern = ', ',
                to_lower = FALSE) %>% 
  select(patientId, evidence_raw) %>% 
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

prep2 <- prep1 %>% 
  mutate(
    evidence_value = suppressWarnings(as.numeric(evidence_value)),
    evidence_value = ifelse(is.na(evidence_value), 1, evidence_value)
    ) %>% 
  pivot_wider(names_from = evidence_raw, 
              values_from = evidence_value, 
              values_fill = 0)

cook_empty_tibble(prep2)


#' @export
unnest_evidences <-
  function(data){
    
    # data <- train_sample
    
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
    
    prep2 <- prep1 %>% 
      mutate(
        evidence_value = suppressWarnings(as.numeric(evidence_value)),
        evidence_value = ifelse(is.na(evidence_value), 1, evidence_value)
      ) %>% 
      pivot_wider(names_from = evidence_raw, 
                  values_from = evidence_value, 
                  values_fill = 0)
    
    # Re-add missing evidence columns
    emp <- cook_empty_tibble(data = prep2)
    
    out <- bind_cols(prep2, emp)
    
    out
    
  }


prepped <- unnest_evidences(train_sample)
map_chr(prepped, class)
# ert