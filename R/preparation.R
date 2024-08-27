

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
  dtplyr[...],
  R/pull,
  R/utils
)

# box::reload(pull)

# Functions ---------------------------------------------------------------


# 1. Unnest Evidences --------------------------------------------------------

#' Get top 10 most frequent values from multiclass evidences
#' 
#' #' @description
#' Get top 10 most frequent values from multiclass evidences
#' 
#' @export
get_top10_multiclass_evidences <-
  function(train){
    
    to_get_top10 <- c("E_55","E_57","E_133","E_152")
    out <- map(to_get_top10,
         .f = function(x){
           
           sums <- train |> 
             select(contains(x)) |> 
             summarize(across(everything(),
                              .fns = ~sum(., na.rm = TRUE))) |> 
             collect()
           
           top10 <- sums |> 
             pivot_longer(everything()) |> 
             arrange(desc(value)) |> 
             slice(1:10) |> # get top 10
             pull(name)
           
           top10
          
         }, .progress = TRUE)
    
    multiclass_top10 <- out
    names(multiclass_top10) <- c("E_55",
                                 "E_57","E_133","E_152")
    readr::write_rds(multiclass_top10,
                     here("data","processed",
                          "multiclass_top10.rds"))
    message("saved in data/processed/multiclass_top10.rds")
  }


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
    multiclass_top10 <- 
      readr::read_rds(here("data","processed",
                           "multiclass_top10.rds"))
    
    # In the following, we reduce the possible values
    # for the multiclass evidences. Convert infrequent
    # values to "Others"
    evidence_classes_df <- evidence_classes_df |> 
      mutate(
        name = case_when(
          str_detect(name, "E_55") ~ ifelse(name %in% 
                                              multiclass_top10$E_55,
                                            name, "E_55_@_Others"),
          str_detect(name, "E_57") ~ ifelse(name %in% 
                                              multiclass_top10$E_57,
                                            name, "E_57_@_Others"),
          str_detect(name, "E_133") ~ ifelse(name %in% 
                                               multiclass_top10$E_133,
                                             name, "E_133_@_Others"),
          str_detect(name, "E_152") ~ ifelse(name %in% 
                                               multiclass_top10$E_152,
                                             name, "E_152_@_Others"),
          .default = name
        )
      ) |> 
      distinct()
    
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
    # data <- train
    # data <- d[[1]]
    
    evidence_classes_df <-
      arrow::read_feather(
        file = here("configs","evidence_classes.arrow")
      )
    evidences <- pull$pull_evidences(pull_from_raw = FALSE)
    
    multiclass_top10 <- 
      readr::read_rds(here("data","processed",
                           "multiclass_top10.rds"))
    
    
    prep1 <- data %>% 
      unnest_tokens(token = "regex", 
                    output = "evidence_raw",
                    input = EVIDENCES,
                    pattern = ', ',
                    to_lower = FALSE)
    
    prep1_lazy <- lazy_dt(prep1)
    # prep1_lazy <- prep1 # FOR TESTING
    
    prep1.1 <- prep1_lazy %>% 
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
                                as.character(1),
                                evidence_value)) 
    
    # Get data type of all evidences. 
    # For categorical variables, only keep 
    # evidence name. So during pivot wider,
    # only evidence name is kept as the column.
    # For multiclass variables, keep the value
    # within the evidence name so that each value 
    # will be given its own column
    
    
    to_join <- evidences %>% 
      select(name, data_type) %>% 
      rename(evidence_name = name)
    
    prep1.5 <- prep1.1 %>% 
      mutate(evidence_name = str_extract(evidence_raw,
                                         pattern = "E_\\d+")) %>% 
      as_tibble() %>% 
      left_join(y = to_join,
                by = "evidence_name") %>% 
      lazy_dt() %>%
      mutate(
        evidence_raw = ifelse(data_type == "C",
                              evidence_name,
                              evidence_raw),
        evidence_value = ifelse(data_type == "M",
                                as.character(1),
                                evidence_value)
      )
    
    
    multiclass_top10
    # TEST
    multiclass_rows <- prep1.5 |> 
      filter(evidence_name %in% names(multiclass_top10)) |> 
      mutate(
        evidence_raw = case_when(
          evidence_name == "E_55" ~ ifelse(evidence_raw %in% 
                                             multiclass_top10$E_55,
                                           evidence_raw, "E_55_@_Others"),
          evidence_name == "E_57" ~ ifelse(evidence_raw %in% 
                                             multiclass_top10$E_57,
                                           evidence_raw, "E_57_@_Others"),
          evidence_name == "E_133" ~ ifelse(evidence_raw %in% 
                                              multiclass_top10$E_133,
                                            evidence_raw, "E_133_@_Others"),
          evidence_name == "E_152" ~ ifelse(evidence_raw %in% 
                                              multiclass_top10$E_152,
                                            evidence_raw, "E_152_@_Others")
        )
      ) |> 
      group_by(across(-evidence_value)) |> 
      summarize(evidence_value = sum(as.numeric(evidence_value),
                                     na.rm = TRUE)) |> 
      ungroup() |> 
      mutate(evidence_value = as.character(evidence_value))
    
    prep1.6 <- prep1.5 |> 
      filter(! evidence_name %in% names(multiclass_top10)) |> 
      as_tibble() |> 
      bind_rows(as_tibble(multiclass_rows)) |> 
      lazy_dt() |> 
      select(-evidence_name, -data_type) |> 
      arrange(patientId)
    
    prep_test <- prep1.6 %>% 
      pivot_wider(names_from = evidence_raw, 
                  values_from = evidence_value)
    
    prep_test_tibble <- as_tibble(prep_test)
    
    profile_cols <- str_subset(colnames(prep_test_tibble),
                               pattern = "E_",
                               negate = TRUE)
    prep_profile <- prep_test_tibble[profile_cols]
    
    # Convert all numeric NA to 0
    prep_numeric_col <- prep_test %>% 
      select(any_of(
        with(evidence_classes_df, 
             name[data_type %in% c("M","B") | class == "numeric"]))) %>% 
      mutate(across(everything(), as.numeric))
    # prep_numeric_col[is.na(prep_numeric_col)] <- 0
    prep_numeric_col <-
      prep_numeric_col %>% 
      as_tibble() |> 
      mutate(across(everything(), ~replace(., is.na(.), 0)))
    
    prep_char_col <-
      prep_test %>% 
      select(any_of(
        with(evidence_classes_df, 
             name[data_type == "C" & class == "character"]))) %>% 
      as_tibble() |> 
      mutate(across(everything(), ~replace(., is.na(.), "NA")))
    # prep_char_col[is.na(prep_char_col)] <- "NA"
    
    
    prep2 <- prep_profile %>% 
      bind_cols(as_tibble(prep_numeric_col))
    
    if (ncol(prep_char_col) > 0){
      rlang::inform("Add character evidences...")
      prep2 <- 
        prep2 %>% 
        bind_cols(as_tibble(prep_char_col)) 
    } else {
      
      # matrix(data = 0, 
      #        nrow = nrow(prep2),
      #        ncol = evidence_classes_df %>% 
      #                        filter(data_type == "C") %>% 
      #                        nrow())
      # 
    }

    
    # Re-add missing evidence columns
    emp <- cook_empty_tibble(data = prep2)
    
    out <- bind_cols(prep2, emp)
    
    out_col_classes <- map_chr(names(out), 
                               .f = function(x){
                                 class(out[[x]])
                               })
    
    cols_to_reclass_df <- tibble(
      name = names(out),
      class_created = out_col_classes
    ) %>% 
      left_join(evidence_classes_df,
                by = join_by(name)) %>% 
      mutate(class = if_else(data_type == "M",
                             "numeric",
                             class)) %>% 
      filter(class_created != class)
    
    # suppose column classes
    walk(cols_to_reclass_df$name,
         .f = function(colname){
           class(out[[colname]]) <<- 
             evidence_classes_df$class[evidence_classes_df$name == colname]
         }, .progress = TRUE
    )
    
    out
    
  }

#' @export
divide_data_to_batches <- function(data,rows_per_batch){
  # number of rows per batch
  batch_rows <- rows_per_batch
  # list of dataframes
  batches <- split(data, 
                   (seq(nrow(data))-1) %/% batch_rows + 1) 
  batches
}

#' @export
unnest_evidences_batch <-
  function(data, chunk_size = 10000, result_path){
    
    # data <- test_adaptive # TEST
    
    n <- nrow(data)
    d_size <- nrow(data)
    
    if (chunk_size > d_size){
      stop("Chunk size is more than number of rows of data")
    }
    

    r  <- rep(1:ceiling(n/chunk_size),each=chunk_size)[1:n]
    d <- split(data,r)
    
    d
    
    out <- map(d,
        unnest_evidences,
        .progress = TRUE)
    
    out <- bind_rows(out)
    
    write_feather(out, 
                  sink = result_path)
    message(paste0("saved unnested evidences to '", result_path,"'"))
    
    out
    
  }


# 2. Unnest Differential Diagnosis ----------------------------------------

#> Description:
#> In this section, we unnest the differential diagnosis.
#> This means separating all the diagnosis per patient
#> and giving each one its own row. This will make it 
#> easier later on to manipulate the data.
#> UPDATE: We now separate the differential diagnosis
#> to different columns. This is because we need this
#> for compositional data analysis.

# data <- 
#   read_feather("data/processed/prep1_unnestevidences.arrow")
# data <- 
#   data %>%
#   slice_head(n = 100000)

#' @export
unnest_differential <- function(data){
  separated_diag_df <- data %>%
    select(patientId, DIFFERENTIAL_DIAGNOSIS) %>% 
    # extract the diagnosis name from DIFFERENTIAL_DIAGNOSIS
    mutate(diagnosis_list = str_extract_all(DIFFERENTIAL_DIAGNOSIS, "(?<=\\[\\')([:graph:]|\\s)+?(?=\\'\\,)")) %>%
    # extract the diagnosis score from DIFFERENTIAL_DIAGNOSIS
    mutate(score_list = str_extract_all(DIFFERENTIAL_DIAGNOSIS, "(?<=\\,\\s)[\\d\\.e\\-]+?(?=\\])")) %>%
    select(-DIFFERENTIAL_DIAGNOSIS) %>% 
    rowwise() %>%
    mutate(diagnosis = list(unlist(diagnosis_list))) %>%
    mutate(differential_score = list(unlist(score_list))) %>%
    unnest_longer(col = c(diagnosis, differential_score)) %>%
    select(-diagnosis_list, -score_list) %>%
    mutate(
      differential_score = as.numeric(differential_score)
    ) %>%
    mutate(diagnosis = paste0("diag_",diagnosis)) %>% 
    pivot_wider(names_from = diagnosis,
                values_from = differential_score)
  
  new_df <- data %>%
    select(-DIFFERENTIAL_DIAGNOSIS) %>% 
    right_join(separated_diag_df, by = "patientId") %>%
    arrange(patientId)
  
  return(new_df)
}

# Sys.time()
# check <- unnest_differential(data)
# names(check)
# Sys.time()
# View(head(check))


#' @export
dd_pivot_wider <-
  function(path_unnested_differential,
           save_path){
    
    # path_unnested_differential <-
    #   "data/processed/unnested_differential/"
    # save_path <-
    #   "data/processed/unnested_differential_wide/"
    path_unnested_differential <-
      "data/processed/validation/unnested_differential/"
    save_path <-
      "data/processed/validation/unnested_differential_wide/"
    
    
    dd_feather_files <- 
      list.files(path_unnested_differential)
    # dd_feather_files <- 
    #   dd_feather_files[1:2]
    
    iwalk(
      dd_feather_files,
      .f = function(x, idx){
        
        # x <- dd_feather_files[1]
        
        dat <- 
          read_feather(paste0(path_unnested_differential,
                              x))
        
        
        
        # Remove unecessary fields and 
        # create initial compositional dependent variable
        # by create a column for each unique diagnosis (49 pathologies = 49 columns)
        train_prep1 <-
          dat %>% 
          select(-INITIAL_EVIDENCE) %>% 
          # mutate(diagnosis = paste0("diag_",diagnosis)) %>% 
          pivot_wider(names_from = diagnosis, values_from = differential_score) %>% 
          as_tibble()
        
        rm(dat)
        
        # save results
        write_feather(
          train_prep1,
          sink = paste0(save_path,
                        "unnested_dd_wide_batch",
                        idx,".arrow")
        )
        
        rm(train_prep1)
        gc()
        
      }, .progress = TRUE
    )
  
    
  }

# 3. Data Preparation for CoDa Modeling ---------------------------------


#' @export
coda_prediction_preparation <-
  function(data){
    
    box::use(
      stats[...],
      caret[...],
      compositions[...]
    )
    

    # data <- check # TEST
    # data <- train_prep1 # TEST
    # data <- prot_dd
    
    
    # Extract outcome variables
    outcome <- 
      data %>% 
      select(matches("diag"))
    
    # Convert all NA diagnosis to very small number
    # Since comp regression cannot accept NA
    outcome[is.na(outcome)] <- 0.0000001
    # Close the compositional data
    outcome <- clo(outcome)
    
    id_cols <-
      c("patientId","PATHOLOGY")
    
    id <- 
      data %>% 
      select(all_of(id_cols))
    
    predictors <-
      data %>% 
      select(-matches("diag")) %>% 
      select(-all_of(id_cols)) %>% 
      mutate(across(where(is.character), as.factor))
    
    
    # Prepare numeric predictors
    pred_num <- 
      predictors %>% 
      select(where(is.numeric))
    # Normalize the numeric variables
    pred_num <- 
      pred_num
      # mutate(across(everything(),~scale(.)[,1]))
    # Convert all NA values to 0
    pred_num[is.na(pred_num)] <- 0
    
    
    # Prepare factor predictors
    # pred_fct <- 
    #   predictors %>% 
    #   select(-where(is.numeric))
    # # dummify the data
    # dmy <- dummyVars(" ~ .", data = pred_fct)
    # pred_fct <- data.frame(predict(dmy, newdata = pred_fct))
    
    pred_fct <- 
      utils$get_dummy_categories_from_pred(predictors)
    
    # Scale the variables
    # pred_fct <- 
    #   pred_fct %>% 
    #   mutate(across(everything(),~scale(.)[,1]))
    
    # Combine all dfs
    # predictor_prepped <-
    #   bind_cols(id, pred_num, pred_fct, outcome)
    
    list(predictors = 
           bind_cols(pred_num, pred_fct),
         id = id,
         outcome = as_tibble(outcome))
    
  }





