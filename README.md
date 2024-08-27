# VITAS

This is VITAS - Virtual Intelligent Tool for Assessment and Symptom Check. This tool predicts differential diagnosis given patient symptoms.

## About the Dataset

This repo uses the DDXPlus Dataset from <https://figshare.com/articles/dataset/DDXPlus_Dataset_English_/22687585>

## Project Structure

The project is organized as follows:

```         
- main.R                  # Main script to run the entire pipeline (part 2)
- data 
    - input               # Unprocessed, original data files
    - processed           # Cleaned and processed data ready for analysis
    - output              # Final outputs such as plots and tables
- src                     # Folder containing python scripts (for part 1)
- R                       # R scripts and modules (for part 2)
    - pull.R              # Module to extract data from sources
    - preparation.R       # Module for data processing
    - coda.R              # Module for feature engineering and modeling
    - evaluation.R        # Module for model evaluation
    - utils.R             # Module for utility functions
- renv                    # R environment managed by renv
- renv.lock               # Lockfile to ensure reproducibility of the R environment
- requirements.txt        # Python package dependencies
- README.md               # Project documentation
- modeling_objects        # `pins` storage for models and metrics
```

### Data

-   **`data/input`:** Contains the raw, unprocessed data files. These files should not be modified directly.

-   **`data/processed`:** Contains the cleaned and processed data, ready for analysis.

### Environment and Dependencies

-   **`renv` and `renv.lock`:** These files manage the R environment for this project. Ensure that `renv.lock` is included to allow others to reproduce your environment.

-   **`requirements.txt`:** Lists the Python packages required for this project. These can be installed using `pip install -r requirements.txt`.

## Part 1: Semantic Search \| Adaptive Questionnaire \| Single Diagnosis Models and Explainers

### Set up

1.  Install Python \>= 3.10
2.  Install venv `bash     pip install virtualenv`
3.  Go to your repo directory `bash     cd <repo_directory>`
4.  Create virtual env `bash     py -3.10 -m venv "vitas-env"`
5.  Activate vitas-env `bash     vitas-env\Scripts\activate`
6.  Install requirements in virtual env `bash     pip install -r requirements.txt`
7.  Download input data from <https://figshare.com/articles/dataset/DDXPlus_Dataset_English_/22687585> and put in `data/input` folder
8.  Configure settings in `src/constants.py`.
9.  Download sentence transformers. Run `src/download_tranformers.ipynb`. Use vitas-env kernel.
10. Train pathology models, run `src/train.ipynb`. This repo trains a decision tree, random forest, gradient boosting, and logistic regression model for each of the 49 diseases in the DDXPlus dataset. Use vitas-env kernel.
11. (Optional) Evaluate pathology models, run `src/eval.ipynb`. Use vitas-env kernel.
12. (Optional) Investigate pathology models errors, run `src/error_analysis.ipynb`. Use vitas-env kernel.
13. (Optional) Test explainer, run `src/explainer.ipynb`. Use vitas-env kernel.
14. Train semantic search, run `src/semantic_search.ipynb`. Use vitas-env kernel.
15. Train adaptive questionnaire, run `src/train_questionnaire.ipynb`. Use vitas-env kernel.
16. (Optional) Evaluate adaptive questionnaire together with pathology models, run `src/eval_questionnaire_experiments.ipynb`. Use vitas-env kernel. Save predictions from experiments by running `src/save_predictions.ipynb`
17. (Optional) Investigate adaptive questionnaire together with pathology models errors, run `src/error_analysis_questionnaire.ipynb`. Use vitas-env kernel.
18. Review outputs in `data/output` folder

### Usage

1.  Go to your repo directory `bash     cd "repo_directory"`
2.  Activate vitas-env `bash     vitas-env\Scripts\activate`
3.  To consult with VITAS, run `bash     cd src     python app.py`
4.  Interact with VITAS. Max text input is 100 words.
5.  Outputs are in `data/output/consultation` folder

## Part 2: Feature Engineering \| CoDa \| Differential Diagnosis Model and Explainers

### Setup Environment

**R Environment:** Install R and use renv::restore() to recreate the environment specified in renv.lock.

### Run the Pipeline

Execute `main.R` to run the entire data processing and analysis pipeline.

The documentation for each step is found in `main.R`

To consult documentation for each function used in `main.R`, run `box::help(insert_function)` .
