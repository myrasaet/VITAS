# VITAS 

This is VITAS - Virtual Intelligent Tool for Assessment and Symptom Check. This tool predicts differential diagnosis given patient symptoms.

## About the Dataset
This repo uses the DDXPlus Dataset from https://figshare.com/articles/dataset/DDXPlus_Dataset_English_/22687585

## Set up
1. Install Python >= 3.10
2. Install venv
    ```bash
    pip install virtualenv
    ```
3. Go to your repo directory
    ```bash
    cd <repo_directory>
    ```
4. Create virtual env
    ```bash
    py -3.10 -m venv "vitas-env"
    ```
5. Activate vitas-env
    ```bash
    vitas-env\Scripts\activate
    ```
6. Install requirements in virtual env
    ```bash
    pip install -r requirements.txt
    ```
7. Download input data from https://figshare.com/articles/dataset/DDXPlus_Dataset_English_/22687585 and put in ```data/input``` folder 
8. Configure settings in ```src/constants.py```.
9. Download sentence transformers. Run ```src/download_tranformers.ipynb```. Use vitas-env kernel.
10. Train pathology models, run ```src/train.ipynb```. This repo trains a decision tree, random forest, gradient boosting, and logistic regression model for each of the 49 diseases in the DDXPlus dataset. Use vitas-env kernel.
11. (Optional) Evaluate pathology models, run ```src/eval.ipynb```. Use vitas-env kernel.
12. (Optional) Investigate pathology models errors, run ```src/error_analysis.ipynb```. Use vitas-env kernel.
13. (Optional) Test explainer, run ```src/explainer.ipynb```. Use vitas-env kernel.
14. Train semantic search, run ```src/semantic_search.ipynb```. Use vitas-env kernel.
15. Train adaptive questionnaire, run ```src/train_questionnaire.ipynb```. Use vitas-env kernel.
16. (Optional) Evaluate adaptive questionnaire together with pathology models, run ```src/eval_questionnaire_experiments.ipynb```. Use vitas-env kernel. Save predictions from experiments by running ```src/save_predictions.ipynb```
17. (Optional) Investigate adaptive questionnaire together with pathology models errors, run ```src/error_analysis_questionnaire.ipynb```. Use vitas-env kernel.
18. Review outputs in ```data/output``` folder

## Usage
1. Go to your repo directory
    ```bash
    cd "repo_directory"
    ```
2. Activate vitas-env
    ```bash
    vitas-env\Scripts\activate
    ```
3. To consult with VITAS, run
    ```bash
    cd src
    python app.py
    ```
4. Interact with VITAS. Max text input is 100 words.
5. Outputs are in ```data/output/consultation``` folder
