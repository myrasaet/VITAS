import pandas as pd
import re
import json
import os
import numpy as np
from sentence_transformers import SentenceTransformer
from treeinterpreter import treeinterpreter as ti
import pickle
import matplotlib.pyplot as plt
from datetime import datetime
from constants import app_n_questions, base_path, model_list
import warnings
warnings.filterwarnings("ignore")

print("App is initializing...")
################### Load Data ###################
# get evidence list
with open(f"{base_path}\\input\\release_evidences.json") as f:
  evidences = json.load(f)
evidences_list = []
evidences_code_to_en = {"AGE": "AGE", "SEX": "SEX"}
evidences_en_to_code = {"AGE": "AGE", "SEX": "SEX"}
for e in evidences.keys():
  # only binary symptoms and antecedents
  if (not evidences[e]["possible-values"]):
    evidences_list.append(evidences[e]["question_en"])
    evidences_code_to_en[e] = evidences[e]["question_en"]
    evidences_en_to_code[evidences[e]["question_en"]] = e

# get disease list
with open(f"{base_path}\\input\\release_conditions.json") as f:
  disease_dict = json.load(f)
disease_list = list(disease_dict.keys())

model_list = list(model_list["tree-based"].keys())
model_list.append("logistic_regression")
model_dict = {model_name:{} for model_name in model_list}
feature_importance = {}
for model_name in model_list:
  for disease in disease_list:
      disease_filename = re.sub('[^a-zA-Z0-9 \n\.]', '', disease).replace(" ", "_")
      with open(f'{base_path}\\output\\diseases\\{disease_filename}\\{model_name}\\{disease_filename}_model.pkl', 'rb') as f:
          model_dict[model_name][disease] = pickle.load(f)
      with open(f'{base_path}\\output\\diseases\\{disease_filename}\\logistic_regression\\feature_importance.json', 'rb') as f:
        feature_importance[disease] = json.load(f)
################### Load Data ###################


################### Load Models ###################
feature_embeddings_df = pd.read_pickle(f'{base_path}\\output\\questionnaire\\questionnaire_embeddings.pkl')
with open(f'{base_path}\\output\\questionnaire\\questionnaire.pkl', 'rb') as f:
    questionnaire = pickle.load(f)
with open(f'{base_path}\\output\\semantic_search\\semantic_search.pkl', 'rb') as f:
    semantic_search = pickle.load(f)
transformer_name = "BioBERT-mnli-snli-scinli-scitail-mednli-stsb"
transformer = SentenceTransformer(f"{base_path}\\input\\{transformer_name}")
################### Load Models ###################


################### Define Functions ###################
def get_next_question(evidences):
    centroid = np.array([feature_embeddings_df.loc[e].values for e in evidences]).mean(axis=0)
    _, indices = questionnaire.kneighbors([centroid])
    ask_list = [evidences_list[i] for i in indices[0] if evidences_list[i] not in evidences]
    try:
        return ask_list
    except:
        return []

def vectorize_input(evidences, age, sex):
    df = pd.DataFrame()
    df["AGE"] = [int(age)]
    # to do: add handling of unexpected input values. assumes strict input values of f/m
    df["SEX"] = [1 if sex=="m" else 0]
    for e in evidences_list:
        df[evidences_en_to_code[e]] = [1 if e in evidences else 0]
    return df

def pred_explain(x):
    # create output path per patient
    datetime_now = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
    img_path = f"{base_path}\\output\\consultations\\{datetime_now}"    
    for model_name in model_list:
      model_output_path = f"{img_path}\\{model_name}"
      if not os.path.exists(model_output_path):
        os.makedirs(model_output_path)

      # predict
      pred_list = []
      for target_disease in disease_list:
          clf_model = model_dict[model_name][target_disease]
          prediction = clf_model.predict_proba(x)
          prediction_proba = prediction[0][1]
          if prediction_proba>0:
              pred_list.append({
                  "disease": target_disease,
                  "probability": prediction_proba
              })
      
      # return all predictions, no rank filter
      # but explain top 1 only - rank allows ties
      if pred_list:
          pred_df = pd.DataFrame(pred_list).set_index('disease')
          pred_df['rank'] = pred_df['probability'].rank(method='min', ascending=False)
          pred_df = pred_df.sort_values(by="rank")
          pred_df = pred_df.sort_values(by="probability", ascending=False)
          pred_dict_all = pred_df.to_dict()["probability"]
          pred_df = pred_df[pred_df["rank"]<=1][["probability"]]
          pred_dict = pred_df.to_dict()["probability"]
          diagnosis_prediction = {
              "diagnosis_prediction": pred_dict_all
          }
          with open(f"{model_output_path}\\diagnosis_prediction.json", "w") as outfile: 
              json.dump(diagnosis_prediction, outfile, indent=True)

          # to do: configurize explainer flag rather than a hard-code condition
          if model_name in ["decision_tree", "random_forest", "logistic_regression"]:
            for target_disease in pred_dict:
                clf_model = model_dict[model_name][target_disease]
                symptoms_en = x.columns.map(evidences_code_to_en)
                symptoms_values = [x[f].values[0] for f in x.columns]
                symptoms_df = pd.DataFrame({"symptoms_en": symptoms_en, "symptoms_values": symptoms_values})
                if model_name=="logistic_regression":
                  prediction = clf_model.predict_proba(x)[0][1]
                  # model_coeffs = clf_model.coef_[0]
                  # get standardized coeffs
                  model_coeffs = [feature_importance[target_disease][evidences_code_to_en[f]] for f in x.columns]
                  contributions_values = [model_coeffs[i]*symptoms_values[i] for i in range(len(symptoms_values))]
                  contributions_df = pd.DataFrame({"symptoms_en": symptoms_en, "contributions_values": contributions_values, "contributions_abs_values": [abs(i) for i in contributions_values]})
                else:
                  prediction, _, contributions = ti.predict(clf_model, x)
                  contributions_values = contributions[0][:,1]
                  contributions_df = pd.DataFrame({"symptoms_en": symptoms_en, "contributions_values": contributions_values, "contributions_abs_values": abs(contributions_values)})
                contributions_df.index  = symptoms_df["symptoms_en"] + "=" + symptoms_df["symptoms_values"].astype(str)
                contributions_df["symptoms_values"] = symptoms_values
                contributions_df = contributions_df[contributions_df["contributions_values"]>0]
                contributions_df = contributions_df[(~contributions_df["symptoms_en"].isin(["AGE", "SEX"])) & (contributions_df["symptoms_values"]>0)]
                contributions_df = contributions_df.sort_values(by="contributions_abs_values", ascending=False).head(10).sort_values(by="contributions_abs_values")
                contributions_df["contributions_values"].plot.barh()
                plt.xlabel("Symptom Importance Score")
                plt.title(f"Probability of {target_disease}: {pred_dict[target_disease]:.3f}\n({model_name})")
                plt.figtext(.01, .99, 'Symptoms with higher importance score support a positive diagnosis.')                
                img_filename = re.sub('[^a-zA-Z0-9 \n\.]', '', target_disease).replace(" ", "_")
                plt.savefig(f"{model_output_path}\\{img_filename}.jpg", bbox_inches='tight')
                plt.clf()
    return img_path
################### Define Functions ###################

       
################### Serve ###################
print("Hi! I'm VITAS! I am a diagnosis support tool. I am here to give information about your health concerns.\nPlease keep in mind that while I can provide information and suggestions, I am not a substitute for professional medical advice. Always consult with a qualified healthcare provider for personalized diagnosis and treatment.")
input_text = input("What is your main health concern?  ")
age = input("What is your age?  ")
sex = input("What is your assigned sex? (f/m)  ")
evidences = []

question_counter = 0
asked = evidences.copy()

# get initial evidence from user input text using semantic search
input_embeddings = transformer.encode([input_text])
_, indices = semantic_search.kneighbors(input_embeddings)
ask_list = [evidences_list[i] for i in indices[0]]
ask = True
next_question_idx = 0
while (question_counter < app_n_questions) and ask:
  if ask_list[next_question_idx] not in asked:
    answer = input(f"Q{question_counter}: {ask_list[next_question_idx]} (y/n)  ")
    asked.append(ask_list[next_question_idx])
    question_counter+=1
    # to do: add handling of unexpected input values. assumes strict input values of y/n
    if answer=="y":
      evidences.append(ask_list[next_question_idx])
      ask = False
    else:
      next_question_idx += 1
  else:
    next_question_idx += 1

# ask relevant evidences using adaptive questionnaire
while question_counter < app_n_questions:
    ask =  True
    next_question_idx = 0
    next_question = get_next_question(evidences)
    while ask and (question_counter < app_n_questions):
      if next_question[next_question_idx] not in asked:
        answer = input(f"Q{question_counter}: {next_question[next_question_idx]} (y/n)  ")
        asked.append(next_question[next_question_idx])
        question_counter+=1
        if answer=="y":
          evidences.append(next_question[next_question_idx])
          ask = False
        else:
          next_question_idx += 1
        if next_question_idx > app_n_questions:
          break
      else:
          next_question_idx += 1

print("Analyzing...")
input_vector = vectorize_input(evidences, age, sex)
output = pred_explain(input_vector)
print(f"Done! Please see output in {output}")
################### Serve ###################

