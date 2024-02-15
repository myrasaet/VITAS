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
from constants import n_questions, base_path
import warnings
warnings.filterwarnings("ignore")

print("App is initializing...")
################### Load Data ###################
# load feature importance from global explainer output
with open(f"{base_path}\\output\\feature_importance.json") as f:
  feature_importance_dict = json.load(f)

# get evidence list
with open(f"{base_path}\\input\\release_evidences.json") as f:
  evidences = json.load(f)
evidences_list = []
evidences_code_to_en = {}
evidences_en_to_code = {}
for e in evidences.keys():
  # only binary symptoms and no antecedents
  if (not evidences[e]["possible-values"]) and (not evidences[e]["is_antecedent"]):
    evidences_list.append(evidences[e]["question_en"])
    evidences_code_to_en[e] = evidences[e]["question_en"]
    evidences_en_to_code[evidences[e]["question_en"]] = e

# get disease list
with open(f"{base_path}\\input\\release_conditions.json") as f:
  disease_dict = json.load(f)
disease_list = list(disease_dict.keys())

model_dict = {}
for disease in disease_list:
    disease_filename = re.sub('[^a-zA-Z0-9 \n\.]', '', disease).replace(" ", "_")
    with open(f'{base_path}\\output\\diseases\\{disease_filename}\\{disease_filename}_model.pkl', 'rb') as f:
        model_dict[disease] = pickle.load(f)

# vectorize feature importance
feature_importance_df = pd.DataFrame()
feature_importance_df["evidence"] = evidences_list
for disease in feature_importance_dict:
    feature_importance_df[disease] = [feature_importance_dict[disease]["top10_relevant_symptoms"].get(evidence, 0) for evidence in evidences_list]
feature_importance_df.set_index('evidence', inplace=True)
feature_importance_df
################### Load Data ###################


################### Load Models ###################
with open(f'{base_path}\\output\\questionnaire\\questionairre.pkl', 'rb') as f:
    questionairre = pickle.load(f)
with open(f'{base_path}\\output\\semantic_search\\semantic_search.pkl', 'rb') as f:
    semantic_search = pickle.load(f)
transformer_name = "BioBERT-mnli-snli-scinli-scitail-mednli-stsb"
transformer = SentenceTransformer(f"{base_path}\\input\\{transformer_name}")
################### Load Models ###################


################### Define Functions ###################
def get_next_question(evidences):
    centroid = np.array([feature_importance_df.loc[e].values for e in evidences]).mean(axis=0)
    _, indices = questionairre.kneighbors([centroid])
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

# TO DO - add params for which model to use RF/Logistic
def pred_explain(x, asked):
    # create output path per patient
    datetime_now = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
    img_path = f"{base_path}\\output\\consultation\\{datetime_now}"
    if not os.path.exists(img_path):
        os.makedirs(img_path)
    
    # predict
    pred_list = []
    for target_disease in disease_list:
        rf_model = model_dict[target_disease]
        prediction = rf_model.predict_proba(x)
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
        # pred_dict = pred_dict_all
        diagnosis_prediction = {
            "diagnosis_prediction": pred_dict_all
        }
        with open(f"{img_path}\\diagnosis_prediction.json", "w") as outfile: 
            json.dump(diagnosis_prediction, outfile, indent=True)
        for target_disease in pred_dict:
            rf_model = model_dict[target_disease]
            prediction, bias, contributions = ti.predict(rf_model, x)
            contributions_values = contributions[0][:,1]
            symptoms_en = x.columns.map(evidences_code_to_en)
            symptoms_values = [str(x[f].values[0]) for f in x.columns]
            symptoms_df = pd.DataFrame({"symptoms_en": symptoms_en, "symptoms_values": symptoms_values})
            contributions_df = pd.DataFrame({"symptoms_en": symptoms_en, "contributions_values": contributions_values, "contributions_abs_values": abs(contributions_values)})
            contributions_df.index  = symptoms_df["symptoms_en"] + "=" + symptoms_df["symptoms_values"]
            contributions_df = contributions_df[contributions_df['symptoms_en'].isin(asked)] # why empty????
            contributions_df = contributions_df.sort_values(by="contributions_abs_values", ascending=False).head(10).sort_values(by="contributions_abs_values")
            contributions_df["contributions_values"].plot.barh()
            plt.xlabel("Symptom Importance Score")
            plt.title(f"Probability of {target_disease}: {pred_dict[target_disease]:.3f}")
            plt.figtext(.01, .99, 'Symptoms with bars pointing to the right support a positive diagnosis.\nSymptoms with bars pointing to the left do not support a positive diagnosis.')
            img_filename = re.sub('[^a-zA-Z0-9 \n\.]', '', target_disease).replace(" ", "_")
            plt.savefig(f"{img_path}\\{img_filename}.jpg", bbox_inches='tight')
            plt.clf()
    return pred_dict_all, img_path
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
while (question_counter < n_questions) and ask:
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

# ask relevant evidences using KNN
while question_counter < n_questions:
    ask =  True
    next_question_idx = 0
    next_question = get_next_question(evidences)
    while ask and (question_counter < n_questions):
      if next_question[next_question_idx] not in asked:
        answer = input(f"Q{question_counter}: {next_question[next_question_idx]} (y/n)  ")
        asked.append(next_question[next_question_idx])
        question_counter+=1
        if answer=="y":
          evidences.append(next_question[next_question_idx])
          ask = False
        else:
          next_question_idx += 1
        if next_question_idx > n_questions:
          break
      else:
          next_question_idx += 1

print("Analyzing...")
input_vector = vectorize_input(evidences, age, sex)
_, output = pred_explain(input_vector, asked)
print(f"Done! Please see output in {output}")
################### Serve ###################

