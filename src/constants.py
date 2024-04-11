import os
from sklearn.linear_model import LogisticRegression
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier

model_list = {
    "logistic_regression": LogisticRegression(random_state=51, solver="lbfgs", max_iter=8000, multi_class="multinomial"),
    "tree-based":{
    "decision_tree": DecisionTreeClassifier(random_state=0),
    "random_forest": RandomForestClassifier(random_state=0),
    "gradient_boost": GradientBoostingClassifier(random_state=0)
    }
}

max_n_questions = 208
app_n_questions = 30
base_path = f"{os.path.dirname(os.getcwd())}\\data"