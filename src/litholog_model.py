import pandas as pd
import numpy as np
from xgboost import XGBClassifier
from sklearn.metrics import accuracy_score



def load_litholog(path="data/litholog1.csv"):
    df = pd.read_csv(path)
    return df


def expand_layers(df):
    data = []

    for _, row in df.iterrows():
        top = int(row["Top"])
        bottom = int(row["Bottom"])
        facies = row["Facies"]

        for d in range(top, bottom):
            data.append([d, facies])

    expanded = pd.DataFrame(data, columns=["Depth", "Facies"])
    return expanded


def encode_facies(df):
    mapping = {
        "sand": 0,
        "mud": 1,
        "coal": 2,
        "carbon_mud": 3
    }

    df["Facies"] = df["Facies"].map(mapping)
    return df, mapping

def sample_data(df, step):
    return df.iloc[::step].reset_index(drop=True)

from xgboost import XGBClassifier
from sklearn.metrics import accuracy_score


def train_and_evaluate(full_df, step):

    sampled = sample_data(full_df, step)

    X_train = sampled[["Depth"]]
    y_train = sampled["Facies"]

    X_test = full_df[["Depth"]]
    y_test = full_df["Facies"]

    # Fix label issue
    unique_classes = np.unique(y_train)
    class_mapping = {c: i for i, c in enumerate(unique_classes)}

    y_train_mapped = y_train.map(class_mapping)

    model = XGBClassifier(n_estimators=50, max_depth=3)
    model.fit(X_train, y_train_mapped)

    # Map predictions back
    inv_mapping = {v: k for k, v in class_mapping.items()}
    preds_mapped = model.predict(X_test)
    preds = [inv_mapping[p] for p in preds_mapped]

    preds = model.predict(X_test)

    acc = accuracy_score(y_test, preds)

    print(f"Step {step} | Accuracy: {acc:.4f}")

    return acc