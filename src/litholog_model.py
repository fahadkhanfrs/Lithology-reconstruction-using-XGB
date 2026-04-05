import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
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

def add_context(df):
    df = df.copy()
    df["prev_facies"] = df["Facies"].shift(1)
    df = df.dropna()
    return df


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

def train_and_evaluate(full_df, step):

    sampled = sample_data(full_df, step)

    X_train = sampled[["Depth", "prev_facies"]]
    y_train = sampled["Facies"]

    X_test = full_df[["Depth", "prev_facies"]]
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
def get_uncertainty(model, X):
    probs = model.predict_proba(X)
    uncertainty = 1 - np.max(probs, axis=1)
    return uncertainty

def uniform_sampling_experiment(full_df, budget=30):

    # Pick evenly spaced indices
    indices = np.linspace(0, len(full_df)-1, budget).astype(int)
    sampled = full_df.iloc[indices]

    X_train = sampled[["Depth", "prev_facies"]]
    y_train = sampled["Facies"]

    unique_classes = np.unique(y_train)
    class_mapping = {c: i for i, c in enumerate(unique_classes)}
    y_train_mapped = y_train.map(class_mapping)

    model = XGBClassifier(n_estimators=50, max_depth=3)
    model.fit(X_train, y_train_mapped)

    X_test = full_df[["Depth", "prev_facies"]]
    y_test = full_df["Facies"]

    preds_mapped = model.predict(X_test)
    inv_mapping = {v: k for k, v in class_mapping.items()}
    preds = np.array([inv_mapping[p] for p in preds_mapped])

    acc = accuracy_score(y_test, preds)

    print(f"Uniform Sampling | Accuracy: {acc:.4f} | Samples used: {len(sampled)}")

    return acc
def uncertainty_sampling_experiment(full_df, initial_step=10, budget=20):

    # Start with sparse sampling
    sampled = full_df.iloc[::initial_step].copy()

    remaining = full_df.drop(sampled.index)

    for i in range(budget):

        X_train = sampled[["Depth", "prev_facies"]]
        y_train = sampled["Facies"]

        # Fix class issue
        unique_classes = np.unique(y_train)
        class_mapping = {c: i for i, c in enumerate(unique_classes)}
        y_train_mapped = y_train.map(class_mapping)

        model = XGBClassifier(n_estimators=50, max_depth=3)
        model.fit(X_train, y_train_mapped)

        # Predict on remaining points
        X_remain = remaining[["Depth", "prev_facies"]]

        # Map predictions
        inv_mapping = {v: k for k, v in class_mapping.items()}
        preds_mapped = model.predict(X_remain)
        preds = np.array([inv_mapping[p] for p in preds_mapped])

        # Get uncertainty
        probs = model.predict_proba(X_remain)
        uncertainty = 1 - probs.max(axis=1)

        # Pick most uncertain point
        idx = np.argmax(uncertainty)

        # Add it to sampled
        new_sample = remaining.iloc[[idx]]
        sampled = pd.concat([sampled, new_sample])

        # Remove from remaining
        remaining = remaining.drop(new_sample.index)

    # Final evaluation
    X_train = sampled[["Depth", "prev_facies"]]
    y_train = sampled["Facies"]

    unique_classes = np.unique(y_train)
    class_mapping = {c: i for i, c in enumerate(unique_classes)}
    y_train_mapped = y_train.map(class_mapping)

    model = XGBClassifier(n_estimators=50, max_depth=3)
    model.fit(X_train, y_train_mapped)

    X_test = full_df[["Depth", "prev_facies"]]
    y_test = full_df["Facies"]

    preds_mapped = model.predict(X_test)
    inv_mapping = {v: k for k, v in class_mapping.items()}
    preds = np.array([inv_mapping[p] for p in preds_mapped])

    acc = accuracy_score(y_test, preds)

    print(f"Uncertainty Sampling | Accuracy: {acc:.4f} | Samples used: {len(sampled)}")

    return acc

def hybrid_sampling_experiment(full_df, initial_step=10, budget=20):

    # Step 1: uniform base (coverage)
    sampled = full_df.iloc[::initial_step].copy()
    remaining = full_df.drop(sampled.index)

    for i in range(budget):

        X_train = sampled[["Depth", "prev_facies"]]
        y_train = sampled["Facies"]

        unique_classes = np.unique(y_train)
        class_mapping = {c: i for i, c in enumerate(unique_classes)}
        y_train_mapped = y_train.map(class_mapping)

        model = XGBClassifier(n_estimators=50, max_depth=3)
        model.fit(X_train, y_train_mapped)

        X_remain = remaining[["Depth", "prev_facies"]]

        probs = model.predict_proba(X_remain)
        uncertainty = 1 - probs.max(axis=1)

        # pick top uncertain
        idx = np.argmax(uncertainty)

        new_sample = remaining.iloc[[idx]]
        sampled = pd.concat([sampled, new_sample])
        remaining = remaining.drop(new_sample.index)

    # Final evaluation
    X_train = sampled[["Depth", "prev_facies"]]
    y_train = sampled["Facies"]

    unique_classes = np.unique(y_train)
    class_mapping = {c: i for i, c in enumerate(unique_classes)}
    y_train_mapped = y_train.map(class_mapping)

    model = XGBClassifier(n_estimators=50, max_depth=3)
    model.fit(X_train, y_train_mapped)

    X_test = full_df[["Depth", "prev_facies"]]
    y_test = full_df["Facies"]

    preds_mapped = model.predict(X_test)
    inv_mapping = {v: k for k, v in class_mapping.items()}
    preds = np.array([inv_mapping[p] for p in preds_mapped])

    acc = accuracy_score(y_test, preds)

    print(f"Hybrid Sampling | Accuracy: {acc:.4f} | Samples used: {len(sampled)}")
    plot_reconstruction(full_df, sampled, model, class_mapping)
    return acc

def plot_reconstruction(full_df, sampled, model, class_mapping):

    inv_mapping = {v: k for k, v in class_mapping.items()}

    X_test = full_df[["Depth", "prev_facies"]]
    preds_mapped = model.predict(X_test)
    preds = np.array([inv_mapping[p] for p in preds_mapped])

    depth = full_df["Depth"].values
    true = full_df["Facies"].values

    # 🎨 Colors
    colors = {
        0: "orange",   # sand
        1: "green",    # mud
        2: "black",    # coal
        3: "purple"    # carbon_mud
    }

    fig, axes = plt.subplots(1, 3, figsize=(9, 10), sharey=True)
    # ---- BLOCK PLOT FUNCTION ----
    def draw_blocks(ax, facies_array, title):
        start = 0

        for i in range(1, len(facies_array)):
            if facies_array[i] != facies_array[i - 1]:
                ax.fill_betweenx(
                    depth[start:i],
                    0,
                    1,
                    color=colors[facies_array[i - 1]]
                )
                start = i

        # last block
        ax.fill_betweenx(
            depth[start:],
            0,
            1,
            color=colors[facies_array[-1]]
        )

        ax.set_title(title)
        ax.set_xticks([])
        ax.set_xlim(0, 1)
        ax.invert_yaxis()

    # ---- TRUE ----
    draw_blocks(axes[0], true, "True Lithology (Ground Truth)")

    # ---- PREDICTED ----
    draw_blocks(axes[1], preds, "Predicted Lithology (Hybrid Sampling)")

    # 🔴 Sampled points
    axes[1].scatter(
        np.random.uniform(0.4, 0.6, size=len(sampled)),
        sampled["Depth"],
        color="red",
        s=20,
        label="Sampled"
    )

    # ❌ Error highlighting
    error_mask = preds != true
    axes[1].scatter(
        np.random.uniform(0.7, 0.9, size=error_mask.sum()),
        depth[error_mask],
        color="yellow",
        s=10,
        label="Error"
    )

    # ---- Legend ----
    labels = {
        0: "Sand",
        1: "Mud",
        2: "Coal",
        3: "Carbon Mud"
    }

    legend_elements = [
        plt.Line2D([0], [0], marker='s', color='w',
                   label=labels[k], markerfacecolor=colors[k], markersize=10)
        for k in labels
    ]

    legend_elements += [
        plt.Line2D([0], [0], marker='o', color='w',
                   label='Sampled Points', markerfacecolor='red', markersize=8),
        plt.Line2D([0], [0], marker='o', color='w',
                   label='Errors', markerfacecolor='yellow', markersize=8)
    ]

    axes[1].legend(handles=legend_elements, loc='center left', bbox_to_anchor=(1, 0.5))
    
    error_band = (preds != true).astype(int)

    axes[2].imshow(
        error_band.reshape(-1, 1),
        aspect='auto',
        cmap='Reds', vmin=0, vmax=1
    )

    axes[2].set_title("Error Distribution")
    axes[2].set_xticks([])
    axes[2].invert_yaxis()

    plt.tight_layout()
    plt.show()