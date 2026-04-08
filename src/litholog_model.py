import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from xgboost import XGBClassifier
from sklearn.metrics import accuracy_score

def load_litholog(path="data/litholog9.csv"):
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
    return df.reset_index(drop=True)

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
    sampled = full_df.iloc[indices].reset_index(drop=True)

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

    return acc, full_df, sampled, model, class_mapping

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
        sampled = pd.concat([sampled, new_sample]).reset_index(drop=True)

        # Remove from remaining
        remaining = remaining.drop(new_sample.index).reset_index(drop=True)

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
        sampled = pd.concat([sampled, new_sample]).reset_index(drop=True)
        remaining = remaining.drop(new_sample.index).reset_index(drop=True)

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
    return acc, full_df, sampled, model, class_mapping

def plot_reconstruction(full_df, sampled, model, class_mapping):

    inv_mapping = {v: k for k, v in class_mapping.items()}

    X_test = full_df[["Depth", "prev_facies"]]
    preds_mapped = model.predict(X_test)
    preds = np.array([inv_mapping[p] for p in preds_mapped])

    # Train uniform sampling model
    budget = len(sampled)
    uniform_indices = np.linspace(0, len(full_df)-1, budget).astype(int)
    uniform_sampled = full_df.iloc[uniform_indices].reset_index(drop=True)
    
    X_train_uniform = uniform_sampled[["Depth", "prev_facies"]]
    y_train_uniform = uniform_sampled["Facies"]
    
    unique_classes_uniform = np.unique(y_train_uniform)
    class_mapping_uniform = {c: i for i, c in enumerate(unique_classes_uniform)}
    y_train_mapped_uniform = y_train_uniform.map(class_mapping_uniform)
    
    model_uniform = XGBClassifier(n_estimators=50, max_depth=3)
    model_uniform.fit(X_train_uniform, y_train_mapped_uniform)
    
    preds_uniform_mapped = model_uniform.predict(X_test)
    inv_mapping_uniform = {v: k for k, v in class_mapping_uniform.items()}
    preds_uniform = np.array([inv_mapping_uniform[p] for p in preds_uniform_mapped])

    depth = full_df["Depth"].to_numpy()
    true = full_df["Facies"].to_numpy()


    colors = {
        0: "orange",   # sand
        1: "green",    # mud
        2: "black",    # coal
        3: "purple",    # carbon_mud
    }

    fig, axes = plt.subplots(1, 3, figsize=(14, 10), sharey=True)
    def draw_blocks(ax, facies_array, depth, title, colors):

        start = 0

        for i in range(1, len(facies_array)):
            if facies_array[i] != facies_array[i - 1]:

                ax.fill_betweenx(
                    [depth[start], depth[i]],
                    0, 1,
                    color=colors[facies_array[i - 1]]
                )
                start = i

        ax.fill_betweenx(
            [depth[start], depth[-1]],
            0, 1,
            color=colors[facies_array[-1]]
        )

        ax.set_xlim(0, 1)
        ax.set_xticks([])
        ax.invert_yaxis()
        ax.set_title(title)

    draw_blocks(axes[0], true, depth, "True Lithology (Ground Truth)", colors)

    draw_blocks(axes[1], preds_uniform, depth, "Predicted Lithology (Uniform Sampling)", colors)

    draw_blocks(axes[2], preds, depth, "Predicted Lithology (Hybrid Sampling)", colors)

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

    axes[2].legend(handles=legend_elements, loc='center left', bbox_to_anchor=(1, 0.5))

    plt.tight_layout()
    plt.show()

def plot_error_analysis(full_df, sampled, model, class_mapping, method=""):
    """
    Enhanced error analysis plot showing:
    - Error distribution along depth
    - Per-class accuracy metrics
    - Error regions highlighted
    """
    
    inv_mapping = {v: k for k, v in class_mapping.items()}
    
    X_test = full_df[["Depth", "prev_facies"]]
    preds_mapped = model.predict(X_test)
    preds = np.array([inv_mapping[p] for p in preds_mapped])
    
    depth = full_df["Depth"].to_numpy()
    true = full_df["Facies"].to_numpy()
    
    # Calculate errors
    errors = (preds != true).astype(int)
    correct = 1 - errors
    
    # Per-class metrics
    facies_names = {0: "Sand", 1: "Mud", 2: "Coal", 3: "Carbon Mud"}
    per_class_acc = {}
    per_class_count = {}
    
    for facies_id in class_mapping.values():
        mask = true == facies_id
        if mask.sum() > 0:
            per_class_acc[facies_names[facies_id]] = correct[mask].mean()
            per_class_count[facies_names[facies_id]] = mask.sum()
    
    colors = {
        0: "orange",
        1: "green",
        2: "black",
        3: "purple",
    }
    
    # Create comprehensive error visualization
    fig = plt.figure(figsize=(16, 12))
    gs = fig.add_gridspec(3, 2, hspace=0.45, wspace=0.35, top=0.93, bottom=0.08, left=0.08, right=0.95)
    
    # 1. Error distribution along depth (line plot)
    ax1 = fig.add_subplot(gs[0, :])
    # Calculate rolling error rate (windowed)
    window = max(10, len(depth) // 50)
    rolling_error = []
    window_depths = []
    for i in range(0, len(depth), window):
        end = min(i + window, len(depth))
        error_rate = 1 - correct[i:end].mean()
        rolling_error.append(error_rate)
        window_depths.append(depth[i:end].mean())
    
    ax1.plot(window_depths, rolling_error, linewidth=2, color='red', alpha=0.7)
    ax1.fill_between(window_depths, rolling_error, alpha=0.3, color='red')
    ax1.set_xlabel("Depth", fontsize=10)
    ax1.set_ylabel("Error Rate (windowed)", fontsize=10)
    ax1.set_title("Prediction Error Rate Along Depth", fontsize=11, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    
    # 2. Error regions visualization (similar to reconstruction)
    ax2 = fig.add_subplot(gs[1, 0])
    start = 0
    for i in range(1, len(errors)):
        if errors[i] != errors[i-1]:
            color = 'red' if errors[i-1] == 1 else 'lightgreen'
            ax2.fill_betweenx([depth[start], depth[i]], 0, 1, color=color, alpha=0.6)
            start = i
    color = 'red' if errors[-1] == 1 else 'lightgreen'
    ax2.fill_betweenx([depth[start], depth[-1]], 0, 1, color=color, alpha=0.6)
    
    ax2.set_xlim(0, 1)
    ax2.set_xticks([])
    ax2.invert_yaxis()
    ax2.set_ylabel("Depth", fontsize=10)
    ax2.set_title("Error Regions (Red=Error, Green=Correct)", fontsize=11, fontweight='bold')
    
    # Add legend for error regions
    from matplotlib.patches import Patch
    error_legend = [
        Patch(facecolor='red', alpha=0.6, label='Incorrect Prediction'),
        Patch(facecolor='lightgreen', alpha=0.6, label='Correct Prediction')
    ]
    ax2.legend(handles=error_legend, loc='lower center', bbox_to_anchor=(0.5, -0.15), ncol=2, fontsize=9)
    
    # 3. Per-class accuracy bar chart
    ax3 = fig.add_subplot(gs[1, 1])
    class_names = list(per_class_acc.keys())
    accuracies = list(per_class_acc.values())
    bars = ax3.barh(class_names, accuracies, color=['orange', 'green', 'black', 'purple'][:len(class_names)], alpha=0.7)
    ax3.set_xlabel("Accuracy", fontsize=10)
    ax3.set_title("Per-Class Accuracy", fontsize=11, fontweight='bold')
    ax3.set_xlim(0, 1)
    
    # Add value labels on bars
    for i, (bar, acc) in enumerate(zip(bars, accuracies)):
        ax3.text(acc + 0.02, bar.get_y() + bar.get_height()/2, f'{acc:.2%}', 
                va='center', fontsize=9)
    
    # 4. Error count by facies
    ax4 = fig.add_subplot(gs[2, 0])
    error_counts = {}
    total_counts = {}
    for facies_id in class_mapping.values():
        mask = true == facies_id
        error_counts[facies_names[facies_id]] = errors[mask].sum()
        total_counts[facies_names[facies_id]] = mask.sum()
    
    class_names = list(error_counts.keys())
    error_vals = list(error_counts.values())
    bars = ax4.bar(class_names, error_vals, color=['orange', 'green', 'black', 'purple'][:len(class_names)], alpha=0.7)
    ax4.set_ylabel("Error Count", fontsize=10)
    ax4.set_title("Number of Misclassifications by Facies", fontsize=11, fontweight='bold')
    
    # Add value labels on bars
    for bar, count in zip(bars, error_vals):
        height = bar.get_height()
        ax4.text(bar.get_x() + bar.get_width()/2., height,
                f'{int(count)}', ha='center', va='bottom', fontsize=9)
    
    # 5. Overall statistics
    ax5 = fig.add_subplot(gs[2, 1])
    ax5.axis('off')
    
    overall_acc = correct.mean()
    overall_error_rate = errors.mean()
    
    stats_text = f"""
    OVERALL PERFORMANCE
    ─────────────────────────
    Total Accuracy:     {overall_acc:.2%}
    Error Rate:         {overall_error_rate:.2%}
    Total Samples:      {len(true):,}
    Correct:            {correct.sum():,}
    Incorrect:          {errors.sum():,}
    
    Sampled Points:     {len(sampled):,}
    Sampling %:         {len(sampled)/len(full_df):.1%}
    """
    
    ax5.text(0.1, 0.5, stats_text, fontsize=10, family='monospace',
            verticalalignment='center', bbox=dict(boxstyle='round', 
            facecolor='wheat', alpha=0.5))
    
    plt.suptitle(f"Lithology Prediction Error Analysis ({method})", fontsize=13, fontweight='bold', y=0.98)
    plt.show()