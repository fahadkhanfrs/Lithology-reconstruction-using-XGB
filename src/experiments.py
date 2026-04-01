from xgboost import XGBRegressor
from sklearn.metrics import mean_absolute_error
import matplotlib.pyplot as plt


def plot_feature_importance(model, feature_names, title):
    importance = model.feature_importances_

    plt.figure()
    plt.bar(feature_names, importance)
    plt.title(f"Feature Importance: {title}")
    plt.xlabel("Features")
    plt.ylabel("Importance")
    plt.show()


def run_experiments(df):

    target = "RHOB"

    feature_sets = {
        "GR only": ["GR"],
        "GR + SP": ["GR", "SP"],
        "GR + SP + RILD": ["GR", "SP", "RILD"],
        "GR + SP + RILD + CNLS": ["GR", "SP", "RILD", "CNLS"],
    }

    costs = {
        "GR": 1,
        "SP": 2,
        "RILD": 4,
        "CNLS": 3
    }

    results = []

    for name, features in feature_sets.items():

        cols = ["Depth"] + features + [target]
        data = df[cols].dropna().sort_values(by="Depth")

        # Depth-based split
        split = int(0.8 * len(data))
        train_df = data.iloc[:split]
        test_df = data.iloc[split:]

        X_train = train_df[features]
        y_train = train_df[target]

        X_test = test_df[features]
        y_test = test_df[target]

        # Train model
        model = XGBRegressor(
            n_estimators=100,
            max_depth=5,
            learning_rate=0.1
        )

        model.fit(X_train, y_train)

        # Evaluate
        preds = model.predict(X_test)
        mae = mean_absolute_error(y_test, preds)

        total_cost = sum(costs[f] for f in features)

        print(f"{name} | MAE: {mae:.4f} | Cost: {total_cost}")

        results.append((name, mae, total_cost))

        # 🔥 Plot feature importance ONLY for important cases
        if name in ["GR + SP", "GR + SP + RILD"]:
            plot_feature_importance(model, features, name)

    return results