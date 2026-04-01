from xgboost import XGBRegressor
from sklearn.metrics import mean_absolute_error

def greedy_selection(df):

    target = "RHOB"
    all_features = ["GR", "SP", "RILD", "CNLS"]

    costs = {
        "GR": 1,
        "SP": 2,
        "RILD": 4,
        "CNLS": 3
    }

    selected = []
    remaining = all_features.copy()

    best_mae = float("inf")

    while remaining:
        best_feature = None

        for f in remaining:
            features = selected + [f]

            data = df[["Depth"] + features + [target]].dropna().sort_values(by="Depth")

            split = int(0.8 * len(data))
            train_df = data.iloc[:split]
            test_df = data.iloc[split:]

            X_train = train_df[features]
            y_train = train_df[target]

            X_test = test_df[features]
            y_test = test_df[target]

            model = XGBRegressor(n_estimators=100, max_depth=5, learning_rate=0.1)
            model.fit(X_train, y_train)

            preds = model.predict(X_test)
            mae = mean_absolute_error(y_test, preds)

            if mae < best_mae:
                best_mae = mae
                best_feature = f

        if best_feature is None:
            break

        selected.append(best_feature)
        remaining.remove(best_feature)

        print(f"Selected: {selected} | MAE: {best_mae}")