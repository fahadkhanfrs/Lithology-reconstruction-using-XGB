from xgboost import XGBRegressor

def train_model(X, y):
    model = XGBRegressor(
        n_estimators=100,
        max_depth=5,
        learning_rate=0.1
    )

    model.fit(X, y)
    return model