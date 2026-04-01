from sklearn.metrics import mean_absolute_error

def evaluate(model, X, y):
    preds = model.predict(X)
    mae = mean_absolute_error(y, preds)
    print(f"MAE: {mae}")