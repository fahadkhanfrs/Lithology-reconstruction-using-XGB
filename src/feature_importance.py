import matplotlib.pyplot as plt

def plot_feature_importance(model, feature_names):

    importance = model.feature_importances_

    plt.figure()
    plt.bar(feature_names, importance)
    plt.title("Feature Importance")
    plt.xlabel("Features")
    plt.ylabel("Importance")
    plt.show()