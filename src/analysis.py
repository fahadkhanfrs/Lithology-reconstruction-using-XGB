import seaborn as sns
import matplotlib.pyplot as plt

def plot_correlation(df):
    cols = ["GR", "SP", "RILD", "CNLS", "RHOB"]

    corr = df[cols].dropna().corr()

    plt.figure()
    sns.heatmap(corr, annot=True)
    plt.title("Correlation Matrix")
    plt.show()