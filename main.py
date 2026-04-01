from src.load_data import load_data
from src.preprocess import preprocess
from src.train_model import train_model
from src.evaluate import evaluate
from src.experiments import run_experiments
from src.plot_results import plot_results
from src.greedy_selection import greedy_selection
from src.analysis import plot_correlation

df = load_data()

X_train, X_test, y_train, y_test = preprocess(df)

model = train_model(X_train, y_train)

evaluate(model, X_test, y_test)

results = run_experiments(df)

plot_results(results)

greedy_selection(df)

plot_correlation(df)