from src.load_data import load_data
from src.preprocess import preprocess
from src.train_model import train_model
from src.evaluate import evaluate
from src.experiments import run_experiments
from src.plot_results import plot_results
from src.greedy_selection import greedy_selection
from src.analysis import plot_correlation
from src.litholog_model import (
    load_litholog, expand_layers, encode_facies,
    train_and_evaluate, uncertainty_sampling_experiment, add_context,
    uniform_sampling_experiment, hybrid_sampling_experiment, plot_reconstruction,
    plot_error_analysis
)

RUN_BASELINE = False
RUN_EXPERIMENTS = True

litho = load_litholog()
expanded = expand_layers(litho)
expanded, mapping = encode_facies(expanded)
expanded = add_context(expanded)


if RUN_BASELINE:
    df = load_data()

    X_train, X_test, y_train, y_test = preprocess(df)

    model = train_model(X_train, y_train)

    evaluate(model, X_test, y_test)

    results = run_experiments(df)

    plot_results(results)

    greedy_selection(df)

    plot_correlation(df)

    for step in [1, 2, 5, 10]:
        train_and_evaluate(expanded, step)


if RUN_EXPERIMENTS:
    print("\n--- Uniform ---")
    acc, full_df, sampled, model, class_mapping = uniform_sampling_experiment(expanded, budget=30)

    print("\n--- Error Analysis (Uniform) ---")
    plot_error_analysis(full_df, sampled, model, class_mapping, method="Uniform Sampling")

    print("\n--- Hybrid ---")
    acc, full_df, sampled, model, class_mapping = hybrid_sampling_experiment(expanded, initial_step=10, budget=20)

    print("\n--- Error Analysis (Hybrid) ---")
    plot_error_analysis(full_df, sampled, model, class_mapping, method="Hybrid Sampling")
