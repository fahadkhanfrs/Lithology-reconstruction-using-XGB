import matplotlib.pyplot as plt

def plot_results(results):
    costs = [r[2] for r in results]
    maes = [r[1] for r in results]
    labels = [r[0] for r in results]

    plt.figure()
    plt.scatter(costs, maes)

    for i, label in enumerate(labels):
        plt.text(costs[i], maes[i], label)

    plt.xlabel("Cost")
    plt.ylabel("MAE")
    plt.title("Accuracy vs Cost")

    plt.show()