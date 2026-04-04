import matplotlib.pyplot as plt

def plot_sampling_results():
    steps = [1, 2, 5, 10]
    acc = [0.9892, 0.8817, 0.7957, 0.7742]

    plt.figure()
    plt.plot(steps, acc, marker='o')
    plt.xlabel("Sampling Step (Depth Interval)")
    plt.ylabel("Accuracy")
    plt.title("Lithology Reconstruction vs Sampling Density")
    #plt.gca().invert_xaxis()  # optional: shows dense → sparse left to right
    plt.show()

plot_sampling_results()