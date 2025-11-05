import numpy as np
import matplotlib.pyplot as plt


def prob_no_conflict(n):
    A = n**2
    return (A * (A - 5) * (A-10)*(A-15))/ A**4


def plot_no_conflict(n_vals):


    pvals = [prob_no_conflict(n) for n in n_vals]

    plt.scatter(n_vals, pvals)
    plt.xlabel("L")
    plt.ylabel("P(No Conflicts)")
    plt.yscale("log")
    plt.xscale("log")
    plt.grid(which="both", axis="both", visible=True)
    plt.show()