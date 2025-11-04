"""
Take in a number of files by name which contain lines:
`Convergence for step X: visited = V, energy = E, prop_visited = P, avg_energy_per_spin = A`
produce two plots: P vs step and A vs step, containing all files' data as different series.
"""

import matplotlib.pyplot as plt
import re
import sys
import numpy as np
import os
import argparse

max_step = 200000

expected_energy_convergence = -1.414


def parse_file(filename):
    steps = []
    prop_visited = []
    avg_energy_per_spin = []
    pattern = re.compile(r'Convergence for step (\d+): visited = \d+, energy = [\d\.\-e]+, prop_visited = ([\d\.e\-]+), avg_energy_per_spin = ([\d\.e\-]+)')
    with open(filename, 'r') as f:
        for line in f:
            match = pattern.search(line)
            if match:
                if int(match.group(1)) > max_step:
                    continue
                steps.append(int(match.group(1)))
                prop_visited.append(float(match.group(2)))
                avg_energy_per_spin.append(float(match.group(3)))
    return steps, prop_visited, avg_energy_per_spin

def plot_convergence(filenames):
    plt.figure(figsize=(12, 6))

    # Plot prop_visited
    plt.subplot(1, 2, 1)
    for filename in filenames:
        steps, prop_visited, _ = parse_file(filename)
        plt.plot(steps, prop_visited, label=os.path.basename(filename))
    # Draw horizontal line at 1.0
    plt.axhline(y=1.0, color='r', linestyle='--', label='Full Coverage (1.0)')
    plt.xlabel('Step')
    plt.xscale('log')
    plt.ylabel('Proportion Visited')
    plt.title('Proportion of Lattice Visited Over Time')
    plt.legend()
    plt.grid(True)

    # Plot avg_energy_per_spin
    plt.subplot(1, 2, 2)
    for filename in filenames:
        steps, _, avg_energy_per_spin = parse_file(filename)
        plt.plot(steps, avg_energy_per_spin, label=os.path.basename(filename))
    # Draw horizontal line at expected energy convergence
    plt.axhline(y=expected_energy_convergence, color='r', linestyle='--', label=f'Expected Energy ({expected_energy_convergence})')
    plt.xlabel('Step')
    plt.xscale('log')
    plt.ylabel('Average Energy per Spin')
    plt.title('Average Energy per Spin Over Time')
    plt.legend()
    plt.grid(True)

    plt.tight_layout()
    plt.savefig('convergence_plots.png')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot convergence data from files.')
    parser.add_argument('files', metavar='F', type=str, nargs='+', help='Files containing convergence data')
    args = parser.parse_args()
    plot_convergence(args.files)