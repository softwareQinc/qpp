#!/opt/homebrew/bin/python3
import sys
import matplotlib.pyplot as plt
import numpy as np
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.sans-serif": "Helvetica",
})
data_file_name = sys.argv[1]

file = open(data_file_name, "r")

lines = file.readlines()
number_of_trials = int(lines[0])
sizes = []
data_sets = []

for line in lines[1:]:
    curr = line.strip().split(' ')
    curr_set = np.array([float(n)/number_of_trials for n in curr[1:]])
    sizes.append(int(curr[0]))
    data_sets.append(curr_set)

data_means = np.array([np.mean(np.log2(data_set)) for data_set in data_sets])
data_stds = np.array([np.std(np.log2(data_set)) for data_set in data_sets])
plt.errorbar(np.log2(sizes), data_means, data_stds, marker='o', linestyle='')
plt.title(r"Average apply() run time for Two Qubit Gate on $2^n$ qubits"
          "(10000 trials)")
plt.xlabel(r"$n$")
plt.ylabel(r"$\log_2(t)$ (ms)")
plt.show()
