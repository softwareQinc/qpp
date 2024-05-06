#!/opt/homebrew/bin/python3

import sys
import subprocess

experiment_name = sys.argv[1]
number_of_trials = int(sys.argv[2])
runs_per_trial = int(sys.argv[3])
max_size_exp = int(sys.argv[4])

f = open(experiment_name, 'w')
f.write(str(runs_per_trial))
f.write("\n")

sizes = [2**n for n in range(max_size_exp+1)]

for size in sizes:
    f.write(str(size))
    result_string = " "
    for trial in range(number_of_trials):
        result = subprocess.run(['./build/single_qubit_gate_apply_test',
                                 str(runs_per_trial),
                                 str(size)], capture_output=True, text=True)
        result_string += result.stdout.strip() + " "
    f.write(result_string)
    f.write("\n")
