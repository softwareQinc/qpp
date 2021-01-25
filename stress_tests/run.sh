#!/bin/bash

# $1 - path to stress test executable
# $2 - directory where results will be collected

if [ $# -lt 2 ]; then
    echo "Usage: $0 <path_to_stress_test_executable> <results_output_directory>" >&2
    exit 1
fi

stress_test="$1"
results_dir="$2"

max_cores=8                                  # number of cores
nq_start=2                                   # number of qubits to start with
nq_end=24                                    # number of qubits to end with

current_date=$(date '+%Y-%m-%d_%a_%H-%M-%S') # current date

if [ ! -d "$results_dir" ]; then
    mkdir "$2"
fi
out_file="$results_dir"/results_$current_date.csv # result output file

filename=$(basename -- "$stress_test")
extension="${filename##*.}"
if [ "$extension" == "py" ]; then
    interpreter="python3 $stress_test"
else
    interpreter=$stress_test
fi

echo "Stress test, see $results_dir/results_$current_date.csv"
echo "$stress_test" >"$out_file"
echo "num_cores, num_qubits, time_seconds" | tee -a "$out_file"
for ((num_cores = 1; num_cores <= max_cores; num_cores++)); do
    for ((n = nq_start; n <= nq_end; n++)); do
        $interpreter "$num_cores" "$n" | tee -a "$out_file"
    done
done
