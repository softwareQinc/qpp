#!/bin/bash

NUM_CORES=4                      # number of cores
NQ_START=2                       # number of qubits to start with
NQ_END=24                        # number of qubits to end with
DATE=`date '+%Y-%m-%d_%H-%M-%S'` # current date
OUT_FILE=results_$DATE.csv       # results output file

echo "Stress test, see results_$DATE.csv"
echo "num_cores num_qubits"
echo "numcores, nqubits, time_seconds" > $OUT_FILE
for ((numcores=1; numcores <= NUM_CORES; numcores++));
do
    for ((n=NQ_START; n <= NQ_END; n++));
    do
        echo "$numcores         $n"
        ./build/stress_test $numcores $n >> $OUT_FILE 
    done;
done;
