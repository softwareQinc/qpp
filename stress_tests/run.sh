#!/bin/bash

NUMCORES=8
NQ_START=2
NQ_END=20
DATE=`date '+%Y-%m-%d_%H:%M:%S'`
OUTFILE=results_$DATE.csv

echo "QFT stress test, see results_$DATE.csv"
echo "num_cores num_qubits"
echo "numcores, nqubits, time_seconds" > $OUTFILE
for ((numcores=1; numcores < NUMCORES; numcores++));
do
    for ((n=NQ_START; n < NQ_END; n++));
    do
        echo "$numcores         $n"
        ./build/qft $numcores $n >> $OUTFILE 
    done;
done;
