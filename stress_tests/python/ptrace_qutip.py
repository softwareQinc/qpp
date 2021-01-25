#!/usr/bin/env python3

# QuTiP partial trace stress tests

import os
import qutip
import sys
import timeit

if len(sys.argv) != 3:
    sys.exit("Please specify the maximum number of cores and qubits!")

num_cores = int(sys.argv[1])  # max number of cores
n = int(sys.argv[2])          # number of qubits
D = 2 ** n                    # total dimension

os.environ['OPENBLAS_NUM_THREADS'] = str(num_cores)
os.environ['MKL_NUM_THREADS'] = str(num_cores)
qutip.settings.num_cpus = num_cores


result = qutip.rand_herm(D, dims=[[2] * n, [2] * n])

# start timing
start_time = timeit.default_timer()

# partial trace over the first qubit
result = result.ptrace(range(1, n))

elapsed = timeit.default_timer() - start_time
# end timing

print("{0}, {1}, {2}".format(num_cores, n, elapsed))
