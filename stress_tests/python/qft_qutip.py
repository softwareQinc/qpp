# QuTiP QFT stress tests

import numpy
import os
import qutip
import sys
import timeit

if len(sys.argv) != 3:
    sys.exit("Please specify the number of cores and qubits!")

num_cores = int(sys.argv[1])  # number of cores
n = int(sys.argv[2])  # number of qubits
N = 2 ** n

os.environ['OPENBLAS_NUM_THREADS'] = str(num_cores)
os.environ['MKL_NUM_THREADS'] = str(num_cores)
qutip.settings.num_cpus = num_cores


def qft_gate_sequence(n=1, swapping=True):
    """
    Quantum Fourier Transform operator on N qubits returning the gate sequence.

    Parameters
    ----------
    n: int
        Number of qubits.
    swapping: boolean
        Flag indicating sequence of swap gates to be applied at the end or not.

    Returns
    -------
    qc: instance of QubitCircuit
        Gate sequence of Hadamard and controlled rotation gates implementing
        QFT.
    """

    if n < 1:
        raise ValueError("Minimum value of n can be 1")

    qc = qutip.QubitCircuit(n)
    if n == 1:
        qc.add_gate("SNOT", targets=[0])
    else:
        for i in range(n):
            qc.add_gate("SNOT", targets=[i])
            for j in range(2, n - i + 1):
                qc.add_gate(r"CPHASE", targets=[i], controls=[i + j - 1],
                            arg_label=r"{\pi/2^{%d}}" % (j - 1),
                            arg_value=numpy.pi / (2 ** (j - 1)))
        if swapping is True:
            for i in range(n // 2):
                qc.add_gate(r"SWAP", targets=[i, n - 1 - i])

    return qc


psi = qutip.ket('0' * n)

# start timing
start_time = timeit.default_timer()

qc0 = qft_gate_sequence(n, True)
for gate in qc0.propagators():
    psi = gate * psi

elapsed = timeit.default_timer() - start_time
# end timing

print("{0}, {1}, {2}".format(num_cores, n, elapsed))
