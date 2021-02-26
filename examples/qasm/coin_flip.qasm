// simple coin flipping, run it, e.g., with 'qasm < coin_flip.qasm 1000'
OPENQASM 2.0;
include "qelib1.inc";

qreg q[1];
creg c[1];

h q[0];
measure q[0] -> c[0];
