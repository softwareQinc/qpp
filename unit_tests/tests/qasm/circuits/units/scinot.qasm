OPENQASM 2.0;
include "qelib1.inc";

qreg q[1];

h q[0];
rz(1.0E-3) q[0];
h q[0];
