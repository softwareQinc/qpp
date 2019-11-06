OPENQASM 2.0;
include "qelib1.inc";

qreg a[1];
qreg b[2];
creg c[2];

U(0.1, 0.2, 0.3) a[0];

// proportional to X
U(pi, 0, pi) b[0];
CX b[0],b[1];

measure b[0] -> c[0];
measure b[1] -> c[1];
