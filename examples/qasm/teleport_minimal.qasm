// quantum teleportation minimal example
OPENQASM 2.0;
include "qelib1.inc";

// declarations, 3 qubits, 2 classical bits
qreg q[3];
creg c0[1];
creg c1[1];

// initial state
u3(0.3,0.2,0.1) q[0];

// teleportation circuit
h q[1];
cx q[1],q[2];
barrier q;
cx q[0],q[1];
h q[0];

// measurements (non-destructive)
measure q[0] -> c0[0];
measure q[1] -> c1[0];

// classically controlled corrections 
if(c1==1) x q[2];
if(c0==1) z q[2];

// u3 adjoint, final expected state should be |0>
// u3^\dagger(a,b,c) = u3(-a,-c,-b)
u3(-0.3,-0.1,-0.2) q[2]
