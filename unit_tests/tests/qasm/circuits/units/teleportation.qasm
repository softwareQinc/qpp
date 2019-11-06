OPENQASM 2.0;
include "qelib1.inc";

qreg a[1];
qreg b[2];
creg c[1];
creg d[1];

// prepare bell state
h b[0];
cx b[0],b[1];

// prepare teleported state
u3(0.3,0.2,0.1) a[0];

// perform measurement
cx a[0],b[0];
h a[0];
measure a[0] -> c[0];
measure b[0] -> d[0];

// classically controlled corrections
if (d==1) x b[1];
if (c==1) z b[1];

// u3 adjoint, final expected state should be |0>
// u3^\dagger(a,b,c) = u3(-a,-c,-b)
u3(-0.3,-0.1,-0.2) b[1]
