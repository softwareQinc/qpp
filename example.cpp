/*
 * File:   example.cpp
 * Author: vlad
 * 
 * Created on December 12, 2013, 10:38 PM
 */

#include "qpp.h"

//#include "matlab.h" // support for MATLAB

// TODO: testing
// TODO: write documentation
// TODO: more robust exception parameter checking

using namespace std;
using namespace qpp;
using namespace qpp::types;
using namespace qpp::ct;

cplx pow3(const cplx& z)
{
	return std::pow(z, 3);
}

int main()
{
	cout << "Starting qpp..." << endl;

	// output format
	//cout << std::scientific;
	cout << std::fixed; // use fixed format for nice formatting
	cout << std::setprecision(4); // only for fixed or scientific modes

	 // TESTING

	 // testing channel and gate
	 cout << endl << "Testing channel(...) and gate(...)." << endl;
	 cmat rho = randrho(16);
	 cmat K = kron(gt.Id2, gt.X, gt.Y, gt.Z);
	 vector<size_t> p = randperm(4); // permutation
	 cout << "Permutation: ";
	 displn(p, ", ");
	 vector<size_t> invp = invperm(p); // inverse permutation
	 cout << "Inverse permutation: ";
	 displn(invp, ", ");
	 cmat r1 = channel(rho, { K }, p, { 2, 2, 2, 2 });
	 cmat r2 = syspermute(channel(syspermute(rho, p, { 2, 2, 2, 2 }), { K }, { 0,
	 1, 2, 3 }, { 2, 2, 2, 2 }), invp, { 2, 2, 2, 2 });
	 cout << norm(r1 - r2) << endl << endl;

	 r1 = gate(rho, K, p, { 2, 2, 2, 2 });
	 r2 = syspermute(
	 gate(syspermute(rho, p, { 2, 2, 2, 2 }), K, { 0, 1, 2, 3 }, { 2, 2,
	 2, 2 }), invp, { 2, 2, 2, 2 });
	 cout << norm(r1 - r2) << endl << endl;

	 displn(channel(prj(mket( { 0, 1 })), { gt.CNOTab }, { 1, 0 }, { 2, 2 }));
	 cout << endl;
	 displn(gate(mket( { 0, 0 }), gt.CNOTab, { 0, 1 }, { 2, 2 }));

	 // quantum teleportation
	 cout << endl << "Qudit teleportation." << endl;
	 ket psi = randket(2); // a random state;
	 cout << "|psi><psi|:" << endl;
	 displn(prj(psi));
	 cmat telecircuit = expandout(gt.H, { 0 }, { 2, 2, 2 })
	 * gt.CTRL(gt.X, { 0 }, { 1 }, 3);
	 ket psiin = kron(psi, st.b00); // input state
	 ket psiout = telecircuit * psiin; // output state before measurement
	 // measure Alice's qubits, measurement results are 1 0
	 psiout = kron(prj(st.z1), prj(st.z0), gt.Id2) * psiout;
	 // apply correction
	 psiout = expandout(powm(gt.Z, 1) * powm(gt.X, 0), { 2 }, { 2, 2, 2 })
	 * psiout;
	 // not necessary to normalize, prj() takes care of it below
	 cmat rhoout = ptrace(prj(psiout), { 0, 1 }, { 2, 2, 2 });
	 cout << endl << "Teleported state:" << endl;
	 displn(rhoout);
	 cout << "Difference in norm: " << norm(prj(psi) - rhoout) << endl;

	 // qudit measurements
	 cout << endl << "Qudit measurements." << endl;
	 cout << "Initially in state |0><0|." << endl;
	 ket zd0(3);
	 zd0 << 1, 0, 0;
	 Qudit q(prj(zd0));
	 cout << "Measuring Z operator non-destructively. Results:" << endl;
	 cout << q.measure() << endl;
	 cout << q.measure() << endl;
	 cout << q.measure() << endl;
	 cout << "Measuring X operator non-destructively. Results:" << endl;
	 cout << q.measure(gt.Xd(3)) << endl;
	 cout << q.measure(gt.Xd(3)) << endl;
	 cout << q.measure(gt.Xd(3)) << endl;
	 // von Neumann projective measurement
	 cout << "Measuring X operator destructively (collapse). Results:" << endl;
	 cout << q.measure(gt.Xd(3), true) << endl;
	 cout << q.measure(gt.Xd(3)) << endl;
	 cout << q.measure(gt.Xd(3)) << endl;
	 cout << "Finally measuring Z operator destructively. Results:" << endl;
	 cout << q.measure(true) << endl;
	 cout << q.measure() << endl;
	 cout << q.measure() << endl;
	 cout << "Final state of qudit:" << endl;
	 displn(q.getRho());

	 // Bell state generator
	 cout << endl << "Bell state generator: " << endl;
	 cmat circuit;
	 circuit = gt.CTRL(gt.X, { 0 }, { 1 }, 2) * expandout(gt.H, 0, { 2, 2 });
	 cmat input = kron(st.z0, st.z0);
	 cmat output = circuit * input;
	 cout << "Circuit matrix representation: " << endl;
	 displn(circuit);
	 cout << endl << "Output (|Bell_0> state) of the circuit on |00>: " << endl;
	 displn(output);

	 // 3-qubit repetion code
	 cout << endl << "3-qubit repetition code: " << endl;
	 cmat rep;
	 rep = gt.CTRL(gt.X, { 0 }, { 2 }, 3) * gt.CTRL(gt.X, { 0 }, { 1 }, 3);
	 input = kron(st.z1, st.z0, st.z0);
	 output = rep * input;
	 cout << "Circuit acting on |000> produces |111>. Check: " << endl;
	 displn(output);

	 // functor test
	 cout << endl << "Functor z^3 acting on:" << endl;
	 cmat a(2, 2);
	 a << 1, 2, 3, 4;
	 displn(a);
	 cout << "Result (with lambda):" << endl;
	 // functor z^3 componentwise, specify OutputScalar and Derived for lambdas
	 displn(cwise<cplx, cmat>(a, [](const cplx& z)->cplx
	 {	return z*z*z;}));
	 cout << "Result (with proper function):" << endl;
	 // automatic type deduction for proper functions
	 displn(cwise(a, &pow3));

	 // Gram-Schmidt
	 cout << endl << "Gram-Schmidt on matrix:" << endl;
	 cmat A(3, 3);
	 A << 1, 1, 0, 0, 2, 0, 0, 0, 0;
	 displn(A);
	 cmat Ags = grams(A);
	 cout << endl << "Result:" << endl;
	 displn(Ags);
	 cout << endl << "Projector is:" << endl;
	 displn(Ags * adjoint(Ags));

	 // spectral decomposition test
	 cout << endl << "Spectral decomposition tests." << endl;
	 size_t D = 4;
	 cmat rH = randH(D);
	 dmat evalsH = hevals(rH);
	 cmat evectsH = hevects(rH);
	 cmat spec = cmat::Zero(D, D);
	 for (size_t i = 0; i < D; i++)
	 spec += evalsH(i) * prj((cmat) evectsH.col(i));
	 cout << "Original matrix: " << endl;
	 displn(rH);
	 cout << endl << "Reconstructed from spectral decomposition: " << endl;
	 displn(spec);
	 cout << "Difference in norm: " << norm(spec - rH) << endl;

	 // channel tests
	 cout << endl << "Channel tests." << endl;
	 size_t nk = 10, d = 2; // nk Kraus on d-dimensional system
	 cout << "Generating a random channel with " << nk
	 << " Kraus operators on a " << d << " dimensional space..." << endl;
	 std::vector<cmat> Ks = randkraus(nk, d);

	 cmat rho_in = randrho(d); // input state
	 cmat rho_out = channel(rho_in, Ks); // output state

	 cout << "Computing its Choi matrix..." << endl;
	 cmat choim = choi(Ks);
	 cout << "Choi matrix:" << endl;
	 displn(choim);
	 cout << endl << "The eigenvalues of the Choi matrix are: " << endl;
	 displn(transpose(hevals(choim)));
	 cout << endl << "Their sum is: " << sum(hevals(choim)) << endl;
	 std::vector<cmat> Kperps = choi2kraus(choim);
	 cout << endl << "The Kraus rank of the channel is: " << Kperps.size()
	 << endl;
	 cmat rho_out1 = channel(rho_in, Kperps);
	 cout << endl << "Difference in norm on output states: "
	 << norm(rho_out1 - rho_out) << endl;
	 cout << endl << "Superoperator matrix:" << endl;
	 cmat smat = super(Ks);
	 displn(smat);
	 cout << endl << "The eigenvalues of the superoperator matrix are: " << endl;
	 cmat evalsupop = evals(smat);
	 displn(transpose(evalsupop));
	 cout << endl << "Their absolute values are: " << endl;
	 for (size_t i = 0; i < (size_t) evalsupop.size(); i++)
	 cout << std::abs((cplx) evalsupop(i)) << " ";
	 cout << endl << endl << "Diference in norm for superoperator action: ";
	 cmat rho_out2 = transpose(
	 (cmat) reshape(smat * reshape(transpose(rho_in), d * d, 1), d, d));
	 cout << norm(rho_out - rho_out2) << endl;

	 // statistics tests
	 cout << endl << "Statistics tests." << endl;
	 std::vector<cplx> ampl = { 1. + ii, 1. - ii };
	 cmat va(1, 4);
	 va << 0.1, 1, 1. + ii, 1. + 2. * ii;
	 DiscreteDistributionAbsSquare dc(va);
	 cout << "The probabilities are: ";
	 displn(dc.probabilities(), ", ", "{", "}");

	 // 	// TIMING tests
	 cout << endl << "Timing tests..." << endl;
	 size_t n = 12; // number of qubits
	 size_t N = std::pow(2, n);
	 vector<size_t> dims; // local dimensions
	 for (size_t i = 0; i < n; i++)
	 dims.push_back(2);
	 cout << "n = " << n << " qubits, matrix size " << N << " x " << N << "."
	 << endl;

	 // matrix initialization
	 cout << endl << "Matrix initialization timing." << endl;
	 // start the timer, automatic tic() in the constructor
	 Timer t, total;
	 cmat randcmat = cmat::Random(N, N);
	 t.toc(); // read the time
	 cout << "Took " << t << " seconds." << endl;

	 // lazy matrix product
	 cout << endl << "Lazy matrix product timing." << endl;
	 t.tic();
	 auto lazyprod = randcmat * randcmat; // lazyprod has type GenMatProduct
	 t.toc(); // read the time
	 cout << "Took " << t << " seconds." << endl;

	 // ptrace1 timing
	 cout << endl << "ptrace1 timing." << endl;
	 t.tic(); // reset the chronometer
	 // trace away half of the qubits
	 ptrace1(randcmat, { (size_t) std::sqrt(N), (size_t) std::sqrt(N) });
	 t.toc(); // read the time
	 cout << "Took " << t << " seconds." << endl;

	 // ptrace2 timing
	 cout << endl << "ptrace2 timing." << endl;
	 t.tic(); // reset the chronometer
	 // trace away half of the qubits
	 ptrace2(randcmat, { (size_t) std::sqrt(N), (size_t) std::sqrt(N) });
	 t.toc(); // read the time
	 cout << "Took " << t << " seconds." << endl;

	 // ptrace
	 cout << endl << "ptrace timing." << endl;
	 vector<size_t> subsys_ptrace = { 0 };
	 cout << "Subsytem(s): ";
	 displn(subsys_ptrace, ", ");
	 t.tic();
	 ptrace(randcmat, subsys_ptrace, dims);
	 t.toc();
	 cout << "Took " << t << " seconds." << endl;

	 // ptranspose
	 cout << endl << "ptranspose timing." << endl;
	 vector<size_t> subsys_ptranspose; // partially transpose n-1 subsystems
	 for (size_t i = 0; i < n - 1; i++)
	 subsys_ptranspose.push_back(i);
	 cout << "Subsytem(s): ";
	 displn(subsys_ptranspose, ", ");
	 t.tic();
	 ptranspose(randcmat, subsys_ptranspose, dims);
	 t.toc();
	 cout << "Took " << t << " seconds." << endl;

	 // syspermute
	 cout << endl << "syspermute timing." << endl;
	 vector<size_t> perm; // left-shift all subsystems by 1
	 for (size_t i = 0; i < n; i++)
	 perm.push_back((i + 1) % n);
	 cout << "Subsytem(s): ";
	 displn(perm, ", ");
	 t.tic();
	 syspermute(randcmat, perm, dims);
	 t.toc();
	 cout << "Took " << t << " seconds." << endl;

	 // matrix product
	 cout << endl << "Matrix product timing." << endl;
	 t.tic(); // reset the chronometer
	 cmat prodmat = randcmat * randcmat; // explicit cmat now
	 t.toc(); // read the time
	 cout << "Took " << t << " seconds." << endl;

	 // END TIMING
	 total.toc(); // read the total running time
	 cout << endl << "Total time: " << total.seconds() << " seconds.";

	cout << endl << "Exiting qpp..." << endl;
}
