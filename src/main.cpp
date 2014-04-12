/*
 * File:   main.cpp
 * Author: vlad
 * 
 * Created on December 12, 2013, 10:38 PM
 */

#include "qpp.h"
#include "classes/qudit.h"
#include "internal.h"

#include <functional>
//#include "matlab.h" // support for MATLAB

// TODO: use .data() raw pointer instead of looping
// TODO: optimize syspermute: ?Eigen::Map?
// TODO: IMPORTANT Rewrite partial trace without syspermute
// TODO: further parallelization
// TODO: testing

using namespace std;
using namespace qpp;
using namespace qpp::types;

int main()
{
	cout << "Starting qpp..." << endl;

	// output format
	// cout << std::scientific;
	cout << std::fixed; // use fixed format for nice formatting
	cout << std::setprecision(4); // only for fixed or scientific modes

	ket zd0(3);
	zd0 << 1, 0, 0;
	Qudit q(proj(zd0));
	cmat U = randU(3);

	// von Neumann projective measurement
	cout << q.measure(gt.Zd(3), true) << endl << endl;

	cout << q.measure() << endl;
	cout << q.measure(true) << endl;
	cout << endl;
	cout << q.measure(gt.Xd(3)) << endl;
	cout << q.measure(gt.Xd(3)) << endl;
	cout << q.measure(gt.Xd(3)) << endl;
	cout << q.measure(gt.Xd(3), true) << endl;
	cout << q.measure(gt.Xd(3)) << endl;
	cout << q.measure(gt.Xd(3)) << endl;

	// TESTING

	// Bell state generator
	cout << endl << "Bell state generator: " << endl;
	cmat circuit;
	circuit = gt.CTRL(gt.X, { 0 }, { 1 }, 2) * expandout(gt.H, 0, { 2, 2 });
	cmat z0(2, 1);
	z0 << 1, 0;
	cmat z1(2, 1);
	z1 << 0, 1;
	cmat input = kron(z0, z0);
	cmat output = circuit * input;
	cout << "Circuit matrix representation: " << endl;
	displn(circuit);
	cout << endl << "Output (|Bell_0> state) of the circuit on |00>: " << endl;
	displn(output);

	// 3-qubit repetion code
	cout << endl << "3-qubit repetition code: " << endl;
	cmat rep;
	rep = gt.CTRL(gt.X, { 0 }, { 2 }, 3) * gt.CTRL(gt.X, { 0 }, { 1 }, 3);
	input = kronlist<cmat>( { z1, z0, z0 });
	output = rep * input;
	cout << "Circuit acting on |000> produces |111>. Check: " << endl;
	displn(output);

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
	cmat evalsH = hevals(rH);
	cmat evectsH = hevects(rH);
	cmat spec = cmat::Zero(D, D);
	for (size_t i = 0; i < D; i++)
		spec += evalsH(i) * proj(evectsH.col(i));
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
	std::vector<cmat> Ks = randKraus(nk, d);

	cmat rho_in = randrho(d); // input state
	cmat rho_out = channel(rho_in, Ks); // output state

	cout << "Computing its Choi matrix..." << endl;
	cmat choim = choi(Ks);
	cout << "Choi matrix:" << endl;
	displn(choim);
	cout << endl << "The eigenvalues of the Choi matrix are: " << endl;
	displn(transpose(hevals(choim)));
	cout << endl << "Their sum is: " << sum(hevals(choim)).real() << endl;
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
			reshape(smat * reshape(transpose(rho_in), d * d, 1), d, d));
	cout << norm(rho_out - rho_out2) << endl;

	// statistics tests
	cout << endl << "Statistics tests." << endl;
	std::vector<cplx> ampl = { 1. + ct::ii, 1. - ct::ii };
	cmat va(1, 4);
	va << 0.1, 1, 1. + ct::ii, 1. + 2. * ct::ii;
	stat::DiscreteDistributionFromComplex dc(va);
	cout << "The probabilities are: ";
	displnSTL(dc.probabilities(), ", ", "{", "}");

	// other tests
	cout << endl << "Timing tests..." << endl;
	size_t n = 12; // number of qubits
	size_t N = std::pow(2, n);
	vector<size_t> dims; // local dimensions
	for (size_t i = 0; i < n; i++)
		dims.push_back(2);
	cout << "n = " << n << " qubits, matrix size " << N << " x " << N << "."
			<< endl;

	// TIMING
	Timer t, total; // start the timer, automatic tic() in the constructor

	// matrix initialization
	cout << endl << "Matrix initialization timing." << endl;
	cmat randcmat = cmat::Random(N, N);
	t.toc(); // read the time
	cout << "Took " << t << " seconds." << endl;

	// lazy matrix product
	cout << endl << "Lazy matrix product timing." << endl;
	t.tic();
	auto lazyprod = randcmat * randcmat; // lazyprod has type GenMatProduct
	t.toc(); // read the time
	cout << "Took " << t << " seconds." << endl;

	// matrix product
	cout << endl << "Matrix product timing." << endl;
	t.tic(); // reset the chronometer
	cmat prodmat = randcmat * randcmat; // explicit cmat now
	t.toc(); // read the time
	cout << "Took " << t << " seconds." << endl;

	// ptrace2 timing
	cout << endl << "ptrace2 timing." << endl;
	t.tic(); // reset the chronometer
	// trace away half of the qubits
	ptrace2(randcmat, { (size_t) std::sqrt(N), (size_t) std::sqrt(N) });
	t.toc(); // read the time
	cout << "Took " << t << " seconds." << endl;

	// ptrace SLOW SLOW SLOW (because of syspermute, do it without)
	cout << endl << "ptrace timing." << endl;
	vector<size_t> subsys_ptrace = { 0 }; // trace away the first qubit
	cout << "Subsytem(s): ";
	displnSTL(subsys_ptrace, ", ");
	t.tic();
	ptrace(randcmat, subsys_ptrace, dims);
	t.toc();
	cout << "Took " << t << " seconds." << endl;

	// ptranspose
	cout << endl << "ptranspose timing." << endl;
	vector<size_t> subsys_ptranspose; // partially transpose all subsystems
	for (size_t i = 0; i < n; i++)
		subsys_ptranspose.push_back(i);
	cout << "Subsytem(s): ";
	displnSTL(subsys_ptranspose, ", ");
	t.tic();
	ptranspose(randcmat, subsys_ptranspose, dims);
	t.toc();
	cout << "Took " << t << " seconds." << endl;

	// syspermute SLOW SLOW SLOW
	cout << endl << "syspermute timing." << endl;
	vector<size_t> perm; // left-shift all subsystems by 1
	for (size_t i = 0; i < n; i++)
		perm.push_back((i + 1) % n);
	cout << "Subsytem(s): ";
	displnSTL(perm, ", ");
	t.tic();
	syspermute(randcmat, perm, dims);
	t.toc();
	cout << "Took " << t << " seconds." << endl;

	// END TIMING
	total.toc(); // read the total running time
	cout << endl << "Total time: " << total.seconds() << " seconds.";

	cout << endl << "Exiting qpp..." << endl;
}
