/*
 * File:   main.cpp
 * Author: vlad
 * 
 * Created on December 12, 2013, 10:38 PM
 */

#include "qpp.h"

//#include "matlab.h" // support for MATLAB

// TODO: expandout function
// TODO: ip (inner product) function, make it general to return matrices
// TODO: use .data() raw pointer instead of looping
// TODO: optimize syspermute: ?Eigen::Map?
// TODO: IMPORTANT Rewrite partial trace without syspermute
// TODO: further parallelization

using namespace std;

using namespace qpp;
using namespace qpp::types;

int main()
{
	_init(); // ALWAYS call _init() at the beginning of main()

	cout << "Starting qpp..." << endl;

	// output format
	cout << std::fixed; // use fixed format for nice formatting
	// cout << std::scientific;
	cout << std::setprecision(4); // only for fixed or scientific modes

	// spectral decomposition test
	cout << endl << "Spectral decomposition tests." << endl;
	size_t D = 4;
	cmat rH = randH(D);
	cmat evalsH = hevals(rH);
	cmat evectsH = hevects(rH);
	cmat spec = cmat::Zero(D, D);
	for (size_t i = 0; i < D; i++)
		spec += evalsH(i) * proj((cmat) (evectsH.col(i)));
	cout << "Original matrix: " << endl;
	displn(rH);
	cout << endl << "Reconstructed from spectral decomposition: " << endl;
	displn(spec);
	cout << "Difference in norm: " << norm((cmat) (spec - rH)) << endl;
	cmat tmp(2, 1);
	tmp << 2, 2; // some un-normalized ket
	cout << endl << "Projector constructed from un-normalized ket: " << endl;
	displn(proj(tmp));
	cout << endl << "Now the un-normalized corresponding dyad: " << endl;
	displn(dya(tmp));
	cout << endl;

	// statistics tests
	cout << endl << "Statistics tests." << endl;
	std::vector<cplx> ampl = { 1. + ct::ii, 1. - ct::ii };
	cmat va(1, 4);
	va << 0.1, 1, 1. + ct::ii, 1. + 2. * ct::ii;
	stat::DiscreteDistributionFromComplex dc(va);
	cout << "The probabilities are: [";
	disp(dc.probabilities(), ", ");
	cout << "]" << endl;

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
	Timer t, total;  // start the timer, automatic tic() in the constructor

	// matrix initialization
	cout << endl << "Matrix initialization timing." << endl;
	cmat randcmat = cmat::Random(N, N);
	t.toc(); // read the time
	cout << "Took " << t << " seconds." << endl;

	// lazy matrix product
	cout << endl << "Lazy matrix product timing." << endl;
	auto lazyprod = randcmat * randcmat;
	t.toc(); // read the time
	cout << "Took " << t << " seconds." << endl;

	// matrix product
	cout << endl << "Matrix product timing." << endl;
	t.tic(); // reset the chronometer
	cmat prodmat;
	prodmat = randcmat * randcmat; // need this (otherwise lazy evaluation)
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
	cout << "Subsytem(s): [";
	disp(subsys_ptrace, ", ");
	cout << "]" << endl;
	t.tic();
	ptrace(randcmat, subsys_ptrace, dims);
	t.toc();
	cout << "Took " << t << " seconds." << endl;

	// ptranspose
	cout << endl << "ptranspose timing." << endl;
	vector<size_t> subsys_ptranspose; // partially transpose all subsystems
	for (size_t i = 0; i < n; i++)
		subsys_ptranspose.push_back(i);
	cout << "Subsytem(s): [";
	disp(subsys_ptranspose, ", ");
	cout << "]" << endl;
	t.tic();
	ptranspose(randcmat, subsys_ptranspose, dims);
	t.toc();
	cout << "Took " << t << " seconds." << endl;

	// syspermute SLOW SLOW SLOW
	cout << endl << "syspermute timing." << endl;
	vector<size_t> perm; // left-shift all subsystems by 1
	for (size_t i = 0; i < n; i++)
		perm.push_back((i + 1) % n);
	cout << "Subsytem(s): [";
	disp(perm, ", ");
	cout << "]" << endl;
	t.tic();
	syspermute(randcmat, perm, dims);
	t.toc();
	cout << "Took " << t << " seconds." << endl;

	// END TIMING
	total.toc(); // read the total running time
	cout << endl << "Total time: " << total.seconds() << " seconds.";

	cout << endl << "Exiting qpp..." << endl;
}
