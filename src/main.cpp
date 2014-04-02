/*
 * File:   main.cpp
 * Author: vlad
 * 
 * Created on December 12, 2013, 10:38 PM
 */

#include <iostream>

#include "qpp.h"
#include "internal.h"
//#include "matlab.h" // support for MATLAB

// TODO: expandout function
// TODO: dyad function
// TODO: proj (dya) function
// TODO: ip (inner product function) function, make it general to return matrices
// TODO: Error class
// TODO: change all for(s) to column major order
// TODO: use .data() raw pointer instead of looping
// TODO: optimize syspermute: column major ordering + ?Eigen::Map?
// TODO: Rewrite partial trace without syspermute IMPORTANT!!!!
// TODO: check agains zero-size matrices (i.e. non-initialized)

using namespace std;

using namespace qpp;
using namespace qpp::types;

//int main(int argc, char **argv)
int main()
{
	cout << "Starting qpp..." << endl;
	_init();

	// Display formatting
	cout << std::fixed; // use fixed format for nice formatting
//	cout << std::scientific;
	cout << std::setprecision(4); // only for fixed or scientific modes

	// statistics tests
	cout << endl << "Statistics tests." << endl;
	std::vector<cplx> ampl = { 1. + ct::ii, 1. - ct::ii };
	cmat va(1, 3);
	va << 1, 1. + ct::ii, 1. + 2. * ct::ii;
	stat::DiscreteDistributionFromComplex dc(va);
	cout << "The probabilities are: [";
	internal::_disp_container(dc.probabilities());
	cout << "]" << endl;

	// other tests
	size_t n = 12; // number of qubits
	size_t N = std::pow(2, n);
	vector<size_t> dims; // local dimensions
	for (size_t i = 0; i < n; i++)
		dims.push_back(2);
	cout << "n = " << n << " qubits, matrix size " << N << " x " << N << "."
			<< endl;

	// TIMING
	Timer t, total;  // start the timer, automatic tic() in the constructor

	// Matrix initialization
	cout << endl << "Matrix initialization timing." << endl;
	cmat randcmat = cmat::Random(N, N);
	t.toc(); // read the time
	cout << "Took " << t << " seconds." << endl;

	// Lazy matrix product
	cout << endl << "Lazy matrix product timing." << endl;
	auto b = randcmat * randcmat;
	t.toc(); // read the time
	cout << "Took " << t << " seconds." << endl;

	// Matrix product
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
	internal::_disp_container(subsys_ptrace);
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
	internal::_disp_container(subsys_ptranspose);
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
	internal::_disp_container(perm);
	cout << "]" << endl;
	t.tic();
	syspermute(randcmat, perm, dims);
	t.toc();
	cout << "Took " << t << " seconds." << endl;

	total.toc(); // read the total running time
	cout << endl << "Total time: " << total.seconds() << " seconds." << endl;
	// END TIMING

	cout << "Exiting qpp..." << endl;
}
