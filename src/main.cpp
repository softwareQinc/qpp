/*
 * File:   main.cpp
 * Author: vlad
 * 
 * Created on December 12, 2013, 10:38 PM
 */

#include <iostream>

#include "qpp.h"
//#include "matlab.h" // support for MATLAB

// TODO: expandout function
// TODO: dyad function
// TODO: proj (dya) function
// TODO: ip (inner product function) function, make it general to return matrices
// TODO: Error class
// TODO: change all for(s) to column major order
// TODO: use .data() raw pointer instead of looping

using namespace std;

using namespace qpp;
using namespace qpp::types;

//int main(int argc, char **argv)
int main()
{
	_init();

	// Display formatting
	std::cout << std::fixed; // use fixed format for nice formatting
//	std::cout << std::scientific;
	std::cout << std::setprecision(4); // only for fixed or scientific modes

	cout << "Starting qpp..." << endl;

	// TIMING
	Timer t, total;  // start the timer, automatic tic() in the constructor
	size_t N = 10000;
	cmat testmat = cmat::Random(N, N);
	t.toc();
	cout << "It took me " << t.seconds() << " seconds to initialize a " << N
			<< " x " << N << " random complex matrix." << endl;
	t.tic();

	cmat b;
	b = testmat * testmat;
	t.toc(); // read the time
	cout << "It took me " << t.seconds()
			<< " seconds to perform matrix multiplication (including initialization)."
			<< endl;

	total.toc();
	cout << "Total time: " << total.seconds() << " seconds." << endl;
	// END TIMING

	cout << endl << "Exiting qpp..." << endl;
}
