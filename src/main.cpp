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

	size_t n = 12; // 12 qubits
	size_t N = std::pow(n, 12);
	std::cout << "n=" << n << " qubits, matrix size " << N << " x " << N << " ."
			<< endl;

	// TIMING
	Timer t, total;  // start the timer, automatic tic() in the constructor
	cmat randcmat = cmat::Random(N, N);
	t.toc(); // read the time
	cout << "Took " << t.seconds() << " seconds.";

	t.tic(); // reset the chronometer
	cmat prodmat;
	prodmat = randcmat * randcmat;
	t.toc(); // read the time
	cout << "Took " << t.seconds() << " seconds.";

	total.toc(); // read the total running time
	cout << "Total time: " << total.seconds() << " seconds." << endl;
	// END TIMING

	cout << endl << "Exiting qpp..." << endl;
}
