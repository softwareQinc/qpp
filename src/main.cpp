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
	Timer t,total;  // start the timer
	size_t N = 10000;
	cmat testmat = cmat::Random(N, N);
	t.toc();
	cout << "It took me " << t.ticks() << " ticks (" << t.secs()
				<< " seconds) to initialize a " << N << " x " << N
				<< " random complex matrix."<<endl;
	t.reset();
	cout << "The norm of a " << N << " x " << N
			<< " random complex matrix is: "<<endl;
	cout << norm(testmat) << endl;
	t.toc(); // read the time
	cout << "It took me " << t.ticks() << " ticks (" << t.secs()
			<< " seconds) to compute the norm." << endl;

	cout << "Total time: " << total.ticks() << " ticks (" << total.secs()
				<< " seconds)." << endl;
	// END TIMING

	cout << endl << "Exiting qpp..." << endl;
}
