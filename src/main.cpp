/*
 * File:   main.cpp
 * Author: vlad
 * 
 * Created on December 12, 2013, 10:38 PM
 */

#include <iostream>
#include <ctime>

#include "qpp.h"
//#include "matlab.h" // support for MATLAB

// TODO: expandout function
// TODO: dyad function
// TODO: proj (dya) function
// TODO: ip (inner product function) function, make it general to return matrices
// TODO: Error class
// TODO: change all for(s) to column major order
// TODO: use .data() raw pointer instead of looping
// TODO: write a timer class, tic/toc

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


	/* TIMING */
	clock_t start, end;
	double time;
	start = clock();
	/* process starts here */

	/* process ends here */
	end = clock();
	time = (double) (end - start);
	cout << "It took me " << time << " clicks ("
			<< ((float) time) / CLOCKS_PER_SEC << " seconds).";
	/* END TIMING */

	cout << endl << "Exiting qpp..." << endl;
}
