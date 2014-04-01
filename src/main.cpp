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
	Timer t; // timer starts
	/* process starts here */

	/* process ends here */
	cout << "It took me " << t.ticks() << " ticks ("
			<< t.secs() << " seconds)."<<endl;

	cout << "It took me " << t.ticks() << " ticks ("
				<< t.secs() << " seconds)."<<endl;

	cout << "It took me " << t.ticks() << " ticks ("
				<< t.secs() << " seconds)."<<endl;
	t.reset();
	cout << "It took me " << t.ticks() << " ticks ("
				<< t.secs() << " seconds)."<<endl;
	cout << "It took me " << t.ticks() << " ticks ("
				<< t.secs() << " seconds)."<<endl;
	/* END TIMING */

	cout << endl << "Exiting qpp..." << endl;
}
