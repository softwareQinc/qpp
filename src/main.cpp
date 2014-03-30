/*
 * File:   main.cpp
 * Author: vlad
 * 
 * Created on December 12, 2013, 10:38 PM
 */

#include <iostream>

#include "qpp.h"
#include "matlab.h" // support for MATLAB

// TODO: expandout function
// TODO: dyad function
// TODO: proj (dya) function
// TODO: ip (inner product function) function, make it general to return matrices
// TODO: check that everything works for expressions in ALL FILES!!!!
// TODO: Error class
// TODO: change all for(s) to column major order
// TODO: use .data() raw pointer instead of looping
// TODO: use a Singleton Engine class (with static members) to get rid of qpp.cpp
// TODO: look at unaryExpr for functors!!!!

using namespace std;

using namespace qpp;
using namespace qpp::types;

float myfunc(const cplx &z)
{
	return std::abs(z);
}

//int main(int argc, char **argv)
int main()
{
	_init();

	cout << "Starting qpp..." << endl;

	auto randu = rand_unitary(100);
	saveMATLAB<cplx>(randu, "/Users/vlad/tmp/test.mat",
			"randu", "w");
	auto res = loadMATLAB<cplx>("/Users/vlad/tmp/test.mat", "randu");
	cout << endl;
	displn(norm(randu-res));

	cout << endl << "Exiting qpp..." << endl;
}
