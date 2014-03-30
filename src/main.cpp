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
// TODO: check that everything works for expressions in ALL FILES!!!!
// TODO: Error class
// TODO: change all for(s) to column major order
// TODO: use .data() raw pointer instead of looping
// TODO: use an Engine class (with static members) to get rid of qpp.cpp

using namespace std;

using namespace qpp;
using namespace qpp::types;

//int main(int argc, char **argv)
int main()
{
	_init();

	cout << "Starting qpp..." << endl;

	displn(rand_unitary(3));

	cout << endl << "Exiting qpp..." << endl;
}
