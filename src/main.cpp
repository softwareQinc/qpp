/*
 * File:   main.cpp
 * Author: vlad
 * 
 * Created on December 12, 2013, 10:38 PM
 */

#include <iostream>
#include <cmath>

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
// TODO: test that everything works with GenProducts!
// TODO: implement selfadjoint eigensolver

using namespace std;

using namespace qpp;
using namespace qpp::types;

int myfunc(const cplx &z)
{
	return std::abs(z);
}

//int main(int argc, char **argv)
int main()
{
	_init();

	std::cout << std::fixed; // use fixed format for nice formatting
//	std::cout << std::scientific;
	std::cout << std::setprecision(4); // only for fixed or scientific modes

	cout << "Starting qpp..." << endl;

	// MATLAB interface testing
	cmat randu = rand_unitary(3);
	cout << endl;
	displn(randu);

	saveMATLAB<cplx>(randu, "/Users/vlad/tmp/test.mat", "randu", "w");
	cmat res = loadMATLAB<cplx>("/Users/vlad/tmp/test.mat", "randu");
	cout << endl;
	displn(randu - res);

	// functor testing
	auto lambda = [](const cplx& z)->int
	{	return abs(z);};

	cmat mat1(3, 3);
	mat1 << 1, -2.56, 335.2321, -4, 5.244, -6.1, 7, -8, 9. + 2.78 * ct::ii;

	cout << endl;
	displn(mat1);

	cout << endl;
	displn(fun<cmat, int>(mat1 * mat1, lambda));

	cout << endl;
	displn(fun<cmat>(mat1 * mat1, myfunc));

	// other functions
	cout << endl;
	cmat mat11=mat1*mat1;
	displn(kron(mat11, mat11));

	//mat11=randu;
	// eigenvalue test
	cout << endl << evals(mat11)(0) << endl;
	displn(mat11*evects(mat1*mat1).col(0)-(evals(mat11)(0))*evects(mat11).col(0));




	cout << endl << "Exiting qpp..." << endl;
}
