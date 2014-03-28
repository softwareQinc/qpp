/*
 * File:   main.cpp
 * Author: vlad
 * 
 * Created on December 12, 2013, 10:38 PM
 */

#include <iostream>
#include <fstream>
#include <string>
#include <stdexcept>
#include <cmath>

#include "qpp.h"

using namespace std;

using namespace qpp;
using namespace qpp::types;

double f(const cplx &z)
{
	return std::abs(z);
}

//int main(int argc, char **argv)
int main()
{
	_init(); // this will be done automatically in the framework

	cmat a(2, 2);
	a << 1. + ct::ii, 2. - ct::ii, 3. + 2.2*ct::ii, 4.1;
	disp(a);
	cout << endl << endl;

	// Automatic type deduction for all 3 template parameters
	// FunctionInputType, FunctionOutputType, MatrixInputType
	cout << fun(a, f) << endl << endl;

	// Need to specify FunctionInputType and FunctionOutputType as std::exp is overloaded
	// and has templated output type
	disp(fun<cplx, cplx>(a, std::exp));
	cout << endl << endl;

	// Lambda invocation, cannot deduce FunctionInputType and FunctionOutputType
	cout << fun<cplx, double>(a, [](const cplx & z)->double
	{	return imag(z);}) << endl << endl;


	// Truncating the result using an MatrixXi as output
	imat tmp = fun<cplx, int>(a, [](const cplx & z)->int
	{	return imag(z);});
	cout << tmp << endl << endl;

}

