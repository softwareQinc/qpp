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

	auto ru = rand_unitary(3);
	displn(ru);
	cout << endl;

	auto ruabs = ru.unaryExpr([]( std::complex<double> z)->double
	{	return std::abs(z);});
	cout << endl;
	displn(ruabs.cast<cplx>());

	auto ruabs1 = fun(ru, myfunc);
	cout << endl << typeid(ruabs1).name() << endl;
	displn(ruabs1);

	cout<<endl;
	displn(abs(ruabs));

	cmat a=cmat::Random(3,3);
	a(1,1)=1;
	cout<<endl;
	displn(a);
	displn(a+a);



	cout << endl << "Exiting qpp..." << endl;
}
