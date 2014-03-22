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

#include "qpp.h"

using namespace std;
using namespace qpp;
using namespace qpp::types;

int main()
{
	_init(); // this will be done automatically in the framework

	vector<size_t> dims={2,3,4,5,6,7}; // prod(dims) = 5040
	vector<size_t> subsys={1,2,4};
	size_t dim=5040;

	size_t cnt=0;
	cmat A(dim,dim);
	for(auto i=0;i<dim;i++)
		for(auto j=0; j<dim; j++)
			A(i,j)=cnt++;
	 disp(ptrace(A, dims, subsys));
	 cout<<endl<<norm(ptrace(A, dims, subsys))<<endl;

	// do the same thing with trAB, without syspermute (i.e, traceout the last subsystems)
	/*
	 size_t cdimsAB[]={6,120};
	 cmat result = trace(A2)*trace(A3)*trace(A4)*kron(A0,A1);
	 std::vector<size_t> dimsAB(cdimsAB,cdimsAB+sizeof(cdimsAB)/sizeof(cdimsAB[0]));
	 cout<<norm(ptrace2(kronMat,dimsAB)-result);
	 */

	// TODO: optimize syspermute, that's where time is spent!
// C++11 code
	/*
	 #define CPP11
	 #ifdef CPP11
	 vector<cmat> vmat={gt::X,gt::Y,gt::Z};
	 cmat kronvmat=kron_list(vmat);
	 cout<<endl;
	 disp(kronvmat);
	 #endif
	 */

/*
	stat::UniformRealDistribution a(-2, 2);
	cout << endl << a.sample() << endl;

	stat::NormalDistribution b;
	cout << endl << b.sample() << endl;

	cmat A = cmat::Random(2, 2);
	cout << endl;
	disp(A);

	cout << endl;
	int n = 41;
	cout << "The " << n << " root of unity is: " << ct::omega(n) << endl;
*/
}
