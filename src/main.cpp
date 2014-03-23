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

int main(int argc, char **argv)
{
	_init(); // this will be done automatically in the framework

	// test partial trace speed
//	size_t n = atoi(argv[1]); // number of qubits
//	size_t nout = atoi(argv[2]); // discard the first ndiscarded qubits

	size_t n = 3;
	size_t nout = n-1; // trace down first n-1 qubits

	size_t dim = pow(2, n); // total dimension

	ivect dims(n);
	for (size_t i = 0; i < n; i++)
		dims(i)=2;

	ivect subsys(n);
	for (size_t i = 0; i < nout; i++)
		subsys(i)=i;

	ivect perm(n);
	for(size_t i = 1; i<n; i++)
		perm(i-1)=i;
	perm(n-1)=0;

	// generate a random 2^n x 2^n matrix
	cmat A = randn(dim);

//	// take the partial trace
//	cmat B = ptrace(A, subsys, dims);
//
//	cout << "Partial trace of " << nout << " out of " << n
//			<< " qubits. Total dimension=2^" << n << "=" << dim << "." << endl;
//	cout << endl << "Remaining matrix:"<<endl;
//	if(B.cols()<32)
//		disp(B);
//	else
//		cout<<" Size too large to be useful :)";

	cvect v(3);
	v<<ct::ii,2,3;
	disp(v);
	cout<<endl<<endl;
	cout<<v<<endl<<endl;

	ivect a(2);
	a<<1,2;

	cmat c=randn(3);
	disp(mat_pow(c,1));
	cout<<endl<<endl;
	cmat result=transpose(v)*v;
	cplx r=result(0);
	disp(r);
	//cplx vv=trace(c);

}
