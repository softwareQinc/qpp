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

//int main(int argc, char **argv)
int main()
{
	_init(); // this will be done automatically in the framework

//	// test partial trace speed
////	size_t n = atoi(argv[1]); // number of qubits
////	size_t nout = atoi(argv[2]); // discard the first ndiscarded qubits
//
//	size_t n = 10; // number of qubits
//	size_t nout = n - 2; // trace down first n-1 qubits
//
//	size_t dim = pow(2, n); // total dimension
//
//	std::vector<size_t> dims;
//	std::vector<size_t> subsys;
//
//	for (size_t i = 0; i < n; i++)
//		dims.push_back(2);
//	for (size_t i = 0; i < nout; i++)
//		subsys.push_back(i);
//
//	// generate a random 2^n x 2^n matrix
//	cmat A = randn(dim);
//
//	// take the partial trace
//	cmat B = ptrace(A, subsys, dims);
//
//	cout << "Partial trace of " << nout << " out of " << n
//			<< " qubits. Total dimension=2^" << n << "=" << dim << "." << endl;
//	cout << endl << "Remaining matrix:" << endl;
//	if (B.cols() < 32)
//		disp(B);
//	else
//		cout << " Size too large to be useful :)";
//
//	cout << endl << "The eigenvalues of B are:" << endl;
//	disp(evals(B));
//	cout << endl << "The eigenvectors of B are:" << endl;
//	disp(evects(B));
//
//	cout << endl;
//
//	disp(1.0 + ct::ii);
//	cout << endl;
//	disp(1.0 + ct::ii);
//	cout << endl;
//	disp_container(dims);
//	cout << endl;
//	disp_container(dims);
//	cout << endl;
//	cout<< "The eigenvalues of Fd(9) are: "<<endl;
//	disp(evals(gt::Fd(9)));
//	cout<<endl;
//	disp(gt::Fd(3));
//	cout<<endl<<"Precision: "<<ct::eps<<endl;
//	cout<<endl;
//	disp(ct::eps);
//	cout<<endl;

	// ptranspose test
	size_t n = 11; // number of qubits
	size_t dim = pow(2,n); // total dimension

	cmat mat(dim,dim);
	size_t cnt=0;
	for(size_t i=0;i<dim;i++)
		for(size_t j=0;j<dim;j++)
			mat(i,j)=cnt++;

	// dimension vector
	vector<size_t> dm;
	for(size_t i=0; i<n;i++)
		dm.push_back(2);

	// subsystems
	vector<size_t> sbs={0,1};

	cout<<"Matrix: "<< endl;
	//disp(mat);
	cout<<endl<<endl;

	cout<<"Partial transpose"<< endl;
	//disp(ptranspose(mat, sbs, dm));
	ptranspose(mat, sbs, dm);
    cout<<endl<<endl;

	cout<<"Eigenvalues of partial transpose"<<endl;
	//disp(evals(ptranspose(mat, sbs, dm)));
	cout<<endl<<endl;
}
