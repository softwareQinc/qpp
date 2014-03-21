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

	/*
	 // test random matrices
	 cmat A=cmat::Random(2,2);
	 cmat B=cmat::Random(3,3);
	 cmat C=cmat::Random(4,4);
	 cmat D=cmat::Random(5,5);
	 cmat kronMat=kron(kron(kron(A,B),C),D);

	 size_t cdims[]={2,3,4,5};
	 std::vector<size_t> dims(cdims,cdims+sizeof(cdims)/sizeof(cdims[0]));

	 size_t cperm[]={3,2,1,0}; // permutation
	 std::vector<size_t> perm(cperm,cperm+sizeof(cperm)/sizeof(cperm[0]));
	 cmat permutedMat=kron(kron(kron(D,C),B),A);

	 disp(permutedMat);
	 cout<<"Norm of the difference: "<<norm(permutedMat-syspermute(kronMat,dims,perm))<<endl;

	 size_t cdims1[]={2,2};
	 std::vector<size_t> dims1(cdims1,cdims1+sizeof(cdims1)/sizeof(cdims1[0]));

	 size_t cperm1[]={1,0}; // permutation
	 std::vector<size_t> perm1(cperm1,cperm1+sizeof(cperm1)/sizeof(cperm1[0]));


	 cmat CNOT12 = gt::CNOT;
	 // CNOT12 -> CNOT 21
	 cmat CNOT21 = syspermute(CNOT12, dims1, perm1);

	 cout<<"CNOT12:"<<endl;
	 disp(CNOT12);
	 cout<<endl<<"CNOT21"<<endl;
	 disp(CNOT21);
	 */
	/*	cmat A0=cmat::Random(3,3);
	 cmat A1=cmat::Random(3,3);
	 cmat A2=cmat::Random(4,4);
	 cmat A3=cmat::Random(5,5);
	 cmat A4=cmat::Random(6,6);


	 vector<cmat> vectMat;
	 vectMat.push_back(A0);
	 vectMat.push_back(A1);
	 vectMat.push_back(A2);
	 vectMat.push_back(A3);
	 vectMat.push_back(A4);


	 cmat kronMat=kron_list(vectMat);
	 //disp(kronMat);
	 cout<<endl;

	 size_t csubsys[]={1,4,2,3}; // trace out A2 A3 A4
	 size_t cdims[]={3,3,4,5,6};

	 std::vector<size_t> dims(cdims,cdims+sizeof(cdims)/sizeof(cdims[0]));
	 std::vector<size_t> subsys(csubsys,csubsys+sizeof(csubsys)/sizeof(csubsys[0]));

	 cmat result = trace(A1)*trace(A2)*trace(A3)*trace(A4)*A0;
	 cout<<norm(ptrace(kronMat, dims, subsys)-result);
	 cout<<endl;
	 cout<<"Size of the remaining systems:"<<result.cols()<<endl;
	 disp(result);
	 */

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
}
