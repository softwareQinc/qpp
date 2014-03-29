/*
 * File:   main.cpp
 * Author: vlad
 * 
 * Created on December 12, 2013, 10:38 PM
 */

#include <iostream>

#include "qpp.h"
#include "matlab.h" // support for MATLAB

using namespace std;

using namespace qpp;
using namespace qpp::types;

//int main(int argc, char **argv)
int main()
{
	cout << "Starting qpp..." << endl << endl;
	_init(); // this will be done automatically in the framework

	std::vector<double> weights;
	weights.push_back(1);
	weights.push_back(2);
	weights.push_back(3);
	stat::DiscreteDistribution d(weights.begin(), weights.end());
	size_t statistics[] = { 0, 0, 0 };
	size_t N = 100000;
	for (size_t i = 0; i < N; i++)
		statistics[d.sample()]++;

	cout << (double) statistics[0] / (double) N << endl;
	cout << (double) statistics[1] / (double) N << endl;
	cout << (double) statistics[2] / (double) N << endl;

	cmat a = randn(2, 2).template cast<cplx>()
			+ ct::ii * randn(2, 2).template cast<cplx>();

	cout << endl;
	stat::DiscreteDistributionFromComplex dc { 1. + ct::ii, ct::ii };
	for (auto i : dc._d.probabilities())
		cout << i << " ";

#define MTYPE cmat

//		MTYPE m1 = MTYPE::Random(30, 40);
//
//		cout << endl << endl;
//		disp(m1);
//		cout << endl << endl;
//		saveMATLAB(m1, "/Users/vlad/tmp/test.mat", "m1","w");
//
//		// we load it back
//		MTYPE m2 = loadMATLAB<MTYPE>("/Users/vlad/tmp/test.mat", "m1");
//		disp(m2);
//		cout << endl << endl;
//
//		// and save it again
//		saveMATLAB(m2, "/Users/vlad/tmp/test.mat", "m2","u");
//
//		// now check the norm difference
//		cout << "Norm difference ||m_1-m_2||=" << norm(m1 - m2) << endl;
	MTYPE mtest = loadMATLAB<MTYPE>("/Users/vlad/tmp/test.mat", "mtest");

	cout << endl << endl;
	cout << mtest.rows() << " x " << mtest.cols() << " loaded matrix: " << endl
			<< endl;
	displn(mtest);

	save(mtest, "/Users/vlad/tmp/qpp.out");

	MTYPE mtest1 = load<MTYPE>("/Users/vlad/tmp/qpp.out");
	cout << "My load function: " << mtest1.rows() << " x " << mtest1.cols()
			<< " matrix:" << endl;
	displn(mtest1);

	cout<<endl<<"Norm difference: "<<norm(mtest1-mtest)<<endl;

	cout << "Ending..." << endl;

}
