/*
 * File:   main.cpp
 * Author: vlad
 * 
 * Created on December 12, 2013, 10:38 PM
 */

#include <iostream>
#include <unsupported/Eigen/MatrixFunctions>

#include "qpp.h"

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
	size_t statistics[] =
	{ 0, 0, 0 };
	size_t N = 100000;
	for (int i = 0; i < N; i++)
		statistics[d.sample()]++;

	cout << (double) statistics[0] / (double) N << endl;
	cout << (double) statistics[1] / (double) N << endl;
	cout << (double) statistics[2] / (double) N << endl;

	cmat a = randn(2, 2).template cast<cplx>()
			+ ct::ii * randn(2, 2).template cast<cplx>();
	cout<<a.sqrt()<<endl;

	stat::DiscreteDistributionFromComplex dc{1.+ct::ii, ct::ii };
	for(auto i:dc._d.probabilities())
		cout<<i<<" ";

	//cout<<sqrtm(a)<<endl;

	cout << endl << endl << "Exiting qpp..." << endl;
}
