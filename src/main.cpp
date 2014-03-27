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
#include "internal.h"

using namespace std;
using namespace qpp;
using namespace qpp::types;

//int main(int argc, char **argv)
int main()
{
	_init(); // this will be done automatically in the framework

	Eigen::MatrixXd _a = Eigen::MatrixXd::Random(4, 4);
	save(_a, "/Users/vlad/tmp/_a");

	Eigen::MatrixXd a = load<Eigen::MatrixXd>("/Users/vlad/tmp/_a");

	std::vector<size_t> subsys =
	{ 2, 2 };
	std::vector<size_t> perm =
	{ 1, 0 };

	cout << "Error in norm difference load/save: " << norm(_a - a) << endl;

	disp(ptrace2(a,
	{ 2, 2 }));
	cout << endl << endl;
	disp(ptrace(a,
	{ 0 },
	{ 2, 2 }));
	cout << endl << endl;

	imat kt(3, 1);
	kt << 0, 0, 1;

	imat bt(1, 3);
	bt << 0, 1, 0;

	disp(kron(kt, bt).template cast<double>());
	cout << endl << endl;

	disp(kron(bt, kt));
	cout << endl << endl;

	disp(bt * kt);
	cout << endl << endl;

	size_t dim = 8;
	cout << "generating random unitary..." << endl;
	cmat u = rand_unitary(dim);
	cout << endl << "done generating random unitary";
	disp(u);
	cout << endl;
	double normdiff = norm((cmat) (u * adjoint(u) - cmat::Identity(dim, dim)));
	cout << "Norm difference: " << normdiff << endl;
	disp(normdiff, std::cout, 18);
	if (normdiff > std::numeric_limits<double>::epsilon())
		cout << "YES";
	else
		cout << "NO";
	cout << endl << endl;
	cout << "The eigenvalues of u are: " << endl;
	disp(transpose(evals(u)));
	cout << endl << endl;

	cout << "The absolute values of the eigenvalues of u are: " << endl;
	disp(abs(transpose(evals(u))));
	cout << endl << endl;

	// matrix functional calculus
	cmat mat(2, 2);
	mat << 1. + ct::ii, 2, 3, 4;
	disp(mat);
	cout << endl << endl;
	disp(logm(expm(mat)), std::cout);
	cout << endl << endl;
}
