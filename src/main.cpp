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

//int main(int argc, char **argv)
int main()
{
	_init(); // this will be done automatically in the framework

	cmat psi(4, 1);
	psi << 0.6, 0, 0, 0.8; // A Bell-state
	psi = psi / norm(psi);

	cmat rho = ptrace((cmat) (psi * adjoint(psi)),
	{ 1 },
	{ 2, 2 });

	cout << entropy(rho) << endl << endl;

	cmat vec(6, 1);
	vec << 0.36, 0, 0, 0.64, 0, 0;
	//cout<<expm(vec)<<endl;

	cout << entropy(vec) << endl << endl;

	cmat u = rand_unitary(2);
	disp(u);
	cout << endl;
	disp((cmat) (u * adjoint(u) - cmat::Identity(2, 2)));
	cout << endl << endl;

	cout << u * adjoint(u) - cmat::Identity(2, 2);

	u<<ct::ee,0,0,ct::ee;
	cout << endl << endl;
	cout << "Matrix log: " << endl;
	cout << funm(u, std::log);
}
