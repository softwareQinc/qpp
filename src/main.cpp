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
	_init();

	cout << "Starting qpp..." << endl;

	cmat m(2,2);
	m<<1,2,3,4;

	displn(mpower(m,ct::ii));
	cout<<endl;
	displn(mpower_n(m,2.1));



	cout << endl << "Exiting qpp..." << endl;

}
