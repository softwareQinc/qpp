/*
 * File:   main.cpp
 * Author: vlad
 * 
 * Created on December 12, 2013, 10:38 PM
 */

#include <iostream>

#include "qpp.h"

using namespace std;

using namespace qpp;
using namespace qpp::types;

//int main(int argc, char **argv)
int main()
{
	cout << "Starting qpp..." << endl << endl;
	_init(); // this will be done automatically in the framework

	cout << endl << endl << "Exiting qpp..." << endl;
}
