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

template<typename Derived>
Derived kron_list1(const std::vector<Eigen::MatrixBase<Derived>> &list)
{
	Derived result = list[0];
	std::cout << "First element in list: " << result;
	std::cout << std::endl;
	for (size_t i = 1; i < list.size(); i++)
		result = kron(result, list[i]);
	return result;
}

//int main(int argc, char **argv)
int main()
{
	_init(); // this will be done automatically in the framework

	Eigen::MatrixXd _a = Eigen::MatrixXd::Random(4, 4);
	save(_a, "/Users/vlad/tmp/_a");

	Eigen::MatrixXd a = load<Eigen::MatrixXd>("/Users/vlad/tmp/_a");

	std::vector<size_t> subsys = { 2, 2 };
	std::vector<size_t> perm = { 1, 0 };

	cout << "Error in norm difference load/save: " << norm(_a - a) << endl;

	disp(ptrace2(a, { 2, 2 }));
	cout << endl << endl;
	disp(ptrace(a, { 0 }, { 2, 2 }));
	cout << endl << endl;

}
