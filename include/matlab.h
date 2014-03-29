/*
 * matlab.h
 *
 *  Created on: Mar 28, 2014
 *      Author: vlad
 */

#ifndef MATLAB_H_
#define MATLAB_H_

#define HAVE_MATLAB_INSTALLED // comment if you do not have MATLAB installed

#ifdef HAVE_MATLAB_INSTALLED
#include <Eigen/Dense>
#include <string>
#include "mat.h"
#include "mex.h"

// define the functions here
template<typename Derived>
void saveMATLAB(const Eigen::MatrixBase<Derived> &A,
		const std::string &mat_file, const std::string & var_name)
{
	// cast the input to a double (internal MATLAB format)
	Eigen::MatrixXd tmp = A.template cast<double>();

}

inline void loadMATLAB(Eigen::MatrixXd &A, const std::string &mat_file,
		const std::string & var_name)
{

}

#endif /* HAVE_MATLAB_INSTALLED */

#endif /* MATLAB_H_ */
