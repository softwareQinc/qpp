/*
 * matlab.h
 *
 *  Created on: Mar 28, 2014
 *      Author: vlad
 */

#ifndef MATLAB_H_
#define MATLAB_H_

// add the path to MATLAB/extern/include in include path

#include <Eigen/Dense>
#include <string>
#include <stdexcept>
#include "mat.h"  // path defined in the Makefile
#include "mex.h"  // path defined in the Makefile

// TODO: finish matlab.h
namespace qpp
{

// define the functions here
template<typename Derived>
void saveMATLAB(const Eigen::MatrixBase<Derived> &A, const char* mat_file,
		const char* var_name)
{
	// cast the input to a double (internal MATLAB format)
	Eigen::MatrixXd tmp = A.template cast<double>();

	MATFile *pmat;
	mxArray *pa;
	pmat = matOpen(mat_file, "w");
	if (pmat == NULL)
	{
		throw std::runtime_error(
				"saveMATLAB: can not create MATLAB output file!");
	}
	pa = mxCreateDoubleMatrix(tmp.rows(), tmp.cols(), mxREAL);
	memcpy((void*)(mxGetPr(pa)), (void*)tmp.data(), sizeof(double)*tmp.size() );
	if(matPutVariable(pmat,var_name,pa))
		throw std::runtime_error("saveMATLAB: can not save the variable into MATLAB file!");

	matClose(pmat);

}

inline void loadMATLAB(Eigen::MatrixXd &A, const std::string &mat_file,
		const std::string & var_name)
{

}

}
#endif /* MATLAB_H_ */
