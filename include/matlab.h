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
#include "types.h"

#include "mat.h"  // path to this is defined in the Makefile
#include "mex.h"  // path to this is defined in the Makefile

// TODO: finish matlab.h
namespace qpp
{

// define the functions here
template<typename Derived>
void saveMATLAB(const Eigen::MatrixBase<Derived> &A,
		const std::string & mat_file, const std::string & var_name)
{
	// cast the input to a double (internal MATLAB format)
	Eigen::MatrixXd tmp = A.template cast<double>();

	MATFile *pmat = matOpen(mat_file.c_str(), "w");
	if (pmat == NULL)
	{
		throw std::runtime_error(
				"saveMATLAB: can not create MATLAB output file!");
	}

	mxArray *pa = mxCreateDoubleMatrix(tmp.rows(), tmp.cols(), mxREAL);
	memcpy((void*) (mxGetPr(pa)), (void*) tmp.data(),
			sizeof(double) * tmp.size());

	if (matPutVariable(pmat, var_name.c_str(), pa))
		throw std::runtime_error(
				"saveMATLAB: can not write to MATLAB output file!");

	matClose(pmat);

}

inline Eigen::MatrixXd loadMATLAB(const std::string &mat_file,
		const std::string & var_name)
{
	MATFile *pmat = matOpen(mat_file.c_str(), "r");
	if (pmat == NULL)
	{
		throw std::runtime_error("loadMATLAB: can not open MATLAB input file!");
	}

	mxArray *pa = matGetVariable(pmat, var_name.c_str());
	if (pa == NULL)
		throw std::runtime_error(
				"loadMATLAB: can not load the item from MATLAB input file!");

	size_t rows=mxGetM(pa);
	size_t cols=mxGetN(pa);

	Eigen::MatrixXd result(rows,cols);

	memcpy( (void*) result.data(), (void*) (mxGetPr(pa)),
			sizeof(double) * mxGetNumberOfElements(pa));
	matClose(pmat);

	return result;
}

}
#endif /* MATLAB_H_ */
