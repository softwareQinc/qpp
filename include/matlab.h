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

// save Eigen::MatrixX to MATLAB .mat file
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

// complex specialization (Eigen::MatrixXcd)
template<>
void saveMATLAB(const Eigen::MatrixBase<Eigen::MatrixXcd> &A,
		const std::string & mat_file, const std::string & var_name)
{
	// cast the input to a double (internal MATLAB format)
	Eigen::MatrixXd tmp_re = A.real().template cast<double>();
	Eigen::MatrixXd tmp_im = A.imag().template cast<double>();

	MATFile *pmat = matOpen(mat_file.c_str(), "w");
	if (pmat == NULL)
	{
		throw std::runtime_error(
				"saveMATLAB: can not create MATLAB output file!");
	}

	mxArray *pa = mxCreateDoubleMatrix(tmp_re.rows(), tmp_re.cols(), mxCOMPLEX);
	memcpy((void*) (mxGetPr(pa)), (void*) tmp_re.data(),
			sizeof(double) * tmp_re.size());
	memcpy((void*) (mxGetPr(pa) + mxGetNumberOfElements(pa) + 2),
			(void*) tmp_im.data(), sizeof(double) * tmp_im.size());

	if (matPutVariable(pmat, var_name.c_str(), pa))
		throw std::runtime_error(
				"saveMATLAB: can not write to MATLAB output file!");

	matClose(pmat);

}

// load Eigen::MatrixX from MATLAB .mat file
template<typename Derived>
Derived loadMATLAB(const std::string &mat_file, const std::string & var_name)
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

	size_t rows = mxGetM(pa);
	size_t cols = mxGetN(pa);

	Eigen::MatrixXd result(rows, cols);

	memcpy((void*) result.data(), (void*) (mxGetPr(pa)),
			sizeof(double) * mxGetNumberOfElements(pa));
	matClose(pmat);

	return result;
}

// complex specialization (Eigen::MatrixXcd)
template<>
Eigen::MatrixXcd loadMATLAB(const std::string &mat_file,
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

	size_t rows = mxGetM(pa);
	size_t cols = mxGetN(pa);

	Eigen::MatrixXd result_re(rows, cols);
	Eigen::MatrixXd result_im(rows, cols);

	memcpy((void*) result_re.data(), (void*) (mxGetPr(pa)),
			sizeof(double) * mxGetNumberOfElements(pa));
	memcpy((void*) result_im.data(),
			(void*) (mxGetPr(pa) + mxGetNumberOfElements(pa) + 2),
			sizeof(double) * mxGetNumberOfElements(pa));
	matClose(pmat);

	return (result_re.template cast<types::cplx>())
			+ ct::ii * (result_im.template cast<types::cplx>());
}

}
#endif /* MATLAB_H_ */
