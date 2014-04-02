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

namespace qpp
{

// load Eigen::MatrixX from MATLAB .mat file
template<typename Scalar>
types::DynMat<Scalar> loadMATLABmatrix(const std::string &mat_file,
		const std::string & var_name)
{
	throw std::logic_error("loadMATLABmatrix: not implemented for this type!");
}

// double specialization
// if var_name is a complex matrix, only the real part is loaded
template<>
inline types::DynMat<double> loadMATLABmatrix(const std::string &mat_file,
		const std::string & var_name)
{
	MATFile *pmat = matOpen(mat_file.c_str(), "r");
	if (pmat == NULL)
	{
		throw std::runtime_error(
				"loadMATLABmatrix: can not open MATLAB input file!");
	}

	mxArray *pa = matGetVariable(pmat, var_name.c_str());
	if (pa == NULL)
		throw std::runtime_error("loadMATLABmatrix: can not load the variable "
				"from MATLAB input file!");

	if (mxGetNumberOfDimensions(pa) != 2) // not a matrix
		throw std::runtime_error(
				"loadMATLABmatrix: loaded variable is not 2-dimensional!");

	if (!mxIsDouble(pa))
		throw std::runtime_error("loadMATLABmatrix: loaded variable is not "
				"in double-precision format!");

	size_t rows = mxGetM(pa);
	size_t cols = mxGetN(pa);

	types::dmat result(rows, cols);

	std::memcpy(result.data(), mxGetPr(pa),
			sizeof(double) * mxGetNumberOfElements(pa));

	mxDestroyArray(pa);
	matClose(pmat);

	return result;
}

// complex specialization
template<>
inline types::DynMat<types::cplx> loadMATLABmatrix(const std::string &mat_file,
		const std::string & var_name)
{
	MATFile *pmat = matOpen(mat_file.c_str(), "r");
	if (pmat == NULL)
	{
		throw std::runtime_error(
				"loadMATLABmatrix: can not open MATLAB input file!");
	}

	mxArray *pa = matGetVariable(pmat, var_name.c_str());
	if (pa == NULL)
		throw std::runtime_error(
				"loadMATLABmatrix: can not load variable from MATLAB input file!");

	if (mxGetNumberOfDimensions(pa) != 2) // not a matrix
		throw std::runtime_error(
				"loadMATLABmatrix: loaded variable is not 2-dimensional!");

	if (!mxIsDouble(pa))
		throw std::runtime_error("loadMATLABmatrix: loaded variable is not "
				"in double-precision format!");

	size_t rows = mxGetM(pa);
	size_t cols = mxGetN(pa);

	types::dmat result_re(rows, cols);
	types::dmat result_im(rows, cols);

	// real part and imaginary part pointers
	double *pa_re = nullptr, *pa_im = nullptr;

	// Populate the real part of the created array.
	pa_re = (double *) mxGetPr(pa);
	std::memcpy(result_re.data(), pa_re,
			sizeof(double) * mxGetNumberOfElements(pa));

	if (mxIsComplex(pa)) // populate the imaginary part if exists
	{
		pa_im = (double *) mxGetPi(pa);
		std::memcpy(result_im.data(), pa_im,
				sizeof(double) * mxGetNumberOfElements(pa));
	}
	else // set to zero the imaginary part
	{
		std::memset(result_im.data(), 0,
				sizeof(double) * mxGetNumberOfElements(pa));
	}

	mxDestroyArray(pa);
	matClose(pmat);

	return (result_re.cast<types::cplx>())
			+ ct::ii * (result_im.cast<types::cplx>());
}

// save Eigen::MatrixX to MATLAB .mat file as a double matrix
// see MATLAB's matOpen(...) documentation
template<typename Scalar>
void saveMATLABmatrix(const types::DynMat<Scalar> &A,
		const std::string & mat_file, const std::string & var_name,
		const std::string & mode)
{
	throw std::logic_error("saveMATLABmatrix: not implemented for this type!");
}

template<> // double specialization
void saveMATLABmatrix(const types::DynMat<double> &A,
		const std::string & mat_file, const std::string & var_name,
		const std::string & mode)
{
	MATFile *pmat = matOpen(mat_file.c_str(), mode.c_str());
	if (pmat == NULL)
	{
		throw std::runtime_error(
				"saveMATLABmatrix: can not open/create MATLAB output file!");
	}

	mxArray *pa = mxCreateDoubleMatrix(A.rows(), A.cols(), mxREAL);
	if (pa == NULL)
		throw std::runtime_error(
				"saveMATLABmatrix: mxCreateDoubleMatrix failed!");

	std::memcpy(mxGetPr(pa), A.data(), sizeof(double) * A.size());

	if (matPutVariable(pmat, var_name.c_str(), pa))
		throw std::runtime_error(
				"saveMATLABmatrix: can not write variable to MATLAB output file!");

	mxDestroyArray(pa);
	matClose(pmat);

}

template<> // complex specialization
void saveMATLABmatrix(const types::DynMat<types::cplx> &A,
		const std::string & mat_file, const std::string & var_name,
		const std::string & mode)
{
	// cast the input to a double (internal MATLAB format)
	types::dmat tmp_re = A.real();
	types::dmat tmp_im = A.imag();

	MATFile *pmat = matOpen(mat_file.c_str(), mode.c_str());
	if (pmat == NULL)
	{
		throw std::runtime_error(
				"saveMATLABmatrix: can not open/create MATLAB output file!");
	}

	mxArray *pa = mxCreateDoubleMatrix(tmp_re.rows(), tmp_re.cols(), mxCOMPLEX);
	if (pa == NULL)
		throw std::runtime_error(
				"saveMATLABmatrix: mxCreateDoubleMatrix failed!");

	double *pa_re, *pa_im;

	/* Populate the real part of the created array. */
	pa_re = (double *) mxGetPr(pa);
	std::memcpy(pa_re, tmp_re.data(), sizeof(double) * tmp_re.size());

	/* Populate the imaginary part of the created array. */
	pa_im = (double *) mxGetPi(pa);
	std::memcpy(pa_im, tmp_im.data(), sizeof(double) * tmp_im.size());

	if (matPutVariable(pmat, var_name.c_str(), pa))
		throw std::runtime_error(
				"saveMATLABmatrix: can not write variable to MATLAB output file!");

	mxDestroyArray(pa);
	matClose(pmat);

}

}
#endif /* MATLAB_H_ */
