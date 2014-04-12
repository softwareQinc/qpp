/*
 * matlab.h
 *
 *  Created on: Mar 28, 2014
 *      Author: vlad
 */

#ifndef MATLAB_H_
#define MATLAB_H_

// MATLAB I/O interfacing
// add the path to MATLAB/extern/include in include path

#include <Eigen/Dense>

#include <stdexcept>
#include <string>

#include "internal.h"
#include "types.h"
#include "classes/exception.h"

#include "mat.h"  // path to this is defined in the Makefile
#include "mex.h"  // path to this is defined in the Makefile

namespace qpp
{

// load Eigen::MatrixX from MATLAB .mat file
template<typename Derived>
Derived loadMATLABmatrix(const std::string &mat_file,
		const std::string & var_name)
{
	throw Exception("loadMATLABmatrix", Exception::Type::UNDEFINED_TYPE);
}

// MatrixXd specialization
// if var_name is a complex matrix, only the real part is loaded
template<>
types::dmat loadMATLABmatrix(const std::string &mat_file,
		const std::string & var_name)
{
	MATFile* pmat = matOpen(mat_file.c_str(), "r");
	if (pmat == NULL)
	{
		throw std::runtime_error(
				"loadMATLABmatrix: Can not open MATLAB file " + mat_file + "!");
	}

	mxArray* pa = matGetVariable(pmat, var_name.c_str());
	if (pa == NULL)
		throw std::runtime_error(
				"loadMATLABmatrix: Can not load the variable " + var_name
						+ " from MATLAB file " + mat_file + "!");

	if (mxGetNumberOfDimensions(pa) != 2) // not a matrix
		throw std::runtime_error(
				"loadMATLABmatrix: Loaded variable " + var_name
						+ " is not 2-dimensional!");

	if (!mxIsDouble(pa))
		throw std::runtime_error(
				"loadMATLABmatrix: Loaded variable " + var_name
						+ " is not in double-precision format!");

	size_t rows = mxGetM(pa);
	size_t cols = mxGetN(pa);

	types::dmat result(rows, cols);

	std::memcpy(result.data(), mxGetPr(pa),
			sizeof(double) * mxGetNumberOfElements(pa));

	mxDestroyArray(pa);
	matClose(pmat);

	return result;
}

// cmat specialization
template<>
types::cmat loadMATLABmatrix(const std::string &mat_file,
		const std::string & var_name)
{
	MATFile* pmat = matOpen(mat_file.c_str(), "r");
	if (pmat == NULL)
	{
		throw std::runtime_error(
				"loadMATLABmatrix: Can not open MATLAB file " + mat_file + "!");
	}

	mxArray* pa = matGetVariable(pmat, var_name.c_str());
	if (pa == NULL)
		throw std::runtime_error(
				"loadMATLABmatrix: Can not load the variable " + var_name
						+ " from MATLAB file " + mat_file + "!");

	if (mxGetNumberOfDimensions(pa) != 2) // not a matrix
		throw std::runtime_error(
				"loadMATLABmatrix: Loaded variable " + var_name
						+ " is not 2-dimensional!");

	if (!mxIsDouble(pa))
		throw std::runtime_error(
				"loadMATLABmatrix: Loaded variable " + var_name
						+ " is not in double-precision format!");

	size_t rows = mxGetM(pa);
	size_t cols = mxGetN(pa);

	types::dmat result_re(rows, cols);
	types::dmat result_im(rows, cols);

	// real part and imaginary part pointers
	double* pa_re = nullptr, *pa_im = nullptr;

	// Populate the real part of the created array.
	pa_re = (double*) mxGetPr(pa);
	std::memcpy(result_re.data(), pa_re,
			sizeof(double) * mxGetNumberOfElements(pa));

	if (mxIsComplex(pa)) // populate the imaginary part if exists
	{
		pa_im = (double*) mxGetPi(pa);
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
template<typename Derived>
void saveMATLABmatrix(const Eigen::MatrixBase<Derived> &A,
		const std::string & mat_file, const std::string & var_name,
		const std::string & mode)
{
	throw Exception("saveMATLABmatrix", Exception::Type::UNDEFINED_TYPE);
}

template<> // Eigen::MatrixXd specialization
void saveMATLABmatrix(const Eigen::MatrixBase<typename types::dmat> &A,
		const std::string & mat_file, const std::string & var_name,
		const std::string & mode)
{
	const types::dmat & rA = A;

	// check zero-size
	if (!internal::_check_nonzero_size(rA))
		throw Exception("saveMATLABmatrix", Exception::Type::ZERO_SIZE);

	MATFile* pmat = matOpen(mat_file.c_str(), mode.c_str());
	if (pmat == NULL)
		throw std::runtime_error(
				"saveMATLABmatrix: Can not open/create MATLAB file " + mat_file
						+ "!");

	mxArray* pa = mxCreateDoubleMatrix(rA.rows(), rA.cols(), mxREAL);
	if (pa == NULL)
		throw std::runtime_error(
				"saveMATLABmatrix: mxCreateDoubleMatrix failed!");

	std::memcpy(mxGetPr(pa), rA.data(), sizeof(double) * rA.size());

	if (matPutVariable(pmat, var_name.c_str(), pa))
		throw std::runtime_error(
				"saveMATLABmatrix: Can not write the variable " + var_name
						+ " to MATLAB file " + mat_file + "!");

	mxDestroyArray(pa);
	matClose(pmat);
}

template<> // Eigen::MatrixXcd specialization
void saveMATLABmatrix(const Eigen::MatrixBase<typename types::cmat> &A,
		const std::string & mat_file, const std::string & var_name,
		const std::string & mode)
{
	const types::cmat & rA = A;

	// check zero-size
	if (!internal::_check_nonzero_size(rA))
		throw Exception("saveMATLABmatrix", Exception::Type::ZERO_SIZE);

	// cast the input to a double (internal MATLAB format)
	types::dmat tmp_re = rA.real();
	types::dmat tmp_im = rA.imag();

	MATFile* pmat = matOpen(mat_file.c_str(), mode.c_str());
	if (pmat == NULL)
		throw std::runtime_error(
				"saveMATLABmatrix: Can not open/create MATLAB file " + mat_file
						+ "!");

	mxArray* pa = mxCreateDoubleMatrix(tmp_re.rows(), tmp_re.cols(), mxCOMPLEX);
	if (pa == NULL)
		throw std::runtime_error(
				"saveMATLABmatrix: mxCreateDoubleMatrix failed!");

	double* pa_re, *pa_im;

	/* Populate the real part of the created array. */
	pa_re = (double*) mxGetPr(pa);
	std::memcpy(pa_re, tmp_re.data(), sizeof(double) * tmp_re.size());

	/* Populate the imaginary part of the created array. */
	pa_im = (double*) mxGetPi(pa);
	std::memcpy(pa_im, tmp_im.data(), sizeof(double) * tmp_im.size());

	if (matPutVariable(pmat, var_name.c_str(), pa))
		throw std::runtime_error(
				"saveMATLABmatrix: Can not write the variable " + var_name
						+ " to MATLAB file " + mat_file + "!");

	mxDestroyArray(pa);
	matClose(pmat);
}

} /* namespace qpp */
#endif /* MATLAB_H_ */
