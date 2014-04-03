/*
 * entropy.h
 *
 *  Created on: Mar 27, 2014
 *      Author: vlad
 */

#ifndef ENTROPY_H_
#define ENTROPY_H_

#include <cmath>
#include "types.h"
#include "util.h"
#include "internal.h"
#include "exception.h"

// entropy functions

namespace qpp
{

// Shannon/von-Neumann entropy
template<typename Scalar>
double shannon(const types::DynMat<Scalar> & A) throw (Exception)
{
	// check zero-size
	if (!internal::_check_nonzero_size(A))
		throw Exception("shannon", Exception::Type::MATRIX_ZERO_SIZE);

	// input is a vector
	if (A.rows() == 1 || A.cols() == 1)
	{
		double result = 0;
		// take the absolut values of the entries
		// to get rid of unwanted imaginary parts
		for (size_t i = 0; i < static_cast<size_t>(A.size()); i++)
			if (std::abs(A(i)) != 0)
				result -= std::abs(A(i)) * std::log2(std::abs(A(i)));

		return result;
	}

	// input is a matrix

	// check square matrix
	if (!internal::_check_square_mat(A))
		throw Exception("shannon", Exception::Type::MATRIX_NOT_SQUARE);

	// get the eigenvalues
	types::cmat ev = evals(A);
	double result = 0;
	// take the absolut values of the entries
	// to get rid of unwanted imaginary parts
	for (size_t i = 0; i < static_cast<size_t>(ev.rows()); i++)
		if (std::abs((types::cplx) ev(i)) != 0)
			result -= std::abs((types::cplx) ev(i))
					* std::log2(std::abs((types::cplx) ev(i)));

	return result;
}

// Renyi-alpha entropy (alpha>=0)
template<typename Scalar>
double renyi(const double alpha, const types::DynMat<Scalar> & A)
		throw (Exception)
{
	if (alpha < 0)
		throw Exception("renyi", Exception::Type::OUT_OF_RANGE);

	if (alpha == 1) // Shannon/von Neumann
		return shannon(A);

	// check zero-size
	if (!internal::_check_nonzero_size(A))
		throw Exception("renyi", Exception::Type::MATRIX_ZERO_SIZE);

	// input is a vector
	if (A.rows() == 1 || A.cols() == 1)
	{
		if (alpha == 0)
			return std::log2(A.size());

		double result = 0;
		// take the absolut values of the entries
		// to get rid of unwanted imaginary parts
		for (size_t i = 0; i < static_cast<size_t>(A.size()); i++)
			if (std::abs((types::cplx) A(i)) != 0)
				result += std::pow(std::abs(A(i)), alpha);

		return std::log2(result) / (1 - alpha);
	}

	// input is a matrix

	// check square matrix
	if (!internal::_check_square_mat(A))
		throw Exception("renyi", Exception::Type::MATRIX_NOT_SQUARE);

	if (alpha == 0)
		return std::log2(A.rows());

	// get the eigenvalues
	types::cmat ev = evals(A);
	double result = 0;
	// take the absolut values of the entries
	// to get rid of unwanted imaginary parts
	for (size_t i = 0; i < static_cast<size_t>(ev.rows()); i++)
		if (std::abs((types::cplx) ev(i)) != 0)
			result += std::pow(std::abs((types::cplx) ev(i)), alpha);

	return std::log2(result) / (1 - alpha);
}

// Renyi-infinity entropy (min entropy)
template<typename Scalar>
double renyi_inf(const types::DynMat<Scalar> & A) throw (Exception)
{
	// check zero-size
	if (!internal::_check_nonzero_size(A))
		throw Exception("renyi_inf", Exception::Type::MATRIX_ZERO_SIZE);

	// input is a vector
	if (A.rows() == 1 || A.cols() == 1)
	{
		double max = 0;
		for (size_t i = 0; i < static_cast<size_t>(A.size()); i++)
			if (std::abs(A(i)) > max)
				max = std::abs(A(i));

		return -std::log2(max);
	}

	// input is a matrix

	// check square matrix
	if (!internal::_check_square_mat(A))
		throw Exception("renyi_inf", Exception::Type::MATRIX_NOT_SQUARE);

	// get the eigenvalues
	types::cmat ev = evals(A);
	double max = 0;
	// take the absolut values of the entries
	// to get rid of unwanted imaginary parts
	for (size_t i = 0; i < static_cast<size_t>(ev.size()); i++)
		if (std::abs((types::cplx) ev(i)) > max)
			max = std::abs((types::cplx) ev(i));

	return -std::log2(max);
}

}
#endif /* ENTROPY_H_ */
