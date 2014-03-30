/*
 * entropy.h
 *
 *  Created on: Mar 27, 2014
 *      Author: vlad
 */

#ifndef ENTROPY_H_
#define ENTROPY_H_

#include <cmath>
#include <stdexcept>
#include "types.h"
#include "util.h"

// entropy functions

namespace qpp
{

// Shannon/von-Neumann entropy with log in given base (default = 2)
template<typename MatrixType>
double shannon(const types::EigenExpression<MatrixType> & A, double base = 2)
{
	// vector
	if (A.rows() == 1 || A.cols() == 1)
	{
		double result = 0;
		// take the absolut values of the entries to get rid of unwanted imaginary parts
		for (size_t i = 0; i < A.size(); i++)
			if (std::abs(A(i)) != 0)
				result -= std::abs(A(i)) * std::log2(std::abs(A(i)));

		return result / std::log2(base);
	}

	// matrix
	// check square matrix
	if (!internal::_check_square_mat(A))
		throw std::runtime_error(
				"shannon: Input must be a row/column vector or a square matrix!");
	// get the eigenvalues
	types::cmat ev = evals(A);
	double result = 0;
	// take the absolut values of the entries to get rid of unwanted imaginary parts
	for (size_t i = 0; i < ev.rows(); i++)
		if (std::abs((types::cplx) ev(i)) != 0)
			result -= std::abs((types::cplx) ev(i))
					* std::log2(std::abs((types::cplx) ev(i)));

	return result / std::log2(base);
}

// Renyi-alpha entropy (alpha>=0) with log in given base (default = 2)
template<typename MatrixType>
double renyi(const double alpha, const types::EigenExpression<MatrixType> & A,
		double base = 2)
{
	if (alpha < 0)
		throw std::runtime_error("renyi: alpha can not be negative!");

	if (alpha == 1) // Shannon/von Neumann
		return shannon(A, base);

	// vector
	if (A.rows() == 1 || A.cols() == 1)
	{
		if (alpha == 0)
			return std::log2(A.size()) / std::log2(base);

		double result = 0;
		// take the absolut values of the entries to get rid of unwanted imaginary parts
		for (size_t i = 0; i < A.size(); i++)
			if (std::abs((types::cplx) A(i)) != 0)
				result += std::pow(std::abs(A(i)), alpha);

		return std::log2(result) / ((1 - alpha) * std::log2(base));
	}

	// matrix
	// check square matrix
	if (!internal::_check_square_mat(A))
		throw std::runtime_error(
				"renyi: Input must be a row/column vector or a square matrix!");

	if (alpha == 0)
		return std::log2(A.rows()) / std::log2(base);

	// get the eigenvalues
	types::cmat ev = evals(A);
	double result = 0;
	// take the absolut values of the entries to get rid of unwanted imaginary parts
	for (size_t i = 0; i < ev.rows(); i++)
		if (std::abs((types::cplx) ev(i)) != 0)
			result += std::pow(std::abs((types::cplx) ev(i)), alpha);

	return std::log2(result) / ((1 - alpha) * std::log2(base));
}

// Renyi-infinity entropy (min entropy) with log in given base (default = 2)
template<typename MatrixType>
double renyi_inf(const types::EigenExpression<MatrixType> & A, double base = 2)
{
	// vector
	if (A.rows() == 1 || A.cols() == 1)
	{
		double max = 0;
		for (size_t i = 0; i < A.size(); i++)
			if (std::abs(A(i)) > max)
				max = std::abs(A(i));

		return -std::log2(max) / std::log2(base);
	}

	// matrix
	// check square matrix
	if (!internal::_check_square_mat(A))
		throw std::runtime_error(
				"renyi_inf: Input must be a row/column vector or a square matrix!");

	// get the eigenvalues
	types::cmat ev = evals(A);
	double max = 0;
	// take the absolut values of the entries to get rid of unwanted imaginary parts
	for (size_t i = 0; i < ev.size(); i++)
		if (std::abs((types::cplx) ev(i)) > max)
			max = std::abs((types::cplx) ev(i));

	return -std::log2(max) / std::log2(base);
}

}
#endif /* ENTROPY_H_ */
