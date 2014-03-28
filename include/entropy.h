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

// von-Neumann entropy with log_base
template<typename Derived>
double entropy(const Eigen::MatrixBase<Derived> & A, double base = 2)
{
	// vector
	if (A.rows() == 1 || A.cols() == 1)
	{
		double result = 0;
		// take the absolut values of the entries to get rid of unwanted imaginary parts
		for (size_t i = 0; i < A.size(); i++)
			if (std::abs((types::cplx) A(i)) != 0)
				result -= std::abs((types::cplx) A(i))
						* std::log2(std::abs((types::cplx) A(i)));
		return result / std::log2(base);
	}

	// matrix
	// check square matrix
	if (!internal::_check_square_mat(A))
		throw std::runtime_error("entropy: Input must be a row/column vector or a square matrix!");
	// get the eigenvalues
	Eigen::MatrixXcd ev = evals(A);
	double result = 0;
	// take the absolut values of the entries to get rid of unwanted imaginary parts
	for (size_t i = 0; i < ev.rows(); i++)
		if (std::abs((types::cplx) ev(i)) != 0)
			result -= std::abs((types::cplx) ev(i))
					* std::log2(std::abs((types::cplx) ev(i)));
	return result / std::log2(base);
}


}
#endif /* ENTROPY_H_ */
