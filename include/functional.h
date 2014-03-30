/*
 * functional.h
 *
 *  Created on: Mar 27, 2014
 *      Author: vlad
 */

#ifndef FUNCTIONAL_H_
#define FUNCTIONAL_H_

#include <stdexcept>
#include <cmath>
#include "types.h"
#include "internal.h"

// Matrix functional calculus

namespace qpp
{

// Computes f(A), where (*f) is the function pointer
/**
 *
 * @param A input matrix
 * @param f function pointer
 * @return types::cmat
 */
template<typename MatrixType>
types::cmat funm(const types::EigenExpression<MatrixType> &A,
		types::cplx (*f)(const types::cplx &))
{
	// check square matrix
	if (!internal::_check_square_mat(A))
		throw std::runtime_error("funm: Matrix must be square!");

	Eigen::ComplexEigenSolver<types::cmat> es(A.template cast<types::cplx>());
	types::cmat evects = es.eigenvectors();
	types::cmat evals = es.eigenvalues();
	for (int i = 0; i < evals.rows(); i++)
		evals(i) = (*f)(evals(i)); // apply f(x) to each eigenvalue

	types::cmat evalsdiag = evals.asDiagonal();

	return evects * evalsdiag * evects.inverse();
}

// Matrix absolute value, note the syntax of Lambda invocation
template<typename MatrixType>
types::cmat absm(const types::EigenExpression<MatrixType> &A)
{
	return funm(adjoint(A) * A, [](const types::cplx & x)->types::cplx
	{	return std::sqrt(x);});
}

// Matrix exponential
template<typename MatrixType>
types::cmat expm(const types::EigenExpression<MatrixType> &A)
{
	return funm(A, std::exp);
}

// Matrix logarithm
template<typename MatrixType>
types::cmat logm(const types::EigenExpression<MatrixType> &A)
{
	return funm(A, std::log);
}

// Matrix square root
template<typename MatrixType>
types::cmat sqrtm(const types::EigenExpression<MatrixType> &A)
{
	return funm(A, std::sqrt);
}

// Matrix sin
template<typename MatrixType>
types::cmat sinm(const types::EigenExpression<MatrixType> &A)
{
	return funm(A, std::sin);
}

// Matrix cos
template<typename MatrixType>
types::cmat cosm(const types::EigenExpression<MatrixType> &A)
{
	return funm(A, std::cos);
}

// Matrix power A^z (CHANGES return type to complex matrix)
template<typename MatrixType>
types::cmat powm(const types::EigenExpression<MatrixType> &A,
		const types::cplx z)

{
	// check square matrix
	if (!internal::_check_square_mat(A))
		throw std::runtime_error("mpower: Matrix must be square!");

	// Define A^0 = Id
	if (real(z) == 0 && imag(z) == 0)
	{
		types::cmat result(A.rows(), A.rows());
		result.setIdentity();
		return result;
	}

	Eigen::ComplexEigenSolver<types::cmat> es(A.template cast<types::cplx>());
	types::cmat evects = es.eigenvectors();
	types::cmat evals = es.eigenvalues();
	for (int i = 0; i < evals.rows(); i++)
		evals(i) = std::pow(static_cast<types::cplx>(evals(i)),
				static_cast<types::cplx>(z));

	types::cmat evalsdiag = evals.asDiagonal();

	return evects * evalsdiag * evects.inverse();

}

// Matrix integer power, preserve return type
// Explicitly multiply the matrix with itself n times
template<typename MatrixType>
types::TemplatedEigenMatrix<MatrixType> powm_int(
		const types::EigenExpression<MatrixType> &A, size_t n)
{
// check square matrix
	if (!internal::_check_square_mat(A))
		throw std::runtime_error("mpower_n: Matrix must be square!");

	types::TemplatedEigenMatrix<MatrixType> result = A;

	if (n == 0)
		return result.setIdentity();

	for (size_t i = 1; i < n; i++)
		result *= A;

	return result;
}

}

#endif /* FUNCTIONAL_H_ */
