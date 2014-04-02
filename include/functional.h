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
template<typename Scalar>
types::cmat funm(const types::DynMat<Scalar> &A,
		types::cplx (*f)(const types::cplx &))
{
	// zero-size
	if (!internal::_check_nonzero_size(A))
		throw std::invalid_argument("funm: Zero-sized input!");

	// check square matrix
	if (!internal::_check_square_mat(A))
		throw std::invalid_argument("funm: Matrix must be square!");

	Eigen::ComplexEigenSolver<types::cmat> es(A.template cast<types::cplx>());
	types::cmat evects = es.eigenvectors();
	types::cmat evals = es.eigenvalues();
	for (size_t i = 0; i < static_cast<size_t>(evals.rows()); i++)
		evals(i) = (*f)(evals(i)); // apply f(x) to each eigenvalue

	types::cmat evalsdiag = evals.asDiagonal();

	return evects * evalsdiag * evects.inverse();
}

// TODO: check this
// Matrix absolute value, note the syntax of Lambda invocation
template<typename Scalar>
types::cmat absm(const types::DynMat<Scalar> &A)
{
	// zero-size
	if (!internal::_check_nonzero_size(A))
		throw std::invalid_argument("absm: Zero-sized input!");

	// check square matrix
	if (!internal::_check_square_mat(A))
		throw std::invalid_argument("absm: Matrix must be square!");

	return funm(adjoint(A) * A, [](const types::cplx & x)->types::cplx
	{	return std::sqrt(x);});
}

// Matrix exponential
template<typename Scalar>
types::cmat expm(const types::DynMat<Scalar> &A)
{
	// zero-size
	if (!internal::_check_nonzero_size(A))
		throw std::invalid_argument("expm: Zero-sized input!");

	// check square matrix
	if (!internal::_check_square_mat(A))
		throw std::invalid_argument("expm: Matrix must be square!");

	return funm(A, std::exp);
}

// Matrix logarithm
template<typename Scalar>
types::cmat logm(const types::DynMat<Scalar> &A)
{
	// zero-size
	if (!internal::_check_nonzero_size(A))
		throw std::invalid_argument("logm: Zero-sized input!");

	// check square matrix
	if (!internal::_check_square_mat(A))
		throw std::invalid_argument("logm: Matrix must be square!");

	return funm(A, std::log);
}

// Matrix square root
template<typename Scalar>
types::cmat sqrtm(const types::DynMat<Scalar> &A)
{
	// zero-size
	if (!internal::_check_nonzero_size(A))
		throw std::invalid_argument("sqrtm: Zero-sized input!");

	// check square matrix
	if (!internal::_check_square_mat(A))
		throw std::invalid_argument("sqrtm: Matrix must be square!");

	return funm(A, std::sqrt);
}

// Matrix sin
template<typename Scalar>
types::cmat sinm(const types::DynMat<Scalar> &A)
{
	// zero-size
	if (!internal::_check_nonzero_size(A))
		throw std::invalid_argument("sinm: Zero-sized input!");

	// check square matrix
	if (!internal::_check_square_mat(A))
		throw std::invalid_argument("sinm: Matrix must be square!");

	return funm(A, std::sin);
}

// Matrix cos
template<typename Scalar>
types::cmat cosm(const types::DynMat<Scalar> &A)
{
	// zero-size
	if (!internal::_check_nonzero_size(A))
		throw std::invalid_argument("cosm: Zero-sized input!");

	// check square matrix
	if (!internal::_check_square_mat(A))
		throw std::invalid_argument("cosm: Matrix must be square!");

	return funm(A, std::cos);
}

// Matrix power A^z (CHANGES return type to complex matrix)
template<typename Scalar>
types::cmat powm(const types::DynMat<Scalar> &A, const types::cplx z)

{
	// zero-size
	if (!internal::_check_nonzero_size(A))
		throw std::invalid_argument("powm: Zero-sized input!");

	// check square matrix
	if (!internal::_check_square_mat(A))
		throw std::invalid_argument("powm: Matrix must be square!");

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
	for (size_t i = 0; i < static_cast<size_t>(evals.rows()); i++)
		evals(i) = std::pow(static_cast<types::cplx>(evals(i)),
				static_cast<types::cplx>(z));

	types::cmat evalsdiag = evals.asDiagonal();

	return evects * evalsdiag * evects.inverse();

}

// Matrix integer power, preserve return type
// Explicitly multiply the matrix with itself n times
template<typename Scalar>
types::DynMat<Scalar> powm_int(const types::DynMat<Scalar> &A, size_t n)
{
	// zero-size
	if (!internal::_check_nonzero_size(A))
		throw std::invalid_argument("powm_int: Zero-sized input!");

	// check square matrix
	if (!internal::_check_square_mat(A))
		throw std::invalid_argument("powm_int: Matrix must be square!");

	types::DynMat<Scalar> result = A;

	if (n == 0)
		return result.setIdentity();

	for (size_t i = 1; i < n; i++)
		result *= A;

	return result;
}

}

#endif /* FUNCTIONAL_H_ */
