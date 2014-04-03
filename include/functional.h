/*
 * functional.h
 *
 *  Created on: Mar 27, 2014
 *      Author: vlad
 */

#ifndef FUNCTIONAL_H_
#define FUNCTIONAL_H_

#include <cmath>
#include "types.h"
#include "internal.h"
#include "exception.h"

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
	// check zero-size
	if (!internal::_check_nonzero_size(A))
		throw Exception("funm", Exception::Type::MATRIX_ZERO_SIZE);

	// check square matrix
	if (!internal::_check_square_mat(A))
		throw Exception("funm", Exception::Type::MATRIX_NOT_SQUARE);

	Eigen::ComplexEigenSolver<types::cmat> es(A.template cast<types::cplx>());
	types::cmat evects = es.eigenvectors();
	types::cmat evals = es.eigenvalues();
	for (size_t i = 0; i < static_cast<size_t>(evals.rows()); i++)
		evals(i) = (*f)(evals(i)); // apply f(x) to each eigenvalue

	types::cmat evalsdiag = evals.asDiagonal();

	return evects * evalsdiag * evects.inverse();
}

// Matrix absolute value, note the syntax of Lambda invocation
template<typename Scalar>
types::cmat absm(const types::DynMat<Scalar> &A)
{
	// check zero-size
	if (!internal::_check_nonzero_size(A))
		throw Exception("absm", Exception::Type::MATRIX_ZERO_SIZE);

	// check square matrix
	if (!internal::_check_square_mat(A))
		throw Exception("absm", Exception::Type::MATRIX_NOT_SQUARE);

	return funm(adjoint(A) * A, [](const types::cplx & x)->types::cplx
	{	return std::sqrt(x);});
}

// Matrix exponential
template<typename Scalar>
types::cmat expm(const types::DynMat<Scalar> &A)
{
	// check zero-size
	if (!internal::_check_nonzero_size(A))
		throw Exception("expm", Exception::Type::MATRIX_ZERO_SIZE);

	// check square matrix
	if (!internal::_check_square_mat(A))
		throw Exception("expm", Exception::Type::MATRIX_NOT_SQUARE);

	return funm(A, std::exp);
}

// Matrix logarithm
template<typename Scalar>
types::cmat logm(const types::DynMat<Scalar> &A)
{
	// check zero-size
	if (!internal::_check_nonzero_size(A))
		throw Exception("logm", Exception::Type::MATRIX_ZERO_SIZE);

	// check square matrix
	if (!internal::_check_square_mat(A))
		throw Exception("logm", Exception::Type::MATRIX_NOT_SQUARE);

	return funm(A, std::log);
}

// Matrix square root
template<typename Scalar>
types::cmat sqrtm(const types::DynMat<Scalar> &A)
{
	// check zero-size
	if (!internal::_check_nonzero_size(A))
		throw Exception("sqrtm", Exception::Type::MATRIX_ZERO_SIZE);

	// check square matrix
	if (!internal::_check_square_mat(A))
		throw Exception("sqrtm", Exception::Type::MATRIX_NOT_SQUARE);

	return funm(A, std::sqrt);
}

// Matrix sin
template<typename Scalar>
types::cmat sinm(const types::DynMat<Scalar> &A)
{
	// check zero-size
	if (!internal::_check_nonzero_size(A))
		throw Exception("sinm", Exception::Type::MATRIX_ZERO_SIZE);

	// check square matrix
	if (!internal::_check_square_mat(A))
		throw Exception("sinm", Exception::Type::MATRIX_NOT_SQUARE);

	return funm(A, std::sin);
}

// Matrix cos
template<typename Scalar>
types::cmat cosm(const types::DynMat<Scalar> &A)
{
	// check zero-size
	if (!internal::_check_nonzero_size(A))
		throw Exception("cosm", Exception::Type::MATRIX_ZERO_SIZE);

	// check square matrix
	if (!internal::_check_square_mat(A))
		throw Exception("cosm", Exception::Type::MATRIX_NOT_SQUARE);

	return funm(A, std::cos);
}

// Matrix power A^z (CHANGES return type to complex matrix)
template<typename Scalar>
types::cmat powm(const types::DynMat<Scalar> &A, const types::cplx z)

{
	// check zero-size
	if (!internal::_check_nonzero_size(A))
		throw Exception("powm", Exception::Type::MATRIX_ZERO_SIZE);

	// check square matrix
	if (!internal::_check_square_mat(A))
		throw Exception("powm", Exception::Type::MATRIX_NOT_SQUARE);

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
	// check zero-size
	if (!internal::_check_nonzero_size(A))
		throw Exception("powm_int", Exception::Type::MATRIX_ZERO_SIZE);

	// check square matrix
	if (!internal::_check_square_mat(A))
		throw Exception("powm_int", Exception::Type::MATRIX_NOT_SQUARE);

	types::DynMat<Scalar> result = A;

	if (n == 0)
		return result.setIdentity();

	for (size_t i = 1; i < n; i++)
		result *= A;

	return result;
}

}

#endif /* FUNCTIONAL_H_ */
