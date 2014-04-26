/*
 * functions.h
 *
 *  Created on: Mar 27, 2014
 *      Author: vlad
 */

#ifndef FUNCTIONS_H_
#define FUNCTIONS_H_

#include <cmath>
#include <numeric>
#include <vector>

#include "constants.h"
#include "internal.h"
#include "types.h"
#include "classes/exception.h"

// Collection of quantum computing useful functions
namespace qpp
{
// Eigen function wrappers

// transpose, preserve return type
template<typename Derived>
types::DynMat<typename Derived::Scalar> transpose(
		const Eigen::MatrixBase<Derived>& A)

{
	const types::DynMat<typename Derived::Scalar> & rA = A;

	// check zero-size
	if (!internal::_check_nonzero_size(rA))
		throw Exception("transpose", Exception::Type::ZERO_SIZE);

	return rA.transpose();
}

// conjugate, preserve return type
template<typename Derived>
types::DynMat<typename Derived::Scalar> conjugate(
		const Eigen::MatrixBase<Derived>& A)

{
	const types::DynMat<typename Derived::Scalar> & rA = A;

	// check zero-size
	if (!internal::_check_nonzero_size(rA))
		throw Exception("conjugate", Exception::Type::ZERO_SIZE);

	return rA.conjugate();
}

// adjoint, preserve return type
template<typename Derived>
types::DynMat<typename Derived::Scalar> adjoint(
		const Eigen::MatrixBase<Derived>& A)
{
	const types::DynMat<typename Derived::Scalar> & rA = A;

	// check zero-size
	if (!internal::_check_nonzero_size(rA))
		throw Exception("adjoint", Exception::Type::ZERO_SIZE);

	return rA.adjoint();
}

// inverse, preserve return type
template<typename Derived>
types::DynMat<typename Derived::Scalar> inverse(
		const Eigen::MatrixBase<Derived>& A)
{
	const types::DynMat<typename Derived::Scalar> & rA = A;

	// check zero-size
	if (!internal::_check_nonzero_size(rA))
		throw Exception("inverse", Exception::Type::ZERO_SIZE);

	return rA.inverse();
}

// trace, preserve return type
template<typename Derived>
typename Derived::Scalar trace(const Eigen::MatrixBase<Derived>& A)

{
	const types::DynMat<typename Derived::Scalar> & rA = A;

	// check zero-size
	if (!internal::_check_nonzero_size(rA))
		throw Exception("trace", Exception::Type::ZERO_SIZE);

	return rA.trace();
}

// determinant, preserve return type
template<typename Derived>
typename Derived::Scalar det(const Eigen::MatrixBase<Derived>& A)

{
	const types::DynMat<typename Derived::Scalar> & rA = A;

	// check zero-size
	if (!internal::_check_nonzero_size(rA))
		throw Exception("det", Exception::Type::ZERO_SIZE);

	return rA.determinant();
}

// logarithm of the determinant, preserve return type
// especially useful when determinant overflows/underflows
// returns +/- inf when determinant overflows/underflows
template<typename Derived>
typename Derived::Scalar logdet(const Eigen::MatrixBase<Derived>& A)

{
	const types::DynMat<typename Derived::Scalar> & rA = A;

	// check zero-size
	if (!internal::_check_nonzero_size(rA))
		throw Exception("logdet", Exception::Type::ZERO_SIZE);

	// check square matrix
	if (!internal::_check_square_mat(rA))
		throw Exception("logdet", Exception::Type::MATRIX_NOT_SQUARE);

	Eigen::PartialPivLU<types::DynMat<typename Derived::Scalar>> lu(rA);
	types::DynMat<typename Derived::Scalar> U =
			lu.matrixLU().template triangularView<Eigen::Upper>();
	typename Derived::Scalar result = std::log(U(0, 0));
	for (size_t i = 1; i < static_cast<size_t>(rA.rows()); i++)
		result += std::log(U(i, i));

	return result;

}

// element-wise sum, preserve return type
template<typename Derived>
typename Derived::Scalar sum(const Eigen::MatrixBase<Derived>& A)

{
	const types::DynMat<typename Derived::Scalar> & rA = A;

// check zero-size
	if (!internal::_check_nonzero_size(rA))
		throw Exception("sum", Exception::Type::ZERO_SIZE);

	return rA.sum();
}

// trace-norm (or Frobenius norm) (CHANGES return type to double matrix)
template<typename Derived>
double norm(const Eigen::MatrixBase<Derived>& A)
{
	const types::DynMat<typename Derived::Scalar> & rA = A;

// check zero-size
	if (!internal::_check_nonzero_size(rA))
		throw Exception("norm", Exception::Type::ZERO_SIZE);

// convert matrix to complex then return its norm
	return (rA.template cast<types::cplx>()).norm();
}

// eigenvalues (CHANGES return type to complex matrix)
template<typename Derived>
types::cmat evals(const Eigen::MatrixBase<Derived>& A)
{
	const types::DynMat<typename Derived::Scalar> & rA = A;

// check zero-size
	if (!internal::_check_nonzero_size(rA))
		throw Exception("evals", Exception::Type::ZERO_SIZE);

// check square matrix
	if (!internal::_check_square_mat(rA))
		throw Exception("evals", Exception::Type::MATRIX_NOT_SQUARE);

	Eigen::ComplexEigenSolver<types::cmat> es(rA.template cast<types::cplx>());
	return es.eigenvalues();
}

// eigenvectors (CHANGES return type to complex matrix)
template<typename Derived>
types::cmat evects(const Eigen::MatrixBase<Derived>& A)
{
	const types::DynMat<typename Derived::Scalar> & rA = A;

// check zero-size
	if (!internal::_check_nonzero_size(rA))
		throw Exception("evects", Exception::Type::ZERO_SIZE);

// check square matrix
	if (!internal::_check_square_mat(rA))
		throw Exception("evects", Exception::Type::MATRIX_NOT_SQUARE);

	Eigen::ComplexEigenSolver<types::cmat> es(rA.template cast<types::cplx>());
	return es.eigenvectors();
}

// eigenvalues of Hermitian matrices (CHANGES return type to double matrix)
template<typename Derived>
types::dmat hevals(const Eigen::MatrixBase<Derived>& A)
{
	const types::DynMat<typename Derived::Scalar> & rA = A;

// check zero-size
	if (!internal::_check_nonzero_size(rA))
		throw Exception("hevals", Exception::Type::ZERO_SIZE);

// check square matrix
	if (!internal::_check_square_mat(rA))
		throw Exception("hevals", Exception::Type::MATRIX_NOT_SQUARE);

	Eigen::SelfAdjointEigenSolver<types::cmat> es(
			rA.template cast<types::cplx>());
	return es.eigenvalues();
}

// eigenvectors of Hermitian matrix (CHANGES return type to complex matrix)
template<typename Derived>
types::cmat hevects(const Eigen::MatrixBase<Derived>& A)
{
	const types::DynMat<typename Derived::Scalar> & rA = A;

// check zero-size
	if (!internal::_check_nonzero_size(rA))
		throw Exception("hevects", Exception::Type::ZERO_SIZE);

// check square matrix
	if (!internal::_check_square_mat(rA))
		throw Exception("hevects", Exception::Type::MATRIX_NOT_SQUARE);

	Eigen::SelfAdjointEigenSolver<types::cmat> es(
			rA.template cast<types::cplx>());
	return es.eigenvectors();
}

// Matrix functional calculus

// Computes f(A), where (*f) is the function pointer
// (CHANGES return type to complex matrix)
/**
 *
 * @param A input matrix
 * @param f function pointer
 * @return types::cmat
 */
template<typename Derived>
types::cmat funm(const Eigen::MatrixBase<Derived> &A,
		types::cplx (*f)(const types::cplx &))
{
	const types::DynMat<typename Derived::Scalar> & rA = A;

// check zero-size
	if (!internal::_check_nonzero_size(rA))
		throw Exception("funm", Exception::Type::ZERO_SIZE);

// check square matrix
	if (!internal::_check_square_mat(rA))
		throw Exception("funm", Exception::Type::MATRIX_NOT_SQUARE);

	Eigen::ComplexEigenSolver<types::cmat> es(rA.template cast<types::cplx>());
	types::cmat evects = es.eigenvectors();
	types::cmat evals = es.eigenvalues();
	for (size_t i = 0; i < static_cast<size_t>(evals.rows()); i++)
		evals(i) = (*f)(evals(i)); // apply f(x) to each eigenvalue

	types::cmat evalsdiag = evals.asDiagonal();

	return evects * evalsdiag * evects.inverse();
}

// Matrix square root
// (CHANGES return type to complex matrix)
template<typename Derived>
types::cmat sqrtm(const Eigen::MatrixBase<Derived> &A)
{
	const types::DynMat<typename Derived::Scalar> & rA = A;

// check zero-size
	if (!internal::_check_nonzero_size(rA))
		throw Exception("sqrtm", Exception::Type::ZERO_SIZE);

// check square matrix
	if (!internal::_check_square_mat(rA))
		throw Exception("sqrtm", Exception::Type::MATRIX_NOT_SQUARE);

	return funm(rA, &std::sqrt);
}

// Matrix absolute value, note the syntax of Lambda invocation
// (CHANGES return type to complex matrix)
template<typename Derived>
types::cmat absm(const Eigen::MatrixBase<Derived> &A)
{
	const types::DynMat<typename Derived::Scalar> & rA = A;

// check zero-size
	if (!internal::_check_nonzero_size(rA))
		throw Exception("absm", Exception::Type::ZERO_SIZE);

// check square matrix
	if (!internal::_check_square_mat(rA))
		throw Exception("absm", Exception::Type::MATRIX_NOT_SQUARE);

	return sqrtm(adjoint(rA) * rA);
}

// Matrix exponential
// (CHANGES return type to complex matrix)
template<typename Derived>
types::cmat expm(const Eigen::MatrixBase<Derived> &A)
{
	const types::DynMat<typename Derived::Scalar> & rA = A;

// check zero-size
	if (!internal::_check_nonzero_size(rA))
		throw Exception("expm", Exception::Type::ZERO_SIZE);

// check square matrix
	if (!internal::_check_square_mat(rA))
		throw Exception("expm", Exception::Type::MATRIX_NOT_SQUARE);

	return funm(rA, &std::exp);
}

// Matrix logarithm
// (CHANGES return type to complex matrix)
template<typename Derived>
types::cmat logm(const Eigen::MatrixBase<Derived> &A)
{
	const types::DynMat<typename Derived::Scalar> & rA = A;

// check zero-size
	if (!internal::_check_nonzero_size(rA))
		throw Exception("logm", Exception::Type::ZERO_SIZE);

// check square matrix
	if (!internal::_check_square_mat(rA))
		throw Exception("logm", Exception::Type::MATRIX_NOT_SQUARE);

	return funm(rA, &std::log);
}

// Matrix sin
// (CHANGES return type to complex matrix)
template<typename Derived>
types::cmat sinm(const Eigen::MatrixBase<Derived> &A)
{
	const types::DynMat<typename Derived::Scalar> & rA = A;

// check zero-size
	if (!internal::_check_nonzero_size(rA))
		throw Exception("sinm", Exception::Type::ZERO_SIZE);

// check square matrix
	if (!internal::_check_square_mat(rA))
		throw Exception("sinm", Exception::Type::MATRIX_NOT_SQUARE);

	return funm(rA, &std::sin);
}

// Matrix cos
// (CHANGES return type to complex matrix)
template<typename Derived>
types::cmat cosm(const Eigen::MatrixBase<Derived> &A)
{
	const types::DynMat<typename Derived::Scalar> & rA = A;

// check zero-size
	if (!internal::_check_nonzero_size(rA))
		throw Exception("cosm", Exception::Type::ZERO_SIZE);

// check square matrix
	if (!internal::_check_square_mat(rA))
		throw Exception("cosm", Exception::Type::MATRIX_NOT_SQUARE);

	return funm(rA, &std::cos);
}

// Matrix power A^z
// (CHANGES return type to complex matrix)
template<typename Derived>
types::cmat spectralpowm(const Eigen::MatrixBase<Derived> &A,
		const types::cplx z)

{
	const types::DynMat<typename Derived::Scalar> & rA = A;

// check zero-size
	if (!internal::_check_nonzero_size(rA))
		throw Exception("spectralpowm", Exception::Type::ZERO_SIZE);

// check square matrix
	if (!internal::_check_square_mat(rA))
		throw Exception("spectralpowm", Exception::Type::MATRIX_NOT_SQUARE);

// Define A^0 = Id, for z IDENTICALLY zero
	if (real(z) == 0 && imag(z) == 0)
	{
		types::cmat result(rA.rows(), rA.rows());
		result.setIdentity();
		return result;
	}

	Eigen::ComplexEigenSolver<types::cmat> es(rA.template cast<types::cplx>());
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
template<typename Derived>
types::DynMat<typename Derived::Scalar> powm(
		const Eigen::MatrixBase<Derived> &A, size_t n)

{
	const types::DynMat<typename Derived::Scalar> & rA = A;

// check zero-size
	if (!internal::_check_nonzero_size(rA))
		throw Exception("powm", Exception::Type::ZERO_SIZE);

// check square matrix
	if (!internal::_check_square_mat(rA))
		throw Exception("powm", Exception::Type::MATRIX_NOT_SQUARE);

	types::DynMat<typename Derived::Scalar> result = rA;

	if (n == 0)
		return result.setIdentity();

	for (size_t i = 1; i < n; i++)
		result *= rA;

	return result;
}

// other functions

// functor; apply f(A) component-wise, where (*f) is the function pointer
// returns a matrix of type OutputScalar
template<typename OutputScalar, typename Derived>
types::DynMat<OutputScalar> cwise(const Eigen::MatrixBase<Derived> &A,
		OutputScalar (*f)(const typename Derived::Scalar &))
{
	const types::DynMat<typename Derived::Scalar> & rA = A;

// check zero-size
	if (!internal::_check_nonzero_size(rA))
		throw Exception("cwise", Exception::Type::ZERO_SIZE);

	types::DynMat<OutputScalar> result(rA.rows(), rA.cols());

	for (size_t j = 0; j < static_cast<size_t>(rA.cols()); j++)
#pragma omp parallel for
		for (size_t i = 0; i < static_cast<size_t>(rA.rows()); i++)
			result(i, j) = (*f)(rA(i, j));

	return result;
}

// Kronecker product of multiple matrices, preserve return type
// variadic template
template<typename T>
types::DynMat<typename T::Scalar> kron(const T& head)
{
	return head;
}
template<typename T, typename ... Args>
types::DynMat<typename T::Scalar> kron(const T& head, const Args&... tail)
{
	return internal::_kron(head, kron(tail...));
}

// Kronecker product of a list (std::vector) of matrices, preserve return type
template<typename Derived>
types::DynMat<typename Derived::Scalar> kron(const std::vector<Derived> &As)

{
	if (As.size() == 0)
		throw Exception("kron", Exception::Type::ZERO_SIZE);

	for (auto it : As)
		if (it.size() == 0)
			throw Exception("kron", Exception::Type::ZERO_SIZE);

	types::DynMat<typename Derived::Scalar> result = As[0];
	for (size_t i = 1; i < As.size(); i++)
	{
		result = kron(result, As[i]);
	}
	return result;
}

// Kronecker product of a list of matrices, preserve return type
// deduce the template parameters from initializer_list
template<typename Derived>
types::DynMat<typename Derived::Scalar> kron(
		const std::initializer_list<Derived> &As)
{
	return kron(std::vector<Derived>(As));
}

// Kronecker product of a matrix with itself $n$ times, preserve return type
template<typename Derived>
types::DynMat<typename Derived::Scalar> kronpow(
		const Eigen::MatrixBase<Derived>& A, size_t n)

{
	const types::DynMat<typename Derived::Scalar> & rA = A;

// check zero-size
	if (!internal::_check_nonzero_size(rA))
		throw Exception("kronpow", Exception::Type::ZERO_SIZE);

// check out of range
	if (n == 0)
		throw Exception("kronpow", Exception::Type::OUT_OF_RANGE);

	types::DynMat<typename Derived::Scalar> result = rA;
	for (size_t i = 1; i < n; i++)
		result = kron(result, rA);
	return result;
}

// reshape the columns of A and returns a matrix with m rows and n columns
// use column-major order (same as MATLAB)
template<typename Derived>
types::DynMat<typename Derived::Scalar> reshape(
		const Eigen::MatrixBase<Derived>& A, size_t rows, size_t cols)
{
	const types::DynMat<typename Derived::Scalar> & rA = A;

	size_t Arows = static_cast<size_t>(rA.rows());
	size_t Acols = static_cast<size_t>(rA.cols());

// check zero-size
	if (!internal::_check_nonzero_size(rA))
		throw Exception("reshape", Exception::Type::ZERO_SIZE);

	if (Arows * Acols != rows * cols)
		throw Exception("reshape", Exception::Type::DIMS_MISMATCH_MATRIX);

	return Eigen::Map<types::DynMat<typename Derived::Scalar>>(
			const_cast<typename Derived::Scalar*>(rA.data()), rows, cols);
}

// permutes the subsystems in a matrix
// perm specifies the destination permutation,
// i.e. the qubit perm[i] is permuted to location i, perm[i]->i
template<typename Derived>
types::DynMat<typename Derived::Scalar> syspermute(
		const Eigen::MatrixBase<Derived>& A, const std::vector<size_t>& perm,
		const std::vector<size_t> &dims)

{
	const types::DynMat<typename Derived::Scalar> & rA = A;

// Error checks

	// check zero-size
	if (!internal::_check_nonzero_size(rA))
		throw Exception("syspermute", Exception::Type::ZERO_SIZE);

	// check that dims is a valid dimension vector
	if (!internal::_check_dims(dims))
		throw Exception("syspermute", Exception::Type::DIMS_INVALID);

	// check that we have a valid permutation
	if (!internal::_check_perm(perm))
		throw Exception("syspermute", Exception::Type::PERM_INVALID);

	// check permutation size
	if (perm.size() != dims.size())
		throw Exception("syspermute", Exception::Type::PERM_INVALID);

	size_t D = static_cast<size_t>(rA.rows());
	size_t numdims = dims.size();

	types::DynMat<typename Derived::Scalar> result;

	auto worker = [](size_t i, size_t numdims, const size_t* cdims,
			const size_t* cperm)
	{
		// use static allocation for speed,
		// double the size for matrices reshaped as vectors
			size_t midx[2 * ct::maxn];
			size_t midxtmp[2 * ct::maxn];
			size_t permdims[2 * ct::maxn];

			/* compute the multi-index */
			internal::_n2multiidx(i, numdims, cdims, midx);

			for (size_t k = 0; k < numdims; k++)
			{
				permdims[k] = cdims[cperm[k]]; // permuted dimensions
				midxtmp[k] = midx[cperm[k]];// permuted multi-indexes
			}
			return internal::_multiidx2n(midxtmp, numdims, permdims);
		};

// check column vector
	if (internal::_check_col_vector(rA)) // we have a column vector
	{
		size_t cdims[ct::maxn];
		size_t cperm[ct::maxn];

		// check that dims match the dimension of rA
		if (!internal::_check_dims_match_cvect(dims, rA))
			throw Exception("syspermute",
					Exception::Type::DIMS_MISMATCH_CVECTOR);

		// copy dims in cdims and perm in cperm
		for (size_t i = 0; i < numdims; i++)
		{
			cdims[i] = dims[i];
			cperm[i] = perm[i];
		}
		result.resize(D, 1);

#pragma omp parallel for
		for (size_t i = 0; i < D; i++)
			result(worker(i, numdims, cdims, cperm)) = rA(i);

		return result;
	}

	else if (internal::_check_square_mat(rA)) // we have a square matrix
	{
		size_t cdims[2 * ct::maxn];
		size_t cperm[2 * ct::maxn];

		// check that dims match the dimension of rA
		if (!internal::_check_dims_match_mat(dims, rA))
			throw Exception("syspermute",
					Exception::Type::DIMS_MISMATCH_MATRIX);

		// copy dims in cdims and perm in cperm
		for (size_t i = 0; i < numdims; i++)
		{
			cdims[i] = dims[i];
			cdims[i + numdims] = dims[i];
			cperm[i] = perm[i];
			cperm[i + numdims] = perm[i] + numdims;
		}
		result.resize(D * D, 1);
		// map A to a column vector
		types::DynMat<typename Derived::Scalar> vectA = Eigen::Map<
				types::DynMat<typename Derived::Scalar>>(
				const_cast<typename Derived::Scalar*>(rA.data()), D * D, 1);

#pragma omp parallel for
		for (size_t i = 0; i < D * D; i++)
			result(worker(i, 2 * numdims, cdims, cperm)) = rA(i);

		return reshape(result, D, D);
	}

	else
		throw Exception("syspermute",
				Exception::Type::MATRIX_NOT_SQUARE_OR_CVECTOR);
}

// Partial trace over subsystem A in a D_A x D_B system
template<typename Derived>
types::DynMat<typename Derived::Scalar> ptrace1(
		const Eigen::MatrixBase<Derived>& A, const std::vector<size_t>& dims)
{
	const types::DynMat<typename Derived::Scalar> & rA = A;

// Error checks

// check zero-size
	if (!internal::_check_nonzero_size(rA))
		throw Exception("ptrace1", Exception::Type::ZERO_SIZE);

// check that dims is a valid dimension vector
	if (!internal::_check_dims(dims))
		throw Exception("ptrace1", Exception::Type::DIMS_INVALID);

// check square matrix
	if (!internal::_check_square_mat(rA))
		throw Exception("ptrace1", Exception::Type::MATRIX_NOT_SQUARE);

// check dims has only 2 elements
	if (dims.size() != 2)
		throw Exception("ptrace1", Exception::Type::NOT_BIPARTITE);

// check that dims match the dimension of A
	if (!internal::_check_dims_match_mat(dims, rA))
		throw Exception("ptrace1", Exception::Type::DIMS_MISMATCH_MATRIX);

	size_t DA = dims[0];
	size_t DB = dims[1];

	types::DynMat<typename Derived::Scalar> result = types::DynMat<
			typename Derived::Scalar>::Zero(DB, DB);

	auto worker = [&](size_t i, size_t j)
	{
		typename Derived::Scalar sum = 0;
		for (size_t m = 0; m < DA; m++)
		sum += rA(m * DB + i, m * DB + j);
		return sum;
	};

	for (size_t j = 0; j < DB; j++) // column major order for speed
#pragma omp parallel for
		for (size_t i = 0; i < DB; i++)
			result(i, j) = worker(i, j);

	return result;
}

// Partial trace over subsystem B in a D_A x D_B system
template<typename Derived>
types::DynMat<typename Derived::Scalar> ptrace2(
		const Eigen::MatrixBase<Derived>& A, const std::vector<size_t>& dims)
{
	const types::DynMat<typename Derived::Scalar> & rA = A;

// Error checks

// check zero-size
	if (!internal::_check_nonzero_size(rA))
		throw Exception("ptrace2", Exception::Type::ZERO_SIZE);

// check that dims is a valid dimension vector
	if (!internal::_check_dims(dims))
		throw Exception("ptrace2", Exception::Type::DIMS_INVALID);

// check square matrix
	if (!internal::_check_square_mat(rA))
		throw Exception("ptrace2", Exception::Type::MATRIX_NOT_SQUARE);

// check dims has only 2 elements
	if (dims.size() != 2)
		throw Exception("ptrace2", Exception::Type::NOT_BIPARTITE);

// check that dims match the dimension of A
	if (!internal::_check_dims_match_mat(dims, rA))
		throw Exception("ptrace2", Exception::Type::DIMS_MISMATCH_MATRIX);

	size_t DA = dims[0];
	size_t DB = dims[1];

	types::DynMat<typename Derived::Scalar> result = types::DynMat<
			typename Derived::Scalar>::Zero(DA, DA);

	for (size_t j = 0; j < DA; j++) // column major order for speed
#pragma omp parallel for
		for (size_t i = 0; i < DA; i++)
			result(i, j) = trace(rA.block(i * DB, j * DB, DB, DB));

	return result;
}

// partial trace
template<typename Derived>
types::DynMat<typename Derived::Scalar> ptrace(
		const Eigen::MatrixBase<Derived>& A, const std::vector<size_t>& subsys,
		const std::vector<size_t>& dims)

{
	const types::DynMat<typename Derived::Scalar> & rA = A;

// error checks

// check zero-size
	if (!internal::_check_nonzero_size(rA))
		throw Exception("ptrace", Exception::Type::ZERO_SIZE);

// check that dims is a valid dimension vector
	if (!internal::_check_dims(dims))
		throw Exception("ptrace", Exception::Type::DIMS_INVALID);

// check square matrix
	if (!internal::_check_square_mat(rA))
		throw Exception("ptrace", Exception::Type::MATRIX_NOT_SQUARE);

// check that dims match the dimension of A
	if (!internal::_check_dims_match_mat(dims, rA))
		throw Exception("ptrace", Exception::Type::DIMS_MISMATCH_MATRIX);

	if (subsys.size() == dims.size())
	{
		types::DynMat<typename Derived::Scalar> result = types::DynMat<
				typename Derived::Scalar>(1, 1);
		result(0, 0) = rA.trace();
		return result;
	}
	if (subsys.size() == 0)
		return rA;
	// check that subsys are valid
	if (!internal::_check_subsys_match_dims(subsys, dims))
		throw Exception("ptrace", Exception::Type::SUBSYS_MISMATCH_DIMS);

	size_t D = static_cast<size_t>(rA.rows());
	size_t n = dims.size();
	size_t nsubsys = subsys.size();
	size_t nsubsysbar = n - nsubsys;
	size_t dimsubsys = 1;
	for (size_t i = 0; i < nsubsys; i++)
		dimsubsys *= dims[subsys[i]];
	size_t dimsubsysbar = D / dimsubsys;

	size_t Cdims[ct::maxn];
	size_t Csubsys[ct::maxn];
	size_t Cdimssubsys[ct::maxn];
	size_t Csubsysbar[ct::maxn];
	size_t Cdimssubsysbar[ct::maxn];

	for (size_t i = 0; i < n; i++)
		Cdims[i] = dims[i];
	for (size_t i = 0; i < nsubsys; i++)
	{
		Csubsys[i] = subsys[i];
		Cdimssubsys[i] = dims[subsys[i]];
	}
	// construct the complement of subsys
	size_t cnt = 0;
	for (size_t i = 0; i < n; i++)
	{
		bool found = false;
		for (size_t m = 0; m < nsubsys; m++)
			if (subsys[m] == i)
			{
				found = true;
				break;
			}
		if (!found)
		{
			Csubsysbar[cnt] = i;
			Cdimssubsysbar[cnt] = dims[i];
			cnt++;
		}
	}

	types::DynMat<typename Derived::Scalar> result = types::DynMat<
			typename Derived::Scalar>(dimsubsysbar, dimsubsysbar);

	auto worker = [&](size_t i, size_t j)
	{
		// use static allocation for speed!

			size_t Cmidxrow[ct::maxn];
			size_t Cmidxcol[ct::maxn];
			size_t Cmidxrowsubsysbar[ct::maxn];
			size_t Cmidxcolsubsysbar[ct::maxn];
			size_t Cmidxsubsys[ct::maxn];

			/* get the row/col multi-indexes of the complement */
			internal::_n2multiidx(i, nsubsysbar, Cdimssubsysbar, Cmidxrowsubsysbar);
			internal::_n2multiidx(j, nsubsysbar, Cdimssubsysbar, Cmidxcolsubsysbar);
			/* write them in the global row/col multi-indexes */
			for(size_t k=0;k<nsubsysbar;k++)
			{
				Cmidxrow[Csubsysbar[k]]=Cmidxrowsubsysbar[k];
				Cmidxcol[Csubsysbar[k]]=Cmidxcolsubsysbar[k];
			}
			typename Derived::Scalar sm = 0;
			for(size_t a=0; a<dimsubsys; a++)
			{
				// get the multi-index over which we do the summation
				internal::_n2multiidx(a, nsubsys, Cdimssubsys, Cmidxsubsys);
				// write it into the global row/col multi-indexes
				for(size_t k=0;k<nsubsys;k++)
				Cmidxrow[Csubsys[k]]=Cmidxcol[Csubsys[k]]=Cmidxsubsys[k];

				// now do the sum
				sm+= rA(internal::_multiidx2n(Cmidxrow,n,Cdims),
						internal::_multiidx2n(Cmidxcol,n,Cdims));
			}

			return sm;
		};

	for (size_t i = 0; i < dimsubsysbar; i++)
#pragma omp parallel for
		for (size_t j = 0; j < dimsubsysbar; j++)
			result(i, j) = worker(i, j);

	return result;
}

// partial transpose
template<typename Derived>
types::DynMat<typename Derived::Scalar> ptranspose(
		const Eigen::MatrixBase<Derived>& A, const std::vector<size_t>& subsys,
		const std::vector<size_t>& dims)

{
	const types::DynMat<typename Derived::Scalar> & rA = A;

// error checks

// check zero-size
	if (!internal::_check_nonzero_size(rA))
		throw Exception("ptranspose", Exception::Type::ZERO_SIZE);

// check that dims is a valid dimension vector
	if (!internal::_check_dims(dims))
		throw Exception("ptranspose", Exception::Type::DIMS_INVALID);

// check square matrix
	if (!internal::_check_square_mat(rA))
		throw Exception("ptranspose", Exception::Type::MATRIX_NOT_SQUARE);

// check that dims match the dimension of A
	if (!internal::_check_dims_match_mat(dims, rA))
		throw Exception("ptranspose", Exception::Type::DIMS_MISMATCH_MATRIX);

	if (subsys.size() == dims.size())
		return rA.transpose();
	if (subsys.size() == 0)
		return rA;
// check that subsys are valid
	if (!internal::_check_subsys_match_dims(subsys, dims))
		throw Exception("ptranspose", Exception::Type::SUBSYS_MISMATCH_DIMS);

	size_t D = static_cast<size_t>(rA.rows());
	size_t numdims = dims.size();
	size_t numsubsys = subsys.size();
	size_t cdims[ct::maxn];
	size_t midxcol[ct::maxn];
	size_t csubsys[ct::maxn];

// copy dims in cdims and subsys in csubsys
	for (size_t i = 0; i < numdims; i++)
		cdims[i] = dims[i];
	for (size_t i = 0; i < numsubsys; i++)
		csubsys[i] = subsys[i];

	types::DynMat<typename Derived::Scalar> result(D, D);

	auto worker = [&](size_t i, size_t j)
	{
		// use static allocation for speed!
			size_t midxcoltmp[ct::maxn];
			size_t midxrow[ct::maxn];

			for (size_t k = 0; k < numdims; k++)
			midxcoltmp[k] = midxcol[k];

			/* compute the row multi-index */
			internal::_n2multiidx(i, numdims, cdims, midxrow);

			for (size_t k = 0; k < numsubsys; k++)
			std::swap(midxcoltmp[csubsys[k]], midxrow[csubsys[k]]);

			/* writes the result */
			result(i, j)=rA(internal::_multiidx2n(midxrow, numdims, cdims),
					internal::_multiidx2n(midxcoltmp, numdims, cdims));

		};

	for (size_t j = 0; j < D; j++)
	{
		// compute the column multi-index
		internal::_n2multiidx(j, numdims, cdims, midxcol);
#pragma omp parallel for
		for (size_t i = 0; i < D; i++)
			worker(i, j);
	}

	return result;
}

// commutator
template<typename Derived1, typename Derived2>
types::DynMat<typename Derived1::Scalar> comm(
		const Eigen::MatrixBase<Derived1> &A,
		const Eigen::MatrixBase<Derived2> &B)
{
	const types::DynMat<typename Derived1::Scalar> & rA = A;
	const types::DynMat<typename Derived2::Scalar> & rB = B;

// check same scalar type
	if (typeid(typename Derived1::Scalar) != typeid(typename Derived2::Scalar))
		throw Exception("comm", Exception::Type::TYPE_MISMATCH);

// check zero-size
	if (!internal::_check_nonzero_size(rA)
			|| !internal::_check_nonzero_size(rB))
		throw Exception("comm", Exception::Type::ZERO_SIZE);

// check square matrices
	if (!internal::_check_square_mat(rA) || !internal::_check_square_mat(rB))
		throw Exception("comm", Exception::Type::MATRIX_NOT_SQUARE);

// check equal dimensions
	if (rA.rows() != rB.rows())
		throw Exception("comm", Exception::Type::DIMS_NOT_EQUAL);

	return rA * rB - rB * rA;
}

// anti-commutator
template<typename Derived1, typename Derived2>
types::DynMat<typename Derived1::Scalar> anticomm(
		const Eigen::MatrixBase<Derived1> &A,
		const Eigen::MatrixBase<Derived2> &B)
{
	const types::DynMat<typename Derived1::Scalar>& rA = A;
	const types::DynMat<typename Derived2::Scalar>& rB = B;

// check same scalar type
	if (typeid(typename Derived1::Scalar) != typeid(typename Derived2::Scalar))
		throw Exception("kron", Exception::Type::TYPE_MISMATCH);

// check zero-size
	if (!internal::_check_nonzero_size(rA)
			|| !internal::_check_nonzero_size(rB))
		throw Exception("anticomm", Exception::Type::ZERO_SIZE);

// check square matrices
	if (!internal::_check_square_mat(rA) || !internal::_check_square_mat(rB))
		throw Exception("anticomm", Exception::Type::MATRIX_NOT_SQUARE);

// check equal dimensions
	if (rA.rows() != rB.rows())
		throw Exception("anticomm", Exception::Type::DIMS_NOT_EQUAL);

	return rA * rB + rB * rA;
}

// projector onto |V><V| (normalized)
template<typename Derived>
types::DynMat<typename Derived::Scalar> prj(const Eigen::MatrixBase<Derived>& V)
{
	const types::DynMat<typename Derived::Scalar> & rV = V;

// check zero-size
	if (!internal::_check_nonzero_size(rV))
		throw Exception("prj", Exception::Type::ZERO_SIZE);

// check column vector
	if (!internal::_check_col_vector(rV))
		throw Exception("prj", Exception::Type::MATRIX_NOT_CVECTOR);

	return rV * adjoint(rV) / trace(rV * adjoint(rV));
}

// optimized, faster than kron(As...)
template<typename Derived>
types::DynMat<typename Derived::Scalar> expandout(
		const Eigen::MatrixBase<Derived>& A, size_t pos,
		const std::vector<size_t>& dims)
{
	const types::DynMat<typename Derived::Scalar> & rA = A;

// check zero-size
	if (!internal::_check_nonzero_size(rA))
		throw Exception("expandout", Exception::Type::ZERO_SIZE);

// check that dims is a valid dimension vector
	if (!internal::_check_dims(dims))
		throw Exception("expandout", Exception::Type::DIMS_INVALID);

// check square matrix
	if (!internal::_check_square_mat(rA))
		throw Exception("expandout", Exception::Type::MATRIX_NOT_SQUARE);

// check that position is valid
	if (pos > dims.size() - 1)
		throw Exception("expandout", Exception::Type::OUT_OF_RANGE);

// check that dims[pos] match the dimension of A
	if (static_cast<size_t>(rA.cols()) != dims[pos])
		throw Exception("expandout", Exception::Type::DIMS_MISMATCH_MATRIX);

	auto multiply = [](size_t x, size_t y)->size_t
	{	return x*y;};

	size_t D = std::accumulate(std::begin(dims), std::end(dims), 1u, multiply);
	types::DynMat<typename Derived::Scalar> result = types::DynMat<
			typename Derived::Scalar>::Identity(D, D);

	size_t Cdims[ct::maxn];
	size_t midx_row[ct::maxn];
	size_t midx_col[ct::maxn];

	for (size_t k = 0; k < dims.size(); k++)
	{
		midx_row[k] = midx_col[k] = 0;
		Cdims[k] = dims[k];
	}

// run over the main diagonal multi-indexes
	for (size_t i = 0; i < D; i++)
	{
		// get row multi_index
		internal::_n2multiidx(i, dims.size(), Cdims, midx_row);
		// get column multi_index (same as row)
		internal::_n2multiidx(i, dims.size(), Cdims, midx_col);
		// run over the gate's row multi-index
		for (size_t a = 0; a < static_cast<size_t>(rA.cols()); a++)
		{
			// construct the total row multi-index
			midx_row[pos] = a;

			// run over the gate's column multi-index
			for (size_t b = 0; b < static_cast<size_t>(rA.cols()); b++)
			{
				// construct the total column multi-index
				midx_col[pos] = b;

				// finally write the values
				result(internal::_multiidx2n(midx_row, dims.size(), Cdims),
						internal::_multiidx2n(midx_col, dims.size(), Cdims)) =
						rA(a, b);
			}
		}
	}

	return result;
}

// applies gate A to part of state vector
// or density matrix specified by subsys
template<typename Derived1, typename Derived2>
types::DynMat<typename Derived1::Scalar> gate(
		const Eigen::MatrixBase<Derived1>& state,
		const Eigen::MatrixBase<Derived2>& A, const std::vector<size_t>& subsys,
		const std::vector<size_t>& dims)
{
	const types::DynMat<typename Derived1::Scalar> & rstate = state;
	const types::DynMat<typename Derived2::Scalar> & rA = A;

	// EXCEPTION CHECKS
	// check same scalar type
	if (typeid(typename Derived1::Scalar) != typeid(typename Derived2::Scalar))
		throw Exception("gate", Exception::Type::TYPE_MISMATCH);

	// check zero sizes
	if (!internal::_check_nonzero_size(rA))
		throw Exception("gate", Exception::Type::ZERO_SIZE);

	// check zero sizes
	if (!internal::_check_nonzero_size(rstate))
		throw Exception("gate", Exception::Type::ZERO_SIZE);

	// check square matrix for the subsys
	if (!internal::_check_square_mat(rA))
		throw Exception("gate", Exception::Type::MATRIX_NOT_SQUARE);

	// check that dims is a valid dimension vector
	if (!internal::_check_dims(dims))
		throw Exception("gate", Exception::Type::DIMS_INVALID);

	// check subsys is valid w.r.t. dims
	if (!internal::_check_subsys_match_dims(subsys, dims))
		throw Exception("gate", Exception::Type::SUBSYS_MISMATCH_DIMS);

	// Use static allocation for speed!
	size_t Cdims[ct::maxn];
	size_t midx_row[ct::maxn];
	size_t midx_rho_row[ct::maxn];

	size_t CdimsA[ct::maxn];
	size_t CsubsysA[ct::maxn];
	size_t midxA_row[ct::maxn];
	size_t midxA_rho_row[ct::maxn];

	size_t CdimsA_bar[ct::maxn];
	size_t CsubsysA_bar[ct::maxn];
	size_t midxA_bar_row[ct::maxn];

	size_t n = dims.size();
	size_t nA = subsys.size();
	size_t nA_bar = n - nA;

	size_t D = 1;
	size_t DA_bar = 1;

	for (size_t k = 0, cnt = 0; k < n; k++)
	{
		midx_row[k] = midx_rho_row[k] = 0;
		Cdims[k] = dims[k];
		D *= dims[k];

		// compute the complement of subsys w.r.t. dims
		if (std::find(std::begin(subsys), std::end(subsys), k)
				== std::end(subsys))
		{
			CsubsysA_bar[cnt] = k;
			CdimsA_bar[cnt] = dims[k];
			midxA_bar_row[cnt] = 0;
			DA_bar *= dims[k];
			cnt++;
		}
	}

	size_t DA = 1;
	for (size_t k = 0; k < nA; k++)
	{
		midxA_row[k] = midxA_rho_row[k] = 0;
		CdimsA[k] = dims[subsys[k]];
		CsubsysA[k] = subsys[k];
		DA *= dims[subsys[k]];
	}

	if (internal::_check_col_vector(rstate)) // we have a ket
	{
		// check that dims match state vector
		if (!internal::_check_dims_match_cvect(dims, rstate))
			throw Exception("gate", Exception::Type::DIMS_MISMATCH_CVECTOR);

		// check that state vector matches the dimensions of the subsys
		if (static_cast<size_t>(rA.cols()) != DA)
			throw Exception("gate", Exception::Type::DIMS_MISMATCH_CVECTOR);

		types::DynMat<typename Derived1::Scalar> result(D, 1);

		// run over the subsys's row multi-index
		for (size_t a = 0; a < DA; a++)
		{
			// get the subsys's row multi-index
			internal::_n2multiidx(a, nA, CdimsA, midxA_row);
			// compute subsys part of the result's row multi-index
			for (size_t k = 0; k < nA; k++)
				midx_row[CsubsysA[k]] = midxA_row[k];

			// run over the complement's row multi-index
			for (size_t i = 0; i < DA_bar; i++)
			{
				// get the complement's row multi-index
				internal::_n2multiidx(i, nA_bar, CdimsA_bar, midxA_bar_row);
				// now compute the complement part of the result's row multi-index
				// and the complement part of the state's total row multi-index
				for (size_t k = 0; k < nA_bar; k++)
					midx_row[CsubsysA_bar[k]] = midx_rho_row[CsubsysA_bar[k]] =
							midxA_bar_row[k];
				// compute the results's row index
				size_t result_row_idx = internal::_multiidx2n(midx_row, n,
						Cdims);

				// compute the coefficient
				typename Derived1::Scalar coeff = 0;
				for (size_t c = 0; c < DA; c++)
				{
					// compute the subsys part state's row multi-index
					internal::_n2multiidx(c, nA, CdimsA, midxA_rho_row);
					// now we have the total state's row multi-index
					for (size_t k = 0; k < nA; k++)
						midx_rho_row[CsubsysA[k]] = midxA_rho_row[k];

					coeff += rA(a, c)
							* rstate(
									internal::_multiidx2n(midx_rho_row, n,
											Cdims));
				}
				// write down the result
				result(result_row_idx) = coeff;
			}
		}
		return result;
	}
	else if (internal::_check_square_mat(rstate)) // we have a matrix
	{

		// check that dims match state matrix
		if (!internal::_check_dims_match_mat(dims, rstate))
			throw Exception("gate", Exception::Type::DIMS_MISMATCH_MATRIX);

		// check that state matrix matches the dimensions of the subsys
		if (static_cast<size_t>(rA.cols()) != DA)
			throw Exception("gate", Exception::Type::DIMS_MISMATCH_MATRIX);

		types::DynMat<typename Derived1::Scalar> result(D, D);

		// run over the subsys's row multi-index
		for (size_t a = 0; a < DA; a++)
		{
			// get the subsys's row multi-index
			internal::_n2multiidx(a, nA, CdimsA, midxA_row);
			// compute subsys part of the result's row multi-index
			for (size_t k = 0; k < nA; k++)
				midx_row[CsubsysA[k]] = midxA_row[k];

			// run over the complement's row multi-index
			for (size_t i = 0; i < DA_bar; i++)
			{
				// get the complement's row multi-index
				internal::_n2multiidx(i, nA_bar, CdimsA_bar, midxA_bar_row);
				// now compute the complement part of the result's row multi-index
				// and the complement part of the state's total row multi-index
				for (size_t k = 0; k < nA_bar; k++)
					midx_row[CsubsysA_bar[k]] = midx_rho_row[CsubsysA_bar[k]] =
							midxA_bar_row[k];
				// compute the results's row index
				size_t result_row_idx = internal::_multiidx2n(midx_row, n,
						Cdims);

				// run over the col index
				for (size_t j = 0; j < D; j++)
				{
					// compute the coefficient
					typename Derived1::Scalar coeff = 0;
					for (size_t c = 0; c < DA; c++)
					{
						// compute the subsys part state's row multi-index
						internal::_n2multiidx(c, nA, CdimsA, midxA_rho_row);
						// now we have the total state's row multi-index
						for (size_t k = 0; k < nA; k++)
							midx_rho_row[CsubsysA[k]] = midxA_rho_row[k];

						coeff += rA(a, c)
								* rstate(
										internal::_multiidx2n(midx_rho_row, n,
												Cdims), j);

					}
					// write down the result
					result(result_row_idx, j) = coeff;
				}
			}
		}
		return result;
	}
	else
		throw Exception("gate", Exception::Type::MATRIX_NOT_SQUARE_OR_CVECTOR);
}

// Gram-Schmidt ortogonalization
template<typename Derived>
types::DynMat<typename Derived::Scalar> grams(const std::vector<Derived>& Vs)
{
// check empty list
	if (!internal::_check_nonzero_size(Vs))
		throw Exception("grams", Exception::Type::ZERO_SIZE);

	for (auto it : Vs)
		if (!internal::_check_nonzero_size(it))
			throw Exception("grams", Exception::Type::ZERO_SIZE);

// check that Vs[0] is a column vector
	if (!internal::_check_col_vector(Vs[0]))
		throw Exception("grams", Exception::Type::MATRIX_NOT_CVECTOR);

// now check that all the rest match the size of the first vector
	for (auto it : Vs)
		if (it.rows() != Vs[0].rows() || it.cols() != 1)
			throw Exception("grams", Exception::Type::DIMS_NOT_EQUAL);

// start the process
	std::vector<types::DynMat<typename Derived::Scalar>> outvecs;
	outvecs.push_back(Vs[0] / norm(Vs[0]));

	types::DynMat<typename Derived::Scalar> cut = types::DynMat<
			typename Derived::Scalar>::Identity(Vs[0].rows(), Vs[0].rows());

	types::DynMat<typename Derived::Scalar> vi = types::DynMat<
			typename Derived::Scalar>::Zero(Vs[0].rows(), 1);

	for (size_t i = 1; i < Vs.size(); i++)
	{
		cut -= prj(outvecs[i - 1]);
		vi = cut * Vs[i];

		if (abs(norm(vi)) > ct::eps) // only adds the non-zero vectors
			outvecs.push_back(vi / norm(vi));
	}

	types::DynMat<typename Derived::Scalar> result(Vs[0].rows(),
			outvecs.size());
	size_t cnt = 0;
	for (auto it : outvecs)
	{
		result.col(cnt) = it;
		cnt++;
	}
	return result;
}

// Gram-Schmidt ortogonalization
// deduce the template parameters from initializer_list
template<typename Derived>
types::DynMat<typename Derived::Scalar> grams(
		const std::initializer_list<Derived>& Vs)
{
	return grams(std::vector<Derived>(Vs));
}

// Gram-Schmidt ortogonalization
template<typename Derived>
types::DynMat<typename Derived::Scalar> grams(
		const Eigen::MatrixBase<Derived>& A)
{
	const types::DynMat<typename Derived::Scalar> & rA = A;

	if (!internal::_check_nonzero_size(rA))
		throw Exception("grams", Exception::Type::ZERO_SIZE);

	std::vector<types::DynMat<typename Derived::Scalar>> input;

	for (size_t i = 0; i < static_cast<size_t>(rA.cols()); i++)
		input.push_back(
				static_cast<types::DynMat<typename Derived::Scalar> >(rA.col(i)));

	return grams<types::DynMat<typename Derived::Scalar>>(input);
}

// integer index to multi-index
// standard lexicographical order, e.g. 00, 01, 10, 11
std::vector<size_t> n2multiidx(size_t n, const std::vector<size_t>& dims)
{
	if (!internal::_check_dims(dims))
		throw Exception("n2multiidx", Exception::Type::DIMS_INVALID);

	auto multiply = [](const size_t x, const size_t y)->size_t
	{	return x*y;};

	if (n >= std::accumulate(std::begin(dims), std::end(dims), 1u, multiply))
		throw Exception("n2multiidx", Exception::Type::OUT_OF_RANGE);

	std::vector<size_t> result(dims.size());
	size_t _n = n;
	for (size_t i = 0; i < dims.size(); i++)
	{
		result[dims.size() - i - 1] = _n
				% static_cast<int>(dims[dims.size() - i - 1]);
		_n = _n / static_cast<int>(dims[dims.size() - i - 1]);
	}

	return result;
}

// multi-index to integer index
// standard lexicographical order, e.g. 00->0, 01->1, 10->2, 11->3
size_t multiidx2n(const std::vector<size_t>& midx,
		const std::vector<size_t>& dims)
{
	if (!internal::_check_dims(dims))
		throw Exception("multiidx2n", Exception::Type::DIMS_INVALID);

	for (size_t i = 0; i < dims.size(); i++)
		if (midx[i] >= dims[i])
			throw Exception("multiidx2n", Exception::Type::OUT_OF_RANGE);

	std::vector<size_t> part_prod(dims.size());

	part_prod[dims.size() - 1] = 1;
	for (size_t j = 1; j < dims.size(); j++)
		part_prod[dims.size() - j - 1] = part_prod[dims.size() - j]
				* dims[dims.size() - j];

	size_t result = 0;
	for (size_t i = 0; i < dims.size(); i++)
		result += midx[i] * part_prod[i];

	return result;
}

// constructs a multi-qubit ket
types::ket mket(const std::vector<size_t>& mask)
{
	size_t n = mask.size();
	size_t D = std::pow(2, n);
// check zero size
	if (n == 0)
		throw Exception("mket", Exception::Type::ZERO_SIZE);
// check mask is a valid binary vector
	for (auto it : mask)
		if (it > 1)
			throw Exception("mket", Exception::Type::NOT_QUBIT_SUBSYS);
	std::vector<size_t> dims(n, 2);
	types::ket result = types::ket::Zero(D);
	size_t pos = multiidx2n(mask, dims);
	result(pos) = 1;
	return result;
}

// constructs a multi-qudit ket
types::ket mket(const std::vector<size_t>& mask,
		const std::vector<size_t>& dims)
{
	size_t n = mask.size();
	auto multiply = [](size_t x, size_t y)->size_t
	{	return x*y;};

	size_t D = std::accumulate(std::begin(dims), std::end(dims), 1u, multiply);

// check zero size
	if (n == 0)
		throw Exception("mket", Exception::Type::ZERO_SIZE);
// check valid dims
	if (!internal::_check_dims(dims))
		throw Exception("mket", Exception::Type::DIMS_INVALID);
// check mask and dims have the same size
	if (mask.size() != dims.size())
		throw Exception("mket", Exception::Type::SUBSYS_MISMATCH_DIMS);
// check mask is a valid vector
	for (size_t i = 0; i < n; i++)
		if (mask[i] >= dims[i])
			throw Exception("mket", Exception::Type::SUBSYS_MISMATCH_DIMS);

	types::ket result = types::ket::Zero(D);
	size_t pos = multiidx2n(mask, dims);
	result(pos) = 1;
	return result;
}

// constructs a multi-qudit ket where all subsystems have equal size d
types::ket mket(const std::vector<size_t>& mask, size_t d)
{
	size_t n = mask.size();
	size_t D = std::pow(d, n);

// check zero size
	if (n == 0)
		throw Exception("mket", Exception::Type::ZERO_SIZE);
// check valid dims
	if (d == 0)
		throw Exception("mket", Exception::Type::DIMS_INVALID);
// check mask is a valid vector
	for (size_t i = 0; i < n; i++)
		if (mask[i] >= d)
			throw Exception("mket", Exception::Type::SUBSYS_MISMATCH_DIMS);

	types::ket result = types::ket::Zero(D);
	std::vector<size_t> dims(n, d);
	size_t pos = multiidx2n(mask, dims);
	result(pos) = 1;
	return result;
}

// inverse permutation
std::vector<size_t> invperm(const std::vector<size_t>& perm)
{
	if (!internal::_check_perm(perm))
		throw Exception("invperm", Exception::Type::PERM_INVALID);

	// construct the inverse
	std::vector<size_t> result(perm.size());
	for (size_t i = 0; i < perm.size(); i++)
		result[perm[i]] = i;

	return result;
}

// compose permutations, perm(sigma)
std::vector<size_t> compperm(const std::vector<size_t>& perm,
		const std::vector<size_t>& sigma)
{
	if (!internal::_check_perm(perm))
		throw Exception("compperm", Exception::Type::PERM_INVALID);
	if (!internal::_check_perm(sigma))
		throw Exception("compperm", Exception::Type::PERM_INVALID);
	if (perm.size() != sigma.size())
		throw Exception("compperm", Exception::Type::PERM_INVALID);

	// construct the composition perm(sigma)
	std::vector<size_t> result(perm.size());
	for (size_t i = 0; i < perm.size(); i++)
		result[i] = perm[sigma[i]];

	return result;
}

} /* namespace qpp */

#endif /* FUNCTIONS_H_ */
