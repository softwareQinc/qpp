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

// trace-norm (or Frobenius norm) (CHANGES return type to double)
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

// eigenvalues (CHANGES return type to complex)
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

// eigenvalues of Hermitian matrices (CHANGES return type to complex)
template<typename Derived>
types::cmat hevals(const Eigen::MatrixBase<Derived>& A)
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
	return es.eigenvalues().template cast<types::cplx>();
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

	return funm(rA, std::sqrt);
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

	return funm(rA, std::exp);
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

	return funm(rA, std::log);
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

	return funm(rA, std::sin);
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

	return funm(rA, std::cos);
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
types::DynMat<OutputScalar> fun(const Eigen::MatrixBase<Derived> &A,
		OutputScalar (*f)(const typename Derived::Scalar &))
{
	const types::DynMat<typename Derived::Scalar> & rA = A;

	// check zero-size
	if (!internal::_check_nonzero_size(rA))
		throw Exception("fun", Exception::Type::ZERO_SIZE);

	types::DynMat<OutputScalar> result(rA.rows(), rA.cols());

	for (size_t j = 0; j < static_cast<size_t>(rA.cols()); j++)
#pragma omp parallel for
		for (size_t i = 0; i < static_cast<size_t>(rA.rows()); i++)
			result(i, j) = (*f)(rA(i, j));

	return result;
}

// Kronecker product of 2 matrices, preserve return type
template<typename Derived1, typename Derived2>
types::DynMat<typename Derived1::Scalar> kron(
		const Eigen::MatrixBase<Derived1>& A,
		const Eigen::MatrixBase<Derived2>& B)
{
	const types::DynMat<typename Derived1::Scalar> & rA = A;
	const types::DynMat<typename Derived2::Scalar> & rB = B;

	// check same scalar type
	if (typeid(typename Derived1::Scalar) != typeid(typename Derived2::Scalar))
		throw Exception("kron", Exception::Type::TYPE_MISMATCH);

	// check zero-size
	if (!internal::_check_nonzero_size(rA))
		throw Exception("kron", Exception::Type::ZERO_SIZE);

	// check zero-size
	if (!internal::_check_nonzero_size(rB))
		throw Exception("kron", Exception::Type::ZERO_SIZE);

	size_t Acols = static_cast<size_t>(rA.cols());
	size_t Arows = static_cast<size_t>(rA.rows());
	size_t Bcols = static_cast<size_t>(rB.cols());
	size_t Brows = static_cast<size_t>(rB.rows());

	types::DynMat<typename Derived1::Scalar> result;
	result.resize(Arows * Brows, Acols * Bcols);

	for (size_t j = 0; j < Acols; j++)
#pragma omp parallel for
		for (size_t i = 0; i < Arows; i++)
			result.block(i * Brows, j * Bcols, Brows, Bcols) = rA(i, j) * rB;

	return result;

}

// Kronecker product of a list of matrices, preserve return type
// works only for list of matrices, not list of expressions
template<typename Derived>
types::DynMat<typename Derived::Scalar> kronlist(
		const std::vector<types::DynMat<typename Derived::Scalar> > &As)

{
	if (As.size() == 0)
		throw Exception("kronlist", Exception::Type::ZERO_SIZE);

	for (auto it : As)
		if (it.size() == 0)
			throw Exception("kronlist", Exception::Type::ZERO_SIZE);

	types::DynMat<typename Derived::Scalar> result = As[0];
	for (size_t i = 1; i < As.size(); i++)
	{
		result = kron(result, As[i]);
	}
	return result;
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

// check that the size of the permutation is OK
	if (!internal::_check_perm(perm, dims))
		throw Exception("syspermute", Exception::Type::PERM_MISMATCH_DIMS);

	size_t dim = static_cast<size_t>(rA.rows());
	size_t numdims = dims.size();

	types::DynMat<typename Derived::Scalar> result;

	// check column vector
	if (internal::_check_col_vector(rA)) // we have a column vector
	{
		// check that dims match the dimension of rA
		if (!internal::_check_dims_match_cvect(dims, rA))
			throw Exception("syspermute",
					Exception::Type::DIMS_MISMATCH_CVECTOR);

		size_t* cdims = new size_t[numdims];
		size_t* cperm = new size_t[numdims];
		size_t* midx = new size_t[numdims];

		// copy dims in cdims and perm in cperm
		for (size_t i = 0; i < numdims; i++)
		{
			cdims[i] = dims[i];
			cperm[i] = perm[i];
		}
		result.resize(dim, 1);
		size_t iperm = 0;
#pragma omp parallel for
		for (size_t i = 0; i < dim; i++)
			internal::_syspermute_worker(numdims, cdims, cperm, i, iperm, rA,
					result);
		delete[] cdims;
		delete[] cperm;
		delete[] midx;
		return result;
	}

	else if (internal::_check_square_mat(rA)) // we have a square matrix
	{
		// check that dims match the dimension of rA
		if (!internal::_check_dims_match_mat(dims, rA))
			throw Exception("syspermute",
					Exception::Type::DIMS_MISMATCH_MATRIX);

		size_t* cdims = new size_t[2 * numdims];
		size_t* cperm = new size_t[2 * numdims];
		size_t* midx = new size_t[2 * numdims];

		// copy dims in cdims and perm in cperm
		for (size_t i = 0; i < numdims; i++)
		{
			cdims[i] = cdims[i + numdims] = dims[i];
			cperm[i] = perm[i];
			cperm[i + numdims] = perm[i] + numdims;
		}
		result.resize(dim * dim, 1);
		// map A to a column vector
		types::DynMat<typename Derived::Scalar> vectA = Eigen::Map<
				types::DynMat<typename Derived::Scalar>>(
				const_cast<typename Derived::Scalar*>(rA.data()), dim * dim, 1);
		size_t iperm = 0;
#pragma omp parallel for
		for (size_t i = 0; i < dim * dim; i++)
			internal::_syspermute_worker(2 * numdims, cdims, cperm, i, iperm,
					vectA, result);
		delete[] cdims;
		delete[] cperm;
		delete[] midx;
		return reshape(result, dim, dim);
	}

	else
		throw Exception("syspermute",
				Exception::Type::MATRIX_NOT_SQUARE_OR_CVECTOR);
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

// check that subsys are valid
	if (!internal::_check_subsys(subsys, dims))
		throw Exception("ptrace", Exception::Type::SUBSYS_MISMATCH_DIMS);

	size_t dim = static_cast<size_t>(rA.rows());
	size_t numsubsys = subsys.size(); // number of subsystems we trace out
	size_t numdims = dims.size(); // total number of subsystems;
	std::vector<size_t> perm(numdims, 0); // the permutation vector
	std::vector<size_t> permdims; // the permuted dimensions

	types::DynMat<typename Derived::Scalar> result;

// the total dimension of the traced-out subsystems
	size_t dimsubsys = 1;
	for (size_t i = 0; i < numsubsys; i++)
		dimsubsys *= dims[subsys[i]];

	std::vector<size_t> sizeAB;
	sizeAB.push_back(dim / dimsubsys);
	sizeAB.push_back(dimsubsys);

// construct the permutation that bring the traced-out subsystems to the end
	size_t cnt0 = 0;
	size_t cnt1 = 0;
	for (size_t i = 0; i < numdims; i++)
	{
		// we find that i belongs to the subsystem
		if (std::find(std::begin(subsys), std::end(subsys), i)
				!= std::end(subsys))
		{
			perm[numdims - numsubsys + cnt0] = i;
			cnt0++;
		}
		else
		{
			perm[cnt1] = i;
			cnt1++;
		}
	}

	return ptrace2(syspermute(rA, perm, dims), sizeAB);
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

	// check that subsys are valid
	if (!internal::_check_subsys(subsys, dims))
		throw Exception("ptranspose", Exception::Type::SUBSYS_MISMATCH_DIMS);

	size_t dim = static_cast<size_t>(rA.rows());
	size_t numdims = dims.size();
	size_t numsubsys = subsys.size();
	size_t* cdims = new size_t[numdims];
	size_t* midxcol = new size_t[numdims];
	size_t* csubsys = new size_t[numsubsys];

	types::DynMat<typename Derived::Scalar> result = rA;

// copy dims in cdims and subsys in csubsys
	for (size_t i = 0; i < numdims; i++)
		cdims[i] = dims[i];
	for (size_t i = 0; i < numsubsys; i++)
		csubsys[i] = subsys[i];

	size_t iperm = 0;
	size_t jperm = 0;
	for (size_t j = 0; j < dim; j++)
	{
		// compute the column multi-index
		internal::_n2multiidx(j, numdims, cdims, midxcol);
#pragma omp parallel for
		for (size_t i = 0; i < dim; i++)
			internal::_ptranspose_worker(midxcol, numdims, numsubsys, cdims,
					csubsys, i, j, iperm, jperm, rA, result);
	}

	delete[] midxcol;
	delete[] cdims;
	delete[] csubsys;

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
types::DynMat<typename Derived::Scalar> proj(
		const Eigen::MatrixBase<Derived>& V)
{
	const types::DynMat<typename Derived::Scalar> & rV = V;

	// check zero-size
	if (!internal::_check_nonzero_size(rV))
		throw Exception("proj", Exception::Type::ZERO_SIZE);

	// check column vector
	if (!internal::_check_col_vector(rV))
		throw Exception("proj", Exception::Type::MATRIX_NOT_CVECTOR);

	return rV * adjoint(rV) / trace(rV * adjoint(rV));
}

// optimized, faster than kronlist(I, A, I, ...)
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

	size_t* Cdims = new size_t[dims.size()];
	size_t* midx_row = new size_t[dims.size()];
	size_t* midx_col = new size_t[dims.size()];

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

	delete[] Cdims;
	delete[] midx_row;
	delete[] midx_col;

	return result;
}

// Gram-Schmidt ortogonalization
template<typename Derived>
types::DynMat<typename Derived::Scalar> grams(
		const std::vector<types::DynMat<typename Derived::Scalar> >& Vs)
{
	// check empty list
	if (!internal::_check_nonzero_size(Vs))
		throw Exception("grams", Exception::Type::ZERO_SIZE);

	for (auto it : Vs)
		if (!internal::_check_nonzero_size(it))
			throw Exception("grams", Exception::Type::ZERO_SIZE);

	// check that Vs[0] is indeed a column vector
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
		cut -= proj(outvecs[i - 1]);
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
template<typename Derived>
types::DynMat<typename Derived::Scalar> grams(
		const Eigen::MatrixBase<Derived>& A)
{
	const types::DynMat<typename Derived::Scalar> & rA = A;

	if (!internal::_check_nonzero_size(rA))
		throw Exception("grams", Exception::Type::ZERO_SIZE);

	std::vector<types::DynMat<typename Derived::Scalar>> input;

	for (size_t i = 0; i < static_cast<size_t>(rA.cols()); i++)
	input.push_back(static_cast<types::DynMat<typename Derived::Scalar>>(rA.col(i)));

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

} /* namespace qpp */

#endif /* FUNCTIONS_H_ */
