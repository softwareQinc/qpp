/*
 * functional.h
 *
 *  Created on: Mar 27, 2014
 *      Author: vlad
 */

#ifndef FUNCTIONAL_H_
#define FUNCTIONAL_H_

#include <vector>
#include <cmath>
#include "types.h"
#include "internal.h"
#include "exception.h"
#include "constants.h"
#include "stat.h"

// Collection of quantum computing useful functions

namespace qpp
{
// Eigen function wrappers

// transpose, preserve return type
template<typename Scalar>
types::DynMat<Scalar> transpose(const types::DynMat<Scalar>& A)

{
	// check zero-size
	if (!internal::_check_nonzero_size(A))
		throw Exception("transpose", Exception::Type::MATRIX_ZERO_SIZE);

	return A.transpose();
}

// conjugate, preserve return type
template<typename Scalar>
types::DynMat<Scalar> conjugate(const types::DynMat<Scalar>& A)

{
	// check zero-size
	if (!internal::_check_nonzero_size(A))
		throw Exception("conjugate", Exception::Type::MATRIX_ZERO_SIZE);

	return A.conjugate();
}

// adjoint, preserve return type
template<typename Scalar>
types::DynMat<Scalar> adjoint(const types::DynMat<Scalar>& A)
{
	// check zero-size
	if (!internal::_check_nonzero_size(A))
		throw Exception("adjoint", Exception::Type::MATRIX_ZERO_SIZE);

	return (A).adjoint();
}

// trace, preserve return type
template<typename Scalar>
Scalar trace(const types::DynMat<Scalar>& A)
{
	// check zero-size
	if (!internal::_check_nonzero_size(A))
		throw Exception("trace", Exception::Type::MATRIX_ZERO_SIZE);

	return A.trace();
}

// element-wise sum, preserve return type
template<typename Scalar>
Scalar sum(const types::DynMat<Scalar>& A)
{
	// check zero-size
	if (!internal::_check_nonzero_size(A))
		throw Exception("trace", Exception::Type::MATRIX_ZERO_SIZE);

	return A.sum();
}

// trace-norm (or Frobenius norm) (CHANGES return type to double)
template<typename Scalar>
double norm(const types::DynMat<Scalar>& A)
{
	// check zero-size
	if (!internal::_check_nonzero_size(A))
		throw Exception("norm", Exception::Type::MATRIX_ZERO_SIZE);

	// convert matrix to complex then return its norm
	return (A.template cast<types::cplx>()).norm();
}

// eigenvalues (CHANGES return type to complex)
template<typename Scalar>
types::cmat evals(const types::DynMat<Scalar>& A)
{
	// check zero-size
	if (!internal::_check_nonzero_size(A))
		throw Exception("evals", Exception::Type::MATRIX_ZERO_SIZE);

	// check square matrix
	if (!internal::_check_square_mat(A))
		throw Exception("evals", Exception::Type::MATRIX_NOT_SQUARE);

	Eigen::ComplexEigenSolver<types::cmat> es(A.template cast<types::cplx>());
	return es.eigenvalues();
}

// eigenvectors (CHANGES return type to complex matrix)
template<typename Scalar>
types::cmat evects(const types::DynMat<Scalar>& A)
{
	// check zero-size
	if (!internal::_check_nonzero_size(A))
		throw Exception("evects", Exception::Type::MATRIX_ZERO_SIZE);

	// check square matrix
	if (!internal::_check_square_mat(A))
		throw Exception("evects", Exception::Type::MATRIX_NOT_SQUARE);

	Eigen::ComplexEigenSolver<types::cmat> es(A.template cast<types::cplx>());
	return es.eigenvectors();
}

// eigenvalues of Hermitian matrices (CHANGES return type to complex)
template<typename Scalar>
types::cmat hevals(const types::DynMat<Scalar>& A)
{
	// check zero-size
	if (!internal::_check_nonzero_size(A))
		throw Exception("hevals", Exception::Type::MATRIX_ZERO_SIZE);

	// check square matrix
	if (!internal::_check_square_mat(A))
		throw Exception("hevals", Exception::Type::MATRIX_NOT_SQUARE);

	Eigen::SelfAdjointEigenSolver<types::cmat> es(
			A.template cast<types::cplx>());
	return es.eigenvalues().template cast<types::cplx>();
}

// eigenvectors of Hermitian matrix (CHANGES return type to complex matrix)
template<typename Scalar>
types::cmat hevects(const types::DynMat<Scalar>& A)
{
	// check zero-size
	if (!internal::_check_nonzero_size(A))
		throw Exception("hevects", Exception::Type::MATRIX_ZERO_SIZE);

	// check square matrix
	if (!internal::_check_square_mat(A))
		throw Exception("hevects", Exception::Type::MATRIX_NOT_SQUARE);

	Eigen::SelfAdjointEigenSolver<types::cmat> es(
			A.template cast<types::cplx>());
	return es.eigenvectors();
}

// Matrix functional calculus

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

// other functions

// functor; apply f(A) component-wise, where (*f) is the function pointer
// returns a matrix of type OutputScalar
template<typename InputScalar, typename OutputScalar>
types::DynMat<OutputScalar> fun(const types::DynMat<InputScalar> &A,
		OutputScalar (*f)(const InputScalar &))
{
	// check zero-size
	if (!internal::_check_nonzero_size(A))
		throw Exception("fun", Exception::Type::MATRIX_ZERO_SIZE);

	types::DynMat<OutputScalar> result(A.rows(), A.cols());

	for (size_t j = 0; j < static_cast<size_t>(A.cols()); j++)
#pragma omp parallel for
		for (size_t i = 0; i < static_cast<size_t>(A.rows()); i++)
			result(i, j) = (*f)(A(i, j));

	return result;
}

// Kronecker product of 2 matrices, preserve return type
template<typename Scalar>
types::DynMat<Scalar> kron(const types::DynMat<Scalar> &A,
		const types::DynMat<Scalar> &B)
{
	// check zero-size
	if (!internal::_check_nonzero_size(A))
		throw Exception("kron", Exception::Type::MATRIX_ZERO_SIZE);

	// check zero-size
	if (!internal::_check_nonzero_size(B))
		throw Exception("kron", Exception::Type::MATRIX_ZERO_SIZE);

	size_t Acols = static_cast<size_t>(A.cols());
	size_t Arows = static_cast<size_t>(A.rows());
	size_t Bcols = static_cast<size_t>(B.cols());
	size_t Brows = static_cast<size_t>(B.rows());

	types::DynMat<Scalar> result;
	result.resize(Arows * Brows, Acols * Bcols);

	for (size_t j = 0; j < Acols; j++)
#pragma omp parallel for
		for (size_t i = 0; i < Arows; i++)
			result.block(i * Brows, j * Bcols, Brows, Bcols) = A(i, j) * B;

	return result;

}

// Kronecker product of a list of matrices, preserve return type
// <Expression> is forced to be a matrix by invocation of kron
// inside the function
template<typename Scalar>
types::DynMat<Scalar> kron_list(const std::vector<types::DynMat<Scalar>> &list)

{
	for (auto i : list)
		if (i.size() == 0)
			throw Exception("kron_list", Exception::Type::MATRIX_ZERO_SIZE);

	types::DynMat<Scalar> result = list[0];
	for (size_t i = 1; i < list.size(); i++)
		result = kron(result, list[i]);
	return result;
}
// Kronecker product of a matrix with itself $n$ times, preserve return type
template<typename Scalar>
types::DynMat<Scalar> kron_pow(const types::DynMat<Scalar> &A, size_t n)

{
	// check zero-size
	if (!internal::_check_nonzero_size(A))
		throw Exception("kron_pow", Exception::Type::MATRIX_ZERO_SIZE);

	std::vector<typename types::DynMat<Scalar>> list;
	for (size_t i = 0; i < n; i++)
		list.push_back(A);
	return kron_list(list);
}

// reshape the columns of A and returns a matrix with m rows and n columns
// use column-major order (same as MATLAB)
template<typename Scalar>
types::DynMat<Scalar> reshape(const types::DynMat<Scalar>& A, size_t rows,
		size_t cols)
{
	size_t Arows = static_cast<size_t>(A.rows());
	size_t Acols = static_cast<size_t>(A.cols());

	// check zero-size
	if (!internal::_check_nonzero_size(A))
		throw Exception("reshape", Exception::Type::MATRIX_ZERO_SIZE);

	if (Arows * Acols != rows * cols)
		throw Exception("reshape", Exception::Type::DIMS_MISMATCH_MATRIX);

	return Eigen::Map<types::DynMat<Scalar>>(A.data(), rows, cols);
}

// permutes the subsystems in a matrix
template<typename Scalar>
types::DynMat<Scalar> syspermute(const types::DynMat<Scalar> &A,
		const std::vector<size_t> perm, const std::vector<size_t> &dims)

{
// Error checks

// check zero-size
	if (!internal::_check_nonzero_size(A))
		throw Exception("syspermute", Exception::Type::MATRIX_ZERO_SIZE);

	// check square matrix
	if (!internal::_check_square_mat(A))
		throw Exception("syspermute", Exception::Type::MATRIX_NOT_SQUARE);

// check that dims is a valid dimension vector
	if (!internal::_check_dims(dims))
		throw Exception("syspermute", Exception::Type::DIMS_HAVE_ZERO);

// check that dims match the dimension of A
	if (!internal::_check_dims_match_mat(dims, A))
		throw Exception("syspermute", Exception::Type::DIMS_MISMATCH_MATRIX);

// check that the size of the permutation is OK
	if (!internal::_check_perm(perm, dims))
		throw Exception("syspermute", Exception::Type::PERM_MISMATCH_DIMS);

	size_t dim = static_cast<size_t>(A.rows());
	size_t numdims = dims.size();
	size_t* cdims = new size_t[numdims];
	size_t* cperm = new size_t[numdims];
	size_t* midxcol = new size_t[numdims];
	types::DynMat<Scalar> result(dim, dim);

// copy dims in cdims and perm in cperm
	for (size_t i = 0; i < numdims; i++)
	{
		cdims[i] = dims[i];
		cperm[i] = perm[i];
	}

	size_t iperm = 0;
	size_t jperm = 0;
	for (size_t j = 0; j < dim; j++)
	{
		// compute the col multi-index
		internal::_n2multiidx(j, numdims, cdims, midxcol);
#pragma omp parallel for
		for (size_t i = 0; i < dim; i++)
			internal::_syspermute_worker(midxcol, numdims, cdims, cperm, i, j,
					iperm, jperm, A, result);
	}

	delete[] cdims;
	delete[] cperm;
	delete[] midxcol;

	return result; // the permuted matrix
}

// Partial trace over subsystem B in a D_A x D_B system
template<typename Scalar>
types::DynMat<Scalar> ptrace2(const types::DynMat<Scalar> &A,
		const std::vector<size_t> dims)
{
// Error checks

	// check zero-size
	if (!internal::_check_nonzero_size(A))
		throw Exception("ptrace2", Exception::Type::MATRIX_ZERO_SIZE);

	// check square matrix
	if (!internal::_check_square_mat(A))
		throw Exception("ptrace2", Exception::Type::MATRIX_NOT_SQUARE);

// check dims has only 2 elements
	if (dims.size() != 2)
		throw Exception("ptrace2", Exception::Type::DIMS_MISMATCH_MATRIX);

// check that dim is a valid dimension vector
	if (!internal::_check_dims(dims))
		throw Exception("ptrace2", Exception::Type::DIMS_HAVE_ZERO);

// check that dims match the dimension of A
	if (!internal::_check_dims_match_mat(dims, A))
		throw Exception("ptrace2", Exception::Type::DIMS_MISMATCH_MATRIX);

	size_t DA = dims[0];
	size_t DB = dims[1];

	types::DynMat<Scalar> result = types::DynMat<Scalar>::Zero(DA, DA);

	for (size_t j = 0; j < DA; j++) // column major order for speed
#pragma omp parallel for
		for (size_t i = 0; i < DA; i++)
		{
			result(i, j) = trace<Scalar>(A.block(i * DB, j * DB, DB, DB));
		}
	return result;
}

// partial trace
template<typename Scalar>
types::DynMat<Scalar> ptrace(const types::DynMat<Scalar> &A,
		const std::vector<size_t> &subsys, const std::vector<size_t> &dims)

{
// error checks

	// check zero-size
	if (!internal::_check_nonzero_size(A))
		throw Exception("ptrace", Exception::Type::MATRIX_ZERO_SIZE);

	// check square matrix
	if (!internal::_check_square_mat(A))
		throw Exception("ptrace", Exception::Type::MATRIX_NOT_SQUARE);

// check that dims is a valid dimension vector
	if (!internal::_check_dims(dims))
		throw Exception("ptrace", Exception::Type::DIMS_HAVE_ZERO);

// check that dims match the dimension of A
	if (!internal::_check_dims_match_mat(dims, A))
		throw Exception("ptrace", Exception::Type::DIMS_MISMATCH_MATRIX);

// check that subsys are valid
	if (!internal::_check_subsys(subsys, dims))
		throw Exception("ptrace", Exception::Type::SUBSYS_MISMATCH_DIMS);

	size_t dim = static_cast<size_t>(A.rows());
	size_t numsubsys = subsys.size(); // number of subsystems we trace out
	size_t numdims = dims.size(); // total number of subsystems;
	std::vector<size_t> perm(numdims, 0); // the permutation vector
	std::vector<size_t> permdims; // the permuted dimensions

	types::DynMat<Scalar> result;

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

	return ptrace2(syspermute(A, perm, dims), sizeAB);
}

// partial transpose
template<typename Scalar>
types::DynMat<Scalar> ptranspose(const types::DynMat<Scalar>& A,
		const std::vector<size_t>& subsys, const std::vector<size_t>& dims)

{
// error checks

	// check zero-size
	if (!internal::_check_nonzero_size(A))
		throw Exception("ptranspose", Exception::Type::MATRIX_ZERO_SIZE);

	// check square matrix
	if (!internal::_check_square_mat(A))
		throw Exception("ptranspose", Exception::Type::MATRIX_NOT_SQUARE);

// check that dims is a valid dimension vector
	if (!internal::_check_dims(dims))
		throw Exception("ptranspose", Exception::Type::DIMS_HAVE_ZERO);

// check that dims match the dimension of A
	if (!internal::_check_dims_match_mat(dims, A))
		throw Exception("ptranspose", Exception::Type::DIMS_MISMATCH_MATRIX);

	// check that subsys are valid
	if (!internal::_check_subsys(subsys, dims))
		throw Exception("ptranspose", Exception::Type::SUBSYS_MISMATCH_DIMS);

	size_t dim = static_cast<size_t>(A.rows());
	size_t numdims = dims.size();
	size_t numsubsys = subsys.size();
	size_t* cdims = new size_t[numdims];
	size_t* midxcol = new size_t[numdims];
	size_t* csubsys = new size_t[numsubsys];

	types::DynMat<Scalar> result = A;

// copy dims in cdims and subsys in csubsys
	for (size_t i = 0; i < numdims; i++)
		cdims[i] = dims[i];
	for (size_t i = 0; i < numsubsys; i++)
		csubsys[i] = subsys[i];

	size_t iperm = 0;
	size_t jperm = 0;
	for (size_t j = 0; j < dim; j++)
	{
		// compute the col multi-index
		internal::_n2multiidx(j, numdims, cdims, midxcol);
#pragma omp parallel for
		for (size_t i = 0; i < dim; i++)
			internal::_ptranspose_worker(midxcol, numdims, numsubsys, cdims,
					csubsys, i, j, iperm, jperm, A, result);
	}

	delete[] midxcol;
	delete[] cdims;
	delete[] csubsys;

	return result;
}

// commutator
template<typename Scalar>
types::DynMat<Scalar> comm(const types::DynMat<Scalar> &A,
		const types::DynMat<Scalar> &B)
{
	// check zero-size
	if (!internal::_check_nonzero_size(A) || !internal::_check_nonzero_size(A))
		throw Exception("comm", Exception::Type::MATRIX_ZERO_SIZE);

	// check square matrices
	if (!internal::_check_square_mat(A) || !internal::_check_square_mat(B))
		throw Exception("comm", Exception::Type::MATRIX_NOT_SQUARE);

	// check equal dimensions
	if (A.rows() != B.rows())
		throw Exception("comm", Exception::Type::DIMS_NOT_EQUAL);

	return A * B - B * A;
}

// anti-commutator of 2 matrices
template<typename Scalar>
types::DynMat<Scalar> anticomm(const types::DynMat<Scalar> &A,
		const types::DynMat<Scalar> &B)
{
	// check zero-size
	if (!internal::_check_nonzero_size(A) || !internal::_check_nonzero_size(A))
		throw Exception("anticomm", Exception::Type::MATRIX_ZERO_SIZE);

	// check square matrices
	if (!internal::_check_square_mat(A) || !internal::_check_square_mat(B))
		throw Exception("anticomm", Exception::Type::MATRIX_NOT_SQUARE);

	// check equal dimensions
	if (A.rows() != B.rows())
		throw Exception("anticomm", Exception::Type::DIMS_NOT_EQUAL);

	return A * B + B * A;
}

}

#endif /* FUNCTIONAL_H_ */
