/* 
 * File:   util.h
 * Author: vlad
 *
 * Created on December 12, 2013, 10:41 PM
 */

#ifndef UTIL_H_
#define	UTIL_H_

#include <stdexcept>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "types.h"
#include "util.h"
#include "constants.h"

// utility functions

namespace qpp
{

// Kronecker product of 2 matrices
types::cmat kron(const types::cmat &A, const types::cmat &B);

// Kronecker product of a list of matrices
types::cmat kron_list(const std::vector<types::cmat> & list);

// Kronecker product of a matrix with itself $n$ times
types::cmat kron_n(const types::cmat &A, int n);

// Partial trace over subsystem B in a D_A X D_B system
types::cmat ptrace2(const types::cmat &AB, const std::vector<size_t> dims);

// permutes the subsystems in a cmat
types::cmat syspermute(const types::cmat &A, const std::vector<size_t> &dims,
		const std::vector<size_t> perm);

// partial trace
types::cmat ptrace(const types::cmat &A, const std::vector<size_t> &dims,
		const std::vector<size_t> &subsys);

// TODO: expandout function

// Matrix power A^z
types::cmat mat_pow(const types::cmat &A, const types::cplx z);

// Matrix functional calculus
// Computes f(A), where f() is the function pointer
types::cmat mat_f(const types::cmat &A, types::cplx (*)(const types::cplx &));

// Matrix exponential
types::cmat mat_exp(const types::cmat &A);

// Random matrix with entries in Uniform[0,1]
types::cmat rand(const size_t rows, const size_t cols);

// Random square matrix with entries in Uniform[0,1]
types::cmat rand(const size_t rows);

// Random matrix with entries in Normal(0,1)
types::cmat randn(const size_t rows, const size_t cols);

// Random square matrix with entries in Normal(0,1)
types::cmat randn(const size_t rows);

// Random unitary matrix
types::cmat rand_unitary(const size_t size);

// Displays a complex Eigen::Matrix (types::cmat) in friendly form
void disp(const types::cmat &A, std::ostream& os = std::cout,
		unsigned int precision = 4, double eps = 1e-16);

// save matrix in a text file (text format, lacks precision)
void save_text(const types::cmat & A, const std::string& fname,
		size_t precision = 16);

// load matrix from text file
types::cmat load_text(const std::string& fname);

// save matrix to a binary file in double precision
void save(const types::cmat & A, const std::string& fname);

// load matrix from binary file
types::cmat load(const std::string& fname);

// reshape the columns of A and returns a cmat with m rows and n columns
// use column-major order (same as MATLAB)
types::cmat reshape(const types::cmat& A, size_t rows, size_t cols);

// inline templates

// Displays a standard container that supports STL iterators
template<typename T>
inline void print_container(const T& x)
{
	for (typename T::const_iterator it = x.begin(); it != x.end(); it++)
		std::cout << *it << " ";
	std::cout << std::endl;
}

// integer index to multi-index
inline std::vector<size_t> n2multiidx(const size_t &n,
		const std::vector<size_t> &dims)
{
	size_t numdims = dims.size();
	std::vector<size_t> result(numdims, 0);
	int _n = n;
	size_t maxn = 1;
	for (size_t i = 0; i < numdims; i++)
		maxn *= dims[i];
	if (n > maxn - 1)
		throw std::runtime_error("Number too large, out of bounds!");

	size_t tmp = 0;
	for (size_t i = 0; i < numdims; i++)
	{
		tmp = _n % static_cast<int>(dims[numdims - i - 1]);
		result[numdims - i - 1] = tmp;
		_n = _n / static_cast<int>(dims[numdims - i - 1]);
	}
	return result;
}

// multi index to integer index
inline size_t multiidx2n(const std::vector<size_t> &midx,
		const std::vector<size_t> &dims)
{
	size_t numdims = dims.size();
	std::vector<size_t> part_prod(numdims, 1); // partial products
	size_t result = 0;
	for (size_t i = 0; i < numdims; i++)
		if (midx[i] >= dims[i])
			throw std::runtime_error(
					"Sub-index exceeds corresponding dimension!");

	for (size_t j = 0; j < numdims; j++)
		if (j == 0)
			part_prod[numdims - 1] = 1;
		else
			part_prod[numdims - j - 1] = part_prod[numdims - j]
					* dims[numdims - j];

	for (size_t i = 0; i < numdims; i++)
		result += midx[i] * part_prod[i];

	return result;
}

// Eigen function wrappers

// transpose
inline types::cmat transpose(const types::cmat& A)
{
	return A.transpose();
}

// conjugate
inline types::cmat conjugate(const types::cmat& A)
{
	return A.conjugate();
}

// adjoint
inline types::cmat adjoint(const types::cmat& A)
{
	return A.adjoint();
}

// trace
inline types::cplx trace(const types::cmat& A)
{
	return A.trace();
}

// trace-norm (or Frobenius norm)
inline double norm(const types::cmat& A)
{
	return A.norm();
}

}

#endif /* UTIL_H_ */
