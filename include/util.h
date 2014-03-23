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

// eigenvalues
inline types::cmat evals(const types::cmat &A)
{
	Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es(A);
	return es.eigenvalues();
}

// eigenvectors
inline types::cmat evects(const types::cmat &A)
{
	Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es(A);
	return es.eigenvectors();
}

// Kronecker product of 2 matrices
types::cmat kron(const types::cmat &A, const types::cmat &B);

// Kronecker product of a list of matrices
types::cmat kron_list(const std::vector<types::cmat> & list);

// Kronecker product of a matrix with itself $n$ times
types::cmat kron_pow(const types::cmat &A, int n);

// integer index to multi-index, use C-style array for speed
void _n2multiidx(const size_t n, const size_t numdims, const size_t *dims,
		size_t *result);

// multi index to integer index, use C-style array for speed
size_t _multiidx2n(const size_t *midx, const size_t numdims,
		const size_t *dims);

// Partial trace over subsystem B in a D_A X D_B system
types::cmat ptrace2(const types::cmat &AB, const std::vector<size_t> dims);

// permutes the subsystems in a cmat
types::cmat syspermute(const types::cmat &A, const std::vector<size_t> perm,
		const std::vector<size_t> &dims);

// partial trace
types::cmat ptrace(const types::cmat &A, const std::vector<size_t> &subsys,
		const std::vector<size_t> &dims);

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

// Displays a complex number in friendly form
void disp(const types::cplx &c, std::ostream& os = std::cout,
		unsigned int precision = 4, double eps = 1e-16);

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
inline void disp_container(const T& x)
{
	std::cout<<"[ ";
	for (typename T::const_iterator it = x.begin(); it != x.end(); it++)
		std::cout << *it << " ";
	std::cout<<"]";
}

// used inside the #pragma omp parallel for in syspermute
void _syspermute_worker(const size_t numdims, const size_t *cdims,
		const size_t *cperm, const size_t i, const size_t j, size_t &iperm,
		size_t &jperm, const types::cmat &A, types::cmat &result);

}

#endif /* UTIL_H_ */
