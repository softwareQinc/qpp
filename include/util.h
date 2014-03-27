/* 
 * File:   util.h
 * Author: vlad
 *
 * Created on December 12, 2013, 10:41 PM
 */

#ifndef UTIL_H_
#define	UTIL_H_

#include <Eigen/QR>
#include <stdexcept>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "types.h"
#include "util.h"
#include "constants.h"
#include "internal.h"
#include "stat.h"

// utility functions

// TODO: expandout function
// TODO: dyad function
// TODO: proj (dya) function
// TODO: ip (inner product function) function, make it general to return matrices
// TODO: entropy functions etc

namespace qpp
{

// Eigen function wrappers (inlines)

// transpose, preserve return type
template<typename Derived>
Derived transpose(const Eigen::MatrixBase<Derived>& A)
{
	return A.transpose();
}

// conjugate, preserve return type
template<typename Derived>
Derived conjugate(const Eigen::MatrixBase<Derived>& A)
{
	return A.conjugate();
}

// adjoint, preserve return type
template<typename Derived>
Derived adjoint(const Eigen::MatrixBase<Derived>& A)
{
	return A.adjoint();
}

// trace, preserve return type
template<typename Derived>
typename Derived::Scalar trace(const Eigen::MatrixBase<Derived>& A)
{
	return A.trace();
}

// absolute values component-wise, returns double
template<typename Derived>
Derived abs(const Eigen::MatrixBase<Derived>& A)
{
	Derived result = Derived::Zero(A.rows(), A.cols());
	for (size_t i = 0; i < A.rows(); i++)
		for (size_t j = 0; j < A.cols(); j++)
			result(i, j) = std::abs(A(i, j));
	return result;
}

// trace-norm (or Frobenius norm) (CHANGES return type to double)
template<typename Derived>
typename Eigen::MatrixXd::Scalar norm(const Eigen::MatrixBase<Derived>& A)
{
	// convert matrix to complex then return its norm
	return (A.template cast<types::cplx>()).norm();
}

// eigenvalues (CHANGES return type to complex)
template<typename Derived>
Eigen::MatrixXcd evals(const Eigen::MatrixBase<Derived>& A)
{
	return (A.template cast<types::cplx>()).eigenvalues();
}

// eigenvectors (CHANGES return type to complex matrix)
template<typename Derived>
Eigen::MatrixXcd evects(const Eigen::MatrixBase<Derived>& A)
{
	Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es(
			A.template cast<types::cplx>());
	return es.eigenvectors();
}

// Kronecker product of 2 matrices, preserve return type
template<typename Derived>
Derived kron(const Eigen::MatrixBase<Derived> &A,
		const Eigen::MatrixBase<Derived> &B)
{
	int Acols = A.cols();
	int Arows = A.rows();
	int Bcols = B.cols();
	int Brows = B.rows();

	Derived result;
	result.resize(Arows * Brows, Acols * Bcols);

	for (int i = 0; i < Arows; i++)
		for (int j = 0; j < Acols; j++)
			result.block(i * Brows, j * Bcols, Brows, Bcols) = A(i, j) * B;
	return result;
}

// Kronecker product of a list of matrices, preserve return type
// <Derived> is forced to be a matrix by invocation of kron inside the function
template<typename Derived>
Derived kron_list(const std::vector<Derived> &list)
{
	Derived result = list[0];
	for (size_t i = 1; i < list.size(); i++)
		result = kron(result, list[i]);
	return result;
}
// Kronecker product of a matrix with itself $n$ times, preserve return type
template<typename Derived>
Derived kron_pow(const Eigen::MatrixBase<Derived> &A, size_t n)
{
	std::vector<Derived> list;
	for (size_t i = 0; i < n; i++)
		list.push_back(A);
	return kron_list(list);
}

// Matrix power A^z (CHANGES return type to complex matrix)
template<typename Derived>
Eigen::MatrixXcd mpower(const Eigen::MatrixBase<Derived> &A,
		const types::cplx z)

{
	// check square matrix
	if (!internal::_check_square_mat(A))
		throw std::runtime_error("mpower: Matrix must be square!");

	// Define A^0 = Id
	if (real(z) == 0 && imag(z) == 0)
	{
		Eigen::MatrixXcd result(A.rows(), A.rows());
		result.setIdentity();
		return result;
	}

	Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es(
			A.template cast<types::cplx>());
	Eigen::MatrixXcd evects = es.eigenvectors();
	Eigen::MatrixXcd evals = es.eigenvalues();
	for (int i = 0; i < evals.rows(); i++)
		evals(i) = std::pow(static_cast<types::cplx>(evals(i)),
				static_cast<types::cplx>(z));

	Eigen::MatrixXcd evalsdiag = evals.asDiagonal();

	return evects * evalsdiag * evects.inverse();
}

// Integer matrix power, preserve return type
// Explicitly multiply the matrix with itself n times
template<typename Derived>
Derived mpower(const Eigen::MatrixBase<Derived> &A, const size_t n)
{
	// check square matrix
	if (!internal::_check_square_mat(A))
		throw std::runtime_error("mpower: Matrix must be square!");

	Derived result = A;

	if (n == 0)
		return result.setIdentity();

	for (size_t i = 1; i < n; i++)
		result *= A;

	return result;
}



// Random double matrix with entries in Normal(0,1)
inline Eigen::MatrixXd randn(const size_t rows, const size_t cols)
{
	stat::NormalDistribution nd; // N(0,1)
	Eigen::MatrixXd A(rows, cols);

	for (size_t i = 0; i < rows; i++)
		for (size_t j = 0; j < cols; j++)
			A(i, j) = nd.sample();
	return A;
}

// Random square matrix with entries in Normal(0,1)
inline Eigen::MatrixXd randn(const size_t rows)
{
	return randn(rows, rows);
}

// Random unitary matrix
// TODO: Use QR decomposition followed by multiplication (-1 1 1 -1 ...) as in Toby's example
inline Eigen::MatrixXcd rand_unitary(const size_t size)
{
	Eigen::MatrixXcd X(size, size);

	X.real() = 1. / sqrt(2) * randn(size);
	X.imag() = 1. / sqrt(2) * randn(size);
	Eigen::HouseholderQR<Eigen::MatrixXcd> qr(X);

	Eigen::MatrixXcd Q = qr.householderQ();

	//= randn(size)+ct:ii*randn(size);

	//H = (H + adjoint(H)) / 2;

	return Q;
}

// reshape the columns of A and returns a matrix with m rows and n columns
// use column-major order (same as MATLAB)
template<typename Derived>
Derived reshape(const Eigen::MatrixBase<Derived>& A, size_t rows, size_t cols)
{
	size_t rowsA = A.rows();
	size_t colsA = A.cols();

	if (rowsA * colsA != rows * cols)
		throw std::runtime_error("reshape: Dimension mismatch!");

	return Eigen::Map<Derived>(static_cast<Derived>(A).data(), rows, cols);
}

// permutes the subsystems in a matrix
template<typename Derived>
Derived syspermute(const Eigen::MatrixBase<Derived> &A,
		const std::vector<size_t> perm, const std::vector<size_t> &dims)
{
	// Error checks

	// check square matrix
	if (!internal::_check_square_mat(A))
		throw std::runtime_error("syspermute: Matrix must be square!");

	// check that dims is a valid dimension vector
	if (!internal::_check_dims(dims))
		throw std::runtime_error("syspermute: Invalid dimensions vector!");

	// check that dims match the dimension of A
	if (!internal::_check_dims_match_mat(dims, A))
		throw std::runtime_error(
				"syspermute: Dimenisons vector does not match the dimension of the matrix!");

	// check that the size of the permutation is OK
	if (!internal::_check_perm(perm, dims))
		throw std::runtime_error("syspermute: Invalid permutation size!");

	size_t dim = static_cast<size_t>(A.rows());
	const size_t numdims = dims.size();
	size_t *cdims = new size_t[numdims];
	size_t *cperm = new size_t[numdims];
	Derived result(dim, dim);

	// copy dims in cdims and perm in cperm
	for (size_t i = 0; i < numdims; i++)
	{
		cdims[i] = dims[i];
		cperm[i] = perm[i];
	}

	size_t iperm = 0;
	size_t jperm = 0;
	for (size_t i = 0; i < dim; i++)
#pragma omp parallel for
		for (size_t j = 0; j < dim; j++)
			internal::_syspermute_worker(numdims, cdims, cperm, i, j, iperm,
					jperm, A, result);

	delete[] cdims;
	delete[] cperm;

	return result; // the permuted matrix
}

// Partial trace over subsystem B in a D_A x D_B system
template<typename Derived>
Derived ptrace2(const Eigen::MatrixBase<Derived> &A,
		const std::vector<size_t> dims)
{
	// Error checks
	// error checks

	// check square matrix
	if (!internal::_check_square_mat(A))
		throw std::runtime_error("ptrace2: Matrix must be square!");

	// check dims has only 2 elements
	if (dims.size() != 2)
		throw std::runtime_error("ptrace2: Must have only 2 dimensions!");

	// check that dim is a valid dimension vector
	if (!internal::_check_dims(dims))
		throw std::runtime_error("ptrace2: Invalid dimensions vector!");

	// check that dims match the dimension of A
	if (!internal::_check_dims_match_mat(dims, A))
		throw std::runtime_error(
				"ptrace2: Dimenisons vector does not match the dimension of the matrix!");

	size_t DA = dims[0];
	size_t DB = dims[1];

	Derived result = Derived::Zero(DA, DA);

	for (size_t i = 0; i < DA; i++)
#pragma omp parallel for
		for (size_t j = 0; j < DA; j++)
		{
			result(i, j) = trace(
					static_cast<Derived>(A.block(i * DB, j * DB, DB, DB)));
		}
	return result;
}

// partial trace
template<typename Derived>
Derived ptrace(const Eigen::MatrixBase<Derived> &A,
		const std::vector<size_t> &subsys, const std::vector<size_t> &dims)
{
	// error checks

	// check square matrix
	if (!internal::_check_square_mat(A))
		throw std::runtime_error("ptrace: Matrix must be square!");

	// check that dims is a valid dimension vector
	if (!internal::_check_dims(dims))
		throw std::runtime_error("ptrace: Invalid dimensions vector!");

	// check that dims match the dimension of A
	if (!internal::_check_dims_match_mat(dims, A))
		throw std::runtime_error(
				"ptrace: Dimenisons vector does not match the dimension of the matrix!");

	if (!internal::_check_subsys(subsys, dims))
		throw std::runtime_error("ptrace: Invalid subsystems!");

	size_t dim = static_cast<size_t>(A.rows());
	size_t numsubsys = subsys.size(); // number of subsystems we trace out
	size_t numdims = dims.size(); // total number of subsystems;
	std::vector<size_t> perm(numdims, 0); // the permutation vector
	std::vector<size_t> permdims; // the permuted dimensions

	Derived result;

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
		if (std::find(subsys.begin(), subsys.end(), i) != subsys.end())
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
template<typename Derived>
Derived ptranspose(const Eigen::MatrixBase<Derived>& A,
		const std::vector<size_t>& subsys, const std::vector<size_t>& dims)
{
	// error checks

	// check square matrix
	if (!internal::_check_square_mat(A))
		throw std::runtime_error("ptranspose: Matrix must be square!");

	// check that dims is a valid dimension vector
	if (!internal::_check_dims(dims))
		throw std::runtime_error("ptranspose: Invalid dimensions vector!");

	// check that dims match the dimension of A
	if (!internal::_check_dims_match_mat(dims, A))
		throw std::runtime_error(
				"ptranspose: Dimenisons vector does not match the dimension of the matrix!");

	if (!internal::_check_subsys(subsys, dims))
		throw std::runtime_error("ptranspose: Invalid subsystems!");

	size_t dim = static_cast<size_t>(A.rows());
	const size_t numdims = dims.size();
	const size_t numsubsys = subsys.size();
	size_t *cdims = new size_t[numdims];
	size_t *midxrow = new size_t[numdims];
	size_t *csubsys = new size_t[numsubsys];

	Derived result = A;

	// copy dims in cdims and subsys in csubsys
	for (size_t i = 0; i < numdims; i++)
		cdims[i] = dims[i];
	for (size_t i = 0; i < numsubsys; i++)
		csubsys[i] = subsys[i];

	size_t iperm = 0;
	size_t jperm = 0;
	for (size_t i = 0; i < dim; i++)
	{
		// compute the row multi-index
		internal::_n2multiidx(i, numdims, cdims, midxrow);
#pragma omp parallel for
		for (size_t j = 0; j < dim; j++) // paralelize this code
			internal::_ptranspose_worker(midxrow, numdims, numsubsys, cdims,
					csubsys, i, j, iperm, jperm, A, result);
	}

	delete[] midxrow;
	delete[] cdims;
	delete[] csubsys;

	return result;
}

// Random matrix with entries in Uniform[0,1]
template<typename Derived>
Derived rand(const size_t rows, const size_t cols)
{
	return Derived::Random(rows, cols);
}

// Random square matrix with entries in Uniform[0,1]
template<typename Derived>
Derived rand(const size_t rows)
{
	return rand<Derived>(rows, rows);
}

// save matrix to a binary file in double precision
template<typename Derived>
void save(const Eigen::MatrixBase<Derived> & A, const std::string& fname)
{
	std::fstream fout;
	fout.open(fname.c_str(), std::ios::out | std::ios::binary);

	if (fout.fail())
	{
		throw std::runtime_error(
				"save: Error writing output file \"" + std::string(fname)
						+ "\"!");
	}
	// we write a header to the file
	const char _header[] = "TYPE::Eigen::Matrix";
	fout.write(_header, sizeof(_header));

	size_t rows = A.rows();
	size_t cols = A.cols();
	fout.write((char*) &rows, sizeof(rows));
	fout.write((char*) &cols, sizeof(cols));
	for (size_t i = 0; i < rows; i++)
		for (size_t j = 0; j < cols; j++)
			fout.write((char*) &A(i, j), sizeof(A(i, j)));
	fout.close();
}

// load matrix from binary file
template<typename Derived>
Derived load(const std::string& fname)
{
	std::fstream fin;
	fin.open(fname.c_str(), std::ios::in | std::ios::binary);

	if (fin.fail())
	{
		throw std::runtime_error(
				"load: Error opening input file \"" + std::string(fname)
						+ "\"!");
	}

	const char _header[] = "TYPE::Eigen::Matrix"; // we write the header to the file
	char *_fheader = new char[sizeof(_header)]; // include the zero-termiation string

	// read the header
	fin.read(_fheader, sizeof(_header));
	if (strcmp(_fheader, _header))
	{
		delete[] _fheader;
		throw std::runtime_error(
				"load: Input file \"" + std::string(fname)
						+ "\" is corrupted!");
	}
	delete[] _fheader;

	size_t rows, cols;
	fin.read((char*) &rows, sizeof(rows));
	fin.read((char*) &cols, sizeof(cols));

	Derived A(rows, cols);
	for (size_t i = 0; i < rows; i++)
		for (size_t j = 0; j < cols; j++)
			fin.read((char*) &A(i, j), sizeof(A(i, j)));
	fin.close();
	return A;
}

//TODO: use 1.0e+05 notation

// Displays an Eigen::MatrixX in friendly form
template<typename Derived>
void disp(const Eigen::MatrixBase<Derived> &A, std::ostream& os = std::cout,
		unsigned int precision = 4)
{
	//std::cout << "typeid: " << typeid(Derived).name() << std::endl;
	os << std::setprecision(precision) << std::fixed << A;
}

template<> // complex matrix specialization
inline void disp(const Eigen::MatrixBase<Eigen::MatrixXcd> &A, std::ostream& os,
		unsigned int precision)

{
	if (A.rows() * A.cols() == 0)
	{
		os << "Empty [" << A.rows() << " x " << A.cols() << "] matrix.";
		os << std::endl;
		return;
	};

	// everything smaller than eps is displayed as 0
	double eps = std::numeric_limits<double>::epsilon();

	std::ostringstream ostr;
	std::vector<std::string> vstr;
	std::string strtmp;

	for (int i = 0; i < A.rows(); i++)
	{
		for (int j = 0; j < A.cols(); j++)
		{
			strtmp.clear(); // clear the temporary string
			ostr.clear();
			ostr.str(std::string()); // clear the ostringstream,

			double re = static_cast<types::cplx>(A(i, j)).real();
			double im = static_cast<types::cplx>(A(i, j)).imag();

			if (std::abs(re) < eps && std::abs(im) < eps)
			{
				vstr.push_back("0 ");
			}
			else if (std::abs(re) < eps)
			{
				ostr << std::setprecision(precision) << std::fixed << im;
				vstr.push_back(ostr.str() + "i");
			}
			else if (std::abs(im) < eps)
			{
				ostr << std::setprecision(precision) << std::fixed << re;
				vstr.push_back(ostr.str() + " ");
			}
			else
			{
				ostr << std::setprecision(precision) << std::fixed << re;
				strtmp = ostr.str();

				strtmp += (im > 0 ? " + " : " - ");
				ostr.clear();
				ostr.str(std::string()); // clear
				ostr << std::setprecision(precision) << std::fixed
						<< std::abs(im);
				strtmp += ostr.str();
				strtmp += "i";
				vstr.push_back(strtmp);
			}
		}
	}

	// determine the maximum lenght of the entries in each column
	std::vector<size_t> maxlengthcols(A.cols(), 0);

	for (int i = 0; i < A.rows(); i++)
		for (int j = 0; j < A.cols(); j++)
			if (vstr[i * A.cols() + j].size() > maxlengthcols[j])
				maxlengthcols[j] = vstr[i * A.cols() + j].size();

	// finally display it!
	for (int i = 0; i < A.rows(); i++)
	{
		os << std::setw(static_cast<int>(maxlengthcols[0])) << std::right
				<< vstr[i * A.cols()]; // display first column
		for (int j = 1; j < A.cols(); j++) // then the rest
			os << std::setw(static_cast<int>(maxlengthcols[j] + 2))
					<< std::right << vstr[i * A.cols() + j];
		if (i < A.rows() - 1)
			os << std::endl;
	}
}

// Displays a complex number in friendly form
inline void disp(const types::cplx c, std::ostream& os = std::cout,
		unsigned int precision = 4)
{
	// put the complex number inside an Eigen matrix
	Eigen::MatrixXcd tmp(1, 1);
	tmp(0, 0) = c;
	disp(tmp, os, precision);
}

}

#endif /* UTIL_H_ */
