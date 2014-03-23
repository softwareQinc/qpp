/* 
 * File:   util.cpp
 * Author: vlad
 *
 * Created on December 12, 2013, 10:43 PM
 */

#include <Eigen/Eigenvalues>
#include <Eigen/Dense>
#include <iostream>
#include <cstring>
#include <stdexcept>
#include <vector>
#include <algorithm>
#include "util.h"
#include "types.h"
#include "stat.h"
#include "constants.h"

namespace qpp
{

// Kronecker product of 2 matrices
types::cmat kron(const types::cmat &A, const types::cmat &B)
{
	int Acols = A.cols();
	int Arows = A.rows();
	int Bcols = B.cols();
	int Brows = B.rows();

	types::cmat result(Arows * Brows, Acols * Bcols);

	for (int i = 0; i < Arows; i++)
		for (int j = 0; j < Acols; j++)
			result.block(i * Brows, j * Bcols, Brows, Bcols) = A(i, j) * B;
	return result;
}

// Kronecker product of a list of matrices
types::cmat kron_list(const std::vector<types::cmat> & list)
{
	types::cmat result;
	result = list[0];
	for (unsigned int i = 1; i < list.size(); i++)
		result = kron(result, (types::cmat) list[i]);
	return result;
}

// Kronecker product of a matrix with itself $n$ times
types::cmat kron_pow(const types::cmat &A, size_t n)
{
	std::vector<types::cmat> list;
	for (size_t i = 0; i < n; i++)
		list.push_back(A);
	return kron_list(list);
}

// integer index to multi-index, use C-style array for speed
void _n2multiidx(const size_t n, const size_t numdims, const size_t *dims,
		size_t *result)
{
	size_t maxn = 1;
	for (size_t i = 0; i < numdims; i++)
		maxn *= dims[i];
	if (n > maxn - 1)
		throw std::runtime_error("Number too large, out of bounds!");

	size_t _n = n;
	for (size_t i = 0; i < numdims; i++)
	{
		result[numdims - i - 1] = _n % static_cast<int>(dims[numdims - i - 1]);
		_n = _n / static_cast<int>(dims[numdims - i - 1]);
	}
}

// multi index to integer index, use C-style array for speed
size_t _multiidx2n(const size_t *midx, const size_t numdims, const size_t *dims)
{
	for (size_t i = 0; i < numdims; i++)
		if (midx[i] >= dims[i])
			throw std::runtime_error(
					"Sub-index exceeds corresponding dimension!");

	size_t *part_prod = new size_t[numdims];

	part_prod[numdims - 1] = 1;
	for (size_t j = 1; j < numdims; j++)
		part_prod[numdims - j - 1] = part_prod[numdims - j] * dims[numdims - j];

	size_t result = 0;
	for (size_t i = 0; i < numdims; i++)
		result += midx[i] * part_prod[i];

	delete[] part_prod;

	return result;
}

// Partial trace over subsystem B in a D_A x D_B system
types::cmat ptrace2(const types::cmat &AB, const std::vector<size_t> dims)
{
	// Error checks

	size_t D = static_cast<size_t>(AB.rows());
	if (D != static_cast<size_t>(AB.cols()))
		throw std::runtime_error("Matrix must be square!");
	size_t DA = dims[0];
	size_t DB = dims[1];
	if (DA * DB != D)
		throw std::runtime_error(
				"Product of partial dimensions must equal the dimension of the matrix!");

	types::cmat result = types::cmat::Zero(DA, DA);

	for (size_t i = 0; i < DA; i++)
#pragma omp parallel for
		for (size_t j = 0; j < DA; j++)
		{
			result(i, j) = trace(
					static_cast<types::cmat>(AB.block(i * DB, j * DB, DB, DB)));
		}
	return result;
}

// permutes the subsystems in a cmat
types::cmat syspermute(const types::cmat &A, const std::vector<size_t> perm,
		const std::vector<size_t> &dims)
{
	// Error checks

	// check square matrix
	size_t dimA = static_cast<size_t>(A.rows());
	if (dimA != static_cast<size_t>(A.cols()))
		throw std::runtime_error("Matrix must be square!");

	// check that the size of the permutation is OK
	const size_t numdims = dims.size();
	if (numdims != perm.size())
		throw std::runtime_error("Invalid permutation size!");

	// check that we have consistent dimensions
	size_t dimtotal = 1;
	for (size_t i = 0; i < numdims; i++)
		dimtotal *= dims[i];
	if (dimtotal != dimA)
		throw std::runtime_error(
				"Dimension list does not match the dimension of the matrix!");

	// check that the permutation is valid
	std::vector<size_t> sort_perm = perm;
	std::sort(sort_perm.begin(), sort_perm.end());
	for (size_t i = 0; i < numdims; i++)
	{
		//std::cout<<sort_perm[i]<<" "<<std::endl;
		if (sort_perm[i] != i)
		{
			throw std::runtime_error("Not a valid permutation!");
		}
	}
	// check dimensions
	size_t tmp = 1;
	for (size_t i = 0; i < numdims; i++)
		tmp *= dims[i];
	if (tmp != dimA)
		throw std::runtime_error("Dimension mismatch!");

	types::cmat result(dimA, dimA);

	size_t *cdims = new size_t[numdims];
	size_t *cperm = new size_t[numdims];

	// copy dims in cdims and perm in cperm
	for (size_t i = 0; i < numdims; i++)
	{
		cdims[i] = dims[i];
		cperm[i] = perm[i];
	}

	size_t iperm = 0;
	size_t jperm = 0;
	for (size_t i = 0; i < dimA; i++)
#pragma omp parallel for
		for (size_t j = 0; j < dimA; j++)
			_syspermute_worker(numdims, cdims, cperm, i, j, iperm, jperm, A,
					result);

	delete[] cdims;
	delete[] cperm;

	return result; // the permuted matrix
}

void _syspermute_worker(const size_t numdims, const size_t *cdims,
		const size_t *cperm, const size_t i, const size_t j, size_t &iperm,
		size_t &jperm, const types::cmat &A, types::cmat &result)
{
	size_t *midxrow = new size_t[numdims];
	size_t *midxcol = new size_t[numdims];
	size_t *midxrowtmp = new size_t[numdims];
	size_t *midxcoltmp = new size_t[numdims];
	size_t *permdims = new size_t[numdims];

	for (size_t i = 0; i < numdims; i++)
		permdims[i] = cdims[cperm[i]]; // permuted dimensions

	// compute the row and col multi-indexes
	_n2multiidx(i, numdims, cdims, midxrow);
	_n2multiidx(j, numdims, cdims, midxcol);

	for (size_t k = 0; k < numdims; k++)
	{
		midxrowtmp[k] = midxrow[cperm[k]]; // permuted multi-indexes
		midxcoltmp[k] = midxcol[cperm[k]]; // permuted multi-indexes
	}

	// move back to integer indexes
	iperm = _multiidx2n(midxrowtmp, numdims, permdims);
	jperm = _multiidx2n(midxcoltmp, numdims, permdims);
	result(iperm, jperm) = A(i, j);

	delete[] midxrow;
	delete[] midxcol;
	delete[] midxrowtmp;
	delete[] midxcoltmp;
	delete[] permdims;
}

// partial trace
types::cmat ptrace(const types::cmat &A, const std::vector<size_t> &subsys,
		const std::vector<size_t> &dims)
{
	types::cmat result;
	std::vector<size_t> permdims;

	std::vector<size_t> subsyssort = subsys;

	// sort the subsystems
	std::sort(subsyssort.begin(), subsyssort.end());
	//print_container(subsyssort);

	size_t numsubsys = subsyssort.size(); // number of subsystems we trace out
	size_t numdims = dims.size(); // total number of subsystems;
	std::vector<size_t> perm(numdims, 0); // the permutation vector

	// the total dimension of the traced-out subsystems
	size_t dimsubsys = 1;
	for (size_t i = 0; i < numsubsys; i++)
		dimsubsys *= dims[subsyssort[i]];

	// total dimension of A, must be the same as A.cols() && A.rows()
	size_t dimtotal = 1;
	for (size_t i = 0; i < numdims; i++)
		dimtotal *= dims[i];

	std::vector<size_t> sizeAB;
	sizeAB.push_back(dimtotal / dimsubsys);
	sizeAB.push_back(dimsubsys);

	// construct the permutation that bring the traced-out subsystems to the end
	size_t cnt0 = 0;
	size_t cnt1 = 0;
	for (size_t i = 0; i < numdims; i++)
	{
		// we find that i belongs to the subsystem
		if (std::find(subsyssort.begin(), subsyssort.end(), i)
				!= subsyssort.end())
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

// Matrix power A^x
types::cmat mat_pow(const types::cmat &A, const types::cplx z)
{
	Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es(A);
	types::cmat evects = es.eigenvectors();
	types::cmat evals = es.eigenvalues();
	for (int i = 0; i < evals.rows(); i++)
		evals(i) = std::pow(static_cast<types::cplx>(evals(i)),
				static_cast<types::cplx>(z));

	types::cmat evalsdiag = evals.asDiagonal();

	return evects * evalsdiag * evects.inverse();
}

// Matrix functional calculus
// Computes f(A), where (*f) is the function pointer
types::cmat mat_f(const types::cmat &A, types::cplx (*f)(const types::cplx &))
{
	Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es(A);
	types::cmat evects = es.eigenvectors();
	types::cmat evals = es.eigenvalues();
	for (int i = 0; i < evals.rows(); i++)
		evals(i) = (*f)(evals(i)); // apply f(x) to each eigenvalue

	types::cmat evalsdiag = evals.asDiagonal();

	return evects * evalsdiag * evects.inverse();
}

// Matrix exponential
types::cmat mat_exp(const types::cmat &A)
{
	return mat_f(A, std::exp);
}

// Random matrix with entries in Uniform[0,1]
types::cmat rand(const size_t rows, const size_t cols)
{
	return Eigen::MatrixXcd::Random(rows, cols);
}

// Random square matrix with entries in Uniform[0,1]
types::cmat rand(const size_t rows)
{
	return rand(rows, rows);
}

// Random matrix with entries in Normal(0,1)
types::cmat randn(const size_t rows, const size_t cols)
{
	stat::NormalDistribution nd; // N(0,1)
	types::cmat A(rows, cols);
	double re, im;

	for (size_t i = 0; i < rows; i++)
		for (size_t j = 0; j < cols; j++)
		{
			re = nd.sample();
			im = nd.sample();
			A(i, j) = re + ct::ii * im;
		}
	return A;
}

// Random square matrix with entries in Normal(0,1)
types::cmat randn(const size_t rows)
{
	return randn(rows, rows);
}

// Random unitary matrix
types::cmat rand_unitary(const size_t size)
{
	types::cmat H = randn(size);
	H = (H + adjoint(H)) / 2;

	return mat_exp(static_cast<types::cmat>(ct::ii * H));
}

// Displays a complex Eigen::Matrix (types::cmat) in friendly form
void disp(const types::cmat &A, std::ostream& os, unsigned int precision,
		double eps)
{
	if (A.rows() * A.cols() == 0)
	{
		os << "Empty [" << A.rows() << " x " << A.cols() << "] matrix."
				<< std::endl;
		return;
	};

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
		os << std::endl;
	}
}

// Displays a complex number in friendly form
void disp(const types::cplx &c, std::ostream& os, unsigned int precision,
		double eps)
{
	// put the complex number inside an Eigen matrix
	types::cmat tmp(1, 1);
	tmp(0, 0) = c;
	disp(tmp, os, precision, eps);
}

// save matrix to a binary file in double precision
void save(const types::cmat & A, const std::string& fname)
{
	std::fstream fout;
	fout.open(fname.c_str(), std::ios::out | std::ios::binary);

	if (fout.fail())
	{
		throw std::runtime_error(
				"Error writing output file \"" + std::string(fname) + "\"!");
	}
	// we write a header to the file
	const char _header[] = "TYPES::CMAT";
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
types::cmat load(const std::string& fname)
{
	std::fstream fin;
	fin.open(fname.c_str(), std::ios::in | std::ios::binary);

	if (fin.fail())
	{
		throw std::runtime_error(
				"Error opening input file \"" + std::string(fname) + "\"!");
	}

	const char _header[] = "TYPES::CMAT"; // we write the header to the file
	char *_fheader = new char[sizeof(_header)]; // include the zero-termiation string

	// read the header
	fin.read(_fheader, sizeof(_header));
	if (strcmp(_fheader, _header))
	{
		delete[] _fheader;
		throw std::runtime_error(
				"Input file \"" + std::string(fname) + "\" is corrupted!");
	}
	delete[] _fheader;

	size_t rows, cols;
	fin.read((char*) &rows, sizeof(rows));
	fin.read((char*) &cols, sizeof(cols));

	types::cmat A(rows, cols);
	for (size_t i = 0; i < rows; i++)
		for (size_t j = 0; j < cols; j++)
			fin.read((char*) &A(i, j), sizeof(A(i, j)));
	fin.close();
	return A;
}

// reshape the columns of A and returns a cmat with m rows and n columns
// use column-major order (same as MATLAB)
types::cmat reshape(const types::cmat& A, size_t rows, size_t cols)
{
	size_t rowsA = A.rows();
	size_t colsA = A.cols();

	if (rowsA * colsA != rows * cols)
		throw std::runtime_error("Dimension mismatch, cannot reshape!");

	Eigen::MatrixXd realA = A.real();
	Eigen::MatrixXd imagA = A.imag();

	realA = Eigen::Map<Eigen::MatrixXd>(realA.data(), rows, cols);
	imagA = Eigen::Map<Eigen::MatrixXd>(imagA.data(), rows, cols);

	return realA.cast<types::cplx>() + ct::ii * imagA.cast<types::cplx>();
}

}
