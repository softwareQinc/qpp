/*
 * io.h
 *
 *  Created on: Mar 27, 2014
 *      Author: vlad
 */

#ifndef IO_H_
#define IO_H_

// input/output

#include <stdexcept>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "types.h"
#include "exception.h"

namespace qpp
{

// Displays a standard container that supports std::begin and std::end
template<typename T>
void disp(const T& x, const std::string & separator = " ",
		const std::string& start = "[", const std::string& end = "]",
		std::ostream& os = std::cout)
{
	os << start;
	auto it = std::begin(x);
	for (; it != std::end(x) - 1; it++)
		os << *it << separator;
	os << *(it++);
	os << end;
}

// Displays a standard container that supports std::begin and std::end
// and adds a new line
template<typename T>
void displn(const T& x, const std::string & separator = " ",
		const std::string& start = "[", const std::string& end = "]",
		std::ostream& os = std::cout)
{
	disp(x, separator, start, end, os);
	os << std::endl;
}

// Displays a C-style array
template<typename T>
void disp(const T* x, const size_t n, const std::string & separator = " ",
		const std::string& start = "[", const std::string& end = "]",
		std::ostream& os = std::cout)
{
	os << start;
	for (size_t i = 0; i < n - 1; i++)
		os << x[i] << separator;
	os << x[n - 1];
	os << end;
}

// Displays a C-style array
// and adds a new line
template<typename T>
void displn(const T* x, const size_t n, const std::string & separator = " ",
		const std::string& start = "[", const std::string& end = "]",
		std::ostream& os = std::cout)
{
	disp(x, n, separator, start, end, os);
	os << std::endl;
}

// Displays an Eigen::MatrixX in friendly form
template<typename Scalar>
void disp(const types::DynMat<Scalar> &A, double chop = ct::chop,
		std::ostream& os = std::cout)
{
	if (A.size() == 0)
	{
		os << "Empty [" << A.rows() << " x " << A.cols() << "] matrix";
		return;
	};

	std::ostringstream ostr;
	ostr.flags(os.flags()); // get the formatting flags
	ostr.precision(os.precision()); // set precision

	std::vector<std::string> vstr;
	std::string strA;

	for (size_t i = 0; i < static_cast<size_t>(A.rows()); i++)
	{
		for (size_t j = 0; j < static_cast<size_t>(A.cols()); j++)
		{
			strA.clear(); // clear the temporary string
			ostr.clear();
			ostr.str(std::string()); // clear the ostringstream,

			// convert to complex
			double re = static_cast<types::cplx>(A(i, j)).real();
			double im = static_cast<types::cplx>(A(i, j)).imag();

			if (std::abs(re) < chop && std::abs(im) < chop)
			{
				vstr.push_back("0 ");
			}
			else if (std::abs(re) < chop)
			{
				ostr << im;
				vstr.push_back(ostr.str() + "i");
			}
			else if (std::abs(im) < chop)
			{
				ostr << re;
				vstr.push_back(ostr.str() + " ");
			}
			else
			{
				ostr << re;
				strA = ostr.str();

				strA += (im > 0 ? " + " : " - ");
				ostr.clear();
				ostr.str(std::string()); // clear
				ostr << std::abs(im);
				strA += ostr.str();
				strA += "i";
				vstr.push_back(strA);
			}
		}
	}

// determine the maximum lenght of the entries in each column
	std::vector<size_t> maxlengthcols(A.cols(), 0);

	for (size_t i = 0; i < static_cast<size_t>(A.rows()); i++)
		for (size_t j = 0; j < static_cast<size_t>(A.cols()); j++)
			if (vstr[i * A.cols() + j].size() > maxlengthcols[j])
				maxlengthcols[j] = vstr[i * A.cols() + j].size();

// finally display it!
	for (size_t i = 0; i < static_cast<size_t>(A.rows()); i++)
	{
		os << std::setw(static_cast<int>(maxlengthcols[0])) << std::right
				<< vstr[i * A.cols()]; // display first column
		for (size_t j = 1; j < static_cast<size_t>(A.cols()); j++) // then the rest
			os << std::setw(static_cast<int>(maxlengthcols[j] + 2))
					<< std::right << vstr[i * A.cols() + j];

		if (i < static_cast<size_t>(A.rows()) - 1)
			os << std::endl;
	}
}

// Displays an Eigen::MatrixX in friendly form
// and adds a new line
template<typename Scalar>
void displn(const types::DynMat<Scalar> &A, double chop = ct::chop,
		std::ostream& os = std::cout)
{
	disp(A, chop, os);
	os << std::endl;
}

// Displays a complex number in friendly form
void disp(const types::cplx c, double chop = ct::chop, std::ostream& os =
		std::cout)
{
// put the complex number inside an Eigen matrix
	types::cmat A(1, 1);
	A(0, 0) = c;
	disp(A, chop, os);
}

// Displays a complex number in friendly form
// and adds a new line
void displn(const types::cplx c, double chop = ct::chop, std::ostream& os =
		std::cout)
{
	disp(c, chop, os);
	os << std::endl;
}

// save matrix to a binary file in double precision
template<typename Scalar>
void save(const types::DynMat<Scalar> & A, const std::string& fname)

{
	// check zero-size
	if (!internal::_check_nonzero_size(A))
		throw Exception("save", Exception::Type::ZERO_SIZE);

	std::fstream fout;
	fout.open(fname.c_str(), std::ios::out | std::ios::binary);

	if (fout.fail())
	{
		throw std::runtime_error(
				"save: Error writing output file \"" + std::string(fname)
						+ "\"!");
	}

	// write the header to file
	const char _header[] = "TYPE::Eigen::Matrix";
	fout.write(_header, sizeof(_header));

	size_t rows = static_cast<size_t>(A.rows());
	size_t cols = static_cast<size_t>(A.cols());
	fout.write((char*) &rows, sizeof(rows));
	fout.write((char*) &cols, sizeof(cols));

	fout.write((char*) A.data(), sizeof(Scalar) * rows * cols);

	fout.close();
}

// load matrix from binary file
template<typename Scalar>
types::DynMat<Scalar> load(const std::string& fname)
{
	std::fstream fin;
	fin.open(fname.c_str(), std::ios::in | std::ios::binary);

	if (fin.fail())
	{
		throw std::runtime_error(
				"load: Error opening input file \"" + std::string(fname)
						+ "\"!");
	}

	const char _header[] = "TYPE::Eigen::Matrix";
	char* _fheader = new char[sizeof(_header)];

	// read the header from file
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

	types::DynMat<Scalar> A(rows, cols);

	fin.read((char*) A.data(), sizeof(Scalar) * rows * cols);

	fin.close();
	return A;
}

} /* namespace qpp */

#endif /* IO_H_ */
