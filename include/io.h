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

namespace qpp
{

// Displays an Eigen::MatrixX in friendly form
template<typename Scalar>
void disp(const types::DynMat<Scalar> &A, double chop = ct::chop,
		std::ostream& os = std::cout)
{
	types::DynMat<Scalar> tmp = A;

	if (tmp.size() == 0)
	{
		os << "Empty [" << tmp.rows() << " x " << tmp.cols() << "] matrix";
		return;
	};

	std::ostringstream ostr;
	ostr.flags(os.flags()); // get the formatting flags
	ostr.precision(os.precision()); // set precision

	std::vector<std::string> vstr;
	std::string strtmp;

	for (size_t i = 0; i < static_cast<size_t>(tmp.rows()); i++)
	{
		for (size_t j = 0; j < static_cast<size_t>(tmp.cols()); j++)
		{
			strtmp.clear(); // clear the temporary string
			ostr.clear();
			ostr.str(std::string()); // clear the ostringstream,

			double re = static_cast<types::cplx>(tmp(i, j)).real();
			double im = static_cast<types::cplx>(tmp(i, j)).imag();

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
				strtmp = ostr.str();

				strtmp += (im > 0 ? " + " : " - ");
				ostr.clear();
				ostr.str(std::string()); // clear
				ostr << std::abs(im);
				strtmp += ostr.str();
				strtmp += "i";
				vstr.push_back(strtmp);
			}
		}
	}

// determine the maximum lenght of the entries in each column
	std::vector<size_t> maxlengthcols(tmp.cols(), 0);

	for (size_t i = 0; i < static_cast<size_t>(tmp.rows()); i++)
		for (size_t j = 0; j < static_cast<size_t>(tmp.cols()); j++)
			if (vstr[i * tmp.cols() + j].size() > maxlengthcols[j])
				maxlengthcols[j] = vstr[i * tmp.cols() + j].size();

// finally display it!
	for (size_t i = 0; i < static_cast<size_t>(tmp.rows()); i++)
	{
		os << std::setw(static_cast<int>(maxlengthcols[0])) << std::right
				<< vstr[i * tmp.cols()]; // display first column
		for (size_t j = 1; j < static_cast<size_t>(tmp.cols()); j++) // then the rest
			os << std::setw(static_cast<int>(maxlengthcols[j] + 2))
					<< std::right << vstr[i * tmp.cols() + j];

		if (i < static_cast<size_t>(tmp.rows()) - 1)
			os << std::endl;
	}
}

// Displays an Eigen::MatrixX in friendly form
// Adds new line after display
template<typename Scalar>
void displn(const types::DynMat<Scalar> &A, double chop = ct::chop,
		std::ostream& os = std::cout)
{
	disp(A, chop, os);
	os << std::endl;
}

// Displays a complex number in friendly form
inline void disp(const types::cplx c, double chop = ct::chop, std::ostream& os =
		std::cout)
{
// put the complex number inside an Eigen matrix
	types::cmat tmp(1, 1);
	tmp(0, 0) = c;
	disp(tmp, chop, os);
}

// Displays a complex number in friendly form
// Adds new line after display
inline void displn(const types::cplx c, double chop = ct::chop,
		std::ostream& os = std::cout)
{
	disp(c, chop, os);
	os << std::endl;
}

// save matrix to a binary file in double precision
template<typename Scalar>
void save(const types::DynMat<Scalar> & A, const std::string& fname)
{
	// zero-size
	if (!internal::_check_nonzero_size(A))
		throw std::invalid_argument("save: Zero-sized input!");

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
	fout.write((char *) &rows, sizeof(rows));
	fout.write((char *) &cols, sizeof(cols));

	fout.write((char *) static_cast<Scalar>(A).data(),
			sizeof(types::DynMat<Scalar>) * rows * cols);

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
	char *_fheader = new char[sizeof(_header)];

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
	fin.read((char *) &rows, sizeof(rows));
	fin.read((char *) &cols, sizeof(cols));

	types::DynMat<Scalar> A(rows, cols);

	fin.read((char *) A.data(), sizeof(Scalar) * rows * cols);

	fin.close();
	return A;
}

}

#endif /* IO_H_ */
