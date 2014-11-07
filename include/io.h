/*
 * Quantum++
 *
 * Copyright (c) 2013 - 2014 Vlad Gheorghiu (vgheorgh@gmail.com)
 *
 * This file is part of Quantum++.
 *
 * Quantum++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Quantum++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Quantum++.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef INCLUDE_IO_H_
#define INCLUDE_IO_H_

// input/output

namespace qpp
{


/**
 * \brief Eigen expression or complex number ostream manipulator.
 *
 * \param A Eigen expression
 * \param chop Set to zero the elements smaller in absolute value
 * than \a chop
 * \return qpp::internal::IOManip
 */
template<typename Derived>
IOManip<char, std::vector<char>::iterator> disp(
		const Eigen::MatrixBase<Derived>& A, double chop = qpp::chop)
{
	return IOManip<char, std::vector<char>::iterator>(A, chop);
}

/**
 * \brief Range ostream manipulator
 *
 * \param first Iterator to the first element of the range
 * \param last  Iterator to the last element of the range
 * \param separator Separator
 * \param start Left marking
 * \param end Right marking
 * \return qpp::internal::IOManip
 */
template<typename InputIterator>
IOManip<char, InputIterator> disp(const InputIterator& first,
		const InputIterator& last, const std::string & separator,
		const std::string& start = "[", const std::string& end = "]")
{
	return IOManip<char, InputIterator>(first, last, separator, start, end);
}

/**
 * \brief Standard container ostream manipulator. The container must support
 * std::begin, std::end and forward iteration.
 *
 * \param x Container
 * \param separator Separator
 * \param start Left marking
 * \param end Right marking
 * \return qpp::internal::IOManip
 */
template<typename T>
IOManip<char, typename T::const_iterator> disp(const T& x,
		const std::string & separator, const std::string& start = "[",
		const std::string& end = "]")
{
	return IOManip<char, typename T::const_iterator>(x.begin(), x.end(),
			separator, start, end);
}

/**
 * \brief C-style string ostream manipulator
 *
 * \see \a qpp::displn()
 *
 * \param x Pointer to the first element
 * \param n Number of elements to be displayed
 * \param separator Separator
 * \param start Left marking
 * \param end Right marking
 * \return qpp::internal::IOManip
 */
template<typename T>
IOManip<T, std::vector<char>::iterator> disp(const T* p, std::size_t n,
		const std::string & separator, const std::string& start = "[",
		const std::string& end = "]")
{
	return IOManip<T, std::vector<char>::iterator>(p, n, separator, start, end);
}

/**
 * \brief Saves Eigen expression to a binary file (internal format) in double
 * precision
 *
 * \see qpp::saveMATLABmatrix()
 *
 * \param A Eigen expression
 * \param fname Output file name
 */
template<typename Derived>
void save(const Eigen::MatrixBase<Derived>& A, const std::string& fname)

{
	const DynMat<typename Derived::Scalar> & rA = A;

	// check zero-size
	if (rA.size() == 0)
		throw Exception("save", Exception::Type::ZERO_SIZE);

	std::fstream fout;
	fout.open(fname, std::ios::out | std::ios::binary);

	if (fout.fail())
	{
		throw std::runtime_error(
				"save: Error writing output file \"" + std::string(fname)
						+ "\"!");
	}

	// write the header to file
	const char _header[] = "TYPE::Eigen::Matrix";
	fout.write(_header, sizeof(_header));

	std::size_t rows = static_cast<std::size_t>(rA.rows());
	std::size_t cols = static_cast<std::size_t>(rA.cols());
	fout.write((char*) &rows, sizeof(rows));
	fout.write((char*) &cols, sizeof(cols));

	fout.write((char*) rA.data(), sizeof(Derived::Scalar) * rows * cols);

	fout.close();
}

/**
 * \brief Loads Eigen matrix from a binary file (internal format) in double
 * precision
 *
 * The template parameter cannot be automatically deduced and
 * must be explicitly provided, depending on the scalar field of the matrix
 * that is being loaded.
 *
 * Example:
 * \code
 * // loads a previously saved Eigen dynamic complex matrix from "input.bin"
 * auto mat = load<cmat>("input.bin");
 * \endcode
 *
 * \see qpp::loadMATLABmatrix()
 *
 * \param A Eigen expression
 * \param fname Output file name
 */
template<typename Derived>
DynMat<typename Derived::Scalar> load(const std::string& fname)
{
	std::fstream fin;
	fin.open(fname, std::ios::in | std::ios::binary);

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

	std::size_t rows, cols;
	fin.read((char*) &rows, sizeof(rows));
	fin.read((char*) &cols, sizeof(cols));

	DynMat<typename Derived::Scalar> A(rows, cols);

	fin.read((char*) A.data(), sizeof(typename Derived::Scalar) * rows * cols);

	fin.close();
	return A;
}

} /* namespace qpp */

#endif /* INCLUDE_IO_H_ */
