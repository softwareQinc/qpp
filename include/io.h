/*
 * io.h
 *
 *  Created on: Mar 27, 2014
 *      Author: vlad
 */

#ifndef IO_H_
#define IO_H_

// input/output

namespace qpp
{

/**
 * \brief Displays a standard container that supports std::begin, std::end
 * and forward iteration. Does not add a newline.
 *
 * \see \a qpp::displn()
 *
 * \param x Container
 * \param separator Separator
 * \param start Left marking
 * \param end Right marking
 * \param os Output stream
 */
template<typename T>
void disp(const T& x, const std::string & separator, const std::string& start =
		"[", const std::string& end = "]", std::ostream& os = std::cout)
{
	os << start;

	auto it = std::begin(x);
	auto it_end = std::end(x);

	if (it != it_end)
	{
		// the iterator just before the end, need this for containers
		// that do not have backwards iterators
		decltype(it_end) it_before_end = it;
		while (it_before_end = it, ++it != it_end)
			;

		it = std::begin(x);
		for (; it != it_before_end; ++it)
			os << *it << separator;
		os << *it;
	}

	os << end;
}

/**
 * \brief Displays a standard container that supports std::begin, std::end
 * and forward iteration. Adds a newline.
 *
 * \see \a qpp::disp()
 *
 * \param x Container
 * \param separator Separator
 * \param start Left marking
 * \param end Right marking
 * \param os Output stream
 */
template<typename T>
void displn(const T& x, const std::string & separator,
		const std::string& start = "[", const std::string& end = "]",
		std::ostream& os = std::cout)
{
	disp(x, separator, start, end, os);
	os << std::endl;
}

/**
 * \brief Displays a C-style array. Does not add a newline.
 *
 * \see \a qpp::displn()
 *
 * \param x Pointer to the first element
 * \param n Number of elements to be displayed
 * \param separator Separator
 * \param start Left marking
 * \param end Right marking
 * \param os Output stream
 */
template<typename T>
void disp(const T* x, const std::size_t n, const std::string& separator,
		const std::string& start = "[", const std::string& end = "]",
		std::ostream& os = std::cout)
{
	os << start;

	for (std::size_t i = 0; i < n - 1; i++)
		os << x[i] << separator;
	if (n > 0)
		os << x[n - 1];

	os << end;
}

/**
 * \brief Displays a C-style array. Adds a newline.
 *
 * \see \a qpp::disp()
 *
 * \param x Pointer to the first element
 * \param n Number of elements to be displayed
 * \param separator Separator
 * \param start Left marking
 * \param end Right marking
 * \param os Output stream
 */
template<typename T>
void displn(const T* x, const std::size_t n, const std::string & separator,
		const std::string& start = "[", const std::string& end = "]",
		std::ostream& os = std::cout)
{
	disp(x, n, separator, start, end, os);
	os << std::endl;
}

/**
 * \brief Displays an Eigen expression in matrix friendly form. Does not add a
 * new line.
 *
 * \see \a qpp::displn()
 *
 * \param A Eigen expression
 * \param chop Set to zero the elements smaller in absolute value
 * than \a chop
 * \param os Output stream
 */
template<typename Derived>
void disp(const Eigen::MatrixBase<Derived>& A, double chop = qpp::chop,
		std::ostream& os = std::cout)
{
	const DynMat<typename Derived::Scalar> & rA = A;

	if (rA.size() == 0)
	{
		os << "Empty [" << rA.rows() << " x " << rA.cols() << "] matrix";
		return;
	};

	std::ostringstream ostr{};
	ostr.flags(os.flags()); // get the formatting flags
	ostr.precision(os.precision()); // set precision

	std::vector<std::string> vstr{};
	std::string strA{};

	for (std::size_t i = 0; i < static_cast<std::size_t>(rA.rows()); i++)
	{
		for (std::size_t j = 0; j < static_cast<std::size_t>(rA.cols()); j++)
		{
			strA.clear(); // clear the temporary string
			ostr.clear();
			ostr.str(std::string{}); // clear the ostringstream

			// convert to complex
			double re = static_cast<cplx>(rA(i, j)).real();
			double im = static_cast<cplx>(rA(i, j)).imag();

			if (std::abs(re) < chop && std::abs(im) < chop)
			{
				ostr << "0 "; // otherwise segfault on destruction
							  // if using only vstr.push_back("0 ");
							  // bug in MATLAB's libmx
				vstr.push_back(ostr.str());
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
	std::vector < std::size_t > maxlengthcols(rA.cols(), 0);

	for (std::size_t i = 0; i < static_cast<std::size_t>(rA.rows()); i++)
		for (std::size_t j = 0; j < static_cast<std::size_t>(rA.cols()); j++)
			if (vstr[i * rA.cols() + j].size() > maxlengthcols[j])
				maxlengthcols[j] = vstr[i * rA.cols() + j].size();

	// finally display it!
	for (std::size_t i = 0; i < static_cast<std::size_t>(rA.rows()); i++)
	{
		os << std::setw(static_cast<int>(maxlengthcols[0])) << std::right
				<< vstr[i * rA.cols()]; // display first column
		// then the rest
		for (std::size_t j = 1; j < static_cast<std::size_t>(rA.cols()); j++)
			os << std::setw(static_cast<int>(maxlengthcols[j] + 2))
					<< std::right << vstr[i * rA.cols() + j];

		if (i < static_cast<std::size_t>(rA.rows()) - 1)
			os << std::endl;
	}
}

/**
 * \brief Displays an Eigen expression in matrix friendly form. Adds a newline.
 *
 * \see \a qpp::disp()
 *
 * \param A Eigen expression
 * \param chop Set to zero the elements smaller in absolute value
 * than \a chop
 * \param os Output stream
 */
template<typename Derived>
void displn(const Eigen::MatrixBase<Derived>& A, double chop = qpp::chop,
		std::ostream& os = std::cout)
{
	disp(A, chop, os);
	os << std::endl;
}

/**
 * \brief Displays a number (implicitly converted to std::complex<double>)
 * in friendly form. Does not add a new line.
 *
 * \see \a qpp::displn()
 *
 * \param z Real/complex number
 * \param chop Set to zero the elements smaller in absolute value
 * than \a chop
 * \param os Output stream
 */
void disp(const cplx z, double chop = qpp::chop, std::ostream& os =
		std::cout)
{
// put the complex number inside an Eigen matrix
	cmat A(1, 1);
	A(0, 0) = z;
	disp(A, chop, os);
}


/**
 * \brief Displays a number (implicitly converted to std::complex<double>)
 * in friendly form. Adds a new line.
 *
 * \see \a qpp::disp()
 *
 * \param z Real/complex number
 * \param chop Set to zero the elements smaller in absolute value
 * than \a chop
 * \param os Output stream
 */
void displn(const cplx z, double chop = qpp::chop, std::ostream& os =
		std::cout)
{
	disp(z, chop, os);
	os << std::endl;
}

/**
 * \brief Saves Eigen expression to a binary file (internal format) in double
 * precission
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

	std::size_t rows = static_cast<std::size_t>(rA.rows());
	std::size_t cols = static_cast<std::size_t>(rA.cols());
	fout.write((char*) &rows, sizeof(rows));
	fout.write((char*) &cols, sizeof(cols));

	fout.write((char*) rA.data(), sizeof(Derived::Scalar) * rows * cols);

	fout.close();
}

/**
 * \brief Loads Eigen matrix from a binary file (internal format) in double
 * precission
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

	std::size_t rows, cols;
	fin.read((char*) &rows, sizeof(rows));
	fin.read((char*) &cols, sizeof(cols));

	DynMat<typename Derived::Scalar> A(rows, cols);

	fin.read((char*) A.data(), sizeof(typename Derived::Scalar) * rows * cols);

	fin.close();
	return A;
}

} /* namespace qpp */

#endif /* IO_H_ */
