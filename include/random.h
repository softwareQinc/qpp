/*
 * random.h
 *
 *  Created on: Mar 27, 2014
 *      Author: vlad
 */

#ifndef RANDOM_H_
#define RANDOM_H_

// random matrices/states

namespace qpp
{

/**
 * \brief Generates an Eigen random dynamic matrix with entries uniformly
 * distributed in the interval [a,b)
 *
 * If complex, then both real and imaginary parts are uniformly distributed
 * in [a,b)
 *
 * This is the generic version that always throws
 * \a qpp::Exception::Type::UNDEFINED_TYPE. It is specialized only for
 * \a qpp::dmat and \a qpp::cmat
 */
template<typename Derived>
Derived rand(std::size_t rows, std::size_t cols, double a = 0, double b = 1)
{
	throw Exception("rand", Exception::Type::UNDEFINED_TYPE);
}

/**
 * \brief Generates an Eigen random dynamic matrix with entries uniformly
 * distributed in the interval [a,b),
 * specialization for double matrices (\a qpp::dmat)
 *
 * The template parameter cannot be automatically deduced and
 * must be explicitly provided
 *
 * Example:
 * \code
 * // generates a 3 x 3 Eigen random dynamic double matrix, with entries uniformly distributed in [-1,1)
 * auto mat = rand<dmat>(3, 3, -1, 1);
 * \endcode
 *
 * \param rows Number of rows of the random generated matrix
 * \param cols Number of columns of the random generated matrix
 * \param a
 * \param b
 * \return Eigen random double dynamic matrix (\a qpp::dmat)
 */
template<>
dmat rand(std::size_t rows, std::size_t cols, double a, double b)
{
	if (rows == 0 || cols == 0)
		throw Exception("rand", Exception::Type::ZERO_SIZE);

	return (0.5 * (b - a)
			* (dmat::Random(rows, cols) + dmat::Ones(rows, cols))
			+ a * dmat::Ones(rows, cols));
}

/**
 * \brief Generates an Eigen random dynamic matrix with entries uniformly
 * distributed in the interval [a,b),
 * specialization for complex matrices (\a qpp::cmat)
 *
 * The template parameter cannot be automatically deduced and
 * must be explicitly provided
 *
* Example:
 * \code
 * // generates a 3 x 3 Eigen random dynamic complex matrix, with entries (both real and imaginary) uniformly distributed in [-1,1)
 * auto mat = rand<cmat>(3, 3, -1, 1);
 * \endcode
 *
 * Both the real part and imaginary part of the entries of the resulting random
 * matrix are uniformly distributed in the interval [a,b).
 *
 * \param rows Number of rows of the random generated matrix
 * \param cols Number of columns of the random generated matrix
 * \param a
 * \param b
 * \return Eigen random double dynamic matrix (\a qpp::dmat)
 */
template<>
cmat rand(std::size_t rows, std::size_t cols, double a, double b)
{
	if (rows == 0 || cols == 0)
		throw Exception("rand", Exception::Type::ZERO_SIZE);

	return rand < dmat
			> (rows, cols, a, b).cast<cplx>() + 1_i * rand < dmat
			> (rows, cols, a, b).cast<cplx>();
}

/**
 * \brief Generates a random double uniformly distributed in
 * the interval [a,b)
 * \param a
 * \param b
 * \return Random real number uniformly distributed in
 * the interval [a,b)
 */
double rand(double a = 0, double b = 1)
{
	UniformRealDistribution<> ud(a, b);
	return ud.sample();
}

/**
 * \brief Generates a random long long integer uniformly distributed in
 * the interval [a,b]
 * \param a
 * \param b
 * \return Random long long integer uniformly distributed in
 * the interval [a,b]
 */
long long int randint(long long a = std::numeric_limits<long long int>::min(),\
		long long b = std::numeric_limits<long long int>::max())
{
	UniformIntegerDistribution<long long int> ud(a, b);

	return ud.sample();
}

// random matrix with entries in Normal(mean, sigma)
template<typename Derived>
Derived randn(std::size_t rows, std::size_t cols, double mean = 0,
		double sigma = 1)
{
	throw Exception("randn", Exception::Type::UNDEFINED_TYPE);
}

// random double matrix with entries in Normal(mean, sigma)
template<>
dmat randn(std::size_t rows, std::size_t cols, double mean, double sigma)
{
	if (rows == 0 || cols == 0)
		throw Exception("randn", Exception::Type::ZERO_SIZE);

	NormalDistribution<> nd(mean, sigma);

	return dmat::Zero(rows, cols).unaryExpr([&nd](double)
	{	return nd.sample();});

}

// random complex matrix with entries in Normal(mean, sigma)
template<>
cmat randn(std::size_t rows, std::size_t cols, double mean, double sigma)
{
	if (rows == 0 || cols == 0)
		throw Exception("randn", Exception::Type::ZERO_SIZE);

	NormalDistribution<> nd(mean, sigma);
	return randn < dmat
			> (rows, cols, mean, sigma).cast<cplx>() + 1_i * randn
			< dmat > (rows, cols, mean, sigma).cast<cplx>();
}

// random number in Normal(mean, sigma)
double randn(double mean = 0, double sigma = 1)
{
	NormalDistribution<> nd(mean, sigma);
	return nd.sample();
}

// Random unitary matrix
// ~3 times slower than Toby Cubitt's MATLAB's,
// because Eigen's QR algorithm is not parallelized
cmat randU(std::size_t D)
{
	if (D == 0)
		throw Exception("randU", Exception::Type::DIMS_INVALID);

	cmat X(D, D);

	X = 1 / std::sqrt(2.) * randn<cmat>(D, D);
	Eigen::HouseholderQR<cmat> qr(X);

	cmat Q = qr.householderQ();
	// phase correction so that the resultant matrix is
	// uniformly distributed according to the Haar measure

	Eigen::VectorXcd phases = (rand<dmat>(D, 1)).cast<cplx>();
	for (std::size_t i = 0; i < static_cast<std::size_t>(phases.rows()); i++)
		phases(i) = std::exp((cplx) (2 * pi * 1_i * phases(i)));

	Q = Q * phases.asDiagonal();

	return Q;
}

// Random isometry
cmat randV(std::size_t Din, std::size_t Dout)
{
	if (Din == 0 || Dout == 0 || Din > Dout)
		throw Exception("randV", Exception::Type::DIMS_INVALID);
	return randU(Dout).block(0, 0, Dout, Din);
}

// Random Kraus operators
std::vector<cmat> randkraus(std::size_t n, std::size_t D)
{
	if (n == 0)
		throw Exception("randkraus", Exception::Type::OUT_OF_RANGE);
	if (D == 0)
		throw Exception("randkraus", Exception::Type::DIMS_INVALID);

	std::vector<cmat> result;
	cmat Fk(D, D);
	cmat U = randU(n * D);
	std::size_t dims[2];
	dims[0] = D;
	dims[1] = n;
	std::size_t midx_row[2] = { 0, 0 };
	std::size_t midx_col[2] = { 0, 0 };

	for (std::size_t k = 0; k < n; k++)
	{
		for (std::size_t a = 0; a < D; a++)
			for (std::size_t b = 0; b < D; b++)
			{
				midx_row[0] = a;
				midx_row[1] = k;
				midx_col[0] = b;
				midx_col[1] = 0;
				Fk(a, b) = U(internal::_multiidx2n(midx_row, 2, dims),
						internal::_multiidx2n(midx_col, 2, dims));
			}
		result.push_back(Fk);
	}

	return result;
}

// Random Hermitian matrix
cmat randH(std::size_t D)
{
	if (D == 0)
		throw Exception("randH", Exception::Type::DIMS_INVALID);

	cmat H = 2 * rand<cmat>(D, D)
			- (1. + 1_i) * cmat::Ones(D, D);

	return H + adjoint(H);
}

// random ket of dimension D why randU() ? and not N(0,1)?
ket randket(std::size_t D)
{
	if (D == 0)
		throw Exception("randket", Exception::Type::DIMS_INVALID);

	ket kt = ket::Ones(D);
	ket result = static_cast<ket>(randU(D) * kt);
	return result / norm(result);
}

// random density matrix
cmat randrho(std::size_t D)
{
	if (D == 0)
		throw Exception("randrho", Exception::Type::DIMS_INVALID);
	cmat result = 10 * randH(D);
	result = result * adjoint(result);
	return result / trace(result);
}

// random permutation (using Knuth shuffle method)
std::vector<std::size_t> randperm(std::size_t n)
{
	if (n == 0)
		throw Exception("randperm", Exception::Type::PERM_INVALID);

	std::vector<std::size_t> result(n);

	// fill in increasing order
	std::iota(std::begin(result), std::end(result), 0);
	// shuffle
	std::shuffle(std::begin(result), std::end(result),
			RandomDevices::get_instance()._rng);

	return result;
}

} /* namespace qpp */

#endif /* RANDOM_H_ */
