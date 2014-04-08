/*
 * random.h
 *
 *  Created on: Mar 27, 2014
 *      Author: vlad
 */

#ifndef RANDOM_H_
#define RANDOM_H_

#include "types.h"
#include "stat.h"
#include "constants.h"
#include "exception.h"

// random matrices/states

namespace qpp
{

// random matrix with entries in Uniform(a,b)
template<typename Scalar>
types::DynMat<Scalar> rand(size_t rows, size_t cols, double a = 0, double b = 1)
{
	throw Exception("rand", Exception::Type::UNDEFINED_TYPE);
}

// random double matrix with entries in Uniform(a,b)
template<>
types::DynMat<double> rand(size_t rows, size_t cols, double a, double b)
{
	if (rows == 0 || cols == 0)
		throw Exception("rand", Exception::Type::ZERO_SIZE);

	return (0.5 * (b - a)
			* (types::dmat::Random(rows, cols) + types::dmat::Ones(rows, cols))
			+ a * types::dmat::Ones(rows, cols));
}

// random complex matrix with entries in Uniform(a,b)
template<>
types::DynMat<types::cplx> rand(size_t rows, size_t cols, double a, double b)
{
	if (rows == 0 || cols == 0)
		throw Exception("rand", Exception::Type::ZERO_SIZE);

	return (0.5 * (b - a)
			* (types::cmat::Random(rows, cols)
					+ (1. + ct::ii) * types::cmat::Ones(rows, cols))
			+ a * types::cmat::Ones(rows, cols));
}

// random number in Uniform(a, b)
double rand(double a = 0, double b = 1)
{
	stat::UniformRealDistribution ud(a, b);
	return ud.sample();
}

// random matrix with entries in Normal(mean, sigma)
template<typename Scalar>
types::DynMat<Scalar> randn(size_t rows, size_t cols, double mean = 0,
		double sigma = 1)
{
	throw Exception("randn", Exception::Type::UNDEFINED_TYPE);
}

// random double matrix with entries in Normal(mean, sigma)
template<>
types::DynMat<double> randn(size_t rows, size_t cols, double mean, double sigma)
{
	if (rows == 0 || cols == 0)
		throw Exception("randn", Exception::Type::ZERO_SIZE);

	stat::NormalDistribution nd(mean, sigma);
	auto lambda = [&nd](double)
	{	return nd.sample();};

	//types::dmat result = types::dmat::Zero(rows, cols);
	return types::dmat::Zero(rows, cols).unaryExpr(lambda);

}

// random complex matrix with entries in Normal(mean, sigma)
template<>
types::DynMat<types::cplx> randn(size_t rows, size_t cols, double mean,
		double sigma)
{
	if (rows == 0 || cols == 0)
		throw Exception("randn", Exception::Type::ZERO_SIZE);

	stat::NormalDistribution nd(mean, sigma);
	auto lambda = [&nd](double)
	{	return nd.sample();};

	return (types::dmat::Zero(rows, cols).unaryExpr(lambda)).cast<types::cplx>()
			+ ct::ii
					* (types::dmat::Zero(rows, cols).unaryExpr(lambda)).cast<
							types::cplx>();
}

// random number in Normal(mean, sigma)
double randn(double mean = 0, double sigma = 1)
{
	stat::NormalDistribution nd(mean, sigma);
	return nd.sample();
}

// Random unitary matrix
//(~3 times slower than Toby Cubitt's MATLAB's, because of the QR algorithm)
types::cmat randU(size_t D)
{
	if (D == 0)
		throw Exception("randU", Exception::Type::DIMS_INVALID);

	types::cmat X(D, D);

	X = 1 / std::sqrt(2) * randn<types::cplx>(D, D);
	Eigen::HouseholderQR<types::cmat> qr(X);

	types::cmat Q = qr.householderQ();
	// phase correction so that the resultant matrix is
	// uniformly distributed according to the Haar measure

	Eigen::VectorXcd phases = (rand<double>(D, 1)).cast<types::cplx>();
	for (size_t i = 0; i < static_cast<size_t>(phases.rows()); i++)
		phases(i) = std::exp((types::cplx) (2 * ct::pi * ct::ii * phases(i)));

	Q = Q * phases.asDiagonal();

	return Q;
}

// Random isometry
types::cmat randV(size_t Din, size_t Dout)
{
	if (Din == 0 || Dout == 0 || Din > Dout)
		throw Exception("randV", Exception::Type::DIMS_INVALID);
	return randU(Dout).block(0, 0, Dout, Din);
}

// Random Kraus operators
std::vector<types::cmat> randKraus(size_t n, size_t D)
{
	if (n == 0)
		throw Exception("randKraus", Exception::Type::OUT_OF_RANGE);
	if (D == 0)
		throw Exception("randKraus", Exception::Type::DIMS_INVALID);

	std::vector<types::cmat> result;
	types::cmat Fk(D, D);
	types::cmat U = randU(n * D);
	size_t dims[2];
	dims[0] = D;
	dims[1] = n;
	size_t midx_row[2] = { 0, 0 };
	size_t midx_col[2] = { 0, 0 };

	for (size_t k = 0; k < n; k++)
	{
		for (size_t a = 0; a < D; a++)
			for (size_t b = 0; b < D; b++)
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
types::cmat randH(size_t D)
{
	if (D == 0)
		throw Exception("randH", Exception::Type::DIMS_INVALID);

	types::cmat H = 2 * rand<types::cplx>(D, D)
			- (1. + ct::ii) * types::cmat::Ones(D, D);

	return H + adjoint(H);
}

// random ket
types::cmat randket(size_t D)
{
	if (D == 0)
		throw Exception("randket", Exception::Type::DIMS_INVALID);

	types::cmat kt = types::cmat::Ones(D, 1);
	types::cmat result = randU(D) * kt;
	return result / norm(result);
}

// random density matrix
types::cmat randrho(size_t D)
{
	if (D == 0)
		throw Exception("randrho", Exception::Type::DIMS_INVALID);
	types::cmat result = 10 * randH(D);
	result = result * adjoint(result);
	return result / trace(result);
}

} /* namespace qpp */

#endif /* RANDOM_H_ */
