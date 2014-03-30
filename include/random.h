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

namespace qpp
{

// Random double matrix with entries in Uniform[0,1]
inline types::dmat rand(size_t rows, size_t cols)
{
	return types::dmat::Random(rows, cols);
}

// Random double square matrix with entries in Uniform[0,1]
inline types::dmat rand(size_t rows)
{
	return rand(rows, rows);
}

// Random double matrix with entries in Normal(0,1)
inline types::dmat randn(size_t rows, size_t cols)
{
	stat::NormalDistribution nd; // N(0,1)
	types::dmat A(rows, cols);

	for (size_t i = 0; i < rows; i++)
		for (size_t j = 0; j < cols; j++)
			A(i, j) = nd.sample();
	return A;
}

// Random double square matrix with entries in Normal(0,1)
inline types::dmat randn(size_t rows)
{
	return randn(rows, rows);
}

// Random unitary matrix
inline types::cmat rand_unitary(size_t size)
{
	types::cmat X(size, size);

	X.real() = 1. / sqrt(2) * randn(size);
	X.imag() = 1. / sqrt(2) * randn(size);
	Eigen::HouseholderQR<types::cmat> qr(X);

	types::cmat Q = qr.householderQ();
	// phase correction so that the resultant matrix is
	// uniformly distributed according to the Haar measure

	Eigen::VectorXcd phases = (rand(size, 1)).cast<types::cplx>();
	for(size_t i=0; i<(size_t)phases.rows(); i++)
		phases(i)=std::exp((types::cplx)(2*ct::pi*ct::ii*phases(i)));

	Q = Q * phases.asDiagonal();

	return Q;
}

}

#endif /* RANDOM_H_ */
