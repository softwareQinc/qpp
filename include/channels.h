/*
 * channels.h
 *
 *  Created on: Apr 6, 2014
 *      Author: vlad
 */

#ifndef CHANNELS_H_
#define CHANNELS_H_

#include <vector>
#include <cmath>
#include "types.h"
#include "internal.h"
#include "exception.h"
#include "constants.h"

namespace qpp
{
// channel specified by Kraus operators {Ks} acting on rho
types::cmat channel(const types::cmat& rho, const std::vector<types::cmat>& Ks)
{
	if (!internal::_check_nonzero_size(rho))
		throw Exception("channel", Exception::Type::ZERO_SIZE);
	if (!internal::_check_square_mat(rho))
		throw Exception("channel", Exception::Type::MATRIX_NOT_SQUARE);
	if (!internal::_check_nonzero_size(Ks))
		throw Exception("channel", Exception::Type::ZERO_SIZE);
	if (!internal::_check_square_mat(Ks[0]))
		throw Exception("channel", Exception::Type::MATRIX_NOT_SQUARE);
	if (Ks[0].rows() != rho.rows())
		throw Exception("channel", Exception::Type::DIMS_NOT_EQUAL);
	for (auto it : Ks)
		if (it.rows() != Ks[0].rows() || it.cols() != Ks[0].rows())
			throw Exception("channel", Exception::Type::DIMS_NOT_EQUAL);

	types::cmat result = types::cmat::Zero(rho.rows(), rho.cols());
	for (auto it : Ks)
		result += it * rho * adjoint(it);

	return result;
}

// returns the Choi matrix of the channel specified by Kraus operators {Ks}
types::cmat choi(const std::vector<types::cmat>& Ks)
{
	if (!internal::_check_nonzero_size(Ks))
		throw Exception("choi", Exception::Type::ZERO_SIZE);
	if (!internal::_check_nonzero_size(Ks[0]))
		throw Exception("choi", Exception::Type::ZERO_SIZE);
	if (!internal::_check_square_mat(Ks[0]))
		throw Exception("choi", Exception::Type::MATRIX_NOT_SQUARE);
	for (auto it : Ks)
		if (it.rows() != Ks[0].rows() || it.cols() != Ks[0].rows())
			throw Exception("channel", Exception::Type::DIMS_NOT_EQUAL);
	size_t D = static_cast<size_t>(Ks[0].rows());

	// construct the D x D \sum |jj> vector
	// (un-normalized maximally entangled state)
	types::cmat MES = types::cmat::Zero(D * D, 1);
	for (size_t a = 0; a < D; a++)
		MES(a * D + a) = 1;

	types::cmat Omega = MES * adjoint(MES);

	types::cmat result = types::cmat::Zero(D * D, D * D);
	for (auto it : Ks)
		result += kron((types::cmat) (types::cmat::Identity(D, D)), it) * Omega
				* adjoint(
						kron((types::cmat) (types::cmat::Identity(D, D)), it));

	return result;
}

// extract orthogonal Kraus operators from Choi matrix
std::vector<types::cmat> choi2kraus(const types::cmat& A)
{
	if (!internal::_check_nonzero_size(A))
		throw Exception("choi2kraus", Exception::Type::ZERO_SIZE);
	if (!internal::_check_square_mat(A))
		throw Exception("choi2kraus", Exception::Type::MATRIX_NOT_SQUARE);
	size_t D = static_cast<size_t>(std::sqrt(A.rows()));
	if (D * D != static_cast<size_t>(A.rows()))
		throw Exception("choi2kraus", Exception::Type::DIMS_MISMATCH_MATRIX);

	types::cmat eval = hevals(A);
	types::cmat evec = hevects(A);
	std::vector<types::cmat> result;

	types::cmat vi = types::cmat::Zero(D * D, 1);
	double ei = 0;
	types::cmat Ki = types::cmat::Zero(D, D);

	for (size_t i = 0; i < D * D; i++)
	{
		ei = eval.real()(i);
		if (std::abs(ei) > ct::eps)
		{
			vi = reshape(static_cast<types::cmat>(evec.col(i)), D, D);
			Ki = std::sqrt(ei) * vi;
			result.push_back(Ki);
		}
	}
	return result;
}
}

#endif /* CHANNELS_H_ */
