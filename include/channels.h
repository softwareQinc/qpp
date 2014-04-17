/*
 * channels.h
 *
 *  Created on: Apr 6, 2014
 *      Author: vlad
 */

#ifndef CHANNELS_H_
#define CHANNELS_H_

#include <cmath>
#include <vector>
#include "constants.h"
#include "functions.h"
#include "internal.h"
#include "types.h"
#include "classes/exception.h"

namespace qpp
{
// the output of a channel specified by Kraus operators {Ks} acting on rho
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

// constructs the superoperator matrix in the standard operator basis |i><j|,
// ordered in lexicographical order, e.g. |0><0|, |0><1|, |1><0|, |1><1|
types::cmat super(const std::vector<types::cmat>& Ks)
{
	if (!internal::_check_nonzero_size(Ks))
		throw Exception("super", Exception::Type::ZERO_SIZE);
	if (!internal::_check_nonzero_size(Ks[0]))
		throw Exception("super", Exception::Type::ZERO_SIZE);
	if (!internal::_check_square_mat(Ks[0]))
		throw Exception("super", Exception::Type::MATRIX_NOT_SQUARE);
	for (auto it : Ks)
		if (it.rows() != Ks[0].rows() || it.cols() != Ks[0].rows())
			throw Exception("super", Exception::Type::DIMS_NOT_EQUAL);
	size_t D = static_cast<size_t>(Ks[0].rows());

	size_t midx_row[2] = { 0, 0 };
	size_t midx_col[2] = { 0, 0 };
	size_t dims[2];
	dims[0] = dims[1] = D;

	types::cmat result(D * D, D * D);
	types::cmat MN = types::cmat::Zero(D, D);
	types::cmat BA = types::cmat::Zero(D, D);
	types::cmat EMN = types::cmat::Zero(D, D);

	for (size_t m = 0; m < D; m++)
	{
		midx_col[0] = m;
		for (size_t n = 0; n < D; n++)
		{
			midx_col[1] = n;
			MN(m, n) = 1;
			// compute E(|m><n|)
			for (size_t i = 0; i < Ks.size(); i++)
				EMN += Ks[i] * MN * adjoint(Ks[i]);
			MN(m, n) = 0;
			for (size_t a = 0; a < D; a++)
			{
				midx_row[0] = a;
				for (size_t b = 0; b < D; b++)
				{
					midx_row[1] = b;
					BA(b, a) = 1;

					// compute result(ab,mn)=<a|E(|m><n)|b>

					result(internal::_multiidx2n(midx_row, 2, dims),
							internal::_multiidx2n(midx_col, 2, dims)) = (EMN
							* BA).trace();
					BA(b, a) = 0;
				}
			}
			EMN = types::cmat::Zero(D, D);
		}
	}
	return result;
}

// returns the Choi matrix (or dynamical matrix) of the channel
// specified by Kraus operators {Ks}
// The superoperator L and Choi matrix C are related by
// L_{ab,mn} = C_{ma,nb}
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
			throw Exception("choi", Exception::Type::DIMS_NOT_EQUAL);
	size_t D = static_cast<size_t>(Ks[0].rows());

	// construct the D x D \sum |jj> vector
	// (un-normalized maximally entangled state)
	types::cmat MES = types::cmat::Zero(D * D, 1);
	for (size_t a = 0; a < D; a++)
		MES(a * D + a) = 1;

	types::cmat Omega = MES * adjoint(MES);

	types::cmat result = types::cmat::Zero(D * D, D * D);
	for (auto it : Ks)
		result += kron(types::cmat::Identity(D, D), it) * Omega
				* adjoint(kron(types::cmat::Identity(D, D), it));

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
		throw Exception("choi2kraus", Exception::Type::DIMS_INVALID);

	types::cmat eval = hevals(A);
	types::cmat evec = hevects(A);
	std::vector<types::cmat> result;

	for (size_t i = 0; i < D * D; i++)
	{
		if (std::abs(eval.real()(i)) > ct::eps)
			result.push_back(
					std::sqrt(eval.real()(i)) * reshape(evec.col(i), D, D));
	}
	return result;
}

} /* namespace qpp */

#endif /* CHANNELS_H_ */
