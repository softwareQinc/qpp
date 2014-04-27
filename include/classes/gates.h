/*
 * gates.h
 *
 *  Created on: Apr 7, 2014
 *      Author: vlad
 */

#ifndef GATES_H_
#define GATES_H_

#include <algorithm>

#include "constants.h"
#include "functions.h"
#include "exception.h"
#include "internal.h"
#include "types.h"

namespace qpp
{

class Gates
{
public:
	// one qubit gates
	types::cmat Id2; // Identity matrix
	types::cmat H; // Hadamard matrix
	types::cmat X; // X matrix
	types::cmat Y; // Y matrix
	types::cmat Z; // Z matrix
	types::cmat S; // S gate
	types::cmat T; // T gate

	// two qubit gates
	types::cmat CNOTab; // CNOT ctrl1 target2
	types::cmat CZ; // Controlled-Phase (Controlled-Z)
	types::cmat CNOTba; // CNOT ctrl2 target1
	types::cmat SWAP; // SWAP gate

	// three qubit gates
	types::cmat TOF; // Toffoli
	types::cmat FRED; // Fredkin
private:
	Gates() :
			Id2(types::cmat::Identity(2, 2)), //
			H(types::cmat::Zero(2, 2)), //
			X(types::cmat::Zero(2, 2)), //
			Y(types::cmat::Zero(2, 2)), //
			Z(types::cmat::Zero(2, 2)), //
			S(types::cmat::Zero(2, 2)), //
			T(types::cmat::Zero(2, 2)), //
			CNOTab(types::cmat::Identity(4, 4)), //
			CZ(types::cmat::Identity(4, 4)), //
			CNOTba(types::cmat::Zero(4, 4)), //
			SWAP(types::cmat::Identity(4, 4)), //
			TOF(types::cmat::Identity(8, 8)), //
			FRED(types::cmat::Identity(8, 8)) //
	{
		// initialize the constants and gates
		H << 1 / std::sqrt(2), 1 / std::sqrt(2), 1 / std::sqrt(2), -1
				/ std::sqrt(2);
		X << 0, 1, 1, 0;
		Z << 1, 0, 0, -1;
		Y << 0, -ct::ii, ct::ii, 0;
		S << 1, 0, 0, ct::ii;
		T << 1, 0, 0, std::exp(ct::ii * ct::pi / 4.0);
		CNOTab.block(2, 2, 2, 2) = X;
		CNOTba(0, 0) = 1;
		CNOTba(1, 3) = 1;
		CNOTba(2, 2) = 1;
		CNOTba(3, 1) = 1;
		CZ(3, 3) = -1;

		SWAP.block(1, 1, 2, 2) = X;
		TOF.block(6, 6, 2, 2) = X;
		FRED.block(4, 4, 4, 4) = SWAP;
	}
public:
	Gates(const Gates&) = delete;
	Gates& operator=(const Gates&) = delete;

	static const Gates& getInstance() // const singleton
	{
		static Gates instance; // Guaranteed to be destroyed.
							   // Instantiated on first use.
							   // Thread safe in C++11
		return instance;
	}
	virtual ~Gates() = default;
public:
	// gates with variable dimension

	// one qubit gates
	types::cmat Rtheta(double theta) const
	{
		types::cmat result(2, 2);
		result << 1, 0, 0, exp(ct::ii * theta);
		return result;
	}

	// one quDit gates

	types::cmat Id(size_t D) const
	{
		if (D == 0)
			throw Exception("Id", Exception::Type::DIMS_INVALID);
		return types::cmat::Identity(D, D);
	}

	types::cmat Zd(size_t D) const
	{
		if (D == 0)
			throw Exception("Zd", Exception::Type::DIMS_INVALID);

		types::cmat result(D, D);
		result = types::cmat::Zero(D, D);
		for (size_t i = 0; i < D; i++)
			result(i, i) = std::pow(ct::omega(D), i);
		return result;
	}

	types::cmat Fd(size_t D) const
	{
		if (D == 0)
			throw Exception("Fd", Exception::Type::DIMS_INVALID);

		types::cmat result(D, D);
		result = types::cmat::Zero(D, D);
		for (size_t j = 0; j < D; j++)
			for (size_t i = 0; i < D; i++)
				result(i, j) = 1 / std::sqrt(D) * std::pow(ct::omega(D), i * j);
		return result;
	}

	types::cmat Xd(size_t D) const // X|k>=|k+1>
	{
		if (D == 0)
			throw Exception("Xd", Exception::Type::DIMS_INVALID);

		return Fd(D).inverse() * Zd(D) * Fd(D);
	}

	// -multi-quDit multi-controlled-gate
	types::cmat CTRL(const types::cmat& A, const std::vector<size_t>& ctrl,
			const std::vector<size_t>& subsys, size_t n, size_t d = 2) const
	{
		// EXCEPTION CHECKS
		// check matrix zero size
		if (!internal::_check_nonzero_size(A))
			throw Exception("CTRL", Exception::Type::ZERO_SIZE);

		// check square matrix
		if (!internal::_check_square_mat(A))
			throw Exception("CTRL", Exception::Type::MATRIX_NOT_SQUARE);

		// check lists zero size
		if (ctrl.size() == 0)
			throw Exception("CTRL", Exception::Type::ZERO_SIZE);
		if (subsys.size() == 0)
			throw Exception("CTRL", Exception::Type::ZERO_SIZE);

		// check out of range
		if (n == 0)
			throw Exception("CTRL", Exception::Type::OUT_OF_RANGE);

		// check valid local dimension
		if (d == 0)
			throw Exception("CTRL", Exception::Type::DIMS_INVALID);

		std::vector<size_t> ctrlgate = ctrl; // ctrl + gate subsystem vector
		ctrlgate.insert(std::end(ctrlgate), std::begin(subsys),
				std::end(subsys));

		std::vector<size_t> dims; // local dimensions vector
		dims.insert(std::begin(dims), n, d);

		// check that ctrl + gate subsystem is valid with respect to local dimensions
		if (!internal::_check_subsys_match_dims(ctrlgate, dims))
			throw Exception("CTRL", Exception::Type::SUBSYS_MISMATCH_DIMS);

		// check that subsys list match the dimension of the matrix
		if (A.cols() != std::pow(d, subsys.size()))
			throw Exception("CTRL", Exception::Type::DIMS_MISMATCH_MATRIX);
		// END EXCEPTION CHECKS

		// Use static allocation for speed!
		size_t Cdims[ct::maxn];
		size_t midx_row[ct::maxn];
		size_t midx_col[ct::maxn];

		size_t CdimsA[ct::maxn];
		size_t midxA_row[ct::maxn];
		size_t midxA_col[ct::maxn];

		size_t Cdims_bar[ct::maxn];
		size_t Csubsys_bar[ct::maxn];
		size_t midx_bar[ct::maxn];

		size_t ngate = subsys.size();
		size_t nctrl = ctrl.size();
		size_t nsubsys_bar = n - ctrlgate.size();
		size_t D = std::pow(d, n);
		size_t DA = static_cast<size_t>(A.cols());
		size_t Dsubsys_bar = std::pow(d, nsubsys_bar);

		for (size_t k = 0, cnt = 0; k < n; k++)
		{
			midx_row[k] = midx_col[k] = 0;
			Cdims[k] = d;

			// compute the complementary subsystem of ctrlgate w.r.t. dims
			if (std::find(std::begin(ctrlgate), std::end(ctrlgate), k)
					== std::end(ctrlgate))
			{
				Csubsys_bar[cnt] = k;
				Cdims_bar[cnt] = d;
				midx_bar[cnt] = 0;
				cnt++;
			}
		}

		for (size_t k = 0; k < ngate; k++)
		{
			midxA_row[k] = midxA_col[k] = 0;
			CdimsA[k] = d;
		}

		types::cmat result = types::cmat::Identity(D, D);
		types::cmat Ak;

		// run over the complement indexes
		for (size_t i = 0; i < Dsubsys_bar; i++)
		{
			// get the complement's row multi-index
			internal::_n2multiidx(i, nsubsys_bar, Cdims_bar, midx_bar);
			for (size_t k = 0; k < d; k++)
			{
				Ak = powm(A, k); // compute A^k
				// run over the subsys's row multi-index
				for (size_t a = 0; a < DA; a++)
				{
					// get the subsys's row multi-index
					internal::_n2multiidx(a, ngate, CdimsA, midxA_row);

					// construct the result's row multi-index

					// first the ctrl part (equal for both row and column)
					for (size_t c = 0; c < nctrl; c++)
						midx_row[ctrl[c]] = midx_col[ctrl[c]] = k;

					// then the complement part (equal for column)
					for (size_t c = 0; c < nsubsys_bar; c++)
						midx_row[Csubsys_bar[c]] = midx_col[Csubsys_bar[c]] =
								midx_bar[c];

					// then the subsys part
					for (size_t c = 0; c < ngate; c++)
						midx_row[subsys[c]] = midxA_row[c];

					// run over the subsys's column multi-index
					for (size_t b = 0; b < DA; b++)
					{
						// get the subsys's column multi-index
						internal::_n2multiidx(b, ngate, CdimsA, midxA_col);

						// construct the result's column multi-index
						for (size_t c = 0; c < ngate; c++)
							midx_col[subsys[c]] = midxA_col[c];

						// finally write the values
						result(internal::_multiidx2n(midx_row, n, Cdims),
								internal::_multiidx2n(midx_col, n, Cdims)) = Ak(
								a, b);
					}
				}

			}
		}
		return result;
	}

};
/* class Gates */

} /* namespace qpp */

#endif /* GATES_H_ */
