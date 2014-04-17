/*
 * gates.h
 *
 *  Created on: Apr 7, 2014
 *      Author: vlad
 */

#ifndef GATES_H_
#define GATES_H_

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
	types::cmat C_S; // Controlled-S
	types::cmat CNOTba; // CNOT ctrl2 target1
	types::cmat SWAP; // SWAP gate

	// three qubit gates
	types::cmat TOF; // Toffoli
	types::cmat FRED; // Fredkin

	// Pauli eigen-states
	types::ket x0;
	types::ket x1;
	types::ket y0;
	types::ket y1;
	types::ket z0;
	types::ket z1;

	// projectors onto Pauli eigen-states
	types::cmat px0;
	types::cmat px1;
	types::cmat py0;
	types::cmat py1;
	types::cmat pz0;
	types::cmat pz1;

	// Bell states
	types::ket b00;
	types::ket b01;
	types::ket b10;
	types::ket b11;

	// projectors onto Bell states
	types::cmat pb00;
	types::cmat pb01;
	types::cmat pb10;
	types::cmat pb11;

	// W and GHZ states
	types::ket GHZ;
	types::ket W;

	// projectors onto GHZ and W
	types::cmat pGHZ;
	types::cmat pW;
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
			C_S(types::cmat::Identity(4, 4)), //
			CNOTba(types::cmat::Zero(4, 4)), //
			SWAP(types::cmat::Identity(4, 4)), //
			TOF(types::cmat::Identity(8, 8)), //
			FRED(types::cmat::Identity(8, 8)), //
			x0(types::ket::Zero(2)), //
			x1(types::ket::Zero(2)), //
			y0(types::ket::Zero(2)), //
			y1(types::ket::Zero(2)), //
			z0(types::ket::Zero(2)), //
			z1(types::ket::Zero(2)), //
			px0(types::cmat::Zero(2, 2)), //
			px1(types::cmat::Zero(2, 2)), //
			py0(types::cmat::Zero(2, 2)), //
			py1(types::cmat::Zero(2, 2)), //
			pz0(types::cmat::Zero(2, 2)), //
			pz1(types::cmat::Zero(2, 2)), //
			b00(types::ket::Zero(4)), //
			b01(types::ket::Zero(4)), //
			b10(types::ket::Zero(4)), //
			b11(types::ket::Zero(4)), //
			pb00(types::cmat::Zero(4, 4)), //
			pb01(types::cmat::Zero(4, 4)), //
			pb10(types::cmat::Zero(4, 4)), //
			pb11(types::cmat::Zero(4, 4)), //
			GHZ(types::ket::Zero(8)), //
			W(types::ket::Zero(8)), //
			pGHZ(types::cmat::Zero(8, 8)), //
			pW(types::cmat::Zero(8, 8))
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
		C_S(3, 3) = ct::ii;

		SWAP.block(1, 1, 2, 2) = X;
		TOF.block(6, 6, 2, 2) = X;
		FRED.block(4, 4, 4, 4) = SWAP;

		x0 << 1 / std::sqrt(2), 1 / std::sqrt(2);
		x1 << 1 / std::sqrt(2), -1 / std::sqrt(2);
		y0 << 1 / std::sqrt(2), ct::ii / std::sqrt(2);
		y1 << 1 / std::sqrt(2), -ct::ii / std::sqrt(2);
		z0 << 1, 0;
		z1 << 0, 1;
		px0 = x0 * x0.adjoint();
		px1 = x1 * x1.adjoint();
		py0 = y0 * y0.adjoint();
		py1 = y1 * y1.adjoint();
		pz0 = z0 * z0.adjoint();
		pz1 = z1 * z1.adjoint();

		// Bell states, following convention from Nielsen & Chuang
		// |ij> -> |b_{ij}> by the CNOT*(H x Id) circuit
		b00 << 1 / std::sqrt(2), 0, 0, 1 / std::sqrt(2);// (|00>+|11>)/sqrt(2)
		b01 << 0, 1 / std::sqrt(2), 1 / std::sqrt(2), 0;// (|01>+|10>)/sqrt(2)
		b10 << 1 / std::sqrt(2), 0, 0, -1 / std::sqrt(2);// (|00>-|11>)/sqrt(2)
		b11 << 0, 1 / std::sqrt(2), -1 / std::sqrt(2), 0;// (|01>-|10>)/sqrt(2)

		pb00 = b00 * b00.adjoint();
		pb01 = b01 * b01.adjoint();
		pb10 = b10 * b10.adjoint();
		pb11 = b11 * b11.adjoint();

		GHZ << 1, 0, 0, 0, 0, 0, 0, 1;
		GHZ = GHZ / std::sqrt(2);
		W << 0, 1, 1, 0, 1, 0, 0, 0;
		W = W / std::sqrt(3);

		pGHZ = GHZ * GHZ.adjoint();
		pW = W * W.adjoint();
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
			const std::vector<size_t>& gate, size_t n, size_t D = 2) const
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
		if (gate.size() == 0)
			throw Exception("CTRL", Exception::Type::ZERO_SIZE);

		// check out of range
		if (n == 0)
			throw Exception("CTRL", Exception::Type::OUT_OF_RANGE);

		// check valid local dimension
		if (D == 0)
			throw Exception("CTRL", Exception::Type::DIMS_INVALID);

		std::vector<size_t> ctrlgate = ctrl; // ctrl + gate subsystem vector
		ctrlgate.insert(std::end(ctrlgate), std::begin(gate), std::end(gate));

		std::vector<size_t> dims; // local dimensions vector
		dims.insert(std::begin(dims), n, D);

		// check that ctrl+gate subsystem is valid with respect to local dimensions
		if (!internal::_check_subsys_match_dims(ctrlgate, dims))
			throw Exception("CTRL", Exception::Type::SUBSYS_MISMATCH_DIMS);

		// check that gate list match the dimension of the matrix
		if (A.cols() != std::pow(D, gate.size()))
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

		for (size_t k = 0, cnt = 0; k < n; k++)
		{
			midx_row[k] = midx_col[k] = 0;
			Cdims[k] = D;

			// compute the complementary subsystem
			if (std::find(std::begin(ctrlgate), std::end(ctrlgate), k)
					== std::end(ctrlgate))
			{
				Csubsys_bar[cnt] = k;
				cnt++;
			}
		}

		for (size_t k = 0; k < gate.size(); k++)
		{
			midxA_row[k] = midxA_col[k] = 0;
			CdimsA[k] = D;
		}

		for (size_t k = 0; k < n - ctrlgate.size(); k++)
		{
			midx_bar[k] = 0;
			Cdims_bar[k] = D;
		}

		types::cmat result = types::cmat::Identity(std::pow(D, n),
				std::pow(D, n));
		types::cmat Ak;

		// run over the complement indexes
		for (size_t i = 0; i < std::pow(D, n - ctrlgate.size()); i++)
		{
			// get the complement's row multi-index
			internal::_n2multiidx(i, n - ctrlgate.size(), Cdims_bar, midx_bar);
			for (size_t k = 0; k < D; k++)
			{
				Ak = powm(A, k); // compute A^k
				// run over the gate's row multi-index
				for (size_t a = 0; a < static_cast<size_t>(A.cols()); a++)
				{
					// get the row multi-index of the gate
					internal::_n2multiidx(a, gate.size(), CdimsA, midxA_row);

					// construct the total row multi-index

					// first the ctrl part (equal for both row and column)
					for (size_t c = 0; c < ctrl.size(); c++)
						midx_row[ctrl[c]] = midx_col[ctrl[c]] = k;

					// then the complement part (equal for column)
					for (size_t c = 0; c < n - ctrlgate.size(); c++)
						midx_row[Csubsys_bar[c]] = midx_col[Csubsys_bar[c]] =
								midx_bar[c];

					// then the gate part
					for (size_t c = 0; c < gate.size(); c++)
						midx_row[gate[c]] = midxA_row[c];

					// run over the gate's column multi-index
					for (size_t b = 0; b < static_cast<size_t>(A.cols()); b++)
					{
						// get the column multi-index of the gate
						internal::_n2multiidx(b, gate.size(), CdimsA,
								midxA_col);

						// construct the total column multi-index
						for (size_t c = 0; c < gate.size(); c++)
							midx_col[gate[c]] = midxA_col[c];

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

} /* namespace qpp */

#endif /* GATES_H_ */
