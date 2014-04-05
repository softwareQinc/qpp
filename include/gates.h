/* 
 * File:   gates.h
 * Author: vlad
 *
 * Created on December 12, 2013, 10:42 PM
 */

#ifndef GATES_H_
#define	GATES_H_

#include "types.h"
#include "constants.h"
#include "functions.h"
#include "internal.h"
#include "exception.h"

// quantum gates

// Eigen
// MatrixXcd::Zero (D,D), MatrixXcd::Random(D, D)

namespace qpp
{
namespace gt
{

// one qubit gates
extern types::cmat H; // Hadamard matrix
extern types::cmat Id2; // Identity matrix
extern types::cmat X; // X matrix
extern types::cmat Y; // Y matrix
extern types::cmat Z; // Z matrix
extern types::cmat S; // S gate
extern types::cmat T; // T gate

// two qubit gates
extern types::cmat CNOT; // CNOT
extern types::cmat CP; // Controlled-Phase

// three qubit gates
extern types::cmat TOF; // Toffoli

inline void _init_gates() // Initialize the gates, call it from qpp::_init()
{
	// initialize the constants and gates
	H = Id2 = X = Y = Z = S = T = types::cmat::Zero(2, 2);
	CNOT = CP = types::cmat::Zero(4, 4);
	TOF = types::cmat::Zero(8, 8);

	H << 1 / std::sqrt(2), 1 / std::sqrt(2), 1 / std::sqrt(2), -1
			/ std::sqrt(2);
	Id2 << 1, 0, 0, 1;
	X << 0, 1, 1, 0;
	Z << 1, 0, 0, -1;
	Y(0, 1) = -ct::ii;
	Y(1, 0) = ct::ii;
	S(0, 0) = 1;
	S(1, 1) = ct::ii;
	T(0, 0) = 1;
	T(1, 1) = exp(ct::ii * ct::pi / 4.0);
	CNOT(0, 0) = 1;
	CNOT(1, 1) = 1;
	CNOT(2, 3) = 1;
	CNOT(3, 2) = 1;
	CP(0, 0) = 1;
	CP(1, 1) = 1;
	CP(2, 2) = 1;
	CP(3, 3) = -1;
	TOF(0, 0) = 1;
	TOF(1, 1) = 1;
	TOF(2, 2) = 1;
	TOF(3, 3) = 1;
	TOF(4, 4) = 1;
	TOF(5, 5) = 1;
	TOF(6, 7) = 1;
	TOF(7, 6) = 1;
}

// gates with variable dimension

// one qubit gates
inline types::cmat Rtheta(double theta)
{
	types::cmat result(2, 2);
	result << 1, 0, 0, exp(ct::ii * theta);
	return result;
}

// one quDit gates

inline types::cmat Id(size_t D)
{
	if (D == 0)
		throw Exception("Id", Exception::Type::DIMS_INVALID);
	return types::cmat::Identity(D, D);
}

inline types::cmat Zd(size_t D)
{
	if (D == 0)
		throw Exception("Zd", Exception::Type::DIMS_INVALID);

	types::cmat result(D, D);
	result = types::cmat::Zero(D, D);
	for (size_t i = 0; i < D; i++)
		result(i, i) = pow(ct::omega(D), i);
	return result;
}

inline types::cmat Fd(size_t D)
{
	if (D == 0)
		throw Exception("Fd", Exception::Type::DIMS_INVALID);

	types::cmat result(D, D);
	result = types::cmat::Zero(D, D);
	for (size_t j = 0; j < D; j++)
		for (size_t i = 0; i < D; i++)
			result(i, j) = 1 / std::sqrt(D) * pow(ct::omega(D), i * j);
	return result;
}

inline types::cmat Xd(size_t D) // X|k>=|k+1>
{
	if (D == 0)
		throw Exception("Xd", Exception::Type::DIMS_INVALID);

	return Fd(D).inverse() * Zd(D) * Fd(D);
}

// multi-quDit multi-controlled-gate
// faster than doing sum |j><j| A^j, especially for small A and large n
// for large A relative to D^n, use sum |j><j| A^j (CTRLsum)
inline types::cmat CTRL(const types::cmat& A, const std::vector<size_t>& ctrl,
		const std::vector<size_t>& gate, size_t n, size_t D = 2)
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
	if (!internal::_check_subsys(ctrlgate, dims))
		throw Exception("CTRL", Exception::Type::SUBSYS_MISMATCH_DIMS);

	// check that gate list match the dimension of the matrix
	if (A.cols() != std::pow(D, gate.size()))
		throw Exception("CTRL", Exception::Type::DIMS_MISMATCH_MATRIX);
// END EXCEPTION CHECKS

	size_t* Cdims = new size_t[n];
	size_t* midx_row = new size_t[n];
	size_t* midx_col = new size_t[n];

	size_t* CdimsA = new size_t[gate.size()];
	size_t* midxA_row = new size_t[gate.size()];
	size_t* midxA_col = new size_t[gate.size()];

	size_t* Cdims_ctrl = new size_t[ctrl.size()];
	size_t* midx_ctrl = new size_t[ctrl.size()];

	size_t* Cdims_bar = new size_t[n - ctrlgate.size()];
	size_t* Csubsys_bar = new size_t[n - ctrlgate.size()];
	size_t* midx_bar = new size_t[n - ctrlgate.size()];

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

	for (size_t k = 0; k < ctrl.size(); k++)
	{
		midx_ctrl[k] = 0;
		Cdims_ctrl[k] = D;
	}

	for (size_t k = 0; k < n - ctrlgate.size(); k++)
	{
		midx_bar[k] = 0;
		Cdims_bar[k] = D;
	}

	types::cmat result = types::cmat::Identity(std::pow(D, n), std::pow(D, n));
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

				// first the ctrl part (equal for column)
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
					internal::_n2multiidx(b, gate.size(), CdimsA, midxA_col);

					// construct the total column multi-index
					for (size_t c = 0; c < gate.size(); c++)
						midx_col[gate[c]] = midxA_col[c];

					// finally write the values
					result(internal::_multiidx2n(midx_row, n, Cdims),
							internal::_multiidx2n(midx_col, n, Cdims)) = Ak(a,
							b);
				}
			}

		}
	}

	delete[] Cdims;
	delete[] midx_row;
	delete[] midx_col;

	delete[] CdimsA;
	delete[] midxA_row;
	delete[] midxA_col;

	delete[] Cdims_ctrl;
	delete[] midx_ctrl;

	delete[] Cdims_bar;
	delete[] Csubsys_bar;
	delete[] midx_bar;

	return result;
}

}
}

#endif	/* GATES_H_ */

