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

// Eigen predefined:
// MatrixXcd::Identity(D, D), MatrixXcd::Zero (D,D), MatrixXcd::Random(D, D)

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

	H << 1 / std::sqrt(2.), 1 / std::sqrt(2.), 1 / std::sqrt(2.), -1
			/ std::sqrt(2.);
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

// two qubit gates
inline types::cmat CU(const types::cmat &U)
{
	if (U.cols() != 2 || U.rows() != 2)
		throw Exception("CU", Exception::Type::NOT_QUBIT_GATE);
	types::cmat result = types::cmat::Zero(4, 4);
	result(0, 0) = 1;
	result(1, 1) = 1;
	result.block(2, 2, 2, 2) = U;
	return result;
}

// one quDit gates
inline types::cmat Zd(size_t D)
{
	if (D == 0)
		throw Exception("Zd", Exception::Type::DIMS_HAVE_ZERO);

	types::cmat result(D, D);
	result = types::cmat::Zero(D, D);
	for (size_t i = 0; i < D; i++)
		result(i, i) = pow(ct::omega(D), i);
	return result;
}

inline types::cmat Fd(size_t D)
{
	if (D == 0)
		throw Exception("Fd", Exception::Type::DIMS_HAVE_ZERO);

	types::cmat result(D, D);
	result = types::cmat::Zero(D, D);
	for (size_t j = 0; j < D; j++)
		for (size_t i = 0; i < D; i++)
			result(i, j) = 1 / std::sqrt(D) * pow(ct::omega(D), i * j);
	return result;
}

inline types::cmat Xd(size_t D)
{
	if (D == 0)
		throw Exception("Xd", Exception::Type::DIMS_HAVE_ZERO);

	return Fd(D) * Zd(D) * Fd(D).inverse();
}

// two qudit gates
inline types::cmat CUd(const types::cmat &U)
{
	// check square matrix
	if (!internal::_check_square_mat(U))
		throw Exception("", Exception::Type::MATRIX_NOT_SQUARE);

	size_t D = static_cast<size_t>(U.cols());
	types::cmat result(D * D, D * D);
	result = types::cmat::Zero(D * D, D * D);
	types::cmat tmp(D, D);
	tmp = types::cmat::Zero(D, D); // the dyad |i><i|

	for (size_t i = 0; i < D; i++)
	{
		if (i > 0)
			tmp(i - 1, i - 1) = 0;
		tmp(i, i) = 1;
		result += kron(tmp, powm_int(U, i));
	}
	return result;
}

}
}

#endif	/* GATES_H_ */

