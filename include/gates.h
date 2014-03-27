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
#include "util.h"

// Eigen predefined:
// MatrixXcd::Identity(D, D), MatrixXcd::Zero (D,D), MatrixXcd::Random(D, D)

namespace qpp
{
namespace gt
{

// one qubit gates
extern Eigen::MatrixXcd H; // Hadamard matrix
extern Eigen::MatrixXcd Id2; // Identity matrix
extern Eigen::MatrixXcd X; // X matrix
extern Eigen::MatrixXcd Y; // Y matrix
extern Eigen::MatrixXcd Z; // Z matrix
extern Eigen::MatrixXcd S; // S gate
extern Eigen::MatrixXcd T; // T gate

// two qubit gates
extern Eigen::MatrixXcd CNOT; // CNOT
extern Eigen::MatrixXcd CP; // Controlled-Phase

// three qubit gates
extern Eigen::MatrixXcd TOF; // Toffoli

inline void _init_gates() // Initialize the gates
{
	// initialize the constants and gates
	H = Id2 = X = Y = Z = S = T = Eigen::MatrixXcd::Zero(2,2);
	CNOT = CP = Eigen::MatrixXcd::Zero(4,4);
	TOF = Eigen::MatrixXcd::Zero(8, 8);

	H << 1 / sqrt(2), 1 / sqrt(2), 1 / sqrt(2), -1 / sqrt(2);
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
inline Eigen::MatrixXcd Rtheta(double theta)
{
	Eigen::MatrixXcd result(2,2);
	result << 1, 0, 0, exp(ct::ii * theta);
	return result;
}

// two qubit gates
inline Eigen::MatrixXcd CU(const Eigen::MatrixXcd &U)
{
	Eigen::MatrixXcd result = Eigen::MatrixXcd::Zero(4, 4);
	result(0, 0) = 1;
	result(1, 1) = 1;
	result.block(2, 2, 2, 2) = U;
	return result;
}

// one quDit gates

inline Eigen::MatrixXcd Zd(size_t D)
{
	Eigen::MatrixXcd result(D, D);
	result = Eigen::MatrixXcd::Zero(D, D);
	for (size_t i = 0; i < D; i++)
		result(i, i) = pow(ct::omega(D), i);
	return result;
}

inline Eigen::MatrixXcd Fd(size_t D)
{
	Eigen::MatrixXcd result(D, D);
	result = Eigen::MatrixXcd::Zero(D, D);
	for (size_t i = 0; i < D; i++)
		for (size_t j = 0; j < D; j++)
			result(i, j) = 1 / sqrt(D) * pow(ct::omega(D), i * j);
	return result;
}

inline Eigen::MatrixXcd Xd(size_t D)
{
	return Fd(D) * Zd(D) * Fd(D).inverse();
}

// two qudit gates
inline Eigen::MatrixXcd CUd(const Eigen::MatrixXcd &U)
{
	size_t D = U.cols(); // retrieves the dimension from the size of U
	Eigen::MatrixXcd result(D * D, D * D);
	result = Eigen::MatrixXcd::Zero(D * D, D * D);
	Eigen::MatrixXcd tmp(D, D);
	tmp = Eigen::MatrixXcd::Zero(D, D); // the dyad |i><i|

	for (size_t i = 0; i < D; i++)
	{
		if (i > 0)
			tmp(i - 1, i - 1) = 0;
		tmp(i, i) = 1;
		result += kron(tmp, mpower(U, i));
	}
	return result;
}

}
}

#endif	/* GATES_H_ */

