/* 
 * File:   gates.cpp
 * Author: vlad
 *
 * Created on December 12, 2013, 10:43 PM
 */

#include <iostream>
#include "gates.h"
#include "types.h"
#include "constants.h"
#include "util.h"

// TODO: make everything inline in header file

namespace qpp
{
namespace gt
{

// various matrices (declared extern in "gates.h") ; make the gates visible
types::cmat2 H, Id2, X, Y, Z, S, T;
types::cmat4 CNOT, CP;
types::cmat TOF(8, 8);

types::cmat2 Rtheta(double theta)
{
	types::cmat2 result;
	result << 1, 0, 0, exp(ct::ii * theta);
	return result;
}

types::cmat4 CU(const types::cmat2 &U)
{
	types::cmat4 result;
	result = types::cmat::Zero(4, 4);
	result(0, 0) = 1;
	result(1, 1) = 1;
	result.block(2, 2, 2, 2) = U;
	return result;
}

// one quDit gates

types::cmat Zd(size_t D)
{
	types::cmat result(D, D);
	result = types::cmat::Zero(D, D);
	for (size_t i = 0; i < D; i++)
		result(i, i) = pow(ct::omega(D), i);
	return result;
}

types::cmat Fd(size_t D)
{
	types::cmat result(D, D);
	result = types::cmat::Zero(D, D);
	for (size_t i = 0; i < D; i++)
		for (size_t j = 0; j < D; j++)
			result(i, j) = 1 / sqrt(D) * pow(ct::omega(D), i * j);
	return result;
}

types::cmat Xd(size_t D)
{
	return Fd(D) * Zd(D) * Fd(D).inverse();
}

types::cmat CUd(const types::cmat &U)
{
	int D = U.cols(); // retrieves the dimension from the size of U
	types::cmat result(D * D, D * D);
	result = types::cmat::Zero(D * D, D * D);
	types::cmat tmp(D, D);
	tmp = types::cmat::Zero(D, D); // the dyad |i><i|

	for (int i = 0; i < D; i++)
	{
		if (i > 0)
			tmp(i - 1, i - 1) = 0;
		tmp(i, i) = 1;
		result += kron(tmp, mat_pow(U, i));
	}
	return result;
}

int _init_gates() // Initialize the gates
{
	// initialize the constants and gates
	H = Id2 = X = Y = Z = S = T = types::cmat2::Zero();
	CNOT = types::cmat4::Zero();
	TOF = types::cmat::Zero(8, 8);

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

	return 0;
}

}
}
