/* 
 * File:   gates.h
 * Author: vlad
 *
 * Created on December 12, 2013, 10:42 PM
 */

#ifndef GATES_H_
#define	GATES_H_

#include "types.h"

// Eigen predefined:
// MatrixXcd::Identity(D, D), MatrixXcd::Zero (D,D), MatrixXcd::Random(D, D)

namespace qpp
{
namespace gt
{

// one qubit gates
extern types::cmat2 H; // Hadamard matrix
extern types::cmat2 Id2; // Identity matrix
extern types::cmat2 X; // X matrix
extern types::cmat2 Y; // Y matrix
extern types::cmat2 Z; // Z matrix
extern types::cmat2 S; // S gate
extern types::cmat2 T; // T gate
types::cmat2 Rtheta(double theta); // Rotation of theta around the Z axis

// two qubit gates
extern types::cmat4 CNOT; // CNOT
extern types::cmat4 CP; // Controlled-Phase
// Controlled-U, for arbitrary U
extern types::cmat4 CU(const types::cmat2 &);

// three qubit gates
extern types::cmat TOF; // Toffoli

// one quDit gates
extern types::cmat Zd(size_t); // generalized Z gate
extern types::cmat Fd(size_t); // generalized Fourier gate
extern types::cmat Xd(size_t); // generalized X gate

// two quDit gates
// Controlled-U, for arbitrary U
extern types::cmat CUd(const types::cmat &);

int _init_gates(); // Initialize the gates, internal function

}
}

#endif	/* GATES_H_ */

