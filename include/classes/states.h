/*
 * states.h
 *
 *  Created on: Apr 20, 2014
 *      Author: vlad
 */

#ifndef STATES_H_
#define STATES_H_

namespace qpp
{

class States: public Singleton<const States> // const Singleton
{
	friend class Singleton<const States> ;
public:
	// Pauli eigen-states
	ket x0 { ket::Zero(2) };
	ket x1 { ket::Zero(2) };
	ket y0 { ket::Zero(2) };
	ket y1 { ket::Zero(2) };
	ket z0 { ket::Zero(2) };
	ket z1 { ket::Zero(2) };

	// projectors onto Pauli eigen-states
	cmat px0 { cmat::Zero(2, 2) };
	cmat px1 { cmat::Zero(2, 2) };
	cmat py0 { cmat::Zero(2, 2) };
	cmat py1 { cmat::Zero(2, 2) };
	cmat pz0 { cmat::Zero(2, 2) };
	cmat pz1 { cmat::Zero(2, 2) };

	// Bell states
	ket b00 { ket::Zero(4) };
	ket b01 { ket::Zero(4) };
	ket b10 { ket::Zero(4) };
	ket b11 { ket::Zero(4) };

	// projectors onto Bell states
	cmat pb00 { cmat::Zero(4, 4) };
	cmat pb01 { cmat::Zero(4, 4) };
	cmat pb10 { cmat::Zero(4, 4) };
	cmat pb11 { cmat::Zero(4, 4) };

	// W and GHZ states
	ket GHZ { ket::Zero(8) };
	ket W { ket::Zero(8) };

	// projectors onto GHZ and W
	cmat pGHZ { cmat::Zero(8, 8) };
	cmat pW { cmat::Zero(8, 8) };
private:
	States()
	{
		// initialize
		x0 << 1 / std::sqrt(2.), 1 / std::sqrt(2.);
		x1 << 1 / std::sqrt(2.), -1 / std::sqrt(2.);
		y0 << 1 / std::sqrt(2.), 1_i / std::sqrt(2.);
		y1 << 1 / std::sqrt(2.), -1_i / std::sqrt(2.);
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
		b00 << 1 / std::sqrt(2.), 0, 0, 1 / std::sqrt(2.);// (|00>+|11>)/sqrt(2)
		b01 << 0, 1 / std::sqrt(2.), 1 / std::sqrt(2.), 0;// (|01>+|10>)/sqrt(2)
		b10 << 1 / std::sqrt(2.), 0, 0, -1 / std::sqrt(2.);	// (|00>-|11>)/sqrt(2)
		b11 << 0, 1 / std::sqrt(2.), -1 / std::sqrt(2.), 0;	// (|01>-|10>)/sqrt(2)

		pb00 = b00 * b00.adjoint();
		pb01 = b01 * b01.adjoint();
		pb10 = b10 * b10.adjoint();
		pb11 = b11 * b11.adjoint();

		GHZ << 1, 0, 0, 0, 0, 0, 0, 1;
		GHZ = GHZ / std::sqrt(2.);
		W << 0, 1, 1, 0, 1, 0, 0, 0;
		W = W / std::sqrt(3.);

		pGHZ = GHZ * GHZ.adjoint();
		pW = W * W.adjoint();
	}
};
/* class States */

} /* namespace qpp */

#endif /* STATES_H_ */
