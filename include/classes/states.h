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
	~States()
	{
	}
	;
public:
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
	States() :
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
		// initialize
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
};
/* class States */

} /* namespace qpp */

#endif /* STATES_H_ */
