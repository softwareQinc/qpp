/*
 * Quantum++
 *
 * Copyright (c) 2013 - 2014 Vlad Gheorghiu (vgheorgh@gmail.com)
 *
 * This file is part of Quantum++.
 *
 * Quantum++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Quantum++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Quantum++.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef INCLUDE_CLASSES_STATES_H_
#define INCLUDE_CLASSES_STATES_H_

namespace qpp
{

    /**
    * \class qpp::States
    * \brief const Singleton class that implements most commonly used states
    */
    class States : public internal::Singleton<const States> // const Singleton
    {
        friend class internal::Singleton<const States>;

    public:
        // Pauli eigen-states
        ket x0{ket::Zero(2)}; ///< Pauli Sigma-X 0-eigenstate |+>
        ket x1{ket::Zero(2)}; ///< Pauli Sigma-X 1-eigenstate |->
        ket y0{ket::Zero(2)}; ///< Pauli Sigma-Y 0-eigenstate
        ket y1{ket::Zero(2)}; ///< Pauli Sigma-Y 1-eigenstate
        ket z0{ket::Zero(2)}; ///< Pauli Sigma-Z 0-eigenstate |0>
        ket z1{ket::Zero(2)}; ///< Pauli Sigma-Z 1-eigenstate |1>

        // projectors onto Pauli eigen-states
        cmat px0{cmat::Zero(2, 2)}; ///< Projector onto the Pauli Sigma-X 0-eigenstate |+><+|
        cmat px1{cmat::Zero(2, 2)}; ///< Projector onto the Pauli Sigma-X 1-eigenstate |-><-|
        cmat py0{cmat::Zero(2, 2)}; ///< Projector onto the Pauli Sigma-Y 0-eigenstate
        cmat py1{cmat::Zero(2, 2)}; ///< Projector onto the Pauli Sigma-Y 1-eigenstate
        cmat pz0{cmat::Zero(2, 2)}; ///< Projector onto the Pauli Sigma-Z 0-eigenstate |0><0|
        cmat pz1{cmat::Zero(2, 2)}; ///< Projector onto the Pauli Sigma-Z 1-eigenstate |1><1|

        // Bell states
        ket b00{ket::Zero(4)}; ///< Bell-00 state (following the convention in Nielsen and Chuang)
        ket b01{ket::Zero(4)}; ///< Bell-01 state (following the convention in Nielsen and Chuang)
        ket b10{ket::Zero(4)}; ///< Bell-10 state (following the convention in Nielsen and Chuang)
        ket b11{ket::Zero(4)}; ///< Bell-11 state (following the convention in Nielsen and Chuang)

        // projectors onto Bell states
        cmat pb00{cmat::Zero(4, 4)}; ///< Projector onto the Bell-00 state
        cmat pb01{cmat::Zero(4, 4)}; ///< Projector onto the Bell-01 state
        cmat pb10{cmat::Zero(4, 4)}; ///< Projector onto the Bell-10 state
        cmat pb11{cmat::Zero(4, 4)}; ///< Projector onto the Bell-11 state

        // W and GHZ states
        ket GHZ{ket::Zero(8)}; ///< GHZ state
        ket W{ket::Zero(8)}; ///< W state

        // projectors onto GHZ and W
        cmat pGHZ{cmat::Zero(8, 8)}; ///< Projector onto the GHZ state
        cmat pW{cmat::Zero(8, 8)}; ///< Projector onto the W state
    private:
        /**
        * Initialize the states
        */
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
            b10 << 1 / std::sqrt(2.), 0, 0, -1 / std::sqrt(2.);    // (|00>-|11>)/sqrt(2)
            b11 << 0, 1 / std::sqrt(2.), -1 / std::sqrt(2.), 0;    // (|01>-|10>)/sqrt(2)

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
    }; /* class States */

} /* namespace qpp */

#endif /* INCLUDE_CLASSES_STATES_H_ */
