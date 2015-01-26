/*
 * Quantum++
 *
 * Copyright (c) 2013 - 2015 Vlad Gheorghiu (vgheorgh@gmail.com)
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

/**
* \file classes/gates.h
* \brief Quantum gates
*/

#ifndef CLASSES_GATES_H_
#define CLASSES_GATES_H_

namespace qpp
{

/**
* \class qpp::Gates
* \brief const Singleton class that implements most commonly used gates
*/
class Gates : public internal::Singleton<const Gates> // const Singleton
{
    friend class internal::Singleton<const Gates>;

public:
    // One qubit gates
    cmat Id2{cmat::Identity(2, 2)};     ///< Identity gate
    cmat H{cmat::Zero(2, 2)};           ///< Hadamard gate
    cmat X{cmat::Zero(2, 2)};           ///< Pauli Sigma-X gate
    cmat Y{cmat::Zero(2, 2)};           ///< Pauli Sigma-Y gate
    cmat Z{cmat::Zero(2, 2)};           ///< Pauli Sigma-Z gate
    cmat S{cmat::Zero(2, 2)};           ///< S gate
    cmat T{cmat::Zero(2, 2)};           ///< T gate

    // two qubit gates
    cmat CNOT{cmat::Identity(4, 4)};  ///< Controlled-NOT control target gate
    cmat CZ{cmat::Identity(4, 4)};      ///< Controlled-Phase gate
    cmat CNOTba{cmat::Zero(4, 4)};      ///< Controlled-NOT target control gate
    cmat SWAP{cmat::Identity(4, 4)};    ///< SWAP gate

    // three qubit gates
    cmat TOF{cmat::Identity(8, 8)};     ///< Toffoli gate
    cmat FRED{cmat::Identity(8, 8)};    ///< Fredkin gate
private:
    /**
    * \brief Initializes the gates
    */
    Gates()
    {
        H << 1 / std::sqrt(2.), 1 / std::sqrt(2.), 1 / std::sqrt(2.), -1
                / std::sqrt(2.);
        X << 0, 1, 1, 0;
        Z << 1, 0, 0, -1;
        Y << 0, -1_i, 1_i, 0;
        S << 1, 0, 0, 1_i;
        T << 1, 0, 0, std::exp(1_i * pi / 4.0);
        CNOT.block(2, 2, 2, 2) = X;
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
    // variable gates

    // one qubit gates

    /**
    * \brief Rotation of \a theta about the 3-dimensional real unit vector \a n
    *
    * \param theta Rotation angle
    * \param n 3-dimensional real unit vector
    * \return Rotation gate
    */
    cmat Rn(double theta, std::vector<double> n) const
    {
        if (n.size() != 3) // not a 3-D vector
            throw Exception("qpp::Gates::Rn()", "n is not a 3-D vector!");

        cmat result(2, 2);
        result = std::cos(theta / 2) * Id2
                - 1_i * std::sin(theta / 2) * (n[0] * X + n[1] * Y + n[2] * Z);

        return result;
    }

    // one quDit gates

    /**
    * \brief Generalized Z gate for qudits
    *
    * \note Defined as \f$ Z = \sum_j \exp(2\pi i j/D) |j\rangle\langle j| \f$
    *
    * \param D Dimension of the Hilbert space
    * \return Generalized Z gate for qudits
    */
    cmat Zd(idx D) const
    {
        if (D == 0)
            throw Exception("qpp::Gates::Zd()", Exception::Type::DIMS_INVALID);

        cmat result = cmat::Zero(D, D);
        for (idx i = 0; i < D; ++i)
            result(i, i) = std::pow(omega(D), i);

        return result;
    }

    /**
    * \brief Fourier transform gate for qudits
    *
    * \note Defined as
    * \f$ F = \sum_{jk} \exp(2\pi i jk/D) |j\rangle\langle k| \f$
    *
    * \param D Dimension of the Hilbert space
    * \return Fourier transform gate for qudits
    */
    cmat Fd(idx D) const
    {
        if (D == 0)
            throw Exception("qpp::Gates::Fd()", Exception::Type::DIMS_INVALID);

        cmat result(D, D);
#pragma omp parallel for collapse(2)
        for (idx j = 0; j < D; ++j)
            for (idx i = 0; i < D; ++i)
                result(i, j) = 1 / std::sqrt(static_cast<double>(D))
                        * std::pow(omega(D), i * j);

        return result;
    }

    /**
    * \brief Generalized X gate for qudits
    *
    * \note Defined as \f$ X = \sum_j |j\oplus 1\rangle\langle j| \f$
    *
    * \param D Dimension of the Hilbert space
    * \return Generalized X gate for qudits
    */
    cmat Xd(idx D) const
    {
        if (D == 0)
            throw Exception("qpp::Gates::Xd()", Exception::Type::DIMS_INVALID);

        return Fd(D).inverse() * Zd(D) * Fd(D);
    }

    /**
    * \brief Identity gate
    *
    * \note Can change the return type from complex matrix (default)
    * by explicitly specifying the template parameter
    *
    * \param D Dimension of the Hilbert space
    * \return Identity gate
    */
    template<typename Derived = Eigen::MatrixXcd>
    Derived Id(idx D) const
    {
        if (D == 0)
            throw Exception("qpp::Gates::Id()", Exception::Type::DIMS_INVALID);

        return Derived::Identity(D, D);
    }

    /**
    * \brief Generates the multi-partite multiple-controlled-\a A gate
    * in matrix form
    * \see qpp::applyCTRL()
    *
    * \note The dimension of the gate \a A must match
    * the dimension of \a subsys
    *
    * \param A Eigen expression
    * \param ctrl Control subsystem indexes
    * \param subsys Subsystem indexes where the gate \a A is applied
    * \param n Total number of subsystes
    * \param d Subsystem dimensions
    * \return CTRL-A gate, as a matrix over the same scalar field as \a A
    */
    template<typename Derived>
    dyn_mat<typename Derived::Scalar> CTRL(const Eigen::MatrixBase<Derived>& A,
            const std::vector<idx>& ctrl,
            const std::vector<idx>& subsys,
            idx n, idx d = 2) const
    {
        const dyn_mat<typename Derived::Scalar>& rA = A;

        // EXCEPTION CHECKS
        // check matrix zero size
        if (!internal::_check_nonzero_size(rA))
            throw Exception("qpp::Gates::CTRL()", Exception::Type::ZERO_SIZE);

        // check square matrix
        if (!internal::_check_square_mat(rA))
            throw Exception("qpp::Gates::CTRL()",
                    Exception::Type::MATRIX_NOT_SQUARE);

        // check lists zero size
        if (ctrl.size() == 0)
            throw Exception("qpp::Gates::CTRL()", Exception::Type::ZERO_SIZE);
        if (subsys.size() == 0)
            throw Exception("qpp::Gates::CTRL()", Exception::Type::ZERO_SIZE);

        // check out of range
        if (n == 0)
            throw Exception("qpp::Gates::CTRL()",
                    Exception::Type::OUT_OF_RANGE);

        // check valid local dimension
        if (d == 0)
            throw Exception("qpp::Gates::CTRL()",
                    Exception::Type::DIMS_INVALID);

        // ctrl + gate subsystem vector
        std::vector<idx> ctrlgate = ctrl;
        ctrlgate.insert(std::end(ctrlgate), std::begin(subsys),
                std::end(subsys));
        std::sort(std::begin(ctrlgate), std::end(ctrlgate));

        std::vector<idx> dims(n, d); // local dimensions vector

        // check that ctrl + gate subsystem is valid
        // with respect to local dimensions
        if (!internal::_check_subsys_match_dims(ctrlgate, dims))
            throw Exception("qpp::Gates::CTRL()",
                    Exception::Type::SUBSYS_MISMATCH_DIMS);

        // check that subsys list match the dimension of the matrix
        if (rA.rows() != std::llround(std::pow(d, subsys.size())))
            throw Exception("qpp::Gates::CTRL()",
                    Exception::Type::DIMS_MISMATCH_MATRIX);
        // END EXCEPTION CHECKS

        // Use static allocation for speed!
        idx Cdims[maxn];
        idx midx_row[maxn];
        idx midx_col[maxn];

        idx CdimsA[maxn];
        idx midxA_row[maxn];
        idx midxA_col[maxn];

        idx Cdims_bar[maxn];
        idx Csubsys_bar[maxn];
        idx midx_bar[maxn];

        idx ngate = subsys.size();
        idx nctrl = ctrl.size();
        idx nsubsys_bar = n - ctrlgate.size();
        idx D = static_cast<idx>(std::llround(std::pow(d, n)));
        idx DA = static_cast<idx>(rA.rows());
        idx Dsubsys_bar = static_cast<idx>(
                std::llround(std::pow(d, nsubsys_bar)));

        // compute the complementary subsystem of ctrlgate w.r.t. dims
        std::vector<idx> allsubsys(n); // all subsystems
        std::iota(std::begin(allsubsys), std::end(allsubsys), 0);
        std::set_difference(std::begin(allsubsys), std::end(allsubsys),
                std::begin(ctrlgate), std::end(ctrlgate),
                std::begin(Csubsys_bar));

        for (idx k = 0; k < n; ++k)
        {
            midx_row[k] = midx_col[k] = 0;
            Cdims[k] = d;
        }

        for (idx k = 0; k < nsubsys_bar; ++k)
        {
            Cdims_bar[k] = d;
            midx_bar[k] = 0;
        }

        for (idx k = 0; k < ngate; ++k)
        {
            midxA_row[k] = midxA_col[k] = 0;
            CdimsA[k] = d;
        }

        dyn_mat<typename Derived::Scalar> result = dyn_mat<
                typename Derived::Scalar>::Identity(D, D);
        dyn_mat<typename Derived::Scalar> Ak;

        // run over the complement indexes
        for (idx i = 0; i < Dsubsys_bar; ++i)
        {
            // get the complement row multi-index
            internal::_n2multiidx(i, nsubsys_bar, Cdims_bar, midx_bar);
            for (idx k = 0; k < d; ++k)
            {
                Ak = powm(rA, k); // compute rA^k
                // run over the subsys row multi-index
                for (idx a = 0; a < DA; ++a)
                {
                    // get the subsys row multi-index
                    internal::_n2multiidx(a, ngate, CdimsA, midxA_row);

                    // construct the result row multi-index

                    // first the ctrl part (equal for both row and column)
                    for (idx c = 0; c < nctrl; ++c)
                        midx_row[ctrl[c]] = midx_col[ctrl[c]] = k;

                    // then the complement part (equal for column)
                    for (idx c = 0; c < nsubsys_bar; ++c)
                        midx_row[Csubsys_bar[c]] = midx_col[Csubsys_bar[c]] =
                                midx_bar[c];

                    // then the subsys part
                    for (idx c = 0; c < ngate; ++c)
                        midx_row[subsys[c]] = midxA_row[c];

                    // run over the subsys column multi-index
                    for (idx b = 0; b < DA; ++b)
                    {
                        // get the subsys column multi-index
                        internal::_n2multiidx(b, ngate, CdimsA, midxA_col);

                        // construct the result column multi-index
                        for (idx c = 0; c < ngate; ++c)
                            midx_col[subsys[c]] = midxA_col[c];

                        // finally write the values
                        result(internal::_multiidx2n(midx_row, n, Cdims),
                                internal::_multiidx2n(midx_col, n, Cdims))
                                = Ak(a, b);
                    }
                }

            }
        }

        return result;
    }

    /**
    * \brief Expands out
    * \see qpp::kron()
    *
    *  Expands out \a A as a matrix in a multi-partite system.
    *  Faster than using qpp::kron(I, I, ..., I, A, I, ..., I).
    *
    * \param A Eigen expression
    * \param pos Position
    * \param dims Dimensions of the multi-partite system
    * \return Tensor product
    * \f$ I\otimes\cdots\otimes I\otimes A \otimes I \otimes\cdots\otimes I\f$,
    * with \a A on position \a pos, as a dynamic matrix
    * over the same scalar field as \a A
    */
    template<typename Derived>
    dyn_mat<typename Derived::Scalar> expandout(
            const Eigen::MatrixBase<Derived>& A, idx pos,
            const std::vector<idx>& dims) const
    {
        const dyn_mat<typename Derived::Scalar>& rA = A;

        // check zero-size
        if (!internal::_check_nonzero_size(rA))
            throw Exception("qpp::Gates::expandout()",
                    Exception::Type::ZERO_SIZE);

        // check that dims is a valid dimension vector
        if (!internal::_check_dims(dims))
            throw Exception("qpp::Gates::expandout()",
                    Exception::Type::DIMS_INVALID);

        // check square matrix
        if (!internal::_check_square_mat(rA))
            throw Exception("qpp::Gates::expandout()",
                    Exception::Type::MATRIX_NOT_SQUARE);

        // check that position is valid
        if (pos > dims.size() - 1)
            throw Exception("qpp::Gates::expandout()",
                    Exception::Type::OUT_OF_RANGE);

        // check that dims[pos] match the dimension of A
        if (static_cast<idx>(rA.rows()) != dims[pos])
            throw Exception("qpp::Gates::expandout()",
                    Exception::Type::DIMS_MISMATCH_MATRIX);

        auto multiply = [](idx x, idx y) -> idx
        {
            return x * y;
        };

        idx D = std::accumulate(std::begin(dims), std::end(dims),
                1u, multiply);
        dyn_mat<typename Derived::Scalar> result = dyn_mat<
                typename Derived::Scalar>::Identity(D, D);

        idx Cdims[maxn];
        idx midx_row[maxn];
        idx midx_col[maxn];

        for (idx k = 0; k < dims.size(); ++k)
        {
            midx_row[k] = midx_col[k] = 0;
            Cdims[k] = dims[k];
        }

        // run over the main diagonal multi-indexes
        for (idx i = 0; i < D; ++i)
        {
            // get row multi_index
            internal::_n2multiidx(i, dims.size(), Cdims, midx_row);
            // get column multi_index (same as row)
            internal::_n2multiidx(i, dims.size(), Cdims, midx_col);
            // run over the gate row multi-index
            for (idx a = 0; a < static_cast<idx>(rA.rows());
                 ++a)
            {
                // construct the total row multi-index
                midx_row[pos] = a;

                // run over the gate column multi-index
                for (idx b = 0;
                     b < static_cast<idx>(rA.cols()); ++b)
                {
                    // construct the total column multi-index
                    midx_col[pos] = b;

                    // finally write the values
                    result(internal::_multiidx2n(midx_row, dims.size(), Cdims),
                            internal::_multiidx2n(midx_col, dims.size(), Cdims))
                            = rA(a, b);
                }
            }
        }

        return result;
    }

}; /* class Gates */

} /* namespace qpp */

#endif /* CLASSES_GATES_H_ */
