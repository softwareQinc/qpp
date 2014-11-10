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

#ifndef INCLUDE_INSTRUMENTS_H_
#define INCLUDE_INSTRUMENTS_H_

// measurements
namespace qpp
{

/**
* \brief Measures the state \a A using the set of Kraus operators \a Ks
*
* \param A Eigen expression
* \param Ks Set of Kraus operators
* \return Pair of vector of probabilities and vector of
* post-measurement normalized states
*/
    template<typename Derived>
    std::pair<std::vector<double>, std::vector<cmat>> measure(
            const Eigen::MatrixBase<Derived> &A, const std::vector<cmat> &Ks)
    {
        const DynMat<typename Derived::Scalar> &rA = A;

        // EXCEPTION CHECKS
        // check zero-size
        if (!internal::_check_nonzero_size(rA))
            throw Exception("measure", Exception::Type::ZERO_SIZE);

        // check the Kraus operators
        if (!internal::_check_nonzero_size(Ks))
            throw Exception("measure", Exception::Type::ZERO_SIZE);
        if (!internal::_check_square_mat(Ks[0]))
            throw Exception("measure", Exception::Type::MATRIX_NOT_SQUARE);
        if (Ks[0].rows() != rA.rows())
            throw Exception("measure", Exception::Type::DIMS_MISMATCH_MATRIX);
        for (auto &&it : Ks)
            if (it.rows() != Ks[0].rows() || it.cols() != Ks[0].rows())
                throw Exception("measure", Exception::Type::DIMS_NOT_EQUAL);
        // END EXCEPTION CHECKS

        // probabilities
        std::vector<double> prob(Ks.size());
        // resulting states
        std::vector<cmat> outstates(Ks.size());

        if (internal::_check_square_mat(rA)) // square matrix
        {
            for (std::size_t i = 0; i < Ks.size(); ++i)
            {
                cmat tmp;
                tmp = Ks[i] * rA * adjoint(Ks[i]); // un-normalized
                prob[i] = std::abs(trace(tmp)); // probability
                if (prob[i] > eps)
                    outstates[i] = tmp / prob[i]; // normalized
            }
        }
        else if (internal::_check_col_vector(rA)) // column vector
        {
            for (std::size_t i = 0; i < Ks.size(); ++i)
            {
                outstates[i] = ket::Zero(rA.rows());
                ket tmp;
                tmp = Ks[i] * rA; // un-normalized
                // probability
                prob[i] = std::abs((adjoint(tmp) * tmp).value());
                if (prob[i] > eps)
                    outstates[i] = tmp / std::sqrt(prob[i]); // normalized
            }
        }
        else
            throw Exception("measure",
                    Exception::Type::MATRIX_NOT_SQUARE_OR_CVECTOR);

        return std::make_pair(prob, outstates);

    }

// std::initializer_list overload, avoids ambiguity for 2-element lists, see
// http://stackoverflow.com/questions/26750039/ambiguity-when-using-initializer-list-as-parameter
/**
* \brief Measures the state \a A using the set of Kraus operators \a Ks
* (std::initializer_list overload)
*
* \param A Eigen expression
* \param Ks Set of Kraus operators
* \return Pair of vector of probabilities and vector of
* post-measurement normalized states
*/
    template<typename Derived>
    std::pair<std::vector<double>, std::vector<cmat>> measure(
            const Eigen::MatrixBase<Derived> &A,
            const std::initializer_list<cmat> &Ks)
    {
        return measure(A, std::vector<cmat>(Ks));
    }

/**
* \brief Measures the state \a A in the orthonormal basis
* specified by the eigenvectors of \a M.
*
* \param A Eigen expression
* \param M Normal matrix whose eigenvectors define the measurement basis
* \return Pair of vector of probabilities and vector of
* post-measurement normalized states
*/
    template<typename Derived>
    std::pair<std::vector<double>, std::vector<cmat>> measure(
            const Eigen::MatrixBase<Derived> &A, const cmat &M)
    {
        const DynMat<typename Derived::Scalar> &rA = A;

        // EXCEPTION CHECKS
        // check zero-size
        if (!internal::_check_nonzero_size(rA))
            throw Exception("measure", Exception::Type::ZERO_SIZE);

        // check the gate U
        if (!internal::_check_nonzero_size(M))
            throw Exception("measure", Exception::Type::ZERO_SIZE);
        if (!internal::_check_square_mat(M))
            throw Exception("measure", Exception::Type::MATRIX_NOT_SQUARE);
        if (M.rows() != rA.rows())
            throw Exception("measure", Exception::Type::DIMS_MISMATCH_MATRIX);
        // END EXCEPTION CHECKS

        // probabilities
        std::vector<double> prob(M.rows());
        // resulting states
        std::vector<cmat> outstates(M.rows());

        if (internal::_check_square_mat(rA)) // square matrix
        {
            for (std::size_t i = 0; i < static_cast<std::size_t>(M.rows()); ++i)
            {
                outstates[i] = cmat::Zero(rA.rows(), rA.rows());
                cmat tmp;
                // un-normalized
                tmp = prj((ket) evects(M).col(i)) * rA
                        * prj((ket) evects(M).col(i));
                prob[i] = std::abs(trace(tmp)); // probability
                if (prob[i] > eps)
                    outstates[i] = tmp / prob[i]; // normalized
            }
        }
        else if (internal::_check_col_vector(rA)) // column vector
        {
            for (std::size_t i = 0; i < static_cast<std::size_t>(M.rows()); ++i)
            {
                outstates[i] = ket::Zero(rA.rows());
                ket tmp;
                tmp = prj((ket) evects(M).col(i)) * rA;
                // probability
                prob[i] = std::abs((adjoint(tmp) * tmp).value());
                if (prob[i] > eps)
                    outstates[i] = tmp / std::sqrt(prob[i]); // normalized
            }
        }
        else
            throw Exception("measure",
                    Exception::Type::MATRIX_NOT_SQUARE_OR_CVECTOR);

        return std::make_pair(prob, outstates);
    }

}
/* namespace qpp */

#endif /* INCLUDE_INSTRUMENTS_H_ */
