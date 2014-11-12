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

#ifndef INCLUDE_ENTROPY_H_
#define INCLUDE_ENTROPY_H_

// various entropies, assume as input either
// a normalized Hermitian matrix or a probability vector

namespace qpp
{

    /**
    * \brief Shannon/von-Neumann entropy of the
    * probability distribution/density matrix \a A
    *
    * \param A Eigen expression, representing a probability distribution
    * (real dynamic column vector) or a density matrix (complex dynamic matrix)
    * \return Shannon/von-Neumann entropy, with the logarithm in base 2
    */
    template<typename Derived>
    double shannon(const Eigen::MatrixBase<Derived> &A)
    {
        const DynMat<typename Derived::Scalar> &rA = A;

        // check zero-size
        if (!internal::_check_nonzero_size(rA))
            throw Exception("qpp::shannon()", Exception::Type::ZERO_SIZE);

        // input is a vector
        if (internal::_check_vector(rA))
        {
            double result = 0;
            // take the absolut value to get rid of tiny negatives
            for (std::size_t i = 0; i < static_cast<std::size_t>(rA.size()); ++i)
                if (std::abs(rA(i)) != 0) // not identically zero
                    result -= std::abs(rA(i)) * std::log2(std::abs(rA(i)));

            return result;
        }

        // input is a matrix

        // check square matrix
        if (!internal::_check_square_mat(rA))
            throw Exception("qpp::shannon()", Exception::Type::MATRIX_NOT_SQUARE);

        // get the eigenvalues
        dmat ev = hevals(rA);
        double result = 0;
        // take the absolut value to get rid of tiny negatives
        for (std::size_t i = 0; i < static_cast<std::size_t>(ev.rows()); ++i)
            if (std::abs((cplx) ev(i)) != 0) // not identically zero
                result -= std::abs((cplx) ev(i))
                        * std::log2(std::abs((cplx) ev(i)));

        return result;
    }

    /**
    * \brief Renyi-\f$\alpha\f$ entropy of the
    * probability distribution/density matrix \a A, for \f$ \alpha\geq 0\f$
    *
    * \param A Eigen expression, representing a probability distribution
    * (real dynamic column vector) or a density matrix (complex dynamic matrix)
    * \param alpha Non-negative real number,
    * use qpp::infty for \f$\alpha = \infty\f$
    * \return Renyi-\f$\alpha\f$ entropy, with the logarithm in base 2
    */
    template<typename Derived>
    double renyi(const Eigen::MatrixBase<Derived> &A, double alpha)
    {
        const DynMat<typename Derived::Scalar> &rA = A;

        if (alpha < 0)
            throw Exception("qpp::renyi()", Exception::Type::OUT_OF_RANGE);

        if (alpha == 1) // Shannon/von Neumann
            return shannon(rA);

        // check zero-size
        if (!internal::_check_nonzero_size(rA))
            throw Exception("qpp::renyi()", Exception::Type::ZERO_SIZE);

        // input is a vector
        if (internal::_check_vector(rA))
        {
            if (alpha == 0) // H max
                return std::log2((double) rA.rows());

            if (alpha == infty) // H min
            {
                double max = 0;
                for (std::size_t i = 0; i < static_cast<std::size_t>(rA.size());
                     ++i)
                    if (std::abs(rA(i)) > max)
                        max = std::abs(rA(i));

                return -std::log2(max);
            }

            double result = 0;
            // take the absolut value to get rid of tiny negatives
            for (std::size_t i = 0; i < static_cast<std::size_t>(rA.size()); ++i)
                if (std::abs((cplx) rA(i)) != 0) // not identically zero
                    result += std::pow(std::abs(rA(i)), alpha);

            return std::log2(result) / (1 - alpha);
        }

        // input is a matrix

        // check square matrix
        if (!internal::_check_square_mat(rA))
            throw Exception("qpp::renyi()", Exception::Type::MATRIX_NOT_SQUARE);

        if (alpha == 0) // H max
            return -std::log2((double) rA.rows());

        if (alpha == infty) // H min
        {
            return -std::log2(svals(rA)[0]);
        }

        // get the singular values
        dmat sv = svals(rA);
        double result = 0;
        for (std::size_t i = 0; i < static_cast<std::size_t>(sv.rows()); ++i)
            if (sv(i) != 0) // not identically zero
                result += std::pow((double) sv(i), alpha);

        return std::log2(result) / (1 - alpha);
    }

    /**
    * \brief Tsallis-\f$\alpha\f$ entropy of the
    * probability distribution/density matrix \a A, for \f$ \alpha\geq 0\f$\n
    *
    * When \f$ \alpha\to 1\f$ the Tsallis entropy converges to the
    * Shannon/von-Neumann entropy, with the logarithm in base \f$ e \f$
    *
    * \param A Eigen expression, representing a probability distribution
    * (real dynamic column vector) or a density matrix (complex dynamic matrix)
    * \param alpha Non-negative real number
    *
    * \return Renyi-\f$\alpha\f$ entropy, with the logarithm in base 2
    */
    template<typename Derived>
    double tsallis(const Eigen::MatrixBase<Derived> &A, double alpha)
    {
        const DynMat<typename Derived::Scalar> &rA = A;

        if (alpha < 0)
            throw Exception("qpp::tsallis()", Exception::Type::OUT_OF_RANGE);

        if (alpha == 1) // Shannon/von Neumann with base e logarithm
            return shannon(rA) * std::log(2.);

        // check zero-size
        if (!internal::_check_nonzero_size(rA))
            throw Exception("qpp::tsallis()", Exception::Type::ZERO_SIZE);

        // input is a vector
        if (internal::_check_vector(rA))
        {
            double result = 0;
            // take the absolut value to get rid of tiny negatives
            for (std::size_t i = 0; i < static_cast<std::size_t>(rA.size()); ++i)
                if (std::abs((cplx) rA(i)) != 0) // not identically zero
                    result += std::pow(std::abs(rA(i)), alpha);

            return (result - 1) / (1 - alpha);
        }

        // input is a matrix

        // check square matrix
        if (!internal::_check_square_mat(rA))
            throw Exception("qpp::tsallis()", Exception::Type::MATRIX_NOT_SQUARE);

        // get the eigenvalues
        dmat ev = hevals(rA);
        double result = 0;
        // take the absolut values of the entries
        //of tiny negativesginary parts
        for (std::size_t i = 0; i < static_cast<std::size_t>(ev.rows()); ++i)
            if (std::abs((cplx) ev(i)) != 0) // not identically zero
                result += std::pow(std::abs((cplx) ev(i)), alpha);

        return (result - 1) / (1 - alpha);
    }

    /**
    * \brief Quantum mutual information between 2 subsystems of a composite system
    *
    * \param A Eigen expression
    * \param subsysA Indexes of the first subsystem
    * \param subsysB Indexes of the second subsystem
    * \param dims Dimensions of the multi-partite system
    * \return Mutual information between the 2 subsystems
    */
    template<typename Derived>
    double qmutualinfo(const Eigen::MatrixBase<Derived> &A,
            const std::vector<std::size_t> &subsysA,
            const std::vector<std::size_t> &subsysB,
            const std::vector<std::size_t> &dims)
    {
        const DynMat<typename Derived::Scalar> &rA = A;

        // error checks

        // check zero-size
        if (!internal::_check_nonzero_size(rA))
            throw Exception("qpp::mutualinfo()", Exception::Type::ZERO_SIZE);

        // check that dims is a valid dimension vector
        if (!internal::_check_dims(dims))
            throw Exception("qpp::mutualinfo()", Exception::Type::DIMS_INVALID);

        // check square matrix
        if (!internal::_check_square_mat(rA))
            throw Exception("qpp::mutualinfo()", Exception::Type::MATRIX_NOT_SQUARE);

        // check that dims match the dimension of A
        if (!internal::_check_dims_match_mat(dims, rA))
            throw Exception("qpp::mutualinfo()", Exception::Type::DIMS_MISMATCH_MATRIX);

        // check that subsys are valid
        if (!internal::_check_subsys_match_dims(subsysA, dims)
                || !internal::_check_subsys_match_dims(subsysB, dims))
            throw Exception("qpp::mutualinfo()", Exception::Type::SUBSYS_MISMATCH_DIMS);

        // The full system indexes {0,1,...,n-1}
        std::vector<std::size_t> full_system(dims.size());
        std::iota(std::begin(full_system), std::end(full_system), 0);

        // Sorted input subsystems
        std::vector<std::size_t> subsysAsorted{subsysA};
        std::vector<std::size_t> subsysBsorted{subsysB};

        // sort the input subsystems (as needed by std::set_difference)
        std::sort(std::begin(subsysAsorted), std::end(subsysAsorted));
        std::sort(std::begin(subsysBsorted), std::end(subsysBsorted));

        // construct the complement of subsys
        std::vector<std::size_t> subsysAbar;
        std::vector<std::size_t> subsysBbar;
        std::vector<std::size_t> subsysABbar;
        std::vector<std::size_t> subsysAB;

        std::set_difference(std::begin(full_system), std::end(full_system),
                std::begin(subsysAsorted), std::end(subsysAsorted),
                std::back_inserter(subsysAbar));
        std::set_difference(std::begin(full_system), std::end(full_system),
                std::begin(subsysBsorted), std::end(subsysBsorted),
                std::back_inserter(subsysBbar));
        std::set_union(std::begin(subsysAsorted), std::end(subsysAsorted),
                std::begin(subsysBsorted), std::end(subsysBsorted),
                std::back_inserter(subsysAB));
        std::sort(std::begin(subsysAB), std::end(subsysAB));

        std::set_difference(std::begin(full_system), std::end(full_system),
                std::begin(subsysAB), std::end(subsysAB),
                std::back_inserter(subsysABbar));

        cmat rhoA = ptrace(rA, subsysAbar, dims);
        cmat rhoB = ptrace(rA, subsysBbar, dims);
        cmat rhoAB = ptrace(rA, subsysABbar, dims);

        return shannon(rhoA) + shannon(rhoB) - shannon(rhoAB);
    }

} /* namespace qpp */

#endif /* INCLUDE_ENTROPY_H_ */
