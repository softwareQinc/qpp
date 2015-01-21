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
* \file entropies.h
* \brief Entropy functions
*/

#ifndef ENTROPY_H_
#define ENTROPY_H_

// various entropies, assume as input either
// a normalized Hermitian matrix or a probability vector

namespace qpp
{

/**
* \brief von-Neumann entropy of the density matrix \a A
*
* \param A Eigen expression
* \return von-Neumann entropy, with the logarithm in base 2
*/
template<typename Derived>
double entropy(const Eigen::MatrixBase<Derived>& A)
{
    const dyn_mat<typename Derived::Scalar>& rA = A;

    // check zero-size
    if (!internal::_check_nonzero_size(rA))
        throw Exception("qpp::entropy()", Exception::Type::ZERO_SIZE);

    // check square matrix
    if (!internal::_check_square_mat(rA))
        throw Exception("qpp::entropy()", Exception::Type::MATRIX_NOT_SQUARE);

    dmat ev = svals(rA); // get the singular values
    double result = 0;
    for (idx i = 0; i < static_cast<idx>(ev.rows()); ++i)
        if (ev(i) != 0) // not identically zero
            result -= ev(i) * std::log2(ev(i));

    return result;
}

/**
* \brief Shannon entropy of the probability distribution \a prob
*
* \param prob Real probability vector
* \return Shannon entropy, with the logarithm in base 2
*/
double entropy(const std::vector<double>& prob)
{
    // check zero-size
    if (!internal::_check_nonzero_size(prob))
        throw Exception("qpp::entropy()", Exception::Type::ZERO_SIZE);

    double result = 0;
    for (idx i = 0; i < prob.size(); ++i)
        if (std::abs(prob[i]) != 0) // not identically zero
            result -= std::abs(prob[i]) * std::log2(std::abs(prob[i]));

    return result;
}

/**
* \brief Renyi-\f$\alpha\f$ entropy of the density matrix \a A,
* for \f$\alpha\geq 0\f$
*
*\note When \f$ \alpha\to 1\f$ the Renyi entropy converges to the
* von-Neumann entropy, with the logarithm in base 2
*
* \param A Eigen expression
* \param alpha Non-negative real number,
* use qpp::infty for \f$\alpha = \infty\f$
* \return Renyi-\f$\alpha\f$ entropy, with the logarithm in base 2
*/
template<typename Derived>
double renyi(const Eigen::MatrixBase<Derived>& A, double alpha)
{
    const dyn_mat<typename Derived::Scalar>& rA = A;

    // check zero-size
    if (!internal::_check_nonzero_size(rA))
        throw Exception("qpp::renyi()", Exception::Type::ZERO_SIZE);

    // check square matrix
    if (!internal::_check_square_mat(rA))
        throw Exception("qpp::renyi()", Exception::Type::MATRIX_NOT_SQUARE);

    if (alpha < 0)
        throw Exception("qpp::renyi()", Exception::Type::OUT_OF_RANGE);

    if (alpha == 0) // H max
        return std::log2(static_cast<double>( rA.rows()));

    if (alpha == 1) // Shannon/von-Neumann
        return entropy(rA);

    if (alpha == infty) // H min
        return -std::log2(svals(rA)[0]);

    dmat sv = svals(rA); // get the singular values
    double result = 0;
    for (idx i = 0; i < static_cast<idx>(sv.rows()); ++i)
        if (sv(i) != 0) // not identically zero
            result += std::pow(sv(i), alpha);

    return std::log2(result) / (1 - alpha);
}

/**
* \brief Renyi-\f$\alpha\f$ entropy of the probability distribution \a prob,
* for \f$\alpha\geq 0\f$
*
*\note When \f$ \alpha\to 1\f$ the Renyi entropy converges to the
* Shannon entropy, with the logarithm in base 2
*
* \param prob Real probability vector
* \param alpha Non-negative real number,
* use qpp::infty for \f$\alpha = \infty\f$
* \return Renyi-\f$\alpha\f$ entropy, with the logarithm in base 2
*/
double renyi(const std::vector<double>& prob, double alpha)
{
    // check zero-size
    if (!internal::_check_nonzero_size(prob))
        throw Exception("qpp::renyi()", Exception::Type::ZERO_SIZE);

    if (alpha < 0)
        throw Exception("qpp::renyi()", Exception::Type::OUT_OF_RANGE);

    if (alpha == 0) // H max
        return std::log2(static_cast<double>( prob.size()));

    if (alpha == 1) // Shannon/von-Neumann
        return entropy(prob);

    if (alpha == infty) // H min
    {
        double max = 0;
        for (idx i = 0; i < prob.size(); ++i)
            if (std::abs(prob[i]) > max)
                max = std::abs(prob[i]);

        return -std::log2(max);
    }

    double result = 0;
    for (idx i = 0; i < prob.size(); ++i)
        if (std::abs(prob[i]) != 0) // not identically zero
            result += std::pow(std::abs(prob[i]), alpha);

    return std::log2(result) / (1 - alpha);
}

/**
* \brief Tsallis-\f$q\f$ entropy of the density matrix \a A,
* for \f$q\geq 0\f$
*
* \note When \f$ q\to 1\f$ the Tsallis entropy converges to the
* von-Neumann entropy, with the logarithm in base \f$ e \f$
*
* \param A Eigen expression
* \param q Non-negative real number
*
* \return Tsallis-\f$q\f$ entropy
*/
template<typename Derived>
double tsallis(const Eigen::MatrixBase<Derived>& A, double q)
{
    const dyn_mat<typename Derived::Scalar>& rA = A;

    // check zero-size
    if (!internal::_check_nonzero_size(rA))
        throw Exception("qpp::tsallis()", Exception::Type::ZERO_SIZE);

    // check square matrix
    if (!internal::_check_square_mat(rA))
        throw Exception("qpp::tsallis()", Exception::Type::MATRIX_NOT_SQUARE);

    if (q < 0)
        throw Exception("qpp::tsallis()", Exception::Type::OUT_OF_RANGE);

    if (q == 1) // Shannon/von-Neumann with base e logarithm
        return entropy(rA) * std::log(2.);

    dmat ev = svals(rA); // get the singular values
    double result = 0;
    for (idx i = 0; i < static_cast<idx>(ev.rows()); ++i)
        if (ev(i) != 0) // not identically zero
            result += std::pow(ev(i), q);

    return (result - 1) / (1 - q);
}

/**
* \brief Tsallis-\f$q\f$ entropy of the probability distribution \a prob,
* for \f$q\geq 0\f$
*
* \note When \f$ q\to 1\f$ the Tsallis entropy converges to the
* Shannon entropy, with the logarithm in base \f$ e \f$
*
* \param prob Real probability vector
* \param q Non-negative real number
*
* \return Tsallis-\f$q\f$ entropy
*/
double tsallis(const std::vector<double>& prob, double q)
{
    // check zero-size
    if (!internal::_check_nonzero_size(prob))
        throw Exception("qpp::tsallis()", Exception::Type::ZERO_SIZE);

    if (q < 0)
        throw Exception("qpp::tsallis()", Exception::Type::OUT_OF_RANGE);

    if (q == 1) // Shannon/von-Neumann with base e logarithm
        return entropy(prob) * std::log(2.);

    double result = 0;
    for (idx i = 0; i < prob.size(); ++i)
        if (std::abs(prob[i]) != 0) // not identically zero
            result += std::pow(std::abs(prob[i]), q);

    return (result - 1) / (1 - q);
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
double qmutualinfo(const Eigen::MatrixBase<Derived>& A,
        const std::vector<idx>& subsysA,
        const std::vector<idx>& subsysB,
        const std::vector<idx>& dims)
{
    const dyn_mat<typename Derived::Scalar>& rA = A;

    // error checks

    // check zero-size
    if (!internal::_check_nonzero_size(rA))
        throw Exception("qpp::qmutualinfo()", Exception::Type::ZERO_SIZE);

    // check that dims is a valid dimension vector
    if (!internal::_check_dims(dims))
        throw Exception("qpp::qmutualinfo()", Exception::Type::DIMS_INVALID);

    // check square matrix
    if (!internal::_check_square_mat(rA))
        throw Exception("qpp::qmutualinfo()",
                Exception::Type::MATRIX_NOT_SQUARE);

    // check that dims match the dimension of A
    if (!internal::_check_dims_match_mat(dims, rA))
        throw Exception("qpp::qmutualinfo()",
                Exception::Type::DIMS_MISMATCH_MATRIX);

    // check that subsys are valid
    if (!internal::_check_subsys_match_dims(subsysA, dims)
            || !internal::_check_subsys_match_dims(subsysB, dims))
        throw Exception("qpp::qmutualinfo()",
                Exception::Type::SUBSYS_MISMATCH_DIMS);

    // The full system indexes {0,1,...,n-1}
    std::vector<idx> full_system(dims.size());
    std::iota(std::begin(full_system), std::end(full_system), 0);

    // Sorted input subsystems
    std::vector<idx> subsysAsorted{subsysA};
    std::vector<idx> subsysBsorted{subsysB};

    // sort the input subsystems (as needed by std::set_difference)
    std::sort(std::begin(subsysAsorted), std::end(subsysAsorted));
    std::sort(std::begin(subsysBsorted), std::end(subsysBsorted));

    // construct the complement of subsys
    std::vector<idx> subsysAbar;
    std::vector<idx> subsysBbar;
    std::vector<idx> subsysABbar;
    std::vector<idx> subsysAB;

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

    return entropy(rhoA) + entropy(rhoB) - entropy(rhoAB);
}

/**
* \brief Quantum mutual information between 2 subsystems of a composite system
*
* \param A Eigen expression
* \param subsysA Indexes of the first subsystem
* \param subsysB Indexes of the second subsystem
* \param d Subsystem dimensions
* \return Mutual information between the 2 subsystems
*/
template<typename Derived>
double qmutualinfo(const Eigen::MatrixBase<Derived>& A,
        const std::vector<idx>& subsysA,
        const std::vector<idx>& subsysB,
        idx d = 2)
{
    const dyn_mat<typename Derived::Scalar>& rA = A;

    // error checks

    // check zero-size
    if (!internal::_check_nonzero_size(rA))
        throw Exception("qpp::qmutualinfo()", Exception::Type::ZERO_SIZE);

    idx n =
            static_cast<idx>(std::llround(std::log2(rA.rows()) /
                    std::log2(d)));
    std::vector<idx> dims(n, d); // local dimensions vector

    return qmutualinfo(rA, subsysA, subsysB, dims);
}

} /* namespace qpp */

#endif /* ENTROPY_H_ */
