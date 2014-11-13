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
// partial measurements
/**
* \brief  Measures the part \a subsys of the multi-partite state vector or density matrix \a A
* using the set of Kraus operators \a Ks
*
* \note The dimension of all \a Ks must match the dimension of \a subsys.
*
* \param A Eigen expression
* \param Ks Set of Kraus operators
* \param subsys Subsystem indexes that are measured
* \param dims Dimensions of the multi-partite system
* \return Tuple consisiting of the result of the measurement, the vector of outcome probabilities
* and the vector of post-measurement normalized states
*/
template<typename Derived>
std::tuple<std::size_t, std::vector<double>, std::vector<cmat>>
measure(
        const Eigen::MatrixBase<Derived> &A,
        const std::vector<cmat> &Ks,
        const std::vector<std::size_t> &subsys,
        const std::vector<std::size_t> &dims)
{
    const cmat &rA = A;

    // EXCEPTION CHECKS

    // check zero-size
    if (!internal::_check_nonzero_size(rA))
        throw Exception("qpp::measure()", Exception::Type::ZERO_SIZE);

    // check that dimension is valid
    if (!internal::_check_dims(dims))
        throw Exception("qpp::measure()", Exception::Type::DIMS_INVALID);

    // check that dims match rho matrix
    if (!internal::_check_dims_match_mat(dims, rA))
        throw Exception("qpp::measure()", Exception::Type::DIMS_MISMATCH_MATRIX);

    // check subsys is valid w.r.t. dims
    if (!internal::_check_subsys_match_dims(subsys, dims))
        throw Exception("qpp::measure()", Exception::Type::SUBSYS_MISMATCH_DIMS);

    std::vector<std::size_t> subsys_dims(subsys.size());
    for (std::size_t i = 0; i < subsys.size(); ++i)
        subsys_dims[i] = dims[subsys[i]];

    std::size_t D = 1;
    for (auto &&it: dims)
        D *= it;

    std::size_t Dsubsys = 1;
    for (auto &&it: subsys_dims)
        Dsubsys *= it;

    std::size_t Dbar = D / Dsubsys;

    // check the Kraus operators
    if (!internal::_check_nonzero_size(Ks))
        throw Exception("qpp::destructive_measure()", Exception::Type::ZERO_SIZE);
    if (!internal::_check_square_mat(Ks[0]))
        throw Exception("qpp::destructive_measure()", Exception::Type::MATRIX_NOT_SQUARE);
    if (Dsubsys != static_cast<std::size_t>(Ks[0].rows()))
        throw Exception("qpp::destructive_measure()", Exception::Type::DIMS_MISMATCH_MATRIX);
    for (auto &&it : Ks)
        if (it.rows() != Ks[0].rows() || it.cols() != Ks[0].rows())
            throw Exception("qpp::destructive_measure()", Exception::Type::DIMS_NOT_EQUAL);
    // END EXCEPTION CHECKS

    // probabilities
    std::vector<double> prob(Ks.size());
    // resulting states
    std::vector<cmat> outstates(Ks.size());

    if (internal::_check_square_mat(rA)) // square matrix
    {
        for (std::size_t i = 0; i < Ks.size(); ++i)
        {
            outstates[i] = cmat::Zero(Dbar, Dbar);
            cmat tmp = apply(rA, Ks[i], subsys, dims);
            tmp = ptrace(tmp, subsys, dims);
            prob[i] = std::abs(trace(tmp)); // probability
            if (prob[i] > eps)
                outstates[i] = tmp / prob[i]; // normalized output state corresponding to measurement result i
        }
    }
    else if (internal::_check_col_vector(rA)) // column vector
    {
        for (std::size_t i = 0; i < Ks.size(); ++i)
        {
            outstates[i] = cmat::Zero(Dbar, Dbar);
            cmat tmp = apply(rA, Ks[i], subsys, dims);
            prob[i] = std::pow(norm(tmp), 2);
            if (prob[i] > eps)
                outstates[i] = ptrace(prj(tmp), subsys, dims);
        }
    }
    else
        throw Exception("qpp::measure()",
                Exception::Type::MATRIX_NOT_SQUARE_OR_CVECTOR);

    // sample from the probability distribution
    std::discrete_distribution<std::size_t> dd(std::begin(prob), std::end(prob));
    std::size_t result = dd(RandomDevices::get_instance()._rng);

    return std::make_tuple(result, prob, outstates);
}

// std::initializer_list overload, avoids ambiguity for 2-element lists, see
// http://stackoverflow.com/questions/26750039/ambiguity-when-using-initializer-list-as-parameter
/**
* \brief  Measures the part \a subsys of the multi-partite state vector or density matrix \a A
* using the set of Kraus operators \a Ks
*
* \note The dimension of all \a Ks must match the dimension of \a subsys.
*
* \param A Eigen expression
* \param subsys Subsystem indexes that are measured
* \param dims Dimensions of the multi-partite system
* \param Ks Set of Kraus operators
* \return Tuple consisiting of the result of the measurement, the vector of outcome probabilities
* and the vector of post-measurement normalized states
*/
template<typename Derived>
std::tuple<std::size_t, std::vector<double>, std::vector<cmat>>
measure(
        const Eigen::MatrixBase<Derived> &A,
        const std::initializer_list<cmat> &Ks,
        const std::vector<std::size_t> &subsys,
        const std::vector<std::size_t> &dims)
{
    return measure(A, std::vector<cmat>(Ks), subsys, dims);
}


/**
* \brief  Measures the part \a subsys of the multi-partite state vector or density matrix \a A
* using the set of Kraus operators \a Ks
*
* \note The dimension of all \a Ks must match the dimension of \a subsys.
*
* \param A Eigen expression
* \param subsys Subsystem indexes that are measured
* \param d Subsystem dimensions
* \param Ks Set of Kraus operators
* \return Tuple consisiting of the result of the measurement, the vector of outcome probabilities
* and the vector of post-measurement normalized states
*/
template<typename Derived>
std::tuple<std::size_t, std::vector<double>, std::vector<cmat>>
measure(
        const Eigen::MatrixBase<Derived> &A,
        const std::vector<cmat> &Ks,
        const std::vector<std::size_t> &subsys,
        const std::size_t d = 2)
{
    const cmat &rA = A;

    // check zero size
    if (!internal::_check_nonzero_size(rA))
        throw Exception("qpp::measure()", Exception::Type::ZERO_SIZE);

    std::size_t n = static_cast<std::size_t>(std::log2(rA.rows()) / std::log2(d));
    std::vector<std::size_t> dims(n, d); // local dimensions vector

    return measure(rA, Ks, subsys, dims);
}

// std::initializer_list overload, avoids ambiguity for 2-element lists, see
// http://stackoverflow.com/questions/26750039/ambiguity-when-using-initializer-list-as-parameter
/**
* \brief  Measures the part \a subsys of the multi-partite state vector or density matrix \a A
* using the set of Kraus operators \a Ks
*
* \note The dimension of all \a Ks must match the dimension of \a subsys.
*
* \param A Eigen expression
* \param subsys Subsystem indexes that are measured
* \param d Subsystem dimensions
* \param Ks Set of Kraus operators
* \return Tuple consisiting of the result of the measurement, the vector of outcome probabilities
* and the vector of post-measurement normalized states
*/
template<typename Derived>
std::tuple<std::size_t, std::vector<double>, std::vector<cmat>>
measure(
        const Eigen::MatrixBase<Derived> &A,
        const std::initializer_list<cmat> &Ks,
        const std::vector<std::size_t> &subsys,
        const std::size_t d = 2)
{
    return measure(A, std::vector<cmat>(Ks), subsys, d);
}

/**
* \brief Measures the part \a subsys of the multi-partite state \a A in the orthonormal basis
* specified by the unitary matrix \a U
*
* \note The dimension of \a U must match the dimension of \a subsys.
*
* \param A Eigen expression
* \param subsys Subsystem indexes that are measured
* \param dims Dimensions of the multi-partite system
* \param U Unitary matrix whose columns represent the measurement basis vectors
* \return Tuple consisiting of the result of the measurement, the vector of outcome probabilities
* and the vector of post-measurement normalized states
*/
template<typename Derived>
std::tuple<std::size_t, std::vector<double>, std::vector<cmat>>
measure(
        const Eigen::MatrixBase<Derived> &A,
        const cmat &U,
        const std::vector<std::size_t> &subsys,
        const std::vector<std::size_t> &dims)
{
    const cmat &rA = A;

    // EXCEPTION CHECKS

    // check zero-size
    if (!internal::_check_nonzero_size(rA))
        throw Exception("qpp::measure()", Exception::Type::ZERO_SIZE);

    // check that dimension is valid
    if (!internal::_check_dims(dims))
        throw Exception("qpp::measure()", Exception::Type::DIMS_INVALID);

    // check that dims match rho matrix
    if (!internal::_check_dims_match_mat(dims, rA))
        throw Exception("qpp::measure()", Exception::Type::DIMS_MISMATCH_MATRIX);

    // check subsys is valid w.r.t. dims
    if (!internal::_check_subsys_match_dims(subsys, dims))
        throw Exception("qpp::measure()", Exception::Type::SUBSYS_MISMATCH_DIMS);

    std::vector<std::size_t> subsys_dims(subsys.size());
    for (std::size_t i = 0; i < subsys.size(); ++i)
        subsys_dims[i] = dims[subsys[i]];

    std::size_t Dsubsys = 1;
    for (auto &&it: subsys_dims)
        Dsubsys *= it;

    // check the unitary basis matrix U
    if (!internal::_check_nonzero_size(U))
        throw Exception("qpp::measure()", Exception::Type::ZERO_SIZE);
    if (!internal::_check_square_mat(U))
        throw Exception("qpp::measure()", Exception::Type::MATRIX_NOT_SQUARE);
    if (Dsubsys != static_cast<std::size_t>(U.rows()))
        throw Exception("qpp::measure()", Exception::Type::DIMS_MISMATCH_MATRIX);
    // END EXCEPTION CHECKS

    std::vector<cmat> Ks(U.rows());
    for (std::size_t i = 0; i < static_cast<std::size_t>(U.rows()); i++)
        Ks[i] = U.col(i) * adjoint(U.col(i));

    return measure(rA, Ks, subsys, dims);
}

/**
* \brief Measures the part \a subsys of the multi-partite state \a A in the orthonormal basis
* specified by the unitary matrix \a U
*
* \note The dimension of \a U must match the dimension of \a subsys.
*
* \param A Eigen expression
* \param subsys Subsystem indexes that are measured
* \param d Subsystem dimensions
* \param U Unitary matrix whose columns represent the measurement basis vectors
* \return Tuple consisiting of the result of the measurement, the vector of outcome probabilities
* and the vector of post-measurement normalized states
*/
template<typename Derived>
std::tuple<std::size_t, std::vector<double>, std::vector<cmat>>
measure(
        const Eigen::MatrixBase<Derived> &A,
        const cmat &U,
        const std::vector<std::size_t> &subsys,
        const std::size_t d = 2)
{
    const cmat &rA = A;

    // check zero size
    if (!internal::_check_nonzero_size(rA))
        throw Exception("qpp::measure()", Exception::Type::ZERO_SIZE);

    std::size_t n = static_cast<std::size_t>(std::log2(rA.rows()) / std::log2(d));
    std::vector<std::size_t> dims(n, d); // local dimensions vector

    return measure(rA, U, subsys, dims);
}

// full measurements
/**
* \brief Measures the state \a A using the set of Kraus operators \a Ks
*
* \param A Eigen expression
* \param Ks Set of Kraus operators
* \return Tuple consisiting of the result of the measurement, the vector of outcome probabilities
* and the vector of post-measurement normalized states
*/
template<typename Derived>
std::tuple<std::size_t, std::vector<double>, std::vector<cmat>> measure(
        const Eigen::MatrixBase<Derived> &A, const std::vector<cmat> &Ks)
{
    const DynMat<typename Derived::Scalar> &rA = A;

    // EXCEPTION CHECKS
    // check zero-size
    if (!internal::_check_nonzero_size(rA))
        throw Exception("qpp::measure()", Exception::Type::ZERO_SIZE);

    // check the Kraus operators
    if (!internal::_check_nonzero_size(Ks))
        throw Exception("qpp::measure()", Exception::Type::ZERO_SIZE);
    if (!internal::_check_square_mat(Ks[0]))
        throw Exception("qpp::measure()", Exception::Type::MATRIX_NOT_SQUARE);
    if (Ks[0].rows() != rA.rows())
        throw Exception("qpp::measure()", Exception::Type::DIMS_MISMATCH_MATRIX);
    for (auto &&it : Ks)
        if (it.rows() != Ks[0].rows() || it.cols() != Ks[0].rows())
            throw Exception("qpp::measure()", Exception::Type::DIMS_NOT_EQUAL);
    // END EXCEPTION CHECKS

    // probabilities
    std::vector<double> prob(Ks.size());
    // resulting states
    std::vector<cmat> outstates(Ks.size());

    if (internal::_check_square_mat(rA)) // square matrix
    {
        for (std::size_t i = 0; i < Ks.size(); ++i)
        {
            outstates[i] = cmat::Zero(rA.rows(), rA.rows());
            cmat tmp = Ks[i] * rA * adjoint(Ks[i]); // un-normalized;
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
            ket tmp = Ks[i] * rA; // un-normalized;
            // probability
            prob[i] = std::abs((adjoint(tmp) * tmp).value());
            if (prob[i] > eps)
                outstates[i] = tmp / std::sqrt(prob[i]); // normalized
        }
    }
    else
        throw Exception("qpp::measure()",
                Exception::Type::MATRIX_NOT_SQUARE_OR_CVECTOR);

    // sample from the probability distribution
    std::discrete_distribution<std::size_t> dd(std::begin(prob), std::end(prob));
    std::size_t result = dd(RandomDevices::get_instance()._rng);

    return std::make_tuple(result, prob, outstates);
}

// std::initializer_list overload, avoids ambiguity for 2-element lists, see
// http://stackoverflow.com/questions/26750039/ambiguity-when-using-initializer-list-as-parameter
/**
* \brief Measures the state \a A using the set of Kraus operators \a Ks
*
* \param A Eigen expression
* \param Ks Set of Kraus operators
* \return Tuple consisiting of the result of the measurement, the vector of outcome probabilities
* and the vector of post-measurement normalized states
*/
template<typename Derived>
std::tuple<std::size_t, std::vector<double>, std::vector<cmat>> measure(
        const Eigen::MatrixBase<Derived> &A,
        const std::initializer_list<cmat> &Ks)
{
    return measure(A, std::vector<cmat>(Ks));
}

/**
* \brief Measures the state \a A in the orthonormal basis specified by the unitary matrix \a U
*
* \param A Eigen expression
* \param U Unitary matrix whose columns represent the measurement basis vectors
* \return Tuple consisiting of the result of the measurement, the vector of outcome probabilities
* and the vector of post-measurement normalized states
*/
template<typename Derived>
std::tuple<std::size_t, std::vector<double>, std::vector<cmat>> measure(
        const Eigen::MatrixBase<Derived> &A, const cmat &U)
{
    const DynMat<typename Derived::Scalar> &rA = A;

    // EXCEPTION CHECKS
    // check zero-size
    if (!internal::_check_nonzero_size(rA))
        throw Exception("qpp::measure()", Exception::Type::ZERO_SIZE);

    // check the unitary basis matrix U
    if (!internal::_check_nonzero_size(U))
        throw Exception("qpp::measure()", Exception::Type::ZERO_SIZE);
    if (!internal::_check_square_mat(U))
        throw Exception("qpp::measure()", Exception::Type::MATRIX_NOT_SQUARE);
    if (U.rows() != rA.rows())
        throw Exception("qpp::measure()", Exception::Type::DIMS_MISMATCH_MATRIX);
    // END EXCEPTION CHECKS

    std::vector<cmat> Ks(U.rows());
    for (std::size_t i = 0; i < static_cast<std::size_t>(U.rows()); i++)
        Ks[i] = U.col(i) * adjoint(U.col(i));

    return measure(rA, Ks);
}
}
/* namespace qpp */

#endif /* INCLUDE_INSTRUMENTS_H_ */
