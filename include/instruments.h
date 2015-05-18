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
* \file instruments.h
* \brief Measurement functions
*/

#ifndef INSTRUMENTS_H_
#define INSTRUMENTS_H_

// measurements
namespace qpp
{
// full measurements
/**
* \brief Measures the state \a A using the set of Kraus operators \a Ks
*
* \param A Eigen expression
* \param Ks Set of Kraus operators
* \return Tuple of: 1. Result of the measurement, 2.
* Vector of outcome probabilities, and 3. Vector of post-measurement
* normalized states
*/
template<typename Derived>
std::tuple<idx, std::vector < double>, std::vector<cmat>>

measure(const Eigen::MatrixBase<Derived>& A, const std::vector <cmat>& Ks)
{
    const dyn_mat<typename Derived::Scalar>& rA = A;

    // EXCEPTION CHECKS
    // check zero-size
    if (!internal::_check_nonzero_size(rA))
        throw Exception("qpp::measure()", Exception::Type::ZERO_SIZE);

    // check the Kraus operators
    if (Ks.size() == 0)
        throw Exception("qpp::measure()", Exception::Type::ZERO_SIZE);
    if (!internal::_check_square_mat(Ks[0]))
        throw Exception("qpp::measure()", Exception::Type::MATRIX_NOT_SQUARE);
    if (Ks[0].rows() != rA.rows())
        throw Exception("qpp::measure()",
                        Exception::Type::DIMS_MISMATCH_MATRIX);
    for (auto&& it : Ks)
        if (it.rows() != Ks[0].rows() || it.cols() != Ks[0].rows())
            throw Exception("qpp::measure()", Exception::Type::DIMS_NOT_EQUAL);
    // END EXCEPTION CHECKS

    // probabilities
    std::vector <double> prob(Ks.size());
    // resulting states
    std::vector <cmat> outstates(Ks.size());

    //************ density matrix ************//
    if (internal::_check_square_mat(rA)) // square matrix
    {
        for (idx i = 0; i < Ks.size(); ++i)
        {
            outstates[i] = cmat::Zero(rA.rows(), rA.rows());
            cmat tmp = Ks[i] * rA * adjoint(Ks[i]); // un-normalized;
            prob[i] = std::abs(trace(tmp)); // probability
            if (prob[i] > eps)
                outstates[i] = tmp / prob[i]; // normalized
        }
    }
        //************ ket ************//
    else if (internal::_check_cvector(rA)) // column vector
    {
        for (idx i = 0; i < Ks.size(); ++i)
        {
            outstates[i] = ket::Zero(rA.rows());
            ket tmp = Ks[i] * rA; // un-normalized;
            // probability
            prob[i] = std::pow(norm(tmp), 2);
            if (prob[i] > eps)
                outstates[i] = tmp / std::sqrt(prob[i]); // normalized
        }
    }
    else
        throw Exception("qpp::measure()",
                        Exception::Type::MATRIX_NOT_SQUARE_OR_CVECTOR);

    // sample from the probability distribution
    std::discrete_distribution<idx> dd(std::begin(prob),
                                       std::end(prob));
    idx result = dd(RandomDevices::get_instance()._rng);

    return std::make_tuple(result, prob, outstates);
}

// std::initializer_list overload, avoids ambiguity for 2-element lists, see
// http://stackoverflow.com
// /questions/26750039/ambiguity-when-using-initializer-list-as-parameter
/**
* \brief Measures the state \a A using the set of Kraus operators \a Ks
*
* \param A Eigen expression
* \param Ks Set of Kraus operators
* \return Tuple of: 1. Result of the measurement, 2.
* Vector of outcome probabilities, and 3. Vector of post-measurement
* normalized states
*/
template<typename Derived>
std::tuple<idx, std::vector < double>, std::vector<cmat>>

measure(const Eigen::MatrixBase<Derived>& A,
        const std::initializer_list<cmat>& Ks)
{
    return measure(A, std::vector<cmat>(Ks));
}

/**
* \brief Measures the state \a A in the orthonormal basis
* specified by the unitary matrix \a U
*
* \param A Eigen expression
* \param U Unitary matrix whose columns represent the measurement basis vectors
* \return Tuple of: 1. Result of the measurement, 2.
* Vector of outcome probabilities, and 3. Vector of post-measurement
* normalized states
*/
template<typename Derived>
std::tuple<idx, std::vector < double>, std::vector<cmat>>

measure(const Eigen::MatrixBase<Derived>& A, const cmat& U)
{
    const dyn_mat<typename Derived::Scalar>& rA = A;

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
        throw Exception("qpp::measure()",
                        Exception::Type::DIMS_MISMATCH_MATRIX);
    // END EXCEPTION CHECKS

    std::vector <cmat> Ks(U.rows());
    for (idx i = 0; i < static_cast<idx>(U.rows()); i++)
        Ks[i] = U.col(i) * adjoint(U.col(i));

    return measure(rA, Ks);
}

// partial measurements
/**
* \brief  Measures the part \a subsys of
* the multi-partite state vector or density matrix \a A
* using the set of Kraus operators \a Ks
* \see qpp::measure_seq()
*
* \note The dimension of all \a Ks must match the dimension of \a subsys.
* The measurement is destructive, i.e. the measured subsystems are traced away.
*
* \param A Eigen expression
* \param Ks Set of Kraus operators
* \param subsys Subsystem indexes that are measured
* \param dims Dimensions of the multi-partite system
* \return Tuple of: 1. Result of the measurement, 2.
* Vector of outcome probabilities, and 3. Vector of post-measurement
* normalized states
*/
template<typename Derived>
std::tuple<idx, std::vector < double>, std::vector<cmat>>

measure(const Eigen::MatrixBase<Derived>& A,
        const std::vector <cmat>& Ks,
        const std::vector <idx>& subsys,
        const std::vector <idx>& dims)
{
    const cmat& rA = A;

    // EXCEPTION CHECKS

    // check zero-size
    if (!internal::_check_nonzero_size(rA))
        throw Exception("qpp::measure()", Exception::Type::ZERO_SIZE);

    // check that dimension is valid
    if (!internal::_check_dims(dims))
        throw Exception("qpp::measure()", Exception::Type::DIMS_INVALID);

    // check that dims match rho matrix
    if (!internal::_check_dims_match_mat(dims, rA))
        throw Exception("qpp::measure()",
                        Exception::Type::DIMS_MISMATCH_MATRIX);

    // check subsys is valid w.r.t. dims
    if (!internal::_check_subsys_match_dims(subsys, dims))
        throw Exception("qpp::measure()",
                        Exception::Type::SUBSYS_MISMATCH_DIMS);

    std::vector <idx> subsys_dims(subsys.size());
    for (idx i = 0; i < subsys.size(); ++i)
        subsys_dims[i] = dims[subsys[i]];

    idx D = prod(std::begin(dims), std::end(dims));
    idx Dsubsys = prod(std::begin(subsys_dims), std::end(subsys_dims));
    idx Dbar = D / Dsubsys;

    // check the Kraus operators
    if (Ks.size() == 0)
        throw Exception("qpp::measure()",
                        Exception::Type::ZERO_SIZE);
    if (!internal::_check_square_mat(Ks[0]))
        throw Exception("qpp::measure()",
                        Exception::Type::MATRIX_NOT_SQUARE);
    if (Dsubsys != static_cast<idx>(Ks[0].rows()))
        throw Exception("qpp::measure()",
                        Exception::Type::DIMS_MISMATCH_MATRIX);
    for (auto&& it : Ks)
        if (it.rows() != Ks[0].rows() || it.cols() != Ks[0].rows())
            throw Exception("qpp::measure()",
                            Exception::Type::DIMS_NOT_EQUAL);
    // END EXCEPTION CHECKS

    // probabilities
    std::vector <double> prob(Ks.size());
    // resulting states
    std::vector <cmat> outstates(Ks.size());

    //************ density matrix ************//
    if (internal::_check_square_mat(rA)) // square matrix
    {
        for (idx i = 0; i < Ks.size(); ++i)
        {
            outstates[i] = cmat::Zero(Dbar, Dbar);
            cmat tmp = apply(rA, Ks[i], subsys, dims);
            tmp = ptrace(tmp, subsys, dims);
            prob[i] = std::abs(trace(tmp)); // probability
            if (prob[i] > eps)
            {
                // normalized output state
                // corresponding to measurement result i
                outstates[i] = tmp / prob[i];
            }
        }
    }
        //************ ket ************//
    else if (internal::_check_cvector(rA)) // column vector
    {
        for (idx i = 0; i < Ks.size(); ++i)
        {
            outstates[i] = cmat::Zero(Dbar, Dbar);
            ket tmp = apply(rA, Ks[i], subsys, dims);
            prob[i] = std::pow(norm(tmp), 2);
            if (prob[i] > eps)
            {
                // normalized output state
                // corresponding to measurement result i
                tmp /= std::sqrt(prob[i]);
                outstates[i] = ptrace(tmp, subsys, dims);
            }
        }
    }
    else
        throw Exception("qpp::measure()",
                        Exception::Type::MATRIX_NOT_SQUARE_OR_CVECTOR);

    // sample from the probability distribution
    std::discrete_distribution<idx> dd(std::begin(prob),
                                       std::end(prob));
    idx result = dd(RandomDevices::get_instance()._rng);

    return std::make_tuple(result, prob, outstates);
}

// std::initializer_list overload, avoids ambiguity for 2-element lists, see
// http://stackoverflow.com
// /questions/26750039/ambiguity-when-using-initializer-list-as-parameter
/**
* \brief  Measures the part \a subsys of
* the multi-partite state vector or density matrix \a A
* using the set of Kraus operators \a Ks
* \see qpp::measure_seq()
*
* \note The dimension of all \a Ks must match the dimension of \a subsys.
* The measurement is destructive, i.e. the measured subsystems are traced away.
*
* \param A Eigen expression
* \param subsys Subsystem indexes that are measured
* \param dims Dimensions of the multi-partite system
* \param Ks Set of Kraus operators
* \return Tuple of: 1. Result of the measurement, 2.
* Vector of outcome probabilities, and 3. Vector of post-measurement
* normalized states
*/
template<typename Derived>
std::tuple<idx, std::vector < double>, std::vector<cmat>>

measure(const Eigen::MatrixBase<Derived>& A,
        const std::initializer_list<cmat>& Ks,
        const std::vector <idx>& subsys,
        const std::vector <idx>& dims)
{
    return measure(A, std::vector<cmat>(Ks), subsys, dims);
}

/**
* \brief  Measures the part \a subsys of
* the multi-partite state vector or density matrix \a A
* using the set of Kraus operators \a Ks
* \see qpp::measure_seq()
*
* \note The dimension of all \a Ks must match the dimension of \a subsys.
* The measurement is destructive, i.e. the measured subsystems are traced away.
*
* \param A Eigen expression
* \param subsys Subsystem indexes that are measured
* \param d Subsystem dimensions
* \param Ks Set of Kraus operators
* \return Tuple of: 1. Result of the measurement, 2.
* Vector of outcome probabilities, and 3. Vector of post-measurement
* normalized states
*/
template<typename Derived>
std::tuple<idx, std::vector < double>, std::vector<cmat>>

measure(const Eigen::MatrixBase<Derived>& A,
        const std::vector <cmat>& Ks,
        const std::vector <idx>& subsys,
        const idx d = 2)
{
    const cmat& rA = A;

    // check zero size
    if (!internal::_check_nonzero_size(rA))
        throw Exception("qpp::measure()", Exception::Type::ZERO_SIZE);

    idx n =
            static_cast<idx>(std::llround(std::log2(rA.rows()) /
                                          std::log2(d)));
    std::vector <idx> dims(n, d); // local dimensions vector

    return measure(rA, Ks, subsys, dims);
}

// std::initializer_list overload, avoids ambiguity for 2-element lists, see
// http://stackoverflow.com
// /questions/26750039/ambiguity-when-using-initializer-list-as-parameter
/**
* \brief  Measures the part \a subsys of
* the multi-partite state vector or density matrix \a A
* using the set of Kraus operators \a Ks
* \see qpp::measure_seq()
*
* \note The dimension of all \a Ks must match the dimension of \a subsys.
* The measurement is destructive, i.e. the measured subsystems are traced away.
*
* \param A Eigen expression
* \param subsys Subsystem indexes that are measured
* \param d Subsystem dimensions
* \param Ks Set of Kraus operators
* \return Tuple of: 1. Result of the measurement, 2.
* Vector of outcome probabilities, and 3. Vector of post-measurement
* normalized states
*/
template<typename Derived>
std::tuple<idx, std::vector < double>, std::vector<cmat>>

measure(const Eigen::MatrixBase<Derived>& A,
        const std::initializer_list<cmat>& Ks,
        const std::vector <idx>& subsys,
        const idx d = 2)
{
    return measure(A, std::vector<cmat>(Ks), subsys, d);
}

/**
* \brief Measures the part \a subsys of
* the multi-partite state vector or density matrix \a A
* in the orthonormal basis specified by the unitary matrix \a U
* \see qpp::measure_seq()
*
* \note The dimension of \a U must match the dimension of \a subsys.
* The measurement is destructive, i.e. the measured subsystems are traced away.
*
* \param A Eigen expression
* \param subsys Subsystem indexes that are measured
* \param dims Dimensions of the multi-partite system
* \param U Unitary matrix whose columns represent the measurement basis vectors
* \return Tuple of: 1. Result of the measurement, 2.
* Vector of outcome probabilities, and 3. Vector of post-measurement
* normalized states
*/
template<typename Derived>
std::tuple<idx, std::vector < double>, std::vector<cmat>>

measure(const Eigen::MatrixBase<Derived>& A,
        const cmat& U,
        const std::vector <idx>& subsys,
        const std::vector <idx>& dims)
{
    const cmat& rA = A;

    // EXCEPTION CHECKS

    // check zero-size
    if (!internal::_check_nonzero_size(rA))
        throw Exception("qpp::measure()", Exception::Type::ZERO_SIZE);

    // check that dimension is valid
    if (!internal::_check_dims(dims))
        throw Exception("qpp::measure()", Exception::Type::DIMS_INVALID);

    // check that dims match rho matrix
    if (!internal::_check_dims_match_mat(dims, rA))
        throw Exception("qpp::measure()",
                        Exception::Type::DIMS_MISMATCH_MATRIX);

    // check subsys is valid w.r.t. dims
    if (!internal::_check_subsys_match_dims(subsys, dims))
        throw Exception("qpp::measure()",
                        Exception::Type::SUBSYS_MISMATCH_DIMS);

    std::vector <idx> subsys_dims(subsys.size());
    for (idx i = 0; i < subsys.size(); ++i)
        subsys_dims[i] = dims[subsys[i]];

    idx Dsubsys = prod(std::begin(subsys_dims), std::end(subsys_dims));

    // check the unitary basis matrix U
    if (!internal::_check_nonzero_size(U))
        throw Exception("qpp::measure()", Exception::Type::ZERO_SIZE);
    if (!internal::_check_square_mat(U))
        throw Exception("qpp::measure()", Exception::Type::MATRIX_NOT_SQUARE);
    if (Dsubsys != static_cast<idx>(U.rows()))
        throw Exception("qpp::measure()",
                        Exception::Type::DIMS_MISMATCH_MATRIX);
    // END EXCEPTION CHECKS

    std::vector <cmat> Ks(U.rows());
    for (idx i = 0; i < static_cast<idx>(U.rows()); i++)
        Ks[i] = U.col(i) * adjoint(U.col(i));

    return measure(rA, Ks, subsys, dims);
}

/**
* \brief Measures the part \a subsys of
* the multi-partite state vector or density matrix \a A\ in the orthonormal
* basis specified by the unitary matrix \a U
* \see qpp::measure_seq()
*
* \note The dimension of \a U must match the dimension of \a subsys.
* The measurement is destructive, i.e. the measured subsystems are traced away.
*
* \param A Eigen expression
* \param subsys Subsystem indexes that are measured
* \param d Subsystem dimensions
* \param U Unitary matrix whose columns represent the measurement basis vectors
* \return Tuple of: 1. Result of the measurement, 2.
* Vector of outcome probabilities, and 3. Vector of post-measurement
* normalized states
*/
template<typename Derived>
std::tuple<idx, std::vector < double>, std::vector<cmat>>

measure(const Eigen::MatrixBase<Derived>& A,
        const cmat& U,
        const std::vector <idx>& subsys,
        const idx d = 2)
{
    const cmat& rA = A;

    // check zero size
    if (!internal::_check_nonzero_size(rA))
        throw Exception("qpp::measure()", Exception::Type::ZERO_SIZE);

    idx n =
            static_cast<idx>(std::llround(std::log2(rA.rows()) /
                                          std::log2(d)));
    std::vector <idx> dims(n, d); // local dimensions vector

    return measure(rA, U, subsys, dims);
}

/**
* \brief Sequentially measures the part \a subsys
* of the multi-partite state vector or density matrix \a A
* in the computational basis
* \see qpp::measure()
*
* \param A Eigen expression
* \param subsys Subsystem indexes that are measured
* \param dims Dimensions of the multi-partite system
* \return Tuple of: 1. Vector of outcome results of the
* measurement (ordered in increasing order with respect to \a subsys, i.e. first
* measurement result corresponds to the subsystem with the smallest index), 2.
* Outcome probability, and 3. Post-measurement normalized state
*/
template<typename Derived>
std::tuple<std::vector < idx>, double, cmat>

measure_seq(const Eigen::MatrixBase<Derived>& A,
            std::vector <idx> subsys,
            std::vector <idx> dims)
{
    dyn_mat<typename Derived::Scalar> cA = A;

    // EXCEPTION CHECKS

    // check zero-size
    if (!internal::_check_nonzero_size(cA))
        throw Exception("qpp::measure_seq()", Exception::Type::ZERO_SIZE);

    // check that dimension is valid
    if (!internal::_check_dims(dims))
        throw Exception("qpp::measure_seq()", Exception::Type::DIMS_INVALID);

    // check that dims match rho matrix
    if (!internal::_check_dims_match_mat(dims, cA))
        throw Exception("qpp::measure_seq()",
                        Exception::Type::DIMS_MISMATCH_MATRIX);

    // check subsys is valid w.r.t. dims
    if (!internal::_check_subsys_match_dims(subsys, dims))
        throw Exception("qpp::measure_seq()",
                        Exception::Type::SUBSYS_MISMATCH_DIMS);

    // END EXCEPTION CHECKS

    std::vector <idx> result;
    double prob = 1;

    // sort subsys in decreasing order,
    // the order of measurements does not matter
    std::sort(std::begin(subsys), std::end(subsys), std::greater<idx>{});

    //************ density matrix or column vector ************//
    if (internal::_check_square_mat(cA) || internal::_check_cvector(cA))
    {
        while (subsys.size() > 0)
        {
            auto tmp = measure(cA, Gates::get_instance().Id(dims[subsys[0]]),
                               {subsys[0]}, dims);
            result.push_back(std::get<0>(tmp));
            prob *= std::get<1>(tmp)[std::get<0>(tmp)];
            cA = std::get<2>(tmp)[std::get<0>(tmp)];

            // remove the subsystem
            dims.erase(std::next(std::begin(dims), subsys[0]));
            subsys.erase(std::begin(subsys));
        }
        // order result in increasing order with respect to subsys
        std::reverse(std::begin(result), std::end(result));
    }
    else
        throw Exception("qpp::measure_seq()",
                        Exception::Type::MATRIX_NOT_SQUARE_OR_CVECTOR);

    return std::make_tuple(result, prob, cA);
}

/**
* \brief Sequentially measures the part \a subsys
* of the multi-partite state vector or density matrix \a A
* in the computational basis
* \see qpp::measure()
*
* \param A Eigen expression
* \param subsys Subsystem indexes that are measured
* \param d Subsystem dimensions
* \return Tuple of: 1. Vector of outcome results of the
* measurement (ordered in increasing order with respect to \a subsys, i.e. first
* measurement result corresponds to the subsystem with the smallest index), 2.
* Outcome probability, and 3. Post-measurement normalized state
*/
template<typename Derived>
std::tuple<std::vector < idx>, double, cmat>

measure_seq(const Eigen::MatrixBase<Derived>& A,
            std::vector <idx> subsys, idx d = 2)
{
    const cmat& rA = A;

    // check zero size
    if (!internal::_check_nonzero_size(rA))
        throw Exception("qpp::measure_seq()", Exception::Type::ZERO_SIZE);

    idx n =
            static_cast<idx>(std::llround(std::log2(rA.rows()) /
                                          std::log2(d)));
    std::vector <idx> dims(n, d); // local dimensions vector

    return measure_seq(rA, subsys, dims);
}

} /* namespace qpp */

#endif /* INSTRUMENTS_H_ */
