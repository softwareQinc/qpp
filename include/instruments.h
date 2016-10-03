/*
 * Quantum++
 *
 * Copyright (c) 2013 - 2016 Vlad Gheorghiu (vgheorgh@gmail.com)
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

// instruments
namespace qpp
{
/**
* \brief Generalized inner product
*
* \param phi Column vector Eigen expression
* \param psi Column vector Eigen expression
* \param subsys Subsystem indexes over which \a phi is defined
* \param dims Dimensions of the multi-partite system
* \return The inner product \f$\langle \phi_{subsys}|\psi\rangle\f$, as a scalar
* or column vector over the remaining Hilbert space
*/
template<typename Derived>
dyn_col_vect<typename Derived::Scalar> ip(
        const Eigen::MatrixBase<Derived>& phi,
        const Eigen::MatrixBase<Derived>& psi,
        const std::vector<idx>& subsys,
        const std::vector<idx>& dims)
{
    const dyn_col_vect<typename Derived::Scalar>& rphi = phi.derived();
    const dyn_col_vect<typename Derived::Scalar>& rpsi = psi.derived();

    // EXCEPTION CHECKS

    // check zero-size
    if ( !internal::_check_nonzero_size(rphi))
        throw Exception("qpp::ip()", Exception::Type::ZERO_SIZE);

    // check zero-size
    if ( !internal::_check_nonzero_size(rpsi))
        throw Exception("qpp::ip()", Exception::Type::ZERO_SIZE);

    // check column vector
    if ( !internal::_check_cvector(rphi))
        throw Exception("qpp::schmidtcoeffs()",
                        Exception::Type::MATRIX_NOT_CVECTOR);

    // check column vector
    if ( !internal::_check_cvector(rpsi))
        throw Exception("qpp::schmidtcoeffs()",
                        Exception::Type::MATRIX_NOT_CVECTOR);

    // check that dims is a valid dimension vector
    if ( !internal::_check_dims(dims))
        throw Exception("qpp::ip()", Exception::Type::DIMS_INVALID);

    // check that subsys are valid
    if ( !internal::_check_subsys_match_dims(subsys, dims))
        throw Exception("qpp::ip()",
                        Exception::Type::SUBSYS_MISMATCH_DIMS);

    // check that dims match state vector psi
    if ( !internal::_check_dims_match_cvect(dims, rpsi))
        throw Exception("qpp::ip()",
                        Exception::Type::DIMS_MISMATCH_CVECTOR);

    // check that subsys match state vector phi
    std::vector<idx> subsys_dims(subsys.size());
    for ( idx i = 0; i < subsys.size(); ++i )
        subsys_dims[i] = dims[subsys[i]];
    if ( !internal::_check_dims_match_cvect(subsys_dims, rphi))
        throw Exception("qpp::ip()",
                        Exception::Type::DIMS_MISMATCH_CVECTOR);
    // END EXCEPTION CHECKS

    idx Dsubsys = prod(std::begin(subsys_dims), std::end(subsys_dims));

    idx D = static_cast<idx>(rpsi.rows());
    idx Dsubsysbar = D / Dsubsys;

    idx n = dims.size();
    idx nsubsys = subsys.size();
    idx nsubsysbar = n - nsubsys;

    idx Cdims[maxn];
    idx Csubsys[maxn];
    idx Cdimssubsys[maxn];
    idx Csubsysbar[maxn];
    idx Cdimssubsysbar[maxn];

    std::vector<idx> subsys_bar = complement(subsys, n);
    std::copy(std::begin(subsys_bar), std::end(subsys_bar),
              std::begin(Csubsysbar));

    for ( idx i = 0; i < n; ++i )
    {
        Cdims[i] = dims[i];
    }
    for ( idx i = 0; i < nsubsys; ++i )
    {
        Csubsys[i] = subsys[i];
        Cdimssubsys[i] = dims[subsys[i]];
    }
    for ( idx i = 0; i < nsubsysbar; ++i )
    {
        Cdimssubsysbar[i] = dims[subsys_bar[i]];
    }

    auto worker = [=](idx b) noexcept
            -> typename Derived::Scalar
    {
        idx Cmidxrow[maxn];
        idx Cmidxrowsubsys[maxn];
        idx Cmidxcolsubsysbar[maxn];

        /* get the col multi-indexes of the complement */
        internal::_n2multiidx(b, nsubsysbar,
                              Cdimssubsysbar, Cmidxcolsubsysbar);
        /* write it in the global row multi-index */
        for ( idx k = 0; k < nsubsysbar; ++k )
        {
            Cmidxrow[Csubsysbar[k]] = Cmidxcolsubsysbar[k];
        }

        typename Derived::Scalar result = 0;
        for ( idx a = 0; a < Dsubsys; ++a )
        {
            /* get the row multi-indexes of the subsys */
            internal::_n2multiidx(a, nsubsys,
                                  Cdimssubsys, Cmidxrowsubsys);
            /* write it in the global row multi-index */
            for ( idx k = 0; k < nsubsys; ++k )
            {
                Cmidxrow[Csubsys[k]] = Cmidxrowsubsys[k];
            }
            // compute the row index
            idx i = internal::_multiidx2n(Cmidxrow, n, Cdims);

            result += std::conj(rphi(a)) * rpsi(i);
        }

        return result;
    }; /* end worker */

    dyn_col_vect<typename Derived::Scalar> result(Dsubsysbar);
#ifdef _WITH_OPENMP_
#pragma omp parallel for
#endif
    for ( idx m = 0; m < Dsubsysbar; ++m )
        result(m) = worker(m);

    return result;
}

/**
* \brief Generalized inner product
*
* \param phi Column vector Eigen expression
* \param psi Column vector Eigen expression
* \param subsys Subsystem indexes over which \a phi is defined
* \param d Subsystem dimensions
* \return The inner product \f$\langle \phi_{subsys}|\psi\rangle\f$, as a scalar
* or column vector over the remaining Hilbert space
*/
template<typename Derived>
dyn_col_vect<typename Derived::Scalar> ip(
        const Eigen::MatrixBase<Derived>& phi,
        const Eigen::MatrixBase<Derived>& psi,
        const std::vector<idx>& subsys,
        idx d = 2)
{
    const dyn_col_vect<typename Derived::Scalar>& rphi = phi.derived();
    const dyn_col_vect<typename Derived::Scalar>& rpsi = psi.derived();

    // EXCEPTION CHECKS

    if ( !internal::_check_nonzero_size(rpsi))
        throw Exception("qpp::ip()", Exception::Type::ZERO_SIZE);

    // check valid dims
    if ( d == 0 )
        throw Exception("qpp::ip()", Exception::Type::DIMS_INVALID);
    // END EXCEPTION CHECKS

    idx n =
            static_cast<idx>(std::llround(std::log2(rpsi.rows()) /
                                          std::log2(d)));
    std::vector<idx> dims(n, d); // local dimensions vector
    return ip(phi, psi, subsys, dims);
}

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
std::tuple<idx, std::vector<double>, std::vector<cmat>>
measure(const Eigen::MatrixBase<Derived>& A, const std::vector<cmat>& Ks)
{
    const dyn_mat<typename Derived::Scalar>& rA = A.derived();

    // EXCEPTION CHECKS

    // check zero-size
    if ( !internal::_check_nonzero_size(rA))
        throw Exception("qpp::measure()", Exception::Type::ZERO_SIZE);

    // check the Kraus operators
    if ( Ks.size() == 0 )
        throw Exception("qpp::measure()", Exception::Type::ZERO_SIZE);
    if ( !internal::_check_square_mat(Ks[0]))
        throw Exception("qpp::measure()", Exception::Type::MATRIX_NOT_SQUARE);
    if ( Ks[0].rows() != rA.rows())
        throw Exception("qpp::measure()",
                        Exception::Type::DIMS_MISMATCH_MATRIX);
    for ( auto&& it : Ks )
        if ( it.rows() != Ks[0].rows() || it.cols() != Ks[0].rows())
            throw Exception("qpp::measure()", Exception::Type::DIMS_NOT_EQUAL);
    // END EXCEPTION CHECKS

    // probabilities
    std::vector<double> prob(Ks.size());
    // resulting states
    std::vector<cmat> outstates(Ks.size());

    //************ density matrix ************//
    if ( internal::_check_square_mat(rA)) // square matrix
    {
        for ( idx i = 0; i < Ks.size(); ++i )
        {
            outstates[i] = cmat::Zero(rA.rows(), rA.rows());
            cmat tmp = Ks[i] * rA * adjoint(Ks[i]); // un-normalized;
            prob[i] = std::abs(trace(tmp)); // probability
            if ( prob[i] > eps )
                outstates[i] = tmp / prob[i]; // normalized
        }
    }
        //************ ket ************//
    else if ( internal::_check_cvector(rA)) // column vector
    {
        for ( idx i = 0; i < Ks.size(); ++i )
        {
            outstates[i] = ket::Zero(rA.rows());
            ket tmp = Ks[i] * rA; // un-normalized;
            // probability
            prob[i] = std::pow(norm(tmp), 2);
            if ( prob[i] > eps )
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
std::tuple<idx, std::vector<double>, std::vector<cmat>>
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
std::tuple<idx, std::vector<double>, std::vector<cmat>>
measure(const Eigen::MatrixBase<Derived>& A, const cmat& U)
{
    const dyn_mat<typename Derived::Scalar>& rA = A.derived();

    // EXCEPTION CHECKS

    // check zero-size
    if ( !internal::_check_nonzero_size(rA))
        throw Exception("qpp::measure()", Exception::Type::ZERO_SIZE);

    // check the unitary basis matrix U
    if ( !internal::_check_nonzero_size(U))
        throw Exception("qpp::measure()", Exception::Type::ZERO_SIZE);
    if ( !internal::_check_square_mat(U))
        throw Exception("qpp::measure()", Exception::Type::MATRIX_NOT_SQUARE);
    if ( U.rows() != rA.rows())
        throw Exception("qpp::measure()",
                        Exception::Type::DIMS_MISMATCH_MATRIX);
    // END EXCEPTION CHECKS

    std::vector<cmat> Ks(U.rows());
    for ( idx i = 0; i < static_cast<idx>(U.rows()); ++i )
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
std::tuple<idx, std::vector<double>, std::vector<cmat>>
measure(const Eigen::MatrixBase<Derived>& A,
        const std::vector<cmat>& Ks,
        const std::vector<idx>& subsys,
        const std::vector<idx>& dims)
{
    const typename Eigen::MatrixBase<Derived>::EvalReturnType& rA
            = A.derived();

    // EXCEPTION CHECKS

    // check zero-size
    if ( !internal::_check_nonzero_size(rA))
        throw Exception("qpp::measure()", Exception::Type::ZERO_SIZE);

    // check that dimension is valid
    if ( !internal::_check_dims(dims))
        throw Exception("qpp::measure()", Exception::Type::DIMS_INVALID);

    // check that dims match rho matrix
    if ( !internal::_check_dims_match_mat(dims, rA))
        throw Exception("qpp::measure()",
                        Exception::Type::DIMS_MISMATCH_MATRIX);

    // check subsys is valid w.r.t. dims
    if ( !internal::_check_subsys_match_dims(subsys, dims))
        throw Exception("qpp::measure()",
                        Exception::Type::SUBSYS_MISMATCH_DIMS);

    std::vector<idx> subsys_dims(subsys.size());
    for ( idx i = 0; i < subsys.size(); ++i )
        subsys_dims[i] = dims[subsys[i]];

    idx D = prod(std::begin(dims), std::end(dims));
    idx Dsubsys = prod(std::begin(subsys_dims), std::end(subsys_dims));
    idx Dsubsysbar = D / Dsubsys;

    // check the Kraus operators
    if ( Ks.size() == 0 )
        throw Exception("qpp::measure()",
                        Exception::Type::ZERO_SIZE);
    if ( !internal::_check_square_mat(Ks[0]))
        throw Exception("qpp::measure()",
                        Exception::Type::MATRIX_NOT_SQUARE);
    if ( Dsubsys != static_cast<idx>(Ks[0].rows()))
        throw Exception("qpp::measure()",
                        Exception::Type::DIMS_MISMATCH_MATRIX);
    for ( auto&& it : Ks )
        if ( it.rows() != Ks[0].rows() || it.cols() != Ks[0].rows())
            throw Exception("qpp::measure()",
                            Exception::Type::DIMS_NOT_EQUAL);
    // END EXCEPTION CHECKS

    // probabilities
    std::vector<double> prob(Ks.size());
    // resulting states
    std::vector<cmat> outstates(Ks.size(), cmat::Zero(Dsubsysbar, Dsubsysbar));

    //************ density matrix ************//
    if ( internal::_check_square_mat(rA)) // square matrix
    {
        for ( idx i = 0; i < Ks.size(); ++i )
        {
            cmat tmp = apply(rA, Ks[i], subsys, dims);
            tmp = ptrace(tmp, subsys, dims);
            prob[i] = std::abs(trace(tmp)); // probability
            if ( prob[i] > eps )
            {
                // normalized output state
                // corresponding to measurement result i
                outstates[i] = tmp / prob[i];
            }
        }
    }
        //************ ket ************//
    else if ( internal::_check_cvector(rA)) // column vector
    {
        for ( idx i = 0; i < Ks.size(); ++i )
        {
            ket tmp = apply(rA, Ks[i], subsys, dims);
            prob[i] = std::pow(norm(tmp), 2);
            if ( prob[i] > eps )
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
* \param Ks Set of Kraus operators
* \param subsys Subsystem indexes that are measured
* \param dims Dimensions of the multi-partite system
* \return Tuple of: 1. Result of the measurement, 2.
* Vector of outcome probabilities, and 3. Vector of post-measurement
* normalized states
*/
template<typename Derived>
std::tuple<idx, std::vector<double>, std::vector<cmat>>
measure(const Eigen::MatrixBase<Derived>& A,
        const std::initializer_list<cmat>& Ks,
        const std::vector<idx>& subsys,
        const std::vector<idx>& dims)
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
* \param Ks Set of Kraus operators
* \param subsys Subsystem indexes that are measured
* \param d Subsystem dimensions
* \return Tuple of: 1. Result of the measurement, 2.
* Vector of outcome probabilities, and 3. Vector of post-measurement
* normalized states
*/
template<typename Derived>
std::tuple<idx, std::vector<double>, std::vector<cmat>>
measure(const Eigen::MatrixBase<Derived>& A,
        const std::vector<cmat>& Ks,
        const std::vector<idx>& subsys,
        idx d = 2)
{
    const typename Eigen::MatrixBase<Derived>::EvalReturnType& rA
            = A.derived();

    // EXCEPTION CHECKS

    // check zero size
    if ( !internal::_check_nonzero_size(rA))
        throw Exception("qpp::measure()", Exception::Type::ZERO_SIZE);

    // check valid dims
    if ( d == 0 )
        throw Exception("qpp::measure()", Exception::Type::DIMS_INVALID);
    // END EXCEPTION CHECKS

    idx n =
            static_cast<idx>(std::llround(std::log2(rA.rows()) /
                                          std::log2(d)));
    std::vector<idx> dims(n, d); // local dimensions vector

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
* \param Ks Set of Kraus operators
* \param subsys Subsystem indexes that are measured
* \param d Subsystem dimensions
* \return Tuple of: 1. Result of the measurement, 2.
* Vector of outcome probabilities, and 3. Vector of post-measurement
* normalized states
*/
template<typename Derived>
std::tuple<idx, std::vector<double>, std::vector<cmat>>
measure(const Eigen::MatrixBase<Derived>& A,
        const std::initializer_list<cmat>& Ks,
        const std::vector<idx>& subsys,
        idx d = 2)
{
    return measure(A, std::vector<cmat>(Ks), subsys, d);
}

/**
* \brief Measures the part \a subsys of
* the multi-partite state vector or density matrix \a A
* in the orthonormal basis or rank-1 POVM specified by the matrix \a V
* \see qpp::measure_seq()
*
* \note The dimension of \a V must match the dimension of \a subsys.
* The measurement is destructive, i.e. the measured subsystems are traced away.
*
* \param A Eigen expression
* \param V Matrix whose columns represent the measurement basis vectors or the
* bra parts of the rank-1 POVM
* \param subsys Subsystem indexes that are measured
* \param dims Dimensions of the multi-partite system
* \return Tuple of: 1. Result of the measurement, 2.
* Vector of outcome probabilities, and 3. Vector of post-measurement
* normalized states
*/
template<typename Derived>
std::tuple<idx, std::vector<double>, std::vector<cmat>>

measure(const Eigen::MatrixBase<Derived>& A,
        const cmat& V,
        const std::vector<idx>& subsys,
        const std::vector<idx>& dims)
{
    const typename Eigen::MatrixBase<Derived>::EvalReturnType& rA
            = A.derived();

    // EXCEPTION CHECKS

    // check zero-size
    if ( !internal::_check_nonzero_size(rA))
        throw Exception("qpp::measure()", Exception::Type::ZERO_SIZE);

    // check that dimension is valid
    if ( !internal::_check_dims(dims))
        throw Exception("qpp::measure()", Exception::Type::DIMS_INVALID);

    // check subsys is valid w.r.t. dims
    if ( !internal::_check_subsys_match_dims(subsys, dims))
        throw Exception("qpp::measure()",
                        Exception::Type::SUBSYS_MISMATCH_DIMS);

    std::vector<idx> subsys_dims(subsys.size());
    for ( idx i = 0; i < subsys.size(); ++i )
        subsys_dims[i] = dims[subsys[i]];

    idx Dsubsys = prod(std::begin(subsys_dims), std::end(subsys_dims));

    // check the matrix V
    if ( !internal::_check_nonzero_size(V))
        throw Exception("qpp::measure()", Exception::Type::ZERO_SIZE);
//    if (!internal::_check_square_mat(U))
//        throw Exception("qpp::measure()", Exception::Type::MATRIX_NOT_SQUARE);
    if ( Dsubsys != static_cast<idx>(V.rows()))
        throw Exception("qpp::measure()",
                        Exception::Type::DIMS_MISMATCH_MATRIX);
    // END EXCEPTION CHECKS

    // number of basis (rank-1 POVM) elements
    idx M = static_cast<idx>(V.cols());

    //************ ket ************//
    if ( internal::_check_cvector(rA))
    {
        const ket& rpsi = A.derived();
        // check that dims match state vector
        if ( !internal::_check_dims_match_cvect(dims, rA))
            throw Exception("qpp::measure()",
                            Exception::Type::DIMS_MISMATCH_CVECTOR);

        std::vector<double> prob(M); // probabilities
        std::vector<cmat> outstates(M); // resulting states

#ifdef _WITH_OPENMP_
#pragma omp parallel for
#endif
        for ( idx i = 0; i < M; ++i )
            outstates[i] = ip(static_cast<const ket&>(V.col(i)),
                              rpsi, subsys, dims);

        for ( idx i = 0; i < M; ++i )
        {
            double tmp = norm(outstates[i]);
            prob[i] = tmp * tmp;
            if ( prob[i] > eps )
            {
                // normalized output state
                // corresponding to measurement result m
                outstates[i] /= tmp;
            }
        }

        // sample from the probability distribution
        std::discrete_distribution<idx> dd(std::begin(prob),
                                           std::end(prob));
        idx result = dd(RandomDevices::get_instance()._rng);

        return std::make_tuple(result, prob, outstates);
    }
        //************ density matrix ************//
    else if ( internal::_check_square_mat(rA))
    {
        // check that dims match rho matrix
        if ( !internal::_check_dims_match_mat(dims, rA))
            throw Exception("qpp::measure()",
                            Exception::Type::DIMS_MISMATCH_MATRIX);

        std::vector<cmat> Ks(M);
        for ( idx i = 0; i < M; ++i )
            Ks[i] = V.col(i) * adjoint(V.col(i));

        return measure(rA, Ks, subsys, dims);
    }
    //************ Exception: not ket nor density matrix ************//
    throw Exception("qpp::measure()",
                    Exception::Type::MATRIX_NOT_SQUARE_OR_CVECTOR);
}

/**
* \brief Measures the part \a subsys of
* the multi-partite state vector or density matrix \a A
* in the orthonormal basis or rank-1 POVM specified by the matrix \a V
* \see qpp::measure_seq()
*
* \note The dimension of \a V must match the dimension of \a subsys.
* The measurement is destructive, i.e. the measured subsystems are traced away.
*
* \param A Eigen expression
* \param V Matrix whose columns represent the measurement basis vectors or the
* bra parts of the rank-1 POVM
* \param subsys Subsystem indexes that are measured
* \param d Subsystem dimensions
* \return Tuple of: 1. Result of the measurement, 2.
* Vector of outcome probabilities, and 3. Vector of post-measurement
* normalized states
*/
template<typename Derived>
std::tuple<idx, std::vector<double>, std::vector<cmat>>

measure(const Eigen::MatrixBase<Derived>& A,
        const cmat& V,
        const std::vector<idx>& subsys,
        idx d = 2)
{
    const typename Eigen::MatrixBase<Derived>::EvalReturnType& rA
            = A.derived();

    // EXCEPTION CHECKS

    // check zero size
    if ( !internal::_check_nonzero_size(rA))
        throw Exception("qpp::measure()", Exception::Type::ZERO_SIZE);

    // check valid dims
    if ( d == 0 )
        throw Exception("qpp::measure()", Exception::Type::DIMS_INVALID);
    // END EXCEPTION CHECKS

    idx n =
            static_cast<idx>(std::llround(std::log2(rA.rows()) /
                                          std::log2(d)));
    std::vector<idx> dims(n, d); // local dimensions vector

    return measure(rA, V, subsys, dims);
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
std::tuple<std::vector<idx>, double, cmat>
measure_seq(const Eigen::MatrixBase<Derived>& A,
            std::vector<idx> subsys,
            std::vector<idx> dims)
{
    dyn_mat<typename Derived::Scalar> cA = A.derived();

    // EXCEPTION CHECKS

    // check zero-size
    if ( !internal::_check_nonzero_size(cA))
        throw Exception("qpp::measure_seq()", Exception::Type::ZERO_SIZE);

    // check that dimension is valid
    if ( !internal::_check_dims(dims))
        throw Exception("qpp::measure_seq()", Exception::Type::DIMS_INVALID);

    // check that dims match rho matrix
    if ( !internal::_check_dims_match_mat(dims, cA))
        throw Exception("qpp::measure_seq()",
                        Exception::Type::DIMS_MISMATCH_MATRIX);

    // check subsys is valid w.r.t. dims
    if ( !internal::_check_subsys_match_dims(subsys, dims))
        throw Exception("qpp::measure_seq()",
                        Exception::Type::SUBSYS_MISMATCH_DIMS);
    // END EXCEPTION CHECKS

    std::vector<idx> result;
    double prob = 1;

    // sort subsys in decreasing order,
    // the order of measurements does not matter
    std::sort(std::begin(subsys), std::end(subsys), std::greater<idx>{});

    //************ density matrix or column vector ************//
    if ( internal::_check_square_mat(cA) || internal::_check_cvector(cA))
    {
        while (subsys.size() > 0)
        {
            auto tmp = measure(
                    cA, Gates::get_instance().Id(dims[subsys[0]]),
                    {subsys[0]}, dims
            );
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
std::tuple<std::vector<idx>, double, cmat>
measure_seq(const Eigen::MatrixBase<Derived>& A,
            std::vector<idx> subsys,
            idx d = 2)
{
    const typename Eigen::MatrixBase<Derived>::EvalReturnType& rA
            = A.derived();

    // EXCEPTION CHECKS

    // check zero size
    if ( !internal::_check_nonzero_size(rA))
        throw Exception("qpp::measure_seq()", Exception::Type::ZERO_SIZE);

    // check valid dims
    if ( d == 0 )
        throw Exception("qpp::measure_seq()", Exception::Type::DIMS_INVALID);
    // END EXCEPTION CHECKS

    idx n =
            static_cast<idx>(std::llround(std::log2(rA.rows()) /
                                          std::log2(d)));
    std::vector<idx> dims(n, d); // local dimensions vector

    return measure_seq(rA, subsys, dims);
}

} /* namespace qpp */

#endif /* INSTRUMENTS_H_ */
