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
* \file experimental/experimental.h
* \brief Experimental/test functions/classes
*/

#ifndef EXPERIMENTAL_EXPERIMENTAL_H_
#define EXPERIMENTAL_EXPERIMENTAL_H_

namespace qpp
{
/**
* \namespace qpp::experimental
* \brief Experimental/test functions/classes, do not use or modify
*/
namespace experimental
{

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
* \param subsys Subsystem indexes that are measured
* \param dims Dimensions of the multi-partite system
* \param V Matrix whose columns represent the measurement basis vectors or the
* bra parts of the rank-1 POVM
* \return Tuple of: 1. Result of the measurement, 2.
* Vector of outcome probabilities, and 3. Vector of post-measurement
* normalized states
*/
template<typename Derived>
std::tuple<idx, std::vector<double>, std::vector<cmat>>

_measure(const Eigen::MatrixBase<Derived>& A,
         const cmat& V,
         const std::vector<idx>& subsys,
         const std::vector<idx>& dims)
{
    const cmat& rA = A;

    // EXCEPTION CHECKS

    // check zero-size
    if (!internal::_check_nonzero_size(rA))
        throw Exception("qpp::_measure()", Exception::Type::ZERO_SIZE);

    // check that dimension is valid
    if (!internal::_check_dims(dims))
        throw Exception("qpp::_measure()", Exception::Type::DIMS_INVALID);

    // check subsys is valid w.r.t. dims
    if (!internal::_check_subsys_match_dims(subsys, dims))
        throw Exception("qpp::_measure()",
                        Exception::Type::SUBSYS_MISMATCH_DIMS);

    std::vector<idx> subsys_dims(subsys.size());
    for (idx i = 0; i < subsys.size(); ++i)
        subsys_dims[i] = dims[subsys[i]];

    idx Dsubsys = prod(std::begin(subsys_dims), std::end(subsys_dims));

    // check the matrix V
    if (!internal::_check_nonzero_size(V))
        throw Exception("qpp::_measure()", Exception::Type::ZERO_SIZE);
//    if (!internal::_check_square_mat(U))
//        throw Exception("qpp::measure()", Exception::Type::MATRIX_NOT_SQUARE);
    if (Dsubsys != static_cast<idx>(V.rows()))
        throw Exception("qpp::_measure()",
                        Exception::Type::DIMS_MISMATCH_MATRIX);
    // END EXCEPTION CHECKS

    // number of basis (rank-1 POVM) elements
    idx M = static_cast<idx>(V.cols());

    // TODO: modify this!!! (and document U in case of rank-one POVMs)
    //************ ket ************//
    if (internal::_check_cvector(rA))
    {
        PRINTLN("Ket");
        // check that dims match state vector
        if (!internal::_check_dims_match_cvect(dims, rA))
            throw Exception("qpp::_measure()",
                            Exception::Type::DIMS_MISMATCH_CVECTOR);

        idx D = static_cast<idx>(rA.rows());
        idx Dbar = D / Dsubsys;

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

        for (idx i = 0; i < n; ++i)
        {
            Cdims[i] = dims[i];
        }
        for (idx i = 0; i < nsubsys; ++i)
        {
            Csubsys[i] = subsys[i];
            Cdimssubsys[i] = dims[subsys[i]];
        }
        for (idx i = 0; i < nsubsysbar; ++i)
        {
            Cdimssubsysbar[i] = dims[subsys_bar[i]];
        }

        std::vector<double> prob(M); // probabilities
        std::vector<cmat> outstates(M, ket::Zero(Dbar)); // resulting states

        std::vector<ket> phi(M); // basis (rank-1 POVM) elements
        for (idx i = 0; i < M; ++i)
            phi[i] = V.col(i);

        auto worker = [=](idx m, idx b) noexcept
                -> typename Derived::Scalar
        {
            idx Cmidxrow[maxn];
            idx Cmidxrowsubsys[maxn];
            idx Cmidxcolsubsysbar[maxn];

            /* get the col multi-indexes of the complement */
            internal::_n2multiidx(b, nsubsysbar,
                                  Cdimssubsysbar, Cmidxcolsubsysbar);
            /* write it in the global row multi-index */
            for (idx k = 0; k < nsubsysbar; ++k)
            {
                Cmidxrow[Csubsysbar[k]] = Cmidxcolsubsysbar[k];
            }

            typename Derived::Scalar result{};

            for (idx a = 0; a < Dsubsys; ++a)
            {
                /* get the row multi-indexes of the subsys */
                internal::_n2multiidx(a, nsubsys,
                                      Cdimssubsys, Cmidxrowsubsys);
                /* write it in the global row multi-index */
                for (idx k = 0; k < nsubsys; ++k)
                {
                    Cmidxrow[Csubsys[k]] = Cmidxrowsubsys[k];
                }
                // compute the row index
                idx i = internal::_multiidx2n(Cmidxrow, n, Cdims);

                // TODO: check whether need adjoint(V.col(m)(a))
                result += rA(i) * std::conj(V.col(m)(a));
            }

            return result;
        };

#pragma omp parallel for collapse(2)
        for (idx i = 0; i < M; ++i)
            for (idx b = 0; b < Dbar; ++b)
                outstates[i](b) = worker(i, b);

        for (idx i = 0; i < M; ++i)
        {
            double tmp = norm(outstates[i]);
            prob[i] = tmp * tmp;
            if (prob[i] > eps)
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
    else if (internal::_check_square_mat(rA))
    {
        PRINTLN("Matrix");
        // check that dims match rho matrix
        if (!internal::_check_dims_match_mat(dims, rA))
            throw Exception("qpp::_measure()",
                            Exception::Type::DIMS_MISMATCH_MATRIX);

        std::vector<cmat> Ks(M);
        for (idx i = 0; i < M; ++i)
            Ks[i] = V.col(i) * adjoint(V.col(i));

        return measure(rA, Ks, subsys, dims);
    }
    //************ Exception: not ket nor density matrix ************//
    throw Exception("qpp::_measure()",
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
* \param subsys Subsystem indexes that are measured
* \param d Subsystem dimensions
* \param V Matrix whose columns represent the measurement basis vectors or the
* bra parts of the rank-1 POVM
* \return Tuple of: 1. Result of the measurement, 2.
* Vector of outcome probabilities, and 3. Vector of post-measurement
* normalized states
*/
template<typename Derived>
std::tuple<idx, std::vector<double>, std::vector<cmat>>

_measure(const Eigen::MatrixBase<Derived>& A,
         const cmat& V,
         const std::vector<idx>& subsys,
         const idx d = 2)
{
    const cmat& rA = A;

    // check zero size
    if (!internal::_check_nonzero_size(rA))
        throw Exception("qpp::_measure()", Exception::Type::ZERO_SIZE);

    idx n =
            static_cast<idx>(std::llround(std::log2(rA.rows()) /
                                          std::log2(d)));
    std::vector<idx> dims(n, d); // local dimensions vector

    return _measure(rA, V, subsys, dims);
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

_measure_seq(const Eigen::MatrixBase<Derived>& A,
             std::vector<idx> subsys,
             std::vector<idx> dims)
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

    std::vector<idx> result;
    double prob = 1;

    // sort subsys in decreasing order,
    // the order of measurements does not matter
    std::sort(std::begin(subsys), std::end(subsys), std::greater<idx>{});

    // TODO: modify this so it treats both density matrices and column vectors
    //************ density matrix or column vector ************//
    if (internal::_check_square_mat(cA) || internal::_check_cvector(cA))
    {
        while (subsys.size() > 0)
        {
            auto tmp = experimental::_measure(
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
        throw Exception("qpp::_measure_seq()",
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

_measure_seq(const Eigen::MatrixBase<Derived>& A,
             std::vector<idx> subsys, idx d = 2)
{
    const cmat& rA = A;

    // check zero size
    if (!internal::_check_nonzero_size(rA))
        throw Exception("qpp::_measure_seq()", Exception::Type::ZERO_SIZE);

    idx n =
            static_cast<idx>(std::llround(std::log2(rA.rows()) /
                                          std::log2(d)));
    std::vector<idx> dims(n, d); // local dimensions vector

    return experimental::_measure_seq(rA, subsys, dims);
}

} /* namespace experimental */
} /* namespace qpp */

#endif /* EXPERIMENTAL_EXPERIMENTAL_H_ */
