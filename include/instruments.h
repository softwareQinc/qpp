/*
 * This file is part of Quantum++.
 *
 * MIT License
 *
 * Copyright (c) 2013 - 2019 Vlad Gheorghiu (vgheorgh@gmail.com)
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

/**
 * \file instruments.h
 * \brief Measurement functions
 */

#ifndef INSTRUMENTS_H_
#define INSTRUMENTS_H_

namespace qpp {
/**
 * \brief Generalized inner product
 *
 * \param phi Column vector Eigen expression
 * \param psi Column vector Eigen expression
 * \param subsys Subsystem indexes over which \a phi is defined
 * \param dims Dimensions of the multi-partite system
 * \return Inner product \f$\langle \phi_{subsys}|\psi\rangle\f$, as a scalar or
 * column vector over the remaining Hilbert space
 */
template <typename Derived>
dyn_col_vect<typename Derived::Scalar>
ip(const Eigen::MatrixBase<Derived>& phi, const Eigen::MatrixBase<Derived>& psi,
   const std::vector<idx>& subsys, const std::vector<idx>& dims) {
    const dyn_col_vect<typename Derived::Scalar>& rphi = phi.derived();
    const dyn_col_vect<typename Derived::Scalar>& rpsi = psi.derived();

    // EXCEPTION CHECKS

    // check zero-size
    if (!internal::check_nonzero_size(rphi))
        throw exception::ZeroSize("qpp::ip()");

    // check zero-size
    if (!internal::check_nonzero_size(rpsi))
        throw exception::ZeroSize("qpp::ip()");

    // check column vector
    if (!internal::check_cvector(rphi))
        throw exception::MatrixNotCvector("qpp::ip()");

    // check column vector
    if (!internal::check_cvector(rpsi))
        throw exception::MatrixNotCvector("qpp::ip()");

    // check that dims is a valid dimension vector
    if (!internal::check_dims(dims))
        throw exception::DimsInvalid("qpp::ip()");

    // check that subsys are valid
    if (!internal::check_subsys_match_dims(subsys, dims))
        throw exception::SubsysMismatchDims("qpp::ip()");

    // check that dims match psi column vector
    if (!internal::check_dims_match_cvect(dims, rpsi))
        throw exception::DimsMismatchCvector("qpp::ip()");

    // check that subsys match pho column vector
    std::vector<idx> subsys_dims(subsys.size());
    for (idx i = 0; i < subsys.size(); ++i)
        subsys_dims[i] = dims[subsys[i]];
    if (!internal::check_dims_match_cvect(subsys_dims, rphi))
        throw exception::DimsMismatchCvector("qpp::ip()");
    // END EXCEPTION CHECKS

    idx Dsubsys = prod(std::begin(subsys_dims), std::end(subsys_dims));

    idx D = static_cast<idx>(rpsi.rows());
    idx Dsubsys_bar = D / Dsubsys;

    idx n = dims.size();
    idx n_subsys = subsys.size();
    idx n_subsys_bar = n - n_subsys;

    idx Cdims[maxn];
    idx Csubsys[maxn];
    idx Cdimssubsys[maxn];
    idx Csubsys_bar[maxn];
    idx Cdimssubsys_bar[maxn];

    std::vector<idx> subsys_bar = complement(subsys, n);
    std::copy(std::begin(subsys_bar), std::end(subsys_bar),
              std::begin(Csubsys_bar));

    for (idx i = 0; i < n; ++i) {
        Cdims[i] = dims[i];
    }
    for (idx i = 0; i < n_subsys; ++i) {
        Csubsys[i] = subsys[i];
        Cdimssubsys[i] = dims[subsys[i]];
    }
    for (idx i = 0; i < n_subsys_bar; ++i) {
        Cdimssubsys_bar[i] = dims[subsys_bar[i]];
    }

    auto worker = [&](idx b) noexcept -> typename Derived::Scalar {
        idx Cmidxrow[maxn];
        idx Cmidxrowsubsys[maxn];
        idx Cmidxcolsubsys_bar[maxn];

        /* get the col multi-indexes of the complement */
        internal::n2multiidx(b, n_subsys_bar, Cdimssubsys_bar,
                             Cmidxcolsubsys_bar);
        /* write it in the global row multi-index */
        for (idx k = 0; k < n_subsys_bar; ++k) {
            Cmidxrow[Csubsys_bar[k]] = Cmidxcolsubsys_bar[k];
        }

        typename Derived::Scalar result = 0;
        for (idx a = 0; a < Dsubsys; ++a) {
            /* get the row multi-indexes of the subsys */
            internal::n2multiidx(a, n_subsys, Cdimssubsys, Cmidxrowsubsys);
            /* write it in the global row multi-index */
            for (idx k = 0; k < n_subsys; ++k) {
                Cmidxrow[Csubsys[k]] = Cmidxrowsubsys[k];
            }
            // compute the row index
            idx i = internal::multiidx2n(Cmidxrow, n, Cdims);

            result += std::conj(rphi(a)) * rpsi(i);
        }

        return result;
    }; /* end worker */

    dyn_col_vect<typename Derived::Scalar> result(Dsubsys_bar);
#ifdef WITH_OPENMP_
#pragma omp parallel for
#endif // WITH_OPENMP_
    for (idx m = 0; m < Dsubsys_bar; ++m)
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
 * \return Inner product \f$\langle \phi_{subsys}|\psi\rangle\f$, as a scalar or
 * column vector over the remaining Hilbert space
 */
template <typename Derived>
dyn_col_vect<typename Derived::Scalar>
ip(const Eigen::MatrixBase<Derived>& phi, const Eigen::MatrixBase<Derived>& psi,
   const std::vector<idx>& subsys, idx d = 2) {
    const dyn_col_vect<typename Derived::Scalar>& rphi = phi.derived();
    const dyn_col_vect<typename Derived::Scalar>& rpsi = psi.derived();

    // EXCEPTION CHECKS

    if (!internal::check_nonzero_size(rpsi))
        throw exception::ZeroSize("qpp::ip()");

    // check valid dims
    if (d < 2)
        throw exception::DimsInvalid("qpp::ip()");
    // END EXCEPTION CHECKS

    idx n = internal::get_num_subsys(static_cast<idx>(rpsi.rows()), d);
    std::vector<idx> dims(n, d); // local dimensions vector
    return ip(phi, psi, subsys, dims);
}

// full measurements
/**
 * \brief Measures the state vector or density operator \a A using the set of
 * Kraus operators \a Ks
 *
 * \param A Eigen expression
 * \param Ks Set of Kraus operators
 * \return Tuple of: 1. Result of the measurement, 2. Vector of outcome
 * probabilities, and 3. Vector of post-measurement normalized states
 */
template <typename Derived>
std::tuple<idx, std::vector<double>, std::vector<cmat>>
measure(const Eigen::MatrixBase<Derived>& A, const std::vector<cmat>& Ks) {
    const dyn_mat<typename Derived::Scalar>& rA = A.derived();

    // EXCEPTION CHECKS

    // check zero-size
    if (!internal::check_nonzero_size(rA))
        throw exception::ZeroSize("qpp::measure()");

    // check the Kraus operators
    if (Ks.empty())
        throw exception::ZeroSize("qpp::measure()");
    if (!internal::check_square_mat(Ks[0]))
        throw exception::MatrixNotSquare("qpp::measure()");
    if (Ks[0].rows() != rA.rows())
        throw exception::DimsMismatchMatrix("qpp::measure()");
    for (auto&& elem : Ks)
        if (elem.rows() != Ks[0].rows() || elem.cols() != Ks[0].rows())
            throw exception::DimsNotEqual("qpp::measure()");
    // END EXCEPTION CHECKS

    // probabilities
    std::vector<double> prob(Ks.size());
    // resulting states
    std::vector<cmat> outstates(Ks.size());

    //************ density matrix ************//
    if (internal::check_square_mat(rA)) // square matrix
    {
        for (idx i = 0; i < Ks.size(); ++i) {
            outstates[i] = cmat::Zero(rA.rows(), rA.rows());
            cmat tmp = Ks[i] * rA * adjoint(Ks[i]); // un-normalized;
            prob[i] = std::abs(trace(tmp));         // probability
            if (prob[i] > 0)
                outstates[i] = tmp / prob[i]; // normalized
        }
    }
    //************ ket ************//
    else if (internal::check_cvector(rA)) // column vector
    {
        for (idx i = 0; i < Ks.size(); ++i) {
            outstates[i] = ket::Zero(rA.rows());
            ket tmp = Ks[i] * rA; // un-normalized;
            // probability
            prob[i] = std::pow(norm(tmp), 2);
            if (prob[i] > 0)
                outstates[i] = tmp / std::sqrt(prob[i]); // normalized
        }
    } else
        throw exception::MatrixNotSquareNorCvector("qpp::measure()");

    // sample from the probability distribution
    std::discrete_distribution<idx> dd(std::begin(prob), std::end(prob));
    auto& gen =
#ifdef NO_THREAD_LOCAL_
        RandomDevices::get_instance().get_prng();
#else
        RandomDevices::get_thread_local_instance().get_prng();
#endif
    idx result = dd(gen);

    return std::make_tuple(result, prob, outstates);
}

// std::initializer_list overload, avoids ambiguity for 2-element lists, see
// http://stackoverflow.com
// /questions/26750039/ambiguity-when-using-initializer-list-as-parameter
/**
 * \brief Measures the state vector or density matrix \a A using the set of
 * Kraus operators \a Ks
 *
 * \param A Eigen expression
 * \param Ks Set of Kraus operators
 * \return Tuple of: 1. Result of the measurement, 2. Vector of outcome
 * probabilities, and 3. Vector of post-measurement normalized states
 */
template <typename Derived>
std::tuple<idx, std::vector<double>, std::vector<cmat>>
measure(const Eigen::MatrixBase<Derived>& A,
        const std::initializer_list<cmat>& Ks) {
    return measure(A, std::vector<cmat>(Ks));
}

/**
 * \brief Measures the state vector or density matrix \a A in the orthonormal
 * basis specified by the unitary matrix \a U
 *
 * \param A Eigen expression
 * \param U Unitary matrix whose columns represent the measurement basis vectors
 * \return Tuple of: 1. Result of the measurement, 2. Vector of outcome
 * probabilities, and 3. Vector of post-measurement normalized states
 */
template <typename Derived>
std::tuple<idx, std::vector<double>, std::vector<cmat>>
measure(const Eigen::MatrixBase<Derived>& A, const cmat& U) {
    const dyn_mat<typename Derived::Scalar>& rA = A.derived();

    // EXCEPTION CHECKS

    // check zero-size
    if (!internal::check_nonzero_size(rA))
        throw exception::ZeroSize("qpp::measure()");

    // check the unitary basis matrix U
    if (!internal::check_nonzero_size(U))
        throw exception::ZeroSize("qpp::measure()");
    if (!internal::check_square_mat(U))
        throw exception::MatrixNotSquare("qpp::measure()");
    if (U.rows() != rA.rows())
        throw exception::DimsMismatchMatrix("qpp::measure()");
    // END EXCEPTION CHECKS

    std::vector<cmat> Ks(U.rows());
    for (idx i = 0; i < static_cast<idx>(U.rows()); ++i)
        Ks[i] = U.col(i) * adjoint(U.col(i));

    return measure(rA, Ks);
}

// partial measurements
/**
 * \brief  Measures the part \a subsys of the multi-partite state vector or
 * density matrix \a A using the set of Kraus operators \a Ks
 * \see qpp::measure_seq()
 *
 * \note The dimension of all \a Ks must match the dimension of \a target.
 * If \a destructive is set to true (by default), the measurement is
 * destructive, i.e. the measured subsystems are traced away.
 *
 * \param A Eigen expression
 * \param Ks Set of Kraus operators
 * \param target Subsystem indexes that are measured
 * \param dims Dimensions of the multi-partite system
 * \param destructive Destructive measurement, true by default
 * \return Tuple of: 1. Result of the measurement, 2. Vector of outcome
 * probabilities, and 3. Vector of post-measurement normalized states
 */
template <typename Derived>
std::tuple<idx, std::vector<double>, std::vector<cmat>>
measure(const Eigen::MatrixBase<Derived>& A, const std::vector<cmat>& Ks,
        const std::vector<idx>& target, const std::vector<idx>& dims,
        bool destructive = true) {
    const typename Eigen::MatrixBase<Derived>::EvalReturnType& rA = A.derived();

    // EXCEPTION CHECKS

    // check zero-size
    if (!internal::check_nonzero_size(rA))
        throw exception::ZeroSize("qpp::measure()");

    // check that dimension is valid
    if (!internal::check_dims(dims))
        throw exception::DimsInvalid("qpp::measure()");

    // check that target is valid w.r.t. dims
    if (!internal::check_subsys_match_dims(target, dims))
        throw exception::SubsysMismatchDims("qpp::measure()");

    std::vector<idx> subsys_dims(target.size());
    for (idx i = 0; i < target.size(); ++i)
        subsys_dims[i] = dims[target[i]];

    idx D = prod(std::begin(dims), std::end(dims));
    idx Dsubsys = prod(std::begin(subsys_dims), std::end(subsys_dims));
    idx Dsubsys_bar = D / Dsubsys;

    // check the Kraus operators
    if (Ks.empty())
        throw exception::ZeroSize("qpp::measure()");
    if (!internal::check_square_mat(Ks[0]))
        throw exception::MatrixNotSquare("qpp::measure()");
    if (Dsubsys != static_cast<idx>(Ks[0].rows()))
        throw exception::DimsMismatchMatrix("qpp::measure()");
    for (auto&& elem : Ks)
        if (elem.rows() != Ks[0].rows() || elem.cols() != Ks[0].rows())
            throw exception::DimsNotEqual("qpp::measure()");
    // END EXCEPTION CHECKS

    // probabilities
    std::vector<double> prob(Ks.size());
    // resulting states
    std::vector<cmat> outstates;

    if (destructive)
        outstates.resize(Ks.size(), cmat::Zero(Dsubsys_bar, Dsubsys_bar));
    else
        outstates.resize(Ks.size(), cmat::Zero(D, D));

    //************ density matrix ************//
    if (internal::check_square_mat(rA)) // square matrix
    {
        // check that dims match rho matrix
        if (!internal::check_dims_match_mat(dims, rA))
            throw exception::DimsMismatchMatrix("qpp::measure()");
        for (idx i = 0; i < Ks.size(); ++i) {
            cmat tmp = apply(rA, Ks[i], target, dims);
            if (destructive)
                tmp = ptrace(tmp, target, dims);
            prob[i] = std::abs(trace(tmp)); // probability
            if (prob[i] > 0) {
                // normalized output state
                // corresponding to measurement result i
                outstates[i] = tmp / prob[i];
            }
        }
    }
    //************ ket ************//
    else if (internal::check_cvector(rA)) // column vector
    {
        // check that dims match psi column vector
        if (!internal::check_dims_match_cvect(dims, rA))
            throw exception::DimsMismatchCvector("qpp::measure()");
        for (idx i = 0; i < Ks.size(); ++i) {
            ket tmp = apply(rA, Ks[i], target, dims);
            prob[i] = std::pow(norm(tmp), 2);
            if (prob[i] > 0) {
                // normalized output state
                // corresponding to measurement result i
                tmp /= std::sqrt(prob[i]);
                if (destructive)
                    outstates[i] = ptrace(tmp, target, dims);
                else
                    outstates[i] = tmp;
            }
        }
    } else
        throw exception::MatrixNotSquareNorCvector("qpp::measure()");

    // sample from the probability distribution
    std::discrete_distribution<idx> dd(std::begin(prob), std::end(prob));
    auto& gen =
#ifdef NO_THREAD_LOCAL_
        RandomDevices::get_instance().get_prng();
#else
        RandomDevices::get_thread_local_instance().get_prng();
#endif
    idx result = dd(gen);

    return std::make_tuple(result, prob, outstates);
}

// std::initializer_list overload, avoids ambiguity for 2-element lists, see
// http://stackoverflow.com
// /questions/26750039/ambiguity-when-using-initializer-list-as-parameter
/**
 * \brief  Measures the part \a target of the multi-partite state vector or
 * density matrix \a A using the set of Kraus operators \a Ks
 * \see qpp::measure_seq()
 *
 * \note The dimension of all \a Ks must match the dimension of \a target.
 * If \a destructive is set to true (by default), the measurement is
 * destructive, i.e. the measured subsystems are traced away.
 *
 * \param A Eigen expression
 * \param Ks Set of Kraus operators
 * \param target Subsystem indexes that are measured
 * \param dims Dimensions of the multi-partite system
 * \param destructive Destructive measurement, true by default
 * \return Tuple of: 1. Result of the measurement, 2. Vector of outcome
 * probabilities, and 3. Vector of post-measurement normalized states
 */
template <typename Derived>
std::tuple<idx, std::vector<double>, std::vector<cmat>>
measure(const Eigen::MatrixBase<Derived>& A,
        const std::initializer_list<cmat>& Ks, const std::vector<idx>& target,
        const std::vector<idx>& dims, bool destructive = true) {
    return measure(A, std::vector<cmat>(Ks), target, dims, destructive);
}

/**
 * \brief  Measures the part \a target of the multi-partite state vector or
 * density matrix \a A using the set of Kraus operators \a Ks
 * \see qpp::measure_seq()
 *
 * \note The dimension of all \a Ks must match the dimension of \a target.
 * If \a destructive is set to true (by default), the measurement is
 * destructive, i.e. the measured subsystems are traced away.

 *
 * \param A Eigen expression
 * \param Ks Set of Kraus operators
 * \param target Subsystem indexes that are measured
 * \param d Subsystem dimensions
 * \param destructive Destructive measurement, true by default
 * \return Tuple of: 1. Result of the measurement, 2. Vector of outcome
 * probabilities, and 3. Vector of post-measurement normalized states
 */
template <typename Derived>
std::tuple<idx, std::vector<double>, std::vector<cmat>>
measure(const Eigen::MatrixBase<Derived>& A, const std::vector<cmat>& Ks,
        const std::vector<idx>& target, idx d = 2, bool destructive = true) {
    const typename Eigen::MatrixBase<Derived>::EvalReturnType& rA = A.derived();

    // EXCEPTION CHECKS

    // check zero size
    if (!internal::check_nonzero_size(rA))
        throw exception::ZeroSize("qpp::measure()");

    // check valid dims
    if (d < 2)
        throw exception::DimsInvalid("qpp::measure()");
    // END EXCEPTION CHECKS

    idx n = internal::get_num_subsys(static_cast<idx>(rA.rows()), d);
    std::vector<idx> dims(n, d); // local dimensions vector

    return measure(rA, Ks, target, dims, destructive);
}

// std::initializer_list overload, avoids ambiguity for 2-element lists, see
// http://stackoverflow.com
// /questions/26750039/ambiguity-when-using-initializer-list-as-parameter
/**
 * \brief  Measures the part \a target of the multi-partite state vector or
 * density matrix \a A using the set of Kraus operators \a Ks
 * \see qpp::measure_seq()
 *
 * \note The dimension of all \a Ks must match the dimension of \a target.
 * If \a destructive is set to true (by default), the measurement is
 * destructive, i.e. the measured subsystems are traced away.
 *
 * \param A Eigen expression
 * \param Ks Set of Kraus operators
 * \param target Subsystem indexes that are measured
 * \param d Subsystem dimensions
 * \param destructive Destructive measurement, true by default
 * \return Tuple of: 1. Result of the measurement, 2. Vector of outcome
 * probabilities, and 3. Vector of post-measurement normalized states
 */
template <typename Derived>
std::tuple<idx, std::vector<double>, std::vector<cmat>>
measure(const Eigen::MatrixBase<Derived>& A,
        const std::initializer_list<cmat>& Ks, const std::vector<idx>& target,
        idx d = 2, bool destructive = true) {
    return measure(A, std::vector<cmat>(Ks), target, d, destructive);
}

/**
 * \brief Measures the part \a target of the multi-partite state vector or
 * density matrix \a A in the orthonormal basis or rank-1 projectors specified
 * by the columns of the matrix \a V
 * \see qpp::measure_seq()
 *
 * \note The dimension of \a V must match the dimension of \a target.
 * If \a destructive is set to true (by default), the measurement is
 * destructive, i.e. the measured subsystems are traced away.
 *
 * \param A Eigen expression
 * \param V Matrix whose columns represent the measurement basis vectors or the
 * bra parts of the rank-1 projectors
 * \param target Subsystem indexes that are measured
 * \param dims Dimensions of the multi-partite system
 * \param destructive Destructive measurement, true by default
 * \return Tuple of: 1. Result of the measurement, 2. Vector of outcome
 * probabilities, and 3. Vector of post-measurement normalized states
 */
template <typename Derived>
std::tuple<idx, std::vector<double>, std::vector<cmat>>
measure(const Eigen::MatrixBase<Derived>& A, const cmat& V,
        const std::vector<idx>& target, const std::vector<idx>& dims,
        bool destructive = true) {
    const typename Eigen::MatrixBase<Derived>::EvalReturnType& rA = A.derived();

    // EXCEPTION CHECKS

    // check zero-size
    if (!internal::check_nonzero_size(rA))
        throw exception::ZeroSize("qpp::measure()");

    // check that dimension is valid
    if (!internal::check_dims(dims))
        throw exception::DimsInvalid("qpp::measure()");

    // check that target is valid w.r.t. dims
    if (!internal::check_subsys_match_dims(target, dims))
        throw exception::SubsysMismatchDims("qpp::measure()");

    std::vector<idx> subsys_dims(target.size());
    for (idx i = 0; i < target.size(); ++i)
        subsys_dims[i] = dims[target[i]];

    idx Dsubsys = prod(std::begin(subsys_dims), std::end(subsys_dims));

    // check the matrix V
    if (!internal::check_nonzero_size(V))
        throw exception::ZeroSize("qpp::measure()");
    if (Dsubsys != static_cast<idx>(V.rows()))
        throw exception::DimsMismatchMatrix("qpp::measure()");
    // END EXCEPTION CHECKS

    // number of basis elements or number of rank-1 projectors
    idx M = static_cast<idx>(V.cols());

    //************ ket ************//
    if (internal::check_cvector(rA)) {
        const ket& rpsi = A.derived();
        // check that dims match state vector
        if (!internal::check_dims_match_cvect(dims, rA))
            throw exception::DimsMismatchCvector("qpp::measure()");

        std::vector<double> prob(M);    // probabilities
        std::vector<cmat> outstates(M); // resulting states

#ifdef WITH_OPENMP_
#pragma omp parallel for
#endif // WITH_OPENMP_
        for (idx i = 0; i < M; ++i) {
            if (destructive)
                outstates[i] =
                    ip(static_cast<const ket&>(V.col(i)), rpsi, target, dims);
            else
                outstates[i] = apply(rpsi, prj(V.col(i)), target, dims);
        }

        for (idx i = 0; i < M; ++i) {
            double tmp = norm(outstates[i]);
            prob[i] = tmp * tmp;
            if (prob[i] > 0) {
                // normalized output state
                // corresponding to measurement result m
                outstates[i] /= tmp;
            }
        }

        // sample from the probability distribution
        std::discrete_distribution<idx> dd(std::begin(prob), std::end(prob));
        auto& gen =
#ifdef NO_THREAD_LOCAL_
            RandomDevices::get_instance().get_prng();
#else
            RandomDevices::get_thread_local_instance().get_prng();
#endif
        idx result = dd(gen);

        return std::make_tuple(result, prob, outstates);
    }
    //************ density matrix ************//
    else if (internal::check_square_mat(rA)) {
        // check that dims match rho matrix
        if (!internal::check_dims_match_mat(dims, rA))
            throw exception::DimsMismatchMatrix("qpp::measure()");

        std::vector<cmat> Ks(M);
        for (idx i = 0; i < M; ++i)
            Ks[i] = V.col(i) * adjoint(V.col(i));

        return measure(rA, Ks, target, dims, destructive);
    }
    //************ Exception: not ket nor density matrix ************//
    throw exception::MatrixNotSquareNorCvector("qpp::measure()");
}

/**
 * \brief Measures the part \a target of the multi-partite state vector or
 * density matrix \a A in the orthonormal basis or rank-1 projectors specified
 * by the columns of the matrix \a V
 * \see qpp::measure_seq()
 *
 * \note The dimension of \a V must match the dimension of \a target.
 * If \a destructive is set to true (by default), the measurement is
 * destructive, i.e. the measured subsystems are traced away.
 *
 * \param A Eigen expression
 * \param V Matrix whose columns represent the measurement basis vectors or the
 * bra parts of the rank-1 projectors
 * \param target Subsystem indexes that are measured
 * \param d Subsystem dimensions
 * \param destructive Destructive measurement, true by default
 * \return Tuple of: 1. Result of the measurement, 2.
 * Vector of outcome probabilities, and 3. Vector of post-measurement
 * normalized states
 */
template <typename Derived>
std::tuple<idx, std::vector<double>, std::vector<cmat>>
measure(const Eigen::MatrixBase<Derived>& A, const cmat& V,
        const std::vector<idx>& target, idx d = 2, bool destructive = true) {
    const typename Eigen::MatrixBase<Derived>::EvalReturnType& rA = A.derived();

    // EXCEPTION CHECKS

    // check zero size
    if (!internal::check_nonzero_size(rA))
        throw exception::ZeroSize("qpp::measure()");

    // check valid dims
    if (d < 2)
        throw exception::DimsInvalid("qpp::measure()");
    // END EXCEPTION CHECKS

    idx n = internal::get_num_subsys(static_cast<idx>(rA.rows()), d);
    std::vector<idx> dims(n, d); // local dimensions vector

    return measure(rA, V, target, dims, destructive);
}

/**
 * \brief Sequentially measures the part \a target of the multi-partite state
 * vector or density matrix \a A in the computational basis
 * \see qpp::measure()
 *
 * \note If \a destructive is set to true (by default), the measurement is
 * destructive, i.e. the measured subsystems are traced away.
 *
 * \param A Eigen expression
 * \param target Subsystem indexes that are measured
 * \param dims Dimensions of the multi-partite system
 * \param destructive Destructive measurement, true by default
 * \return Tuple of: 1. Vector of outcome results of the
 * measurement (ordered in increasing order with respect to \a target, i.e.
 * first measurement result corresponds to the subsystem with the smallest
 * index), 2. Outcome probability, and 3. Post-measurement normalized state
 */
template <typename Derived>
std::tuple<std::vector<idx>, double, cmat>
measure_seq(const Eigen::MatrixBase<Derived>& A, std::vector<idx> target,
            std::vector<idx> dims, bool destructive = true) {
    //    typename std::remove_const<
    //            typename Eigen::MatrixBase<Derived>::EvalReturnType
    //    >::type rA = A.derived();

    dyn_mat<typename Derived::Scalar> rA = A.derived();

    // EXCEPTION CHECKS

    // check zero-size
    if (!internal::check_nonzero_size(rA))
        throw exception::ZeroSize("qpp::measure_seq()");
    // check that dimension is valid
    if (!internal::check_dims(dims))
        throw exception::DimsInvalid("qpp::measure_seq()");

    // check square matrix or column vector
    if (internal::check_square_mat(rA)) {
        // check that dims match rho matrix
        if (!internal::check_dims_match_mat(dims, rA))
            throw exception::DimsMismatchMatrix("qpp::measure_seq()");
    } else if (internal::check_cvector(rA)) {
        // check that dims match psi column vector
        if (!internal::check_dims_match_cvect(dims, rA))
            throw exception::DimsMismatchCvector("qpp::measure_seq()");
    } else
        throw exception::MatrixNotSquareNorCvector("qpp::measure_seq()");

    // check that target is valid w.r.t. dims
    if (!internal::check_subsys_match_dims(target, dims))
        throw exception::SubsysMismatchDims("qpp::measure_seq()");
    // END EXCEPTION CHECKS

    std::vector<idx> result;
    double prob = 1;

    // sort target in decreasing order,
    // the order of measurements does not matter
    std::sort(std::begin(target), std::end(target), std::greater<idx>{});

    //************ density matrix or column vector ************//
    while (target.size() > 0) {
        auto tmp = measure(rA, Gates::get_instance().Id(dims[target[0]]),
                           {target[0]}, dims, destructive);
        result.emplace_back(std::get<0>(tmp));
        prob *= std::get<1>(tmp)[std::get<0>(tmp)];
        rA = std::get<2>(tmp)[std::get<0>(tmp)];

        // remove the subsystem
        dims.erase(std::next(std::begin(dims), target[0]));
        target.erase(std::begin(target));
    }
    // order result in increasing order with respect to target
    std::reverse(std::begin(result), std::end(result));

    return std::make_tuple(result, prob, rA);
}

/**
 * \brief Sequentially measures the part \a target of the multi-partite state
 * vector or density matrix \a A in the computational basis
 * \see qpp::measure()
 *
 * \note If \a destructive is set to true (by default), the measurement is
 * destructive, i.e. the measured subsystems are traced away.
 *
 * \param A Eigen expression
 * \param target Subsystem indexes that are measured
 * \param d Subsystem dimensions
 * \param destructive Destructive measurement, true by default
 * \return Tuple of: 1. Vector of outcome results of the
 * measurement (ordered in increasing order with respect to \a target, i.e.
 * first measurement result corresponds to the subsystem with the smallest
 * index), 2. Outcome probability, and 3. Post-measurement normalized state
 */
template <typename Derived>
std::tuple<std::vector<idx>, double, cmat>
measure_seq(const Eigen::MatrixBase<Derived>& A, std::vector<idx> target,
            idx d = 2, bool destructive = true) {
    const typename Eigen::MatrixBase<Derived>::EvalReturnType& rA = A.derived();

    // EXCEPTION CHECKS

    // check zero size
    if (!internal::check_nonzero_size(rA))
        throw exception::ZeroSize("qpp::measure_seq()");

    // check valid dims
    if (d < 2)
        throw exception::DimsInvalid("qpp::measure_seq()");
    // END EXCEPTION CHECKS

    idx n = internal::get_num_subsys(static_cast<idx>(rA.rows()), d);
    std::vector<idx> dims(n, d); // local dimensions vector

    return measure_seq(rA, target, dims, destructive);
}

} /* namespace qpp */

#endif /* INSTRUMENTS_H_ */
