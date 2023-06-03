/*
 * This file is part of Quantum++.
 *
 * Copyright (c) 2013 - 2023 softwareQ Inc. All rights reserved.
 *
 * MIT License
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
 * \file operations.hpp
 * \brief Quantum operation functions
 */

#ifndef QPP_OPERATIONS_HPP_
#define QPP_OPERATIONS_HPP_

namespace qpp {
/**
 * \brief Applies the gate \a A to the part \a target of the multi-partite state
 * vector or density matrix \a state
 *
 * \note The dimension of the gate \a A must match the dimension of \a target
 *
 * \param state Eigen expression
 * \param A Eigen expression
 * \param target Subsystem indexes where the gate \a A is applied
 * \param dims Dimensions of the multi-partite system
 * \return Gate \a A applied to the part \a target of \a state
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] expr_t<Derived1>
apply(const Eigen::MatrixBase<Derived1>& state,
      const Eigen::MatrixBase<Derived2>& A, const std::vector<idx>& target,
      const std::vector<idx>& dims) {
    const expr_t<Derived1>& rstate = state.derived();
    const dyn_mat<typename Derived2::Scalar>& rA = A.derived();

    // EXCEPTION CHECKS

    // check types
    if (!std::is_same<typename Derived1::Scalar,
                      typename Derived2::Scalar>::value)
        throw exception::TypeMismatch("qpp::apply()", "A/state");

    // check zero sizes
    if (!internal::check_nonzero_size(rA))
        throw exception::ZeroSize("qpp::apply()", "A");

    // check zero sizes
    if (!internal::check_nonzero_size(rstate))
        throw exception::ZeroSize("qpp::apply()", "state");

    // check zero sizes
    if (!internal::check_nonzero_size(target))
        throw exception::ZeroSize("qpp::apply()", "target");

    // check square matrix for the gate
    if (!internal::check_square_mat(rA))
        throw exception::MatrixNotSquare("qpp::apply()", "A");

    // check that dimension is valid
    if (!internal::check_dims(dims))
        throw exception::DimsInvalid("qpp::apply()", "dims");

    // check that target is valid w.r.t. dims
    if (!internal::check_subsys_match_dims(target, dims))
        throw exception::SubsysMismatchDims("qpp::apply()", "dims/target");

    // check valid state and matching dimensions
    if (internal::check_cvector(rstate)) {
        if (!internal::check_dims_match_cvect(dims, rstate))
            throw exception::DimsMismatchCvector("qpp::apply()", "dims/state");
    } else if (internal::check_square_mat(rstate)) {
        if (!internal::check_dims_match_mat(dims, rstate))
            throw exception::DimsMismatchMatrix("qpp::apply()", "dims/state");
    } else
        throw exception::MatrixNotSquareNorCvector("qpp::apply()", "state");

    // check that gate matches the dimensions of the target
    idx gate_size = target.size(); // number of subsystems of the target
    std::vector<idx> subsys_dims(gate_size);
    for (idx i = 0; i < gate_size; ++i)
        subsys_dims[i] = dims[target[i]];
    if (!internal::check_dims_match_mat(subsys_dims, rA))
        throw exception::MatrixMismatchSubsys("qpp::apply()", "A/dims/target");
    // END EXCEPTION CHECKS

    idx D = static_cast<idx>(rstate.rows()); // total dimension
    idx n = dims.size();                     // total number of subsystems
    idx DA = static_cast<idx>(rA.rows());    // dimension of gate subsystem

    idx Cdims[internal::maxn];      // local dimensions total
    idx CdimsA[internal::maxn];     // local dimensions gate
    idx CdimsA_bar[internal::maxn]; // local dimensions complement

    // compute the complementary subsystem of gate w.r.t. dims
    std::vector<idx> gate_bar = complement(target, n);
    // number of subsystems that are complementary to the gate
    idx gate_bar_size = gate_bar.size();

    idx DA_bar = 1; // dimension of the complement
    for (idx i = 0; i < gate_bar_size; ++i)
        DA_bar *= dims[gate_bar[i]];

    for (idx k = 0; k < n; ++k)
        Cdims[k] = dims[k];
    for (idx k = 0; k < gate_size; ++k)
        CdimsA[k] = dims[target[k]];
    for (idx k = 0; k < gate_bar_size; ++k)
        CdimsA_bar[k] = dims[gate_bar[k]];

    // worker, computes the coefficient and the index for the ket case, used in
    // #pragma omp parallel for collapse
    auto coeff_idx_ket =
        [&](idx m_,
            idx r_) noexcept -> std::pair<typename Derived1::Scalar, idx> {
        idx indx = 0;
        typename Derived1::Scalar coeff = 0;

        idx Cmidx[internal::maxn];      // the total multi-index
        idx CmidxA[internal::maxn];     // the gate multi-index
        idx CmidxA_bar[internal::maxn]; // the complement multi-index

        // compute the index

        // set the complement multi-index
        internal::n2multiidx(r_, gate_bar_size, CdimsA_bar, CmidxA_bar);
        for (idx k = 0; k < gate_bar_size; ++k) {
            Cmidx[gate_bar[k]] = CmidxA_bar[k];
        }

        // set the gate multi-index
        internal::n2multiidx(m_, gate_size, CdimsA, CmidxA);
        for (idx k = 0; k < gate_size; ++k) {
            Cmidx[target[k]] = CmidxA[k];
        }

        // we now got the total index
        indx = internal::multiidx2n(Cmidx, n, Cdims);

        // compute the coefficient
        for (idx n_ = 0; n_ < DA; ++n_) {
            internal::n2multiidx(n_, gate_size, CdimsA, CmidxA);
            for (idx k = 0; k < gate_size; ++k) {
                Cmidx[target[k]] = CmidxA[k];
            }
            coeff += rA(m_, n_) * rstate(internal::multiidx2n(Cmidx, n, Cdims));
        }

        return std::make_pair(coeff, indx);
    }; /* end coeff_idx_ket */

    // worker, computes the coefficient and the index for the density matrix
    // case, used in #pragma omp parallel for collapse
    auto coeff_idx_rho = [&](idx m1_, idx r1_, idx m2_, idx r2_) noexcept
        -> std::tuple<typename Derived1::Scalar, idx, idx> {
        idx idxrow = 0;
        idx idxcol = 0;
        typename Derived1::Scalar coeff = 0, lhs = 1, rhs = 1;

        idx Cmidxrow[internal::maxn];      // the total row multi-index
        idx Cmidxcol[internal::maxn];      // the total col multi-index
        idx CmidxArow[internal::maxn];     // the gate row multi-index
        idx CmidxAcol[internal::maxn];     // the gate col multi-index
        idx CmidxA_barrow[internal::maxn]; // the complement row multi-index
        idx CmidxA_barcol[internal::maxn]; // the complement col multi-index

        // compute the ket/bra indexes

        // set the complement multi-index
        internal::n2multiidx(r1_, gate_bar_size, CdimsA_bar, CmidxA_barrow);
        internal::n2multiidx(r2_, gate_bar_size, CdimsA_bar, CmidxA_barcol);
        for (idx k = 0; k < gate_bar_size; ++k) {
            Cmidxrow[gate_bar[k]] = CmidxA_barrow[k];
            Cmidxcol[gate_bar[k]] = CmidxA_barcol[k];
        }

        // set the gate multi-index
        internal::n2multiidx(m1_, gate_size, CdimsA, CmidxArow);
        internal::n2multiidx(m2_, gate_size, CdimsA, CmidxAcol);
        for (idx k = 0; k < gate_size; ++k) {
            Cmidxrow[target[k]] = CmidxArow[k];
            Cmidxcol[target[k]] = CmidxAcol[k];
        }

        // we now got the total row/col indexes
        idxrow = internal::multiidx2n(Cmidxrow, n, Cdims);
        idxcol = internal::multiidx2n(Cmidxcol, n, Cdims);

        // compute the coefficient
        for (idx n1_ = 0; n1_ < DA; ++n1_) {
            internal::n2multiidx(n1_, gate_size, CdimsA, CmidxArow);
            for (idx k = 0; k < gate_size; ++k) {
                Cmidxrow[target[k]] = CmidxArow[k];
            }
            idx idxrowtmp = internal::multiidx2n(Cmidxrow, n, Cdims);

            lhs = rA(m1_, n1_);

            for (idx n2_ = 0; n2_ < DA; ++n2_) {
                internal::n2multiidx(n2_, gate_size, CdimsA, CmidxAcol);
                for (idx k = 0; k < gate_size; ++k) {
                    Cmidxcol[target[k]] = CmidxAcol[k];
                }

                rhs = adjoint(rA)(n2_, m2_);

                idx idxcoltmp = internal::multiidx2n(Cmidxcol, n, Cdims);

                coeff += lhs * rstate(idxrowtmp, idxcoltmp) * rhs;
            }
        }

        return std::make_tuple(coeff, idxrow, idxcol);
    }; /* end coeff_idx_rho */

    //************ ket ************//
    if (internal::check_cvector(rstate)) // we have a ket
    {
        // check that dims match state vector
        if (!internal::check_dims_match_cvect(dims, rstate))
            throw exception::DimsMismatchCvector("qpp::apply()", "dims/state");
        if (D == 1)
            return rstate;

        expr_t<Derived1> result = rstate;

#ifdef QPP_OPENMP
// NOLINTNEXTLINE
#pragma omp parallel for collapse(2)
#endif // QPP_OPENMP
        for (idx m = 0; m < DA; ++m)
            for (idx r = 0; r < DA_bar; ++r) {
                result(coeff_idx_ket(m, r).second) = coeff_idx_ket(m, r).first;
            }

        return result;
    }
    //************ density matrix ************//
    else // we have a density operator
    {
        // check that dims match state matrix
        if (!internal::check_dims_match_mat(dims, rstate))
            throw exception::DimsMismatchMatrix("qpp::apply()", "dims/state");

        if (D == 1)
            return rstate;

        expr_t<Derived1> result = rstate;

#ifdef QPP_OPENMP
// NOLINTNEXTLINE
#pragma omp parallel for collapse(4)
#endif // QPP_OPENMP
        for (idx m1 = 0; m1 < DA; ++m1)
            for (idx r1 = 0; r1 < DA_bar; ++r1)
                for (idx m2 = 0; m2 < DA; ++m2)
                    for (idx r2 = 0; r2 < DA_bar; ++r2) {
                        auto coeff_idxes = coeff_idx_rho(m1, r1, m2, r2);
                        result(std::get<1>(coeff_idxes),
                               std::get<2>(coeff_idxes)) =
                            std::get<0>(coeff_idxes);
                    }

        return result;
    }
}

/**
 * \brief Applies the gate \a A to the part \a target of the multi-partite state
 * vector or density matrix \a state
 *
 * \note The dimension of the gate \a A must match the dimension of \a target
 *
 * \param state Eigen expression
 * \param A Eigen expression
 * \param target Subsystem indexes where the gate \a A is applied
 * \param d Subsystem dimensions
 * \return Gate \a A applied to the part \a target of \a state
 */
template <typename Derived1, typename Derived2>
expr_t<Derived1> apply(const Eigen::MatrixBase<Derived1>& state,
                       const Eigen::MatrixBase<Derived2>& A,
                       const std::vector<idx>& target, idx d = 2) {
    const expr_t<Derived1>& rstate = state.derived();
    const dyn_mat<typename Derived1::Scalar>& rA = A.derived();

    // EXCEPTION CHECKS

    // check zero size
    if (!internal::check_nonzero_size(rstate))
        throw exception::ZeroSize("qpp::apply()", "state");

    // check valid dims
    if (d < 2)
        throw exception::DimsInvalid("qpp::apply()", "d");
    // END EXCEPTION CHECKS

    idx n = internal::get_num_subsys(static_cast<idx>(rstate.rows()), d);
    std::vector<idx> dims(n, d); // local dimensions vector

    return apply(rstate, rA, target, dims);
}

/**
 * \brief Applies the channel specified by the set of Kraus operators \a Ks to
 * the density matrix \a A
 *
 * \note The Kraus operators can have their range different from their domain
 * (i.e., they can be rectangular matrices)
 *
 * \param A Eigen expression
 * \param Ks Set of Kraus operators
 * \return Output density matrix after the action of the channel
 */
template <typename Derived>
[[qpp::parallel]] cmat apply(const Eigen::MatrixBase<Derived>& A,
                             const std::vector<cmat>& Ks) {
    const cmat& rA = A.derived();

    // EXCEPTION CHECKS

    if (!internal::check_nonzero_size(rA))
        throw exception::ZeroSize("qpp::apply()", "A");
    if (!internal::check_square_mat(rA))
        throw exception::MatrixNotSquare("qpp::apply()", "A");
    if (Ks.empty())
        throw exception::ZeroSize("qpp::apply()", "Ks");
    if (Ks[0].cols() != rA.rows())
        throw exception::DimsMismatchMatrix("qpp::apply()", "Ks[0]");
    for (auto&& elem : Ks)
        if (elem.rows() != Ks[0].rows() || elem.cols() != Ks[0].cols())
            throw exception::DimsNotEqual("qpp::apply()", "K");
    // END EXCEPTION CHECKS

    idx Dout = Ks[0].rows();
    cmat result = cmat::Zero(Dout, Dout);

#ifdef QPP_OPENMP
// NOLINTNEXTLINE
#pragma omp parallel for
#endif // QPP_OPENMP
    for (const auto& K : Ks) {
#ifdef QPP_OPENMP
#pragma omp critical
#endif // QPP_OPENMP
        { result += K * rA * adjoint(K); }
    }

    return result;
}

/**
 * \brief Applies the channel specified by the set of Kraus operators \a Ks to
 * the part \a target of the multi-partite density matrix \a A
 *
 * \param A Eigen expression
 * \param Ks Set of Kraus operators
 * \param target Subsystem indexes where the Kraus operators \a Ks are applied
 * \param dims Dimensions of the multi-partite system
 * \return Output density matrix after the action of the channel
 */
template <typename Derived>
[[qpp::parallel]] cmat
apply(const Eigen::MatrixBase<Derived>& A, const std::vector<cmat>& Ks,
      const std::vector<idx>& target, const std::vector<idx>& dims) {
    const cmat& rA = A.derived();

    // EXCEPTION CHECKS

    // check zero sizes
    if (!internal::check_nonzero_size(rA))
        throw exception::ZeroSize("qpp::apply()", "A");

    // check zero sizes
    if (!internal::check_nonzero_size(target))
        throw exception::ZeroSize("qpp::apply()", "target");

    // check square matrix for the A
    if (!internal::check_square_mat(rA))
        throw exception::MatrixNotSquare("qpp::apply()", "A");

    // check that dimension is valid
    if (!internal::check_dims(dims))
        throw exception::DimsInvalid("qpp::apply()", "dims");

    // check that target is valid w.r.t. dims
    if (!internal::check_subsys_match_dims(target, dims))
        throw exception::SubsysMismatchDims("qpp::apply()", "dims/target");

    // check valid state and matching dimensions
    if (internal::check_cvector(rA)) {
        if (!internal::check_dims_match_cvect(dims, rA))
            throw exception::DimsMismatchCvector("qpp::apply()", "A/dims");
    } else if (internal::check_square_mat(rA)) {
        if (!internal::check_dims_match_mat(dims, rA))
            throw exception::DimsMismatchMatrix("qpp::apply()", "A/dims");
    } else
        throw exception::MatrixNotSquareNorCvector("qpp::apply()", "A");

    std::vector<idx> subsys_dims(target.size());
    for (idx i = 0; i < static_cast<idx>(target.size()); ++i)
        subsys_dims[i] = dims[target[i]];

    // check the Kraus operators
    if (Ks.empty())
        throw exception::ZeroSize("qpp::apply()", "Ks");
    if (!internal::check_square_mat(Ks[0]))
        throw exception::MatrixNotSquare("qpp::apply()", "Ks[0]");
    if (!internal::check_dims_match_mat(subsys_dims, Ks[0]))
        throw exception::MatrixMismatchSubsys("qpp::apply()",
                                              "dims/Ks[0]/target");
    for (auto&& elem : Ks)
        if (elem.rows() != Ks[0].rows() || elem.cols() != Ks[0].cols())
            throw exception::DimsNotEqual("qpp::apply()", "K");
    // END EXCEPTION CHECKS

    cmat result = cmat::Zero(rA.rows(), rA.cols());

#ifdef QPP_OPENMP
// NOLINTNEXTLINE
#pragma omp parallel for
#endif // QPP_OPENMP
    for (const auto& K : Ks) {
#ifdef QPP_OPENMP
#pragma omp critical
#endif // QPP_OPENMP
        { result += apply(rA, K, target, dims); }
    }

    return result;
}

/**
 * \brief Applies the channel specified by the set of Kraus operators \a Ks to
 * the part \a target of the multi-partite density matrix \a A
 *
 * \param A Eigen expression
 * \param Ks Set of Kraus operators
 * \param target Subsystem indexes where the Kraus operators \a Ks are applied
 * \param d Subsystem dimensions
 * \return Output density matrix after the action of the channel
 */
template <typename Derived>
cmat apply(const Eigen::MatrixBase<Derived>& A, const std::vector<cmat>& Ks,
           const std::vector<idx>& target, idx d = 2) {
    const cmat& rA = A.derived();

    // EXCEPTION CHECKS

    // check zero sizes
    if (!internal::check_nonzero_size(rA))
        throw exception::ZeroSize("qpp::apply()", "A");

    // check valid dims
    if (d < 2)
        throw exception::DimsInvalid("qpp::apply()", "d");
    // END EXCEPTION CHECKS

    idx n = internal::get_num_subsys(static_cast<idx>(rA.rows()), d);
    std::vector<idx> dims(n, d); // local dimensions vector

    return apply(rA, Ks, target, dims);
}

/**
 * \brief Choi matrix
 * \see qpp::choi2kraus()
 *
 * Constructs the Choi matrix of the channel specified by the set of Kraus
 * operators \a Ks in the standard operator basis \f$\{|i\rangle\langle j|\}\f$
 * ordered in lexicographical order, i.e.
 * \f$|0\rangle\langle 0|\f$, \f$|0\rangle\langle 1|\f$ etc.
 *
 * \note The Kraus operators can have their range different from their domain
 * (i.e., they can be rectangular matrices). The superoperator matrix \f$S\f$
 * and the Choi matrix \f$C\f$ are related by \f$S_{ab,mn} = C_{ma,nb}\f$.
 *
 * \param Ks Set of Kraus operators
 * \return Choi matrix
 */
[[qpp::parallel]] inline cmat kraus2choi(const std::vector<cmat>& Ks) {
    // EXCEPTION CHECKS

    if (Ks.empty())
        throw exception::ZeroSize("qpp::kraus2choi()", "Ks");
    if (!internal::check_nonzero_size(Ks[0]))
        throw exception::ZeroSize("qpp::kraus2choi()", "Ks[0]");
    for (auto&& elem : Ks)
        if (elem.rows() != Ks[0].rows() || elem.cols() != Ks[0].cols())
            throw exception::DimsNotEqual("qpp::kraus2choi()", "K");
    // END EXCEPTION CHECKS

    idx Din = static_cast<idx>(Ks[0].cols());
    idx Dout = static_cast<idx>(Ks[0].rows());

    // construct the Din x Din \sum |jj> vector
    // (un-normalized maximally entangled state)
    cmat MES = cmat::Zero(Din * Din, 1);

#ifdef QPP_OPENMP
// NOLINTNEXTLINE
#pragma omp parallel for
#endif // QPP_OPENMP
    for (idx a = 0; a < Din; ++a)
        MES(a * Din + a, 0) = 1;

    cmat Omega = MES * adjoint(MES);

    cmat result = cmat::Zero(Din * Dout, Din * Dout);

#ifdef QPP_OPENMP
// NOLINTNEXTLINE
#pragma omp parallel for
#endif // QPP_OPENMP
    for (const auto& K : Ks) {
#ifdef QPP_OPENMP
#pragma omp critical
#endif // QPP_OPENMP
        {
            result += kron(cmat::Identity(Din, Din), K) * Omega *
                      adjoint(kron(cmat::Identity(Din, Din), K));
        }
    }

    return result;
}

/**
 * \brief Orthogonal Kraus operators from Choi matrix
 * \see qpp::kraus2choi()
 *
 * Extracts a set of orthogonal (under Hilbert-Schmidt operator norm) Kraus
 * operators from the Choi matrix \a A
 *
 * \note The Kraus operators can have their range different from their domain
 * (i.e., they can be rectangular matrices). The Kraus operators satisfy
 * \f$Tr(K_i^\dagger K_j)=\delta_{ij}\f$ for all \f$i\neq j\f$.
 *
 * \param A Choi matrix
 * \param Din Dimension of the Kraus input Hilbert space
 * \param Dout Dimension of the Kraus output Hilbert space
 * \return Set of orthogonal Kraus operators
 */
inline std::vector<cmat> choi2kraus(const cmat& A, idx Din, idx Dout) {
    // EXCEPTION CHECKS

    if (!internal::check_nonzero_size(A))
        throw exception::ZeroSize("qpp::choi2kraus()", "A");
    if (!internal::check_square_mat(A))
        throw exception::MatrixNotSquare("qpp::choi2kraus()", "A");
    // check equal dimensions
    if (Din * Dout != static_cast<idx>(A.rows()))
        throw exception::DimsInvalid("qpp::choi2kraus()", "A/Din/Dout");
    // END EXCEPTION CHECKS

    rmat ev = hevals(A);
    cmat evec = hevects(A);
    std::vector<cmat> result;

    for (idx i = 0; i < Din * Dout; ++i) {
        if (std::abs(ev(i)) > 0)
            result.emplace_back(std::sqrt(std::abs(ev(i))) *
                                reshape(evec.col(i), Dout, Din));
    }

    return result;
}

/**
 * \brief Orthogonal Kraus operators from Choi matrix
 * \see qpp::kraus2choi()
 *
 * Extracts a set of orthogonal (under Hilbert-Schmidt operator norm) Kraus
 * operators from the Choi matrix \a A
 *
 * \note The Kraus operators are assumed to have their range equal to their
 * domain (i.e., they are square matrices). The Kraus operators satisfy
 * \f$Tr(K_i^\dagger K_j)=\delta_{ij}\f$ for all \f$i\neq j\f$.
 *
 * \param A Choi matrix
 * \return Set of orthogonal Kraus operators
 */
inline std::vector<cmat> choi2kraus(const cmat& A) {
    // EXCEPTION CHECKS

    if (!internal::check_nonzero_size(A))
        throw exception::ZeroSize("qpp::choi2kraus()", "A");
    if (!internal::check_square_mat(A))
        throw exception::MatrixNotSquare("qpp::choi2kraus()", "A");
    // check equal dimensions
    idx D = internal::get_num_subsys(static_cast<idx>(A.rows()), 2);
    if (D * D != static_cast<idx>(A.rows()))
        throw exception::DimsInvalid("qpp::choi2kraus()", "A");
    // END EXCEPTION CHECKS

    return choi2kraus(A, D, D);
}

/**
 * \brief Converts Choi matrix to superoperator matrix
 * \see qpp::super2choi()
 *
 * \note The superoperator can have its range different from its domain
 * (i.e., it can be a rectangular matrix)
 *
 * \param A Choi matrix
 * \param Din Dimension of the Kraus input Hilbert space
 * \param Dout Dimension of the Kraus output Hilbert space
 * \return Superoperator matrix
 */
[[qpp::parallel]] inline cmat choi2super(const cmat& A, idx Din, idx Dout) {
    // EXCEPTION CHECKS

    if (!internal::check_nonzero_size(A))
        throw exception::ZeroSize("qpp::choi2super()", "A");
    if (!internal::check_square_mat(A))
        throw exception::MatrixNotSquare("qpp::choi2super()", "A");
    // check equal dimensions
    if (Din * Dout != static_cast<idx>(A.rows()))
        throw exception::DimsInvalid("qpp::choi2super()", "A/Din/Dout");
    // END EXCEPTION CHECKS

    cmat result(Dout * Dout, Din * Din);

#ifdef QPP_OPENMP
// NOLINTNEXTLINE
#pragma omp parallel for collapse(4)
#endif // QPP_OPENMP
    for (idx a = 0; a < Dout; ++a)
        for (idx b = 0; b < Dout; ++b)
            for (idx m = 0; m < Din; ++m)
                for (idx n = 0; n < Din; ++n)
                    result(a * Dout + b, m * Din + n) =
                        A(m * Dout + a, n * Dout + b);

    return result;
}

/**
 * \brief Converts Choi matrix to superoperator matrix
 * \see qpp::super2choi()
 *
 * \note The superoperator is assumed to have the its range equal to its domain
 * (i.e., its Kraus operators are square matrices).
 *
 * \param A Choi matrix
 * \return Superoperator matrix
 */
inline cmat choi2super(const cmat& A) {
    // EXCEPTION CHECKS

    if (!internal::check_nonzero_size(A))
        throw exception::ZeroSize("qpp::choi2super()", "A");
    if (!internal::check_square_mat(A))
        throw exception::MatrixNotSquare("qpp::choi2super()", "A");
    // check equal dimensions
    idx D = internal::get_num_subsys(static_cast<idx>(A.rows()), 2);
    if (D * D != static_cast<idx>(A.rows()))
        throw exception::DimsInvalid("qpp::choi2super()", "A");
    // END EXCEPTION CHECKS

    return choi2super(A, D, D);
}

/**
 * \brief Converts superoperator matrix to Choi matrix
 * \see qpp::choi2super()
 *
 * \note The superoperator can have its range different from its domain
 * (i.e., it can be a rectangular matrix)
 *
 * \param A Superoperator matrix
 * \return Choi matrix
 */
[[qpp::parallel]] inline cmat super2choi(const cmat& A) {
    // EXCEPTION CHECKS

    if (!internal::check_nonzero_size(A))
        throw exception::ZeroSize("qpp::super2choi()", "A");
    idx Din = internal::get_dim_subsys(static_cast<idx>(A.cols()), 2);
    idx Dout = internal::get_dim_subsys(static_cast<idx>(A.rows()), 2);
    // check equal dimensions
    if (Din * Din != static_cast<idx>(A.cols()))
        throw exception::DimsInvalid("qpp::super2choi()", "A");
    if (Dout * Dout != static_cast<idx>(A.rows()))
        throw exception::DimsInvalid("qpp::super2choi()", "A");
    // END EXCEPTION CHECKS

    cmat result(Din * Dout, Din * Dout);

#ifdef QPP_OPENMP
// NOLINTNEXTLINE
#pragma omp parallel for collapse(4)
#endif // QPP_OPENMP
    for (idx a = 0; a < Dout; ++a)
        for (idx b = 0; b < Dout; ++b)
            for (idx m = 0; m < Din; ++m)
                for (idx n = 0; n < Din; ++n)
                    result(m * Dout + a, n * Dout + b) =
                        A(a * Dout + b, m * Din + n);

    return result;
}

/**
 * \brief Superoperator matrix
 * \see qpp::super2kraus()
 *
 * Constructs the superoperator matrix of the channel specified by the set of
 * Kraus operators \a Ks in the standard operator basis
 * \f$\{|i\rangle\langle j|\}\f$ ordered in lexicographical order, i.e.
 * \f$|0\rangle\langle 0|\f$, \f$|0\rangle\langle 1|\f$ etc.
 *
 * \note The Kraus operators can have their range different from their domain
 * (i.e., they can be rectangular matrices)
 *
 * \param Ks Set of Kraus operators
 * \return Superoperator matrix
 */
[[qpp::parallel]] inline cmat kraus2super(const std::vector<cmat>& Ks) {
    // EXCEPTION CHECKS

    if (Ks.empty())
        throw exception::ZeroSize("qpp::kraus2super()", "Ks");
    if (!internal::check_nonzero_size(Ks[0]))
        throw exception::ZeroSize("qpp::kraus2super()", "Ks[0]");
    for (auto&& elem : Ks)
        if (elem.rows() != Ks[0].rows() || elem.cols() != Ks[0].cols())
            throw exception::DimsNotEqual("qpp::kraus2super()", "K");
    // END EXCEPTION CHECKS

    idx Din = static_cast<idx>(Ks[0].cols());
    idx Dout = static_cast<idx>(Ks[0].rows());

    cmat MN = cmat::Zero(Din, Din);
    cmat EMN = cmat::Zero(Dout, Dout);
    bra A = bra::Zero(Dout);
    ket B = ket::Zero(Dout);
    cmat result(Dout * Dout, Din * Din);

#ifdef QPP_OPENMP
// NOLINTNEXTLINE
#pragma omp parallel for collapse(2)
#endif // QPP_OPENMP
    for (idx m = 0; m < Din; ++m) {
        for (idx n = 0; n < Din; ++n) {
#ifdef QPP_OPENMP
#pragma omp critical
#endif        // QPP_OPENMP
            { // DO NOT ERASE THIS CURLY BRACKET!!!! OMP CRITICAL CODE
                // compute E(|m><n|)
                MN(m, n) = 1;
                for (const auto& K : Ks)
                    EMN += K * MN * adjoint(K);
                MN(m, n) = 0;

                for (idx a = 0; a < Dout; ++a) {
                    A(a) = 1;
                    for (idx b = 0; b < Dout; ++b) {
                        // compute result(ab,mn)=<a|E(|m><n)|b>
                        B(b) = 1;
                        result(a * Dout + b, m * Din + n) =
                            static_cast<cmat>(A * EMN * B).value();
                        B(b) = 0;
                    }
                    A(a) = 0;
                }
                EMN = cmat::Zero(Dout, Dout);
            } // DO NOT ERASE THIS CURLY BRACKET!!!! OMP CRITICAL CODE
        }
    }

    return result;
}

/**
 * \brief Orthogonal Kraus operators from superoperator matrix
 * \see qpp::kraus2super()
 *
 * Extracts a set of orthogonal (under the Hilbert-Schmidt operator norm) Kraus
 * operators from the superoperator matrix \a A
 *
 * \note The superoperator can have its range different from its domain
 * (i.e., it can be a rectangular matrix). The Kraus operators satisfy
 * \f$Tr(K_i^\dagger K_j)=\delta_{ij}\f$ for all \f$i\neq j\f$.
 *
 * \param A Superoperator matrix
 * \return Set of orthogonal Kraus operators
 */
inline std::vector<cmat> super2kraus(const cmat& A) {
    // EXCEPTION CHECKS

    if (!internal::check_nonzero_size(A))
        throw exception::ZeroSize("qpp::super2kraus()", "A");
    idx Din = internal::get_dim_subsys(static_cast<idx>(A.cols()), 2);
    idx Dout = internal::get_dim_subsys(static_cast<idx>(A.rows()), 2);
    // check equal dimensions
    if (Din * Din != static_cast<idx>(A.cols()))
        throw exception::DimsInvalid("qpp::super2kraus()", "A");
    if (Dout * Dout != static_cast<idx>(A.rows()))
        throw exception::DimsInvalid("qpp::super2kraus()", "A");
    // END EXCEPTION CHECKS

    return choi2kraus(super2choi(A), Din, Dout);
}

/**
 * \brief Partial trace
 * \see qpp::ptrace2()
 *
 * Partial trace over the first subsystem of bi-partite state vector or density
 * matrix
 *
 * \param A Eigen expression
 * \param dims Dimensions of the bi-partite system
 * \return Partial trace \f$Tr_{A}(\cdot)\f$ over the first subsytem \f$A\f$
 * in a bi-partite system \f$A\otimes B\f$, as a dynamic matrix over the same
 * scalar field as \a A
 */
template <typename Derived>
[[qpp::critical, qpp::parallel]] dyn_mat<typename Derived::Scalar>
ptrace1(const Eigen::MatrixBase<Derived>& A, const std::vector<idx>& dims) {
    const typename Eigen::MatrixBase<Derived>::EvalReturnType& rA = A.derived();

    // EXCEPTION CHECKS

    // check zero-size
    if (!internal::check_nonzero_size(rA))
        throw exception::ZeroSize("qpp::ptrace1()", "A");

    // check that dims is a valid dimension vector
    if (!internal::check_dims(dims))
        throw exception::DimsInvalid("qpp::ptrace1()", "dims");

    // check dims has only 2 elements
    if (dims.size() != 2)
        throw exception::NotBipartite("qpp::ptrace1()", "dims");
    // END EXCEPTION CHECKS

    idx DA = dims[0];
    idx DB = dims[1];

    dyn_mat<typename Derived::Scalar> result =
        dyn_mat<typename Derived::Scalar>::Zero(DB, DB);

    //************ ket ************//
    if (internal::check_cvector(rA)) // we have a ket
    {
        // check that dims match the dimension of A
        if (!internal::check_dims_match_cvect(dims, rA))
            throw exception::DimsMismatchCvector("qpp::ptrace1()", "A/dims");

        auto worker = [&](idx i, idx j) noexcept -> typename Derived::Scalar {
            typename Derived::Scalar sum = 0;
            for (idx m = 0; m < DA; ++m)
                sum += rA(m * DB + i) * std::conj(rA(m * DB + j));

            return sum;
        }; /* end worker */

#ifdef QPP_OPENMP
// NOLINTNEXTLINE
#pragma omp parallel for collapse(2)
#endif // QPP_OPENMP
       // column major order for speed
        for (idx j = 0; j < DB; ++j)
            for (idx i = 0; i < DB; ++i)
                result(i, j) = worker(i, j);
    }
    //************ density matrix ************//
    else if (internal::check_square_mat(rA)) // we have a density operator
    {
        // check that dims match the dimension of A
        if (!internal::check_dims_match_mat(dims, rA))
            throw exception::DimsMismatchMatrix("qpp::ptrace1()", "A/dims");

        auto worker = [&](idx i, idx j) noexcept -> typename Derived::Scalar {
            typename Derived::Scalar sum = 0;
            for (idx m = 0; m < DA; ++m)
                sum += rA(m * DB + i, m * DB + j);

            return sum;
        }; /* end worker */

#ifdef QPP_OPENMP
// NOLINTNEXTLINE
#pragma omp parallel for collapse(2)
#endif // QPP_OPENMP
       // column major order for speed
        for (idx j = 0; j < DB; ++j)
            for (idx i = 0; i < DB; ++i)
                result(i, j) = worker(i, j);
    }
    //************ Exception: not ket nor density matrix ************//
    else
        throw exception::MatrixNotSquareNorCvector("qpp::ptrace1()", "A");

    return result;
}

/**
 * \brief Partial trace
 * \see qpp::ptrace2()
 *
 * Partial trace over the first subsystem of bi-partite state vector or density
 * matrix
 *
 * \param A Eigen expression
 * \param d Subsystem dimensions
 * \return Partial trace \f$Tr_{A}(\cdot)\f$ over the first subsytem \f$A\f$ in
 * a bi-partite system \f$A\otimes B\f$, as a dynamic matrix over the same
 * scalar field as \a A
 */
template <typename Derived>
dyn_mat<typename Derived::Scalar> ptrace1(const Eigen::MatrixBase<Derived>& A,
                                          idx d = 2) {
    const dyn_mat<typename Derived::Scalar>& rA = A.derived();

    // EXCEPTION CHECKS

    // check zero size
    if (!internal::check_nonzero_size(rA))
        throw exception::ZeroSize("qpp::ptrace1()", "A");

    // check valid dims
    if (d == 0)
        throw exception::DimsInvalid("qpp::ptrace1()", "d");
    // END EXCEPTION CHECKS

    std::vector<idx> dims(2, d); // local dimensions vector

    return ptrace1(rA, dims);
}

/**
 * \brief Partial trace
 * \see qpp::ptrace1()
 *
 * Partial trace over the second subsystem of bi-partite state vector or density
 * matrix
 *
 * \param A Eigen expression
 * \param dims Dimensions of the bi-partite system
 * \return Partial trace \f$Tr_{B}(\cdot)\f$ over the second subsytem \f$B\f$ in
 * a bi-partite system \f$A\otimes B\f$, as a dynamic matrix over the same
 * scalar field as \a A
 */
template <typename Derived>
[[qpp::critical, qpp::parallel]] dyn_mat<typename Derived::Scalar>
ptrace2(const Eigen::MatrixBase<Derived>& A, const std::vector<idx>& dims) {
    const typename Eigen::MatrixBase<Derived>::EvalReturnType& rA = A.derived();

    // EXCEPTION CHECKS

    // check zero-size
    if (!internal::check_nonzero_size(rA))
        throw exception::ZeroSize("qpp::ptrace2()", "A");

    // check that dims is a valid dimension vector
    if (!internal::check_dims(dims))
        throw exception::DimsInvalid("qpp::ptrace2()", "dims");

    // check dims has only 2 elements
    if (dims.size() != 2)
        throw exception::NotBipartite("qpp::ptrace2()", "dims");
    // END EXCEPTION CHECKS

    idx DA = dims[0];
    idx DB = dims[1];

    dyn_mat<typename Derived::Scalar> result =
        dyn_mat<typename Derived::Scalar>::Zero(DA, DA);

    //************ ket ************//
    if (internal::check_cvector(rA)) // we have a ket
    {
        // check that dims match the dimension of A
        if (!internal::check_dims_match_cvect(dims, rA))
            throw exception::DimsMismatchCvector("qpp::ptrace2()", "A/dims");

        auto worker = [&](idx i, idx j) noexcept -> typename Derived::Scalar {
            typename Derived::Scalar sum = 0;
            for (idx m = 0; m < DB; ++m)
                sum += rA(i * DB + m) * std::conj(rA(j * DB + m));

            return sum;
        }; /* end worker */

#ifdef QPP_OPENMP
// NOLINTNEXTLINE
#pragma omp parallel for collapse(2)
#endif // QPP_OPENMP
       // column major order for speed
        for (idx j = 0; j < DA; ++j)
            for (idx i = 0; i < DA; ++i)
                result(i, j) = worker(i, j);
    }
    //************ density matrix ************//
    else if (internal::check_square_mat(rA)) // we have a density operator
    {
        // check that dims match the dimension of A
        if (!internal::check_dims_match_mat(dims, rA))
            throw exception::DimsMismatchMatrix("qpp::ptrace2()", "A/dims");

#ifdef QPP_OPENMP
// NOLINTNEXTLINE
#pragma omp parallel for collapse(2)
#endif // QPP_OPENMP
       // column major order for speed
        for (idx j = 0; j < DA; ++j)
            for (idx i = 0; i < DA; ++i)
                result(i, j) = trace(rA.block(i * DB, j * DB, DB, DB));
    }
    //************ Exception: not ket nor density matrix ************//
    else
        throw exception::MatrixNotSquareNorCvector("qpp::ptrace2()", "A");

    return result;
}

/**
 * \brief Partial trace
 * \see qpp::ptrace1()
 *
 * Partial trace over the second subsystem of bi-partite state vector or density
 * matrix
 *
 * \param A Eigen expression
 * \param d Subsystem dimensions
 * \return Partial trace \f$Tr_{B}(\cdot)\f$ over the second subsytem \f$B\f$ in
 * a bi-partite system \f$A\otimes B\f$, as a dynamic matrix over the same
 * scalar field as \a A
 */
template <typename Derived>
dyn_mat<typename Derived::Scalar> ptrace2(const Eigen::MatrixBase<Derived>& A,
                                          idx d = 2) {
    const dyn_mat<typename Derived::Scalar>& rA = A.derived();

    // EXCEPTION CHECKS

    // check zero size
    if (!internal::check_nonzero_size(rA))
        throw exception::ZeroSize("qpp::ptrace2()", "A");

    // check valid dims
    if (d == 0)
        throw exception::DimsInvalid("qpp::ptrace2()", "d");
    // END EXCEPTION CHECKS

    std::vector<idx> dims(2, d); // local dimensions vector

    return ptrace2(rA, dims);
}

/**
 * \brief Partial trace
 * \see qpp::ptrace1(), qpp::ptrace2()
 *
 * Partial trace of the multi-partite state vector or density matrix over the
 * list \a target of subsystems
 *
 * \param A Eigen expression
 * \param target Subsystem indexes
 * \param dims Dimensions of the multi-partite system
 * \return Partial trace \f$Tr_{subsys}(\cdot)\f$ over the subsytems \a target
 * in a multi-partite system, as a dynamic matrix over the same scalar field as
 * \a A
 */
template <typename Derived>
[[qpp::critical, qpp::parallel]] dyn_mat<typename Derived::Scalar>
ptrace(const Eigen::MatrixBase<Derived>& A, const std::vector<idx>& target,
       const std::vector<idx>& dims) {
    const typename Eigen::MatrixBase<Derived>::EvalReturnType& rA = A.derived();

    // EXCEPTION CHECKS

    // check zero-size
    if (!internal::check_nonzero_size(rA))
        throw exception::ZeroSize("qpp::ptrace()", "A");

    // check that dims is a valid dimension vector
    if (!internal::check_dims(dims))
        throw exception::DimsInvalid("qpp::ptrace()", "dims");

    // check that target is valid w.r.t. dims
    if (!internal::check_subsys_match_dims(target, dims))
        throw exception::SubsysMismatchDims("qpp::ptrace()", "dims/target");

    // check valid state and matching dimensions
    if (internal::check_cvector(rA)) {
        if (!internal::check_dims_match_cvect(dims, rA))
            throw exception::DimsMismatchCvector("qpp::ptrace()", "A/dims");
    } else if (internal::check_square_mat(rA)) {
        if (!internal::check_dims_match_mat(dims, rA))
            throw exception::DimsMismatchMatrix("qpp::ptrace()", "A/dims");
    } else
        throw exception::MatrixNotSquareNorCvector("qpp::ptrace()", "A");
    // END EXCEPTION CHECKS

    idx D = static_cast<idx>(rA.rows());
    idx n = dims.size();
    idx n_subsys = target.size();
    idx n_subsys_bar = n - n_subsys;
    idx Dsubsys = 1;
    for (idx i = 0; i < n_subsys; ++i)
        Dsubsys *= dims[target[i]];
    idx Dsubsys_bar = D / Dsubsys;

    idx Cdims[internal::maxn];
    idx Csubsys[internal::maxn];
    idx Cdimssubsys[internal::maxn];
    idx Csubsys_bar[internal::maxn];
    idx Cdimssubsys_bar[internal::maxn];

    idx Cmidxcolsubsys_bar[internal::maxn];

    std::vector<idx> subsys_bar = complement(target, n);
    std::copy(subsys_bar.begin(), subsys_bar.end(), std::begin(Csubsys_bar));

    for (idx i = 0; i < n; ++i) {
        Cdims[i] = dims[i];
    }
    for (idx i = 0; i < n_subsys; ++i) {
        Csubsys[i] = target[i];
        Cdimssubsys[i] = dims[target[i]];
    }
    for (idx i = 0; i < n_subsys_bar; ++i) {
        Cdimssubsys_bar[i] = dims[subsys_bar[i]];
    }

    dyn_mat<typename Derived::Scalar> result =
        dyn_mat<typename Derived::Scalar>(Dsubsys_bar, Dsubsys_bar);

    //************ ket ************//
    if (internal::check_cvector(rA)) // we have a ket
    {
        if (target.size() == dims.size()) {
            result(0, 0) = (adjoint(rA) * rA).value();
            return result;
        }

        if (target.empty())
            return rA * adjoint(rA);

        auto worker = [&](idx i) noexcept -> typename Derived::Scalar {
            // use static allocation for speed!
            idx Cmidxrow[internal::maxn];
            idx Cmidxcol[internal::maxn];
            idx Cmidxrowsubsys_bar[internal::maxn];
            idx Cmidxsubsys[internal::maxn];

            /* get the row multi-indexes of the complement */
            internal::n2multiidx(i, n_subsys_bar, Cdimssubsys_bar,
                                 Cmidxrowsubsys_bar);
            /* write them in the global row/col multi-indexes */
            for (idx k = 0; k < n_subsys_bar; ++k) {
                Cmidxrow[Csubsys_bar[k]] = Cmidxrowsubsys_bar[k];
                Cmidxcol[Csubsys_bar[k]] = Cmidxcolsubsys_bar[k];
            }
            typename Derived::Scalar sm = 0;
            for (idx a = 0; a < Dsubsys; ++a) {
                // get the multi-index over which we do the summation
                internal::n2multiidx(a, n_subsys, Cdimssubsys, Cmidxsubsys);
                // write it into the global row/col multi-indexes
                for (idx k = 0; k < n_subsys; ++k)
                    Cmidxrow[Csubsys[k]] = Cmidxcol[Csubsys[k]] =
                        Cmidxsubsys[k];

                // now do the sum
                sm += rA(internal::multiidx2n(Cmidxrow, n, Cdims)) *
                      std::conj(rA(internal::multiidx2n(Cmidxcol, n, Cdims)));
            }

            return sm;
        }; /* end worker */

        for (idx j = 0; j < Dsubsys_bar; ++j) // column major order for speed
        {
            // compute the column multi-indexes of the complement
            internal::n2multiidx(j, n_subsys_bar, Cdimssubsys_bar,
                                 Cmidxcolsubsys_bar);
#ifdef QPP_OPENMP
// NOLINTNEXTLINE
#pragma omp parallel for
#endif // QPP_OPENMP
            for (idx i = 0; i < Dsubsys_bar; ++i) {
                result(i, j) = worker(i);
            }
        }
    }
    //************ density matrix ************//
    else // we have a density operator
    {
        if (target.size() == dims.size()) {
            result(0, 0) = rA.trace();
            return result;
        }

        if (target.empty())
            return rA;

        auto worker = [&](idx i) noexcept -> typename Derived::Scalar {
            // use static allocation for speed!
            idx Cmidxrow[internal::maxn];
            idx Cmidxcol[internal::maxn];
            idx Cmidxrowsubsys_bar[internal::maxn];
            idx Cmidxsubsys[internal::maxn];

            /* get the row/col multi-indexes of the complement */
            internal::n2multiidx(i, n_subsys_bar, Cdimssubsys_bar,
                                 Cmidxrowsubsys_bar);
            /* write them in the global row/col multi-indexes */
            for (idx k = 0; k < n_subsys_bar; ++k) {
                Cmidxrow[Csubsys_bar[k]] = Cmidxrowsubsys_bar[k];
                Cmidxcol[Csubsys_bar[k]] = Cmidxcolsubsys_bar[k];
            }
            typename Derived::Scalar sm = 0;
            for (idx a = 0; a < Dsubsys; ++a) {
                // get the multi-index over which we do the summation
                internal::n2multiidx(a, n_subsys, Cdimssubsys, Cmidxsubsys);
                // write it into the global row/col multi-indexes
                for (idx k = 0; k < n_subsys; ++k)
                    Cmidxrow[Csubsys[k]] = Cmidxcol[Csubsys[k]] =
                        Cmidxsubsys[k];

                // now do the sum
                sm += rA(internal::multiidx2n(Cmidxrow, n, Cdims),
                         internal::multiidx2n(Cmidxcol, n, Cdims));
            }

            return sm;
        }; /* end worker */

        for (idx j = 0; j < Dsubsys_bar; ++j) // column major order for speed
        {
            // compute the column multi-indexes of the complement
            internal::n2multiidx(j, n_subsys_bar, Cdimssubsys_bar,
                                 Cmidxcolsubsys_bar);
#ifdef QPP_OPENMP
// NOLINTNEXTLINE
#pragma omp parallel for
#endif // QPP_OPENMP
            for (idx i = 0; i < Dsubsys_bar; ++i) {
                result(i, j) = worker(i);
            }
        }
    }

    return result;
}

/**
 * \brief Partial trace
 * \see qpp::ptrace1(), qpp::ptrace2()
 *
 * Partial trace of the multi-partite state vector or density matrix over the
 * list \a target of subsystems
 *
 * \param A Eigen expression
 * \param target Subsystem indexes
 * \param d Subsystem dimensions
 * \return Partial trace \f$Tr_{subsys}(\cdot)\f$ over the subsytems \a target
 * in a multi-partite system, as a dynamic matrix over the same scalar field as
 * \a A
 */
template <typename Derived>
dyn_mat<typename Derived::Scalar> ptrace(const Eigen::MatrixBase<Derived>& A,
                                         const std::vector<idx>& target,
                                         idx d = 2) {
    const typename Eigen::MatrixBase<Derived>::EvalReturnType& rA = A.derived();

    // EXCEPTION CHECKS

    // check zero size
    if (!internal::check_nonzero_size(rA))
        throw exception::ZeroSize("qpp::ptrace()", "A");

    // check valid dims
    if (d < 2)
        throw exception::DimsInvalid("qpp::ptrace()", "d");
    // END EXCEPTION CHECKS

    idx n = internal::get_num_subsys(static_cast<idx>(rA.rows()), d);
    std::vector<idx> dims(n, d); // local dimensions vector

    return ptrace(rA, target, dims);
}

/**
 * \brief Partial transpose
 *
 * Partial transpose of the multi-partite state vector or density matrix over
 * the list \a target of subsystems
 *
 * \param A Eigen expression
 * \param target Subsystem indexes
 * \param dims Dimensions of the multi-partite system
 * \return Partial transpose \f$(\cdot)^{T_{subsys}}\f$ over the subsytems
 * \a target in a multi-partite system, as a dynamic matrix over the same scalar
 * field as \a A
 */
template <typename Derived>
dyn_mat<typename Derived::Scalar> [[qpp::critical, qpp::parallel]] ptranspose(
    const Eigen::MatrixBase<Derived>& A, const std::vector<idx>& target,
    const std::vector<idx>& dims) {
    const typename Eigen::MatrixBase<Derived>::EvalReturnType& rA = A.derived();

    // EXCEPTION CHECKS

    // check zero-size
    if (!internal::check_nonzero_size(rA))
        throw exception::ZeroSize("qpp::ptranspose()", "A");

    // check that dims is a valid dimension vector
    if (!internal::check_dims(dims))
        throw exception::DimsInvalid("qpp::ptranspose()", "dims");

    // check that target is valid w.r.t. dims
    if (!internal::check_subsys_match_dims(target, dims))
        throw exception::SubsysMismatchDims("qpp::ptranspose()", "dims/target");

    // check valid state and matching dimensions
    if (internal::check_cvector(rA)) {
        if (!internal::check_dims_match_cvect(dims, rA))
            throw exception::DimsMismatchCvector("qpp::ptranspose()", "A/dims");
    } else if (internal::check_square_mat(rA)) {
        if (!internal::check_dims_match_mat(dims, rA))
            throw exception::DimsMismatchMatrix("qpp::ptranspose()", "A/dims");
    } else
        throw exception::MatrixNotSquareNorCvector("qpp::ptranspose()", "A");
    // END EXCEPTION CHECKS

    idx D = static_cast<idx>(rA.rows());
    idx n = dims.size();
    idx n_subsys = target.size();
    idx Cdims[internal::maxn];
    idx Cmidxcol[internal::maxn];
    idx Csubsys[internal::maxn];

    // copy dims in Cdims and target in Csubsys
    for (idx i = 0; i < n; ++i)
        Cdims[i] = dims[i];
    for (idx i = 0; i < n_subsys; ++i)
        Csubsys[i] = target[i];

    dyn_mat<typename Derived::Scalar> result(D, D);

    //************ ket ************//
    if (internal::check_cvector(rA)) // we have a ket
    {
        if (target.size() == dims.size())
            return (rA * adjoint(rA)).transpose();

        if (target.empty())
            return rA * adjoint(rA);

        auto worker = [&](idx i) noexcept -> typename Derived::Scalar {
            // use static allocation for speed!
            idx midxcoltmp[internal::maxn];
            idx midxrow[internal::maxn];

            for (idx k = 0; k < n; ++k)
                midxcoltmp[k] = Cmidxcol[k];

            /* compute the row multi-index */
            internal::n2multiidx(i, n, Cdims, midxrow);

            for (idx k = 0; k < n_subsys; ++k)
                std::swap(midxcoltmp[Csubsys[k]], midxrow[Csubsys[k]]);

            /* writes the result */
            return rA(internal::multiidx2n(midxrow, n, Cdims)) *
                   std::conj(rA(internal::multiidx2n(midxcoltmp, n, Cdims)));
        }; /* end worker */

        for (idx j = 0; j < D; ++j) {
            // compute the column multi-index
            internal::n2multiidx(j, n, Cdims, Cmidxcol);

#ifdef QPP_OPENMP
// NOLINTNEXTLINE
#pragma omp parallel for
#endif // QPP_OPENMP
            for (idx i = 0; i < D; ++i)
                result(i, j) = worker(i);
        }
    }
    //************ density matrix ************//
    else // we have a density operator
    {
        if (target.size() == dims.size())
            return rA.transpose();

        if (target.empty())
            return rA;

        auto worker = [&](idx i) noexcept -> typename Derived::Scalar {
            // use static allocation for speed!
            idx midxcoltmp[internal::maxn];
            idx midxrow[internal::maxn];

            for (idx k = 0; k < n; ++k)
                midxcoltmp[k] = Cmidxcol[k];

            /* compute the row multi-index */
            internal::n2multiidx(i, n, Cdims, midxrow);

            for (idx k = 0; k < n_subsys; ++k)
                std::swap(midxcoltmp[Csubsys[k]], midxrow[Csubsys[k]]);

            /* writes the result */
            return rA(internal::multiidx2n(midxrow, n, Cdims),
                      internal::multiidx2n(midxcoltmp, n, Cdims));
        }; /* end worker */

        for (idx j = 0; j < D; ++j) {
            // compute the column multi-index
            internal::n2multiidx(j, n, Cdims, Cmidxcol);

#ifdef QPP_OPENMP
// NOLINTNEXTLINE
#pragma omp parallel for
#endif // QPP_OPENMP
            for (idx i = 0; i < D; ++i)
                result(i, j) = worker(i);
        }
    }

    return result;
}

/**
 * \brief Partial transpose
 *
 * Partial transpose of the multi-partite state vector or density matrix over
 * the list \a target of subsystems
 *
 * \param A Eigen expression
 * \param target Subsystem indexes
 * \param d Subsystem dimensions
 * \return Partial transpose \f$(\cdot)^{T_{subsys}}\f$ over the subsytems
 * \a target in a multi-partite system, as a dynamic matrix over the same scalar
 * field as \a A
 */
template <typename Derived>
dyn_mat<typename Derived::Scalar>
ptranspose(const Eigen::MatrixBase<Derived>& A, const std::vector<idx>& target,
           idx d = 2) {
    const typename Eigen::MatrixBase<Derived>::EvalReturnType& rA = A.derived();

    // EXCEPTION CHECKS

    // check zero size
    if (!internal::check_nonzero_size(rA))
        throw exception::ZeroSize("qpp::ptranspose()", "A");

    // check valid dims
    if (d < 2)
        throw exception::DimsInvalid("qpp::ptranspose()", "d");
    // END EXCEPTION CHECKS

    idx n = internal::get_num_subsys(static_cast<idx>(rA.rows()), d);
    std::vector<idx> dims(n, d); // local dimensions vector

    return ptranspose(rA, target, dims);
}

/**
 * \brief Subsystem permutation
 *
 * Permutes the subsystems of a state vector or density matrix. The qubit
 * \a perm[\a i] is permuted to the location \a i.
 *
 * \param A Eigen expression
 * \param perm Permutation
 * \param dims Dimensions of the multi-partite system
 * \return Permuted system, as a dynamic matrix over the same scalar field as
 * \a A
 */
template <typename Derived>
[[qpp::critical, qpp::parallel]] dyn_mat<typename Derived::Scalar>
syspermute(const Eigen::MatrixBase<Derived>& A, const std::vector<idx>& perm,
           const std::vector<idx>& dims) {
    const typename Eigen::MatrixBase<Derived>::EvalReturnType& rA = A.derived();

    // EXCEPTION CHECKS

    // check zero-size
    if (!internal::check_nonzero_size(rA))
        throw exception::ZeroSize("qpp::syspermute()", "A");

    // check that dims is a valid dimension vector
    if (!internal::check_dims(dims))
        throw exception::DimsInvalid("qpp::syspermute()", "dims");

    // check that we have a valid permutation
    if (!internal::check_perm(perm))
        throw exception::PermInvalid("qpp::syspermute()", "perm");

    // check that permutation match dimensions
    if (perm.size() != dims.size())
        throw exception::PermMismatchDims("qpp::syspermute()", "dims/perm");

    // check valid state and matching dimensions
    if (internal::check_cvector(rA)) {
        if (!internal::check_dims_match_cvect(dims, rA))
            throw exception::DimsMismatchCvector("qpp::syspermute()", "A/dims");
    } else if (internal::check_square_mat(rA)) {
        if (!internal::check_dims_match_mat(dims, rA))
            throw exception::DimsMismatchMatrix("qpp::syspermute()", "A/dims");
    } else
        throw exception::MatrixNotSquareNorCvector("qpp::syspermute()", "A");
    // END EXCEPTION CHECKS

    idx D = static_cast<idx>(rA.rows());
    idx n = dims.size();

    dyn_mat<typename Derived::Scalar> result;

    //************ ket ************//
    if (internal::check_cvector(rA)) // we have a column vector
    {
        idx Cdims[internal::maxn];
        idx Cperm[internal::maxn];

        // copy dims in Cdims and perm in Cperm
        for (idx i = 0; i < n; ++i) {
            Cdims[i] = dims[i];
            Cperm[i] = perm[i];
        }
        result.resize(D, 1);

        auto worker = [&Cdims, &Cperm, n](idx i) noexcept -> idx {
            // use static allocation for speed,
            // allocate twice the size for matrices reshaped as vectors
            idx midx[internal::maxn];
            idx midxtmp[internal::maxn];
            idx permdims[internal::maxn];

            /* compute the multi-index */
            internal::n2multiidx(i, n, Cdims, midx);

            for (idx k = 0; k < n; ++k) {
                permdims[k] = Cdims[Cperm[k]]; // permuted dimensions
                midxtmp[k] = midx[Cperm[k]];   // permuted multi-indexes
            }

            return internal::multiidx2n(midxtmp, n, permdims);
        }; /* end worker */

#ifdef QPP_OPENMP
// NOLINTNEXTLINE
#pragma omp parallel for
#endif // QPP_OPENMP
        for (idx i = 0; i < D; ++i)
            result(worker(i)) = rA(i);
    }
    //************ density matrix ************//
    else // we have a density operator
    {
        idx Cdims[2 * internal::maxn];
        idx Cperm[2 * internal::maxn];

        // copy dims in Cdims and perm in Cperm
        for (idx i = 0; i < n; ++i) {
            Cdims[i] = dims[i];
            Cdims[i + n] = dims[i];
            Cperm[i] = perm[i];
            Cperm[i + n] = perm[i] + n;
        }
        result.resize(D * D, 1);
        // map A to a column vector
        dyn_mat<typename Derived::Scalar> vectA =
            Eigen::Map<dyn_mat<typename Derived::Scalar>>(
                const_cast<typename Derived::Scalar*>(rA.data()), D * D, 1);

        auto worker = [&Cdims, &Cperm, n](idx i) noexcept -> idx {
            // use static allocation for speed,
            // allocate twice the size for matrices reshaped as vectors
            idx midx[2 * internal::maxn];
            idx midxtmp[2 * internal::maxn];
            idx permdims[2 * internal::maxn];

            /* compute the multi-index */
            internal::n2multiidx(i, 2 * n, Cdims, midx);

            for (idx k = 0; k < 2 * n; ++k) {
                permdims[k] = Cdims[Cperm[k]]; // permuted dimensions
                midxtmp[k] = midx[Cperm[k]];   // permuted multi-indexes
            }

            return internal::multiidx2n(midxtmp, 2 * n, permdims);
        }; /* end worker */

#ifdef QPP_OPENMP
// NOLINTNEXTLINE
#pragma omp parallel for
#endif // QPP_OPENMP
        for (idx i = 0; i < D * D; ++i)
            result(worker(i)) = rA(i);

        result = reshape(result, D, D);
    }

    return result;
}

/**
 * \brief Subsystem permutation
 *
 * Permutes the subsystems of a state vector or density matrix. The qubit
 * \a perm[\a i] is permuted to the location \a i.
 *
 * \param A Eigen expression
 * \param perm Permutation
 * \param d Subsystem dimensions
 * \return Permuted system, as a dynamic matrix over the same scalar field as
 * \a A
 */
template <typename Derived>
dyn_mat<typename Derived::Scalar>
syspermute(const Eigen::MatrixBase<Derived>& A, const std::vector<idx>& perm,
           idx d = 2) {
    const typename Eigen::MatrixBase<Derived>::EvalReturnType& rA = A.derived();

    // EXCEPTION CHECKS

    // check zero size
    if (!internal::check_nonzero_size(rA))
        throw exception::ZeroSize("qpp::syspermute()", "A");

    // check valid dims
    if (d < 2)
        throw exception::DimsInvalid("qpp::syspermute()", "d");
    // END EXCEPTION CHECKS

    idx n = internal::get_num_subsys(static_cast<idx>(rA.rows()), d);
    std::vector<idx> dims(n, d); // local dimensions vector

    return syspermute(rA, perm, dims);
}

/**
 * \brief Applies the controlled-gate \a A to the part \a target of the
 * multi-partite state vector or density matrix \a state
 * \see qpp::Gates::CTRL()
 *
 * \note The dimension of the gate \a A must match the dimension of \a target.
 * Also, all control subsystems in \a ctrl must have the same dimension.
 *
 * \param state Eigen expression
 * \param A Eigen expression
 * \param ctrl Control subsystem indexes
 * \param target Subsystem indexes where the gate \a A is applied
 * \param dims Dimensions of the multi-partite system
 * \param shift Optional, performs the control as if the \a ctrl qudits were
 * \f$X\f$-incremented component-wise by \a shift (in the order given by
 * \a ctrl). For example, for two qutrits (D=3), applying a control gate on the
 * initial state \f$|00\rangle\f$ with first qutrit as control, second qutrit as
 * target, and \a shift = {1}, yields the state \f$|01\rangle\f$. If present,
 * the size of \a shift must be the same as the size of \a ctrl.
 * \return CTRL-A gate applied to the part \a target of \a state
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] expr_t<Derived1>
applyCTRL(const Eigen::MatrixBase<Derived1>& state,
          const Eigen::MatrixBase<Derived2>& A, const std::vector<idx>& ctrl,
          const std::vector<idx>& target, const std::vector<idx>& dims,
          std::optional<std::vector<idx>> shift = std::nullopt) {
    const expr_t<Derived1>& rstate = state.derived();
    const dyn_mat<typename Derived2::Scalar>& rA = A.derived();

    // EXCEPTION CHECKS

    // check types
    if (!std::is_same<typename Derived1::Scalar,
                      typename Derived2::Scalar>::value)
        throw exception::TypeMismatch("qpp::applyCTRL()", "A/state");

    // check zero sizes
    if (!internal::check_nonzero_size(rA))
        throw exception::ZeroSize("qpp::applyCTRL()", "A");

    // check zero sizes
    if (!internal::check_nonzero_size(rstate))
        throw exception::ZeroSize("qpp::applyCTRL()", "state");

    // check zero sizes
    if (!internal::check_nonzero_size(ctrl))
        throw exception::ZeroSize("qpp::applyCTRL()", "ctrl");

    // check zero sizes
    if (!internal::check_nonzero_size(target))
        throw exception::ZeroSize("qpp::applyCTRL()", "target");

    // check square matrix for the gate
    if (!internal::check_square_mat(rA))
        throw exception::MatrixNotSquare("qpp::applyCTRL()", "A");

    // check valid state and matching dimensions
    if (internal::check_cvector(rstate)) {
        if (!internal::check_dims_match_cvect(dims, rstate))
            throw exception::DimsMismatchCvector("qpp::applyCTRL()",
                                                 "dims/state");
    } else if (internal::check_square_mat(rstate)) {
        if (!internal::check_dims_match_mat(dims, rstate))
            throw exception::DimsMismatchMatrix("qpp::applyCTRL()",
                                                "dims/state");
    } else
        throw exception::MatrixNotSquareNorCvector("qpp::applyCTRL()", "state");

    // check that ctrl subsystem is valid w.r.t. dims
    if (!internal::check_subsys_match_dims(ctrl, dims))
        throw exception::SubsysMismatchDims("qpp::applyCTRL()", "ctrl/dims");

    // check that all control subsystems have the same dimension
    idx d = dims[ctrl[0]];
    for (idx i = 1; i < static_cast<idx>(ctrl.size()); ++i)
        if (dims[ctrl[i]] != d)
            throw exception::DimsNotEqual("qpp::applyCTRL()", "ctrl");

    // check that dimension is valid
    if (!internal::check_dims(dims))
        throw exception::DimsInvalid("qpp::applyCTRL()", "dims");

    // check that target is valid w.r.t. dims
    if (!internal::check_subsys_match_dims(target, dims))
        throw exception::SubsysMismatchDims("qpp::applyCTRL()", "dims/target");

    // check that ctrl and target don't share common elements
    for (idx elem_ctrl : ctrl)
        for (idx elem_target : target)
            if (elem_ctrl == elem_target)
                throw exception::OutOfRange("qpp::applyCTRL()", "ctrl/target");

    // check that gate matches the dimensions of the target
    std::vector<idx> target_dims(target.size());
    for (idx i = 0; i < static_cast<idx>(target.size()); ++i)
        target_dims[i] = dims[target[i]];
    if (!internal::check_dims_match_mat(target_dims, rA))
        throw exception::MatrixMismatchSubsys("qpp::applyCTRL()", "A/target");

    std::vector<idx> ctrlgate = ctrl; // ctrl + gate subsystem vector
    ctrlgate.insert(ctrlgate.end(), target.begin(), target.end());
    // FIXME if needed
    std::sort(ctrlgate.begin(), ctrlgate.end());

    // check that ctrl + gate subsystem is valid
    // with respect to local dimensions
    if (!internal::check_subsys_match_dims(ctrlgate, dims))
        throw exception::SubsysMismatchDims("qpp::applyCTRL()",
                                            "dims/ctrl/target");

    // check shift
    if (shift.has_value() && (shift.value().size() != ctrl.size()))
        throw exception::SizeMismatch("qpp::applyCTRL()", "ctrl/shift");
    if (shift.has_value())
        for (idx& elem : shift.value()) {
            if (elem >= d)
                throw exception::OutOfRange("qpp::applyCTRL()", "shift");
            elem = d - elem;
        }

    // END EXCEPTION CHECKS

    if (!shift.has_value())
        shift = std::vector<idx>(ctrl.size(), 0);

    // construct the table of A^k
    std::vector<dyn_mat<typename Derived1::Scalar>> Ak;
    for (idx k = 0; k < d; ++k) {
        Ak.emplace_back(powm(rA, k));
    }

    idx D = static_cast<idx>(rstate.rows()); // total dimension

    auto applyCTRL_ket =
        [d, D, &Ak = std::as_const(Ak), &ctrl = std::as_const(ctrl),
         &target = std::as_const(target), &dims = std::as_const(dims),
         &shift = std::as_const(shift)](
            const dyn_col_vect<typename Derived1::Scalar>& state) {
            dyn_col_vect<typename Derived1::Scalar> result =
                dyn_col_vect<typename Derived1::Scalar>::Zero(D);
#ifdef QPP_OPENMP
// NOLINTNEXTLINE
#pragma omp parallel for
#endif // QPP_OPENMP
            for (idx r = 0; r < d; ++r) {
                // copy shift
                std::vector<idx> ctrl_shift = shift.value();
                // increment ctrl_shift by r (mod d)
                std::transform(ctrl_shift.begin(), ctrl_shift.end(),
                               ctrl_shift.begin(),
                               [r, d](idx elem) { return (elem + r) % d; });
                // project on ctrl_shift
                dyn_col_vect<typename Derived1::Scalar> chopped_psi =
                    internal::project_ket_on_dits(state, ctrl_shift, ctrl, dims,
                                                  D);
                // project on complement
                dyn_col_vect<typename Derived1::Scalar> chopped_psi_bar =
                    state / d - chopped_psi;
                // apply gates on target if non-zero
                if (chopped_psi !=
                    dyn_col_vect<typename Derived1::Scalar>::Zero(D)) {
                    chopped_psi = apply(chopped_psi, Ak[r], target, dims);
                }
#ifdef QPP_OPENMP
// NOLINTNEXTLINE
#pragma omp critical
#endif // QPP_OPENMP
                { result += (chopped_psi + chopped_psi_bar); }
            }

            return result;
        }; // end applyCTRL_ket

    //************ ket ************//
    if (internal::check_cvector(rstate)) // we have a ket
    {
        // check that dims match state vector
        if (!internal::check_dims_match_cvect(dims, rstate))
            throw exception::DimsMismatchCvector("qpp::applyCTRL()",
                                                 "dims/state");
        if (D == 1)
            return rstate;

        return applyCTRL_ket(rstate);
    } // end ket

    //************ density matrix ************//
    else // we have a density operator
    {
        // check that dims match state matrix
        if (!internal::check_dims_match_mat(dims, rstate))
            throw exception::DimsMismatchMatrix("qpp::applyCTRL()",
                                                "dims/state");

        if (D == 1)
            return rstate;

        dyn_mat<typename Derived1::Scalar> result =
            dyn_mat<typename Derived1::Scalar>::Zero(D, D);

#ifdef QPP_OPENMP
// NOLINTNEXTLINE
#pragma omp parallel for
#endif // QPP_OPENMP
       // decompose rho as \sum |psi_i><phi_i|
        for (idx i = 0; i < D; ++i) {
            std::vector<idx> midx_i = n2multiidx(i, dims);
            dyn_col_vect<typename Derived1::Scalar> psi_i =
                dyn_col_vect<typename Derived1::Scalar>::Zero(D);
            psi_i(i) = 1;
            dyn_row_vect<typename Derived1::Scalar> phi_i_bra = rstate.row(i);
            dyn_col_vect<typename Derived1::Scalar> phi_i_ket =
                adjoint(phi_i_bra);

            // compute CTRL(|psi_i><phi_i|)CTRL^\dagger
            psi_i = applyCTRL_ket(psi_i);
            phi_i_bra = adjoint(applyCTRL_ket(phi_i_ket));
// add to result
#ifdef QPP_OPENMP
// NOLINTNEXTLINE
#pragma omp critical
#endif // QPP_OPENMP
            { result += kron(psi_i, phi_i_bra); }
        } // end for(i)

        return result;
    } // end density operator
}

/**
 * \brief Applies the controlled-gate \a A to the part \a target of the
 * multi-partite state vector or density matrix \a state
 * \see qpp::Gates::CTRL()
 *
 * \note The dimension of the gate \a A must match the dimension of \a
 * target. Also, all control subsystems in \a ctrl must have the same
 * dimension.
 *
 * \param state Eigen expression
 * \param A Eigen expression
 * \param ctrl Control subsystem indexes
 * \param target Subsystem indexes where the gate \a A is applied
 * \param d Subsystem dimensions
 * \param shift Optional, performs the control as if the \a ctrl qudits were
 * \f$X\f$-incremented component-wise by \a shift. For example, for two qutrits
 * (D=3), applying a control gate on the initial state \f$|00\rangle\f$ with
 * first qutrit as control, second qutrit as target, and \a shift = {1}, yields
 * the state \f$|01\rangle\f$. If present, the size of \a shift must be the same
 * as the size of \a ctrl.
 * \return CTRL-A gate applied to the part \a target of \a state
 */
template <typename Derived1, typename Derived2>
expr_t<Derived1>
applyCTRL(const Eigen::MatrixBase<Derived1>& state,
          const Eigen::MatrixBase<Derived2>& A, const std::vector<idx>& ctrl,
          const std::vector<idx>& target, idx d = 2,
          std::optional<std::vector<idx>> shift = std::nullopt) {
    const expr_t<Derived1>& rstate = state.derived();
    const dyn_mat<typename Derived1::Scalar>& rA = A.derived();

    // EXCEPTION CHECKS

    // check zero size
    if (!internal::check_nonzero_size(rstate))
        throw exception::ZeroSize("qpp::applyCTRL()", "state");

    // check valid dims
    if (d < 2)
        throw exception::DimsInvalid("qpp::applyCTRL()", "d");
    // END EXCEPTION CHECKS

    idx n = internal::get_num_subsys(static_cast<idx>(rstate.rows()), d);
    std::vector<idx> dims(n, d); // local dimensions vector

    return applyCTRL(rstate, rA, ctrl, target, dims, shift);
}

/**
 * \brief Applies the single qudit controlled-gate \a A with multiple
 * control qudits listed in \a ctrl on every qudit listed in \a target,
 * i.e., CTRL-CTRL-...-CTRL-A-A-...-A.
 * \see qpp::Gates::CTRL()
 *
 * \note The dimension of the gate \a A must match the dimension of every qudit
 * in \a target. All control subsystems in \a ctrl must have the same dimension.
 *
 * \param state Eigen expression
 * \param A Eigen expression, single qudit quantum gate
 * \param ctrl Control subsystem indexes
 * \param target Target qudit indexes; the gate \a A is applied on every
 * one of them depending on the values of the control qudits
 * \param dims Dimensions of the multi-partite system
 * \param shift Optional, performs the control as if the \a ctrl qudits were
 * \f$X\f$-incremented component-wise by \a shift (in the order given by
 * \a ctrl). For example, for two qutrits (D=3), applying a control gate on the
 * initial state \f$|00\rangle\f$ with first qutrit as control, second qutrit as
 * target, and \a shift = {1}, yields the state \f$|01\rangle\f$. If present,
 * the size of \a shift must be the same as the size of \a ctrl.
 * \return CTRL-CTRL-...-CTRL-A-A-...-A gate applied to every \a target qudit
 * in \a state
 */
template <typename Derived1, typename Derived2>
[[qpp::critical, qpp::parallel]] expr_t<Derived1>
applyCTRL_fan(const Eigen::MatrixBase<Derived1>& state,
              const Eigen::MatrixBase<Derived2>& A,
              const std::vector<idx>& ctrl, const std::vector<idx>& target,
              const std::vector<idx>& dims,
              std::optional<std::vector<idx>> shift = std::nullopt) {
    const expr_t<Derived1>& rstate = state.derived();
    const dyn_mat<typename Derived2::Scalar>& rA = A.derived();

    // EXCEPTION CHECKS

    // check types
    if (!std::is_same<typename Derived1::Scalar,
                      typename Derived2::Scalar>::value)
        throw exception::TypeMismatch("qpp::applyCTRL_fan()", "A/state");

    // check zero sizes
    if (!internal::check_nonzero_size(rA))
        throw exception::ZeroSize("qpp::applyCTRL_fan()", "A");

    // check zero sizes
    if (!internal::check_nonzero_size(rstate))
        throw exception::ZeroSize("qpp::applyCTRL_fan()", "state");

    // check zero sizes
    if (!internal::check_nonzero_size(ctrl))
        throw exception::ZeroSize("qpp::applyCTRL_fan()", "ctrl");

    // check zero sizes
    if (!internal::check_nonzero_size(target))
        throw exception::ZeroSize("qpp::applyCTRL_fan()", "target");

    // check square matrix for the gate
    if (!internal::check_square_mat(rA))
        throw exception::MatrixNotSquare("qpp::applyCTRL_fan()", "A");

    // check valid state and matching dimensions
    if (internal::check_cvector(rstate)) {
        if (!internal::check_dims_match_cvect(dims, rstate))
            throw exception::DimsMismatchCvector("qpp::applyCTRL_fan()",
                                                 "dims/state");
    } else if (internal::check_square_mat(rstate)) {
        if (!internal::check_dims_match_mat(dims, rstate))
            throw exception::DimsMismatchMatrix("qpp::applyCTRL_fan()",
                                                "dims/state");
    } else
        throw exception::MatrixNotSquareNorCvector("qpp::applyCTRL_fan()",
                                                   "state");

    // check that ctrl subsystem is valid w.r.t. dims
    if (!internal::check_subsys_match_dims(ctrl, dims))
        throw exception::SubsysMismatchDims("qpp::applyCTRL_fan()",
                                            "ctrl/dims");

    // check that all control subsystems have the same dimension
    idx d = dims[ctrl[0]];
    for (idx i = 1; i < static_cast<idx>(ctrl.size()); ++i)
        if (dims[ctrl[i]] != d)
            throw exception::DimsNotEqual("qpp::applyCTRL_fan()", "ctrl");

    // check that dimension is valid
    if (!internal::check_dims(dims))
        throw exception::DimsInvalid("qpp::applyCTRL_fan()", "dims");

    // check that target is valid w.r.t. dims
    if (!internal::check_subsys_match_dims(target, dims))
        throw exception::SubsysMismatchDims("qpp::applyCTRL_fan()",
                                            "dims/target");

    // check that all target subsystems have the same dimension
    for (idx i = 1; i < static_cast<idx>(target.size()); ++i)
        if (dims[target[i]] != d)
            throw exception::DimsNotEqual("qpp::applyCTRL_fan()", "target");

    // check that ctrl and target don't share common elements
    for (idx elem_ctrl : ctrl)
        for (idx elem_target : target)
            if (elem_ctrl == elem_target)
                throw exception::OutOfRange("qpp::applyCTRL_fan()",
                                            "ctrl/target");

    // check that gate matches the dimensions of the target's qudits
    if (rA.rows() != static_cast<Eigen::Index>(d))
        throw exception::MatrixMismatchSubsys("qpp::applyCTRL_fan()",
                                              "A/target");

    std::vector<idx> ctrlgate = ctrl; // ctrl + gate subsystem vector
    ctrlgate.insert(ctrlgate.end(), target.begin(), target.end());
    // FIXME if needed
    std::sort(ctrlgate.begin(), ctrlgate.end());

    // check that ctrl + gate subsystem is valid
    // with respect to local dimensions
    if (!internal::check_subsys_match_dims(ctrlgate, dims))
        throw exception::SubsysMismatchDims("qpp::applyCTRL_fan()",
                                            "dims/ctrl/target");

    // check shift
    if (shift.has_value() && (shift.value().size() != ctrl.size()))
        throw exception::SizeMismatch("qpp::applyCTRL_fan()", "ctrl/shift");
    if (shift.has_value())
        for (idx& elem : shift.value()) {
            if (elem >= d)
                throw exception::OutOfRange("qpp::applyCTRL_fan()", "shift");
            elem = d - elem;
        }

    // END EXCEPTION CHECKS

    if (!shift.has_value())
        shift = std::vector<idx>(ctrl.size(), 0);

    // construct the table of A^k
    std::vector<dyn_mat<typename Derived1::Scalar>> Ak;
    for (idx k = 0; k < d; ++k) {
        Ak.emplace_back(powm(rA, k));
    }

    idx D = static_cast<idx>(rstate.rows()); // total dimension

    auto applyCTRL_fan_ket =
        [d, D, &Ak = std::as_const(Ak), &ctrl = std::as_const(ctrl),
         &target = std::as_const(target), &dims = std::as_const(dims),
         &shift = std::as_const(shift)](
            const dyn_col_vect<typename Derived1::Scalar>& state) {
            dyn_col_vect<typename Derived1::Scalar> result =
                dyn_col_vect<typename Derived1::Scalar>::Zero(D);
#ifdef QPP_OPENMP
// NOLINTNEXTLINE
#pragma omp parallel for
#endif // QPP_OPENMP
            for (idx r = 0; r < d; ++r) {
                // copy shift
                std::vector<idx> ctrl_shift = shift.value();
                // increment ctrl_shift by r (mod d)
                std::transform(ctrl_shift.begin(), ctrl_shift.end(),
                               ctrl_shift.begin(),
                               [r, d](idx elem) { return (elem + r) % d; });
                // project on ctrl_shift
                dyn_col_vect<typename Derived1::Scalar> chopped_psi =
                    internal::project_ket_on_dits(state, ctrl_shift, ctrl, dims,
                                                  D);
                // project on complement
                dyn_col_vect<typename Derived1::Scalar> chopped_psi_bar =
                    state / d - chopped_psi;
                // apply gates on target if non-zero
                if (chopped_psi !=
                    dyn_col_vect<typename Derived1::Scalar>::Zero(D)) {
                    for (idx elem : target) {
                        chopped_psi = apply(chopped_psi, Ak[r], {elem}, dims);
                    }
                }
#ifdef QPP_OPENMP
// NOLINTNEXTLINE
#pragma omp critical
#endif // QPP_OPENMP
                { result += (chopped_psi + chopped_psi_bar); }
            }
            return result;
        }; // end applyCTRL_fan_ket

    //************ ket ************//
    if (internal::check_cvector(rstate)) // we have a ket
    {
        // check that dims match state vector
        if (!internal::check_dims_match_cvect(dims, rstate))
            throw exception::DimsMismatchCvector("qpp::applyCTRL()",
                                                 "dims/state");
        if (D == 1)
            return rstate;

        return applyCTRL_fan_ket(rstate);
    } // end ket

    //************ density matrix ************//
    else // we have a density operator
    {
        // check that dims match state matrix
        if (!internal::check_dims_match_mat(dims, rstate))
            throw exception::DimsMismatchMatrix("qpp::applyCTRL()",
                                                "dims/state");

        if (D == 1)
            return rstate;

        dyn_mat<typename Derived1::Scalar> result =
            dyn_mat<typename Derived1::Scalar>::Zero(D, D);

#ifdef QPP_OPENMP
// NOLINTNEXTLINE
#pragma omp parallel for
#endif // QPP_OPENMP
       // decompose rho as \sum |psi_i><phi_i| and applyCTRL_fan_ket()
       // individually
        for (idx i = 0; i < D; ++i) {
            std::vector<idx> midx_i = n2multiidx(i, dims);
            dyn_col_vect<typename Derived1::Scalar> psi_i =
                dyn_col_vect<typename Derived1::Scalar>::Zero(D);
            psi_i(i) = 1;
            dyn_row_vect<typename Derived1::Scalar> phi_i_bra = rstate.row(i);
            dyn_col_vect<typename Derived1::Scalar> phi_i_ket =
                adjoint(phi_i_bra);

            // compute CTRL_fan(|psi_i><phi_i|)CTRL_fan^\dagger
            psi_i = applyCTRL_fan_ket(psi_i);
            phi_i_bra = adjoint(applyCTRL_fan_ket(phi_i_ket));
// add to result
#ifdef QPP_OPENMP
// NOLINTNEXTLINE
#pragma omp critical
#endif // QPP_OPENMP
            { result += kron(psi_i, phi_i_bra); }
        } // end for(i)

        return result;
    } // end density operator
}

/**
 * \brief Applies the single qudit controlled-gate \a A with multiple
 * control qudits listed in \a ctrl on every qudit listed in \a target,
 * i.e., CTRL-CTRL-...-CTRL-A-A-...-A.
 * \see qpp::Gates::CTRL()
 *
 * \note The dimension of the gate \a A must match the dimension of every qudit
 * in \a target. All control subsystems in \a ctrl must have the same dimension.
 *
 * \param state Eigen expression
 * \param A Eigen expression, single qudit quantum gate
 * \param ctrl Control subsystem indexes
 * \param target Target qudit indexes; the gate \a A is applied on every
 * one of them depending on the values of the control qudits
 * \param d Subsystem dimensions
 * \param shift Optional, performs the control as if the \a ctrl qudits were
 * \f$X\f$-incremented component-wise by \a shift (in the order given by
 * \a ctrl). For example, for two qutrits (D=3), applying a control gate on the
 * initial state \f$|00\rangle\f$ with first qutrit as control, second qutrit as
 * target, and \a shift = {1}, yields the state \f$|01\rangle\f$. If present,
 * the size of \a shift must be the same as the size of \a ctrl.
 * \return CTRL-CTRL-...-CTRL-A-A-...-A gate applied to every \a target qudit
 * in \a state
 */
template <typename Derived1, typename Derived2>
expr_t<Derived1>
applyCTRL_fan(const Eigen::MatrixBase<Derived1>& state,
              const Eigen::MatrixBase<Derived2>& A,
              const std::vector<idx>& ctrl, const std::vector<idx>& target,
              idx d = 2, std::optional<std::vector<idx>> shift = std::nullopt) {
    const expr_t<Derived1>& rstate = state.derived();
    const dyn_mat<typename Derived1::Scalar>& rA = A.derived();

    // EXCEPTION CHECKS

    // check zero size
    if (!internal::check_nonzero_size(rstate))
        throw exception::ZeroSize("qpp::applyCTRL_fan()", "state");

    // check valid dims
    if (d < 2)
        throw exception::DimsInvalid("qpp::applyCTRL_fan()", "d");
    // END EXCEPTION CHECKS

    idx n = internal::get_num_subsys(static_cast<idx>(rstate.rows()), d);
    std::vector<idx> dims(n, d); // local dimensions vector

    return applyCTRL_fan(rstate, rA, ctrl, target, dims, shift);
}

// as in https://arxiv.org/abs/1707.08834
/**
 * \brief Applies the qudit quantum Fourier transform to the part \a target
 * of the multi-partite state vector or density matrix \a A
 * \see qpp::QFT()
 *
 * \param A Eigen expression
 * \param target Subsystem indexes where the QFT is applied
 * \param d Subsystem dimensions
 * \param swap Swaps the qubits/qudits at the end (true by default)
 * \return Qudit Quantum Fourier transform applied to the part \a target of
 * \a A
 */
template <typename Derived>
[[qpp::critical]] expr_t<Derived> applyQFT(const Eigen::MatrixBase<Derived>& A,
                                           const std::vector<idx>& target,
                                           idx d = 2, bool swap = true) {
    const expr_t<Derived>& rA = A.derived();

    // EXCEPTION CHECKS

    // check zero sizes
    if (!internal::check_nonzero_size(rA))
        throw exception::ZeroSize("qpp::applyQFT()", "A");

    // check valid subsystem dimension
    if (d < 2)
        throw exception::DimsInvalid("qpp::applyQFT()", "d");

    // total number of qubits/qudits in the state
    idx n = internal::get_num_subsys(static_cast<idx>(rA.rows()), d);

    std::vector<idx> dims(n, d); // local dimensions vector

    // check that target is valid w.r.t. dims
    if (!internal::check_subsys_match_dims(target, dims))
        throw exception::SubsysMismatchDims("qpp::applyQFT()", "dims/target");

    // check valid state and matching dimensions
    if (internal::check_cvector(rA)) {
        if (!internal::check_dims_match_cvect(dims, rA))
            throw exception::DimsMismatchCvector("qpp::applyQFT()", "A/dims");
    } else if (internal::check_square_mat(rA)) {
        if (!internal::check_dims_match_mat(dims, rA))
            throw exception::DimsMismatchMatrix("qpp::applyQFT()", "A/dims");
    } else
        throw exception::MatrixNotSquareNorCvector("qpp::applyQFT()", "A");
    // END EXCEPTION CHECKS

    expr_t<Derived> result = rA;

    idx n_subsys = target.size();

    if (d == 2) // qubits
    {
        for (idx i = 0; i < n_subsys; ++i) {
            // apply Hadamard on qubit i
            result = apply(result, Gates::get_no_thread_local_instance().H,
                           {target[i]});
            // apply controlled rotations
            for (idx j = 2; j <= n_subsys - i; ++j) {
                // construct Rj
                cmat Rj(2, 2);
                Rj << 1, 0, 0,
                    std::exp(static_cast<cplx::value_type>(2.0 * pi) * 1_i /
                             static_cast<cplx::value_type>(std::pow(2, j)));
                result =
                    applyCTRL(result, Rj, {target[i + j - 1]}, {target[i]});
            }
        }
        if (swap) {
            // we have the qubits in reversed order, we must swap them
            for (idx i = 0; i < n_subsys / 2; ++i) {
                result =
                    apply(result, Gates::get_no_thread_local_instance().SWAP,
                          {target[i], target[n_subsys - i - 1]});
            }
        }

    } else { // qudits
        for (idx i = 0; i < n_subsys; ++i) {
            // apply qudit Fourier on qudit i
            result = apply(result, Gates::get_no_thread_local_instance().Fd(d),
                           {target[i]}, d);
            // apply controlled rotations
            for (idx j = 2; j <= n_subsys - i; ++j) {
                // construct Rj
                cmat Rj = cmat::Zero(d, d);
                for (idx m = 0; m < d; ++m) {
                    Rj(m, m) = std::exp(
                        static_cast<cplx::value_type>(2.0 * pi * m) * 1_i /
                        static_cast<cplx::value_type>(std::pow(d, j)));
                }
                result =
                    applyCTRL(result, Rj, {target[i + j - 1]}, {target[i]}, d);
            }
        }
        if (swap) {
            // we have the qudits in reversed order, we must swap them
            for (idx i = 0; i < n_subsys / 2; ++i) {
                result = apply(result,
                               Gates::get_no_thread_local_instance().SWAPd(d),
                               {target[i], target[n_subsys - i - 1]}, d);
            }
        }
    }

    return result;
}

// as in https://arxiv.org/abs/1707.08834
/**
 * \brief Applies the inverse (adjoint) qudit quantum Fourier transform to
 * the part \a target of the multi-partite state vector or density matrix \a A
 * \see qpp::TFQ()
 *
 * \param A Eigen expression
 * \param target Subsystem indexes where the TFQ is applied
 * \param d Subsystem dimensions
 * \param swap Swaps the qubits/qudits at the end (true by default)
 * \return Inverse (adjoint) qudit Quantum Fourier transform applied to the
 * part \a target of \a A
 */
template <typename Derived>
[[qpp::critical]] expr_t<Derived> applyTFQ(const Eigen::MatrixBase<Derived>& A,
                                           const std::vector<idx>& target,
                                           idx d = 2, bool swap = true) {
    const expr_t<Derived>& rA = A.derived();

    // EXCEPTION CHECKS

    // check zero sizes
    if (!internal::check_nonzero_size(rA))
        throw exception::ZeroSize("qpp::applyTFQ()", "A");

    // check valid subsystem dimension
    if (d < 2)
        throw exception::DimsInvalid("qpp::applyTFQ()", "d");

    // total number of qubits/qudits in the state
    idx n = internal::get_num_subsys(static_cast<idx>(rA.rows()), d);

    std::vector<idx> dims(n, d); // local dimensions vector

    // check that target is valid w.r.t. dims
    if (!internal::check_subsys_match_dims(target, dims))
        throw exception::SubsysMismatchDims("qpp::applyTFQ()", "dims/target");

    // check valid state and matching dimensions
    if (internal::check_cvector(rA)) {
        if (!internal::check_dims_match_cvect(dims, rA))
            throw exception::DimsMismatchCvector("qpp::applyTFQ()", "A/dims");
    } else if (internal::check_square_mat(rA)) {
        if (!internal::check_dims_match_mat(dims, rA))
            throw exception::DimsMismatchMatrix("qpp::applyTFQ()", "A/dims");
    } else
        throw exception::MatrixNotSquareNorCvector("qpp::applyTFQ()", "A");
    // END EXCEPTION CHECKS

    expr_t<Derived> result = rA;

    idx n_subsys = target.size();

    if (d == 2) // qubits
    {
        if (swap) {
            // we have the qubits in reversed order, we must swap them
            for (idx i = n_subsys / 2; i-- > 0;) {
                result =
                    apply(result, Gates::get_no_thread_local_instance().SWAP,
                          {target[i], target[n_subsys - i - 1]});
            }
        }
        for (idx i = n_subsys; i-- > 0;) {
            // apply controlled rotations
            for (idx j = n_subsys - i + 1; j-- > 2;) {
                // construct Rj
                cmat Rj(2, 2);
                Rj << 1, 0, 0,
                    std::exp(static_cast<cplx::value_type>(-2.0 * pi) * 1_i /
                             static_cast<cplx::value_type>(std::pow(2, j)));
                result =
                    applyCTRL(result, Rj, {target[i + j - 1]}, {target[i]});
            }
            // apply Hadamard on qubit i
            result = apply(result, Gates::get_no_thread_local_instance().H,
                           {target[i]});
        }
    } else { // qudits
        if (swap) {
            // we have the qudits in reversed order, we must swap them
            for (idx i = n_subsys / 2; i-- > 0;) {
                result = apply(result,
                               Gates::get_no_thread_local_instance().SWAPd(d),
                               {target[i], target[n_subsys - i - 1]}, d);
            }
        }
        for (idx i = n_subsys; i-- > 0;) {
            // apply controlled rotations
            for (idx j = n_subsys - i + 1; j-- > 2;) {
                // construct Rj
                cmat Rj = cmat::Zero(d, d);
                for (idx m = 0; m < d; ++m) {
                    Rj(m, m) = std::exp(
                        static_cast<cplx::value_type>(-2.0 * pi * m) * 1_i /
                        static_cast<cplx::value_type>(std::pow(d, j)));
                }
                result =
                    applyCTRL(result, Rj, {target[i + j - 1]}, {target[i]}, d);
            }
            // apply qudit Fourier on qudit i
            result = apply(result,
                           adjoint(Gates::get_no_thread_local_instance().Fd(d)),
                           {target[i]}, d);
        }
    }

    return result;
}

// as in https://arxiv.org/abs/1707.08834
/**
 * \brief Qudit quantum Fourier transform
 * \see qpp::applyQFT()
 *
 * \param A Eigen expression
 * \param d Subsystem dimensions
 * \param swap Swaps the qubits/qudits at the end (true by default)
 * \return Qudit quantum Fourier transform applied on \a A
 */
template <typename Derived>
expr_t<Derived> QFT(const Eigen::MatrixBase<Derived>& A, idx d = 2,
                    bool swap = true) {
    const expr_t<Derived>& rA = A.derived();

    // EXCEPTION CHECKS

    // check zero-size
    if (!internal::check_nonzero_size(rA))
        throw exception::ZeroSize("qpp::QFT()", "A");

    // check valid subsystem dimension
    if (d < 2)
        throw exception::DimsInvalid("qpp::QFT()", "d");

    // total number of qubits/qudits in the state
    idx n = internal::get_num_subsys(static_cast<idx>(rA.rows()), d);

    std::vector<idx> dims(n, d); // local dimensions vector

    // check valid state and matching dimensions
    if (internal::check_cvector(rA)) {
        if (!internal::check_dims_match_cvect(dims, rA))
            throw exception::DimsMismatchCvector("qpp::QFT()", "A/dims");
    } else if (internal::check_square_mat(rA)) {
        if (!internal::check_dims_match_mat(dims, rA))
            throw exception::DimsMismatchMatrix("qpp::QFT()", "A/dims");
    } else
        throw exception::MatrixNotSquareNorCvector("qpp::QFT()", "A");
    // END EXCEPTION CHECKS

    std::vector<idx> subsys(n);
    std::iota(subsys.begin(), subsys.end(), 0);
    ket result = applyQFT(rA, subsys, d, swap);

    return result;
}

// as in https://arxiv.org/abs/1707.08834
/**
 * \brief Inverse (adjoint) qudit quantum Fourier transform
 * \see qpp::applyTFQ()
 *
 * \param A Eigen expression
 * \param d Subsystem dimensions
 * \param swap Swaps the qubits/qudits at the end (true by default)
 * \return Inverse (adjoint) qudit quantum Fourier transform applied on \a A
 */
template <typename Derived>
expr_t<Derived> TFQ(const Eigen::MatrixBase<Derived>& A, idx d = 2,
                    bool swap = true) {
    const expr_t<Derived>& rA = A.derived();

    // EXCEPTION CHECKS

    // check zero-size
    if (!internal::check_nonzero_size(rA))
        throw exception::ZeroSize("qpp::TFQ()", "A");

    // check valid subsystem dimension
    if (d < 2)
        throw exception::DimsInvalid("qpp::TFQ()", "d");

    // total number of qubits/qudits in the state
    idx n = internal::get_num_subsys(static_cast<idx>(rA.rows()), d);

    std::vector<idx> dims(n, d); // local dimensions vector

    // check valid state and matching dimensions
    if (internal::check_cvector(rA)) {
        if (!internal::check_dims_match_cvect(dims, rA))
            throw exception::DimsMismatchCvector("qpp::QFT()", "A/dims");
    } else if (internal::check_square_mat(rA)) {
        if (!internal::check_dims_match_mat(dims, rA))
            throw exception::DimsMismatchMatrix("qpp::QFT()", "A/dims");
    } else
        throw exception::MatrixNotSquareNorCvector("qpp::QFT()", "A");
    // END EXCEPTION CHECKS

    std::vector<idx> subsys(n);
    std::iota(subsys.begin(), subsys.end(), 0);
    ket result = applyTFQ(rA, subsys, d, swap);

    return result;
}

/**
 * \brief Quantumly-accessible Random Access Memory (qRAM) over classical
 data,
 * implements \f$\sum_j\alpha_j|j\rangle\stackrel{qRAM}{\longrightarrow}
   \sum_j\alpha_j|j\rangle|m_j\rangle\f$
 *
 * \see
 * <a
 href="https://iopscience.iop.org/article/10.1088/1367-2630/17/12/123010">
 * New Journal of Physics, Vol. 17, No. 12, Pp. 123010 (2015)</a>
 *
 * \param psi Column vector Eigen expression, input amplitudes
 * \param data Vector storing the classical data
 * \param DqRAM qRAM subsystem dimension (must be at least 1 + maximum value
 * stored in the qRAM)
 * \return Superposition over the qRAM values
 */
template <typename Derived>
[[qpp::parallel]] dyn_col_vect<typename Derived::Scalar>
qRAM(const Eigen::MatrixBase<Derived>& psi, const qram& data, idx DqRAM) {
    const dyn_mat<typename Derived::Scalar>& rpsi = psi.derived();

    // EXCEPTION CHECKS

    // check zero-size
    if (!internal::check_nonzero_size(rpsi))
        throw exception::ZeroSize("qpp::qRAM()", "psi");
    if (!internal::check_nonzero_size(data))
        throw exception::ZeroSize("qpp::qRAM()", "data");

    // check column vector
    if (!internal::check_cvector(rpsi))
        throw exception::MatrixNotCvector("qpp::qRAM()", "psi");

    // check equal dimensions
    idx Din = static_cast<idx>(rpsi.rows());
    if (static_cast<idx>(data.size()) < Din)
        throw exception::DimsInvalid("qpp::qRAM()", "data/psi");

    idx max_val_qRAM = *std::max_element(data.begin(), data.end());
    if (DqRAM <= max_val_qRAM)
        throw exception::OutOfRange("qpp::qRAM()", "data/DqRAM");
    // END EXCEPTION CHECKS

    idx Dout = Din * DqRAM;
    ket result(Dout);

#ifdef QPP_OPENMP
// NOLINTNEXTLINE
#pragma omp parallel for
#endif // QPP_OPENMP
    for (idx i = 0; i < Din; ++i) {
        dyn_col_vect<typename Derived::Scalar> ket_i =
            dyn_col_vect<typename Derived::Scalar>::Zero(DqRAM);
        ket_i(data[i]) = 1;
        result.block(i * DqRAM, 0, DqRAM, 1) = rpsi(i) * ket_i;
    }

    return result;
}

/**
 * \brief Quantumly-accessible Random Access Memory (qRAM) over classical
 data,
 * implements \f$\sum_j\alpha_j|j\rangle\stackrel{qRAM}{\longrightarrow}
   \sum_j\alpha_j|j\rangle|m_j\rangle\f$
 *
 * \see
 * <a
 href="https://iopscience.iop.org/article/10.1088/1367-2630/17/12/123010">
 * New Journal of Physics, Vol. 17, No. 12, Pp. 123010 (2015)</a>
 *
 * \note The qRAM subsystem dimension is set to 1 + maximum value stored in
 the
 * qRAM
 *
 * \param psi Column vector Eigen expression, input amplitudes
 * \param data Vector storing the classical data
 * \return Superposition over the qRAM values
 */
template <typename Derived>
dyn_col_vect<typename Derived::Scalar>
qRAM(const Eigen::MatrixBase<Derived>& psi, const qram& data) {
    // EXCEPTION CHECKS

    // check zero-size
    if (!internal::check_nonzero_size(data))
        throw exception::ZeroSize("qpp::qRAM()", "data");
    // END EXCEPTION CHECKS

    idx max_val_qRAM = *std::max_element(data.begin(), data.end());

    return qRAM(psi, data, 1 + max_val_qRAM);
}
} /* namespace qpp */

#endif /* QPP_OPERATIONS_HPP_ */
