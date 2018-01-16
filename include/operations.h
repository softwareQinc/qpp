/*
 * This file is part of Quantum++.
 *
 * MIT License
 *
 * Copyright (c) 2013 - 2018 Vlad Gheorghiu (vgheorgh@gmail.com)
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
* \file operations.h
* \brief Quantum operation functions
*/

#ifndef OPERATIONS_H_
#define OPERATIONS_H_

namespace qpp {

/**
* \brief Applies the controlled-gate \a A to the part \a subsys
* of the multi-partite state vector or density matrix \a state
* \see qpp::Gates::CTRL()
*
* \note The dimension of the gate \a A must match
* the dimension of \a subsys.
* Also, all control subsystems in \a ctrl must have the same dimension.
*
* \param state Eigen expression
* \param A Eigen expression
* \param ctrl Control subsystem indexes
* \param subsys Subsystem indexes where the gate \a A is applied
* \param dims Dimensions of the multi-partite system
* \return CTRL-A gate applied to the part \a subsys of \a state
*/
template <typename Derived1, typename Derived2>
dyn_mat<typename Derived1::Scalar>
applyCTRL(const Eigen::MatrixBase<Derived1>& state,
          const Eigen::MatrixBase<Derived2>& A, const std::vector<idx>& ctrl,
          const std::vector<idx>& subsys, const std::vector<idx>& dims) {
    const typename Eigen::MatrixBase<Derived1>::EvalReturnType& rstate =
        state.derived();
    const dyn_mat<typename Derived2::Scalar>& rA = A.derived();

    // EXCEPTION CHECKS

    // check types
    if (!std::is_same<typename Derived1::Scalar,
                      typename Derived2::Scalar>::value)
        throw exception::TypeMismatch("qpp::applyCTRL()");

    // check zero sizes
    if (!internal::check_nonzero_size(rA))
        throw exception::ZeroSize("qpp::applyCTRL()");

    // check zero sizes
    if (!internal::check_nonzero_size(rstate))
        throw exception::ZeroSize("qpp::applyCTRL()");

    // check square matrix for the gate
    if (!internal::check_square_mat(rA))
        throw exception::MatrixNotSquare("qpp::applyCTRL()");

    // check that all control subsystems have the same dimension
    idx d = ctrl.size() > 0 ? dims[ctrl[0]] : 1;
    for (idx i = 1; i < ctrl.size(); ++i)
        if (dims[ctrl[i]] != d)
            throw exception::DimsNotEqual("qpp::applyCTRL()");

    // check that dimension is valid
    if (!internal::check_dims(dims))
        throw exception::DimsInvalid("qpp::applyCTRL()");

    // check subsys is valid w.r.t. dims
    if (!internal::check_subsys_match_dims(subsys, dims))
        throw exception::SubsysMismatchDims("qpp::applyCTRL()");

    // check that gate matches the dimensions of the subsys
    std::vector<idx> subsys_dims(subsys.size());
    for (idx i = 0; i < subsys.size(); ++i)
        subsys_dims[i] = dims[subsys[i]];
    if (!internal::check_dims_match_mat(subsys_dims, rA))
        throw exception::MatrixMismatchSubsys("qpp::applyCTRL()");

    std::vector<idx> ctrlgate = ctrl; // ctrl + gate subsystem vector
    ctrlgate.insert(std::end(ctrlgate), std::begin(subsys), std::end(subsys));
    std::sort(std::begin(ctrlgate), std::end(ctrlgate));

    // check that ctrl + gate subsystem is valid
    // with respect to local dimensions
    if (!internal::check_subsys_match_dims(ctrlgate, dims))
        throw exception::SubsysMismatchDims("qpp::applyCTRL()");
    // END EXCEPTION CHECKS

    // construct the table of A^i and (A^dagger)^i
    std::vector<dyn_mat<typename Derived1::Scalar>> Ai;
    std::vector<dyn_mat<typename Derived1::Scalar>> Aidagger;
    for (idx i = 0; i < std::max(d, static_cast<idx>(2)); ++i) {
        Ai.push_back(powm(rA, i));
        Aidagger.push_back(powm(adjoint(rA), i));
    }

    idx D = static_cast<idx>(rstate.rows()); // total dimension
    idx N = dims.size();                     // total number of subsystems
    idx ctrlsize = ctrl.size();              // number of ctrl subsystem
    idx ctrlgatesize = ctrlgate.size();      // number of ctrl+gate subsystems
    idx subsyssize = subsys.size(); // number of subsystems of the target
    // dimension of ctrl subsystem
    idx Dctrl = static_cast<idx>(std::llround(std::pow(d, ctrlsize)));
    idx DA = static_cast<idx>(rA.rows()); // dimension of gate subsystem

    idx Cdims[maxn];          // local dimensions
    idx CdimsA[maxn];         // local dimensions
    idx CdimsCTRL[maxn];      // local dimensions
    idx CdimsCTRLA_bar[maxn]; // local dimensions

    // compute the complementary subsystem of ctrlgate w.r.t. dims
    std::vector<idx> ctrlgate_bar = complement(ctrlgate, N);
    // number of subsystems that are complementary to the ctrl+gate
    idx ctrlgate_barsize = ctrlgate_bar.size();

    idx DCTRLA_bar = 1; // dimension of the rest
    for (idx i = 0; i < ctrlgate_barsize; ++i)
        DCTRLA_bar *= dims[ctrlgate_bar[i]];

    for (idx k = 0; k < N; ++k)
        Cdims[k] = dims[k];
    for (idx k = 0; k < subsyssize; ++k)
        CdimsA[k] = dims[subsys[k]];
    for (idx k = 0; k < ctrlsize; ++k)
        CdimsCTRL[k] = d;
    for (idx k = 0; k < ctrlgate_barsize; ++k)
        CdimsCTRLA_bar[k] = dims[ctrlgate_bar[k]];

    // worker, computes the coefficient and the index for the ket case
    // used in #pragma omp parallel for collapse
    auto coeff_idx_ket = [&](idx i_, idx m_, idx r_) noexcept
                             ->std::pair<typename Derived1::Scalar, idx> {
        idx indx = 0;
        typename Derived1::Scalar coeff = 0;

        idx Cmidx[maxn];          // the total multi-index
        idx CmidxA[maxn];         // the gate part multi-index
        idx CmidxCTRLA_bar[maxn]; // the rest multi-index

        // compute the index

        // set the CTRL part
        for (idx k = 0; k < ctrlsize; ++k) {
            Cmidx[ctrl[k]] = i_;
        }

        // set the rest
        internal::n2multiidx(r_, N - ctrlgatesize, CdimsCTRLA_bar,
                             CmidxCTRLA_bar);
        for (idx k = 0; k < N - ctrlgatesize; ++k) {
            Cmidx[ctrlgate_bar[k]] = CmidxCTRLA_bar[k];
        }

        // set the A part
        internal::n2multiidx(m_, subsyssize, CdimsA, CmidxA);
        for (idx k = 0; k < subsyssize; ++k) {
            Cmidx[subsys[k]] = CmidxA[k];
        }

        // we now got the total index
        indx = internal::multiidx2n(Cmidx, N, Cdims);

        // compute the coefficient
        for (idx n_ = 0; n_ < DA; ++n_) {
            internal::n2multiidx(n_, subsyssize, CdimsA, CmidxA);
            for (idx k = 0; k < subsyssize; ++k) {
                Cmidx[subsys[k]] = CmidxA[k];
            }
            coeff +=
                Ai[i_](m_, n_) * rstate(internal::multiidx2n(Cmidx, N, Cdims));
        }

        return std::make_pair(coeff, indx);
    }; /* end coeff_idx_ket */

    // worker, computes the coefficient and the index
    // for the density matrix case
    // used in #pragma omp parallel for collapse
    auto coeff_idx_rho = [&](idx i1_, idx m1_, idx r1_, idx i2_, idx m2_,
                             idx r2_) noexcept
                             ->std::tuple<typename Derived1::Scalar, idx, idx> {
        idx idxrow = 0;
        idx idxcol = 0;
        typename Derived1::Scalar coeff = 0, lhs = 1, rhs = 1;

        idx Cmidxrow[maxn];          // the total row multi-index
        idx Cmidxcol[maxn];          // the total col multi-index
        idx CmidxArow[maxn];         // the gate part row multi-index
        idx CmidxAcol[maxn];         // the gate part col multi-index
        idx CmidxCTRLrow[maxn];      // the control row multi-index
        idx CmidxCTRLcol[maxn];      // the control col multi-index
        idx CmidxCTRLA_barrow[maxn]; // the rest row multi-index
        idx CmidxCTRLA_barcol[maxn]; // the rest col multi-index

        // compute the ket/bra indexes

        // set the CTRL part
        internal::n2multiidx(i1_, ctrlsize, CdimsCTRL, CmidxCTRLrow);
        internal::n2multiidx(i2_, ctrlsize, CdimsCTRL, CmidxCTRLcol);

        for (idx k = 0; k < ctrlsize; ++k) {
            Cmidxrow[ctrl[k]] = CmidxCTRLrow[k];
            Cmidxcol[ctrl[k]] = CmidxCTRLcol[k];
        }

        // set the rest
        internal::n2multiidx(r1_, N - ctrlgatesize, CdimsCTRLA_bar,
                             CmidxCTRLA_barrow);
        internal::n2multiidx(r2_, N - ctrlgatesize, CdimsCTRLA_bar,
                             CmidxCTRLA_barcol);
        for (idx k = 0; k < N - ctrlgatesize; ++k) {
            Cmidxrow[ctrlgate_bar[k]] = CmidxCTRLA_barrow[k];
            Cmidxcol[ctrlgate_bar[k]] = CmidxCTRLA_barcol[k];
        }

        // set the A part
        internal::n2multiidx(m1_, subsyssize, CdimsA, CmidxArow);
        internal::n2multiidx(m2_, subsyssize, CdimsA, CmidxAcol);
        for (idx k = 0; k < subsys.size(); ++k) {
            Cmidxrow[subsys[k]] = CmidxArow[k];
            Cmidxcol[subsys[k]] = CmidxAcol[k];
        }

        // we now got the total row/col indexes
        idxrow = internal::multiidx2n(Cmidxrow, N, Cdims);
        idxcol = internal::multiidx2n(Cmidxcol, N, Cdims);

        // check whether all CTRL row and col multi indexes are equal
        bool all_ctrl_rows_equal = true;
        bool all_ctrl_cols_equal = true;

        idx first_ctrl_row, first_ctrl_col;
        if (ctrlsize > 0) {
            first_ctrl_row = CmidxCTRLrow[0];
            first_ctrl_col = CmidxCTRLcol[0];
        } else {
            first_ctrl_row = first_ctrl_col = 1;
        }

        for (idx k = 1; k < ctrlsize; ++k) {
            if (CmidxCTRLrow[k] != first_ctrl_row) {
                all_ctrl_rows_equal = false;
                break;
            }
        }
        for (idx k = 1; k < ctrlsize; ++k) {
            if (CmidxCTRLcol[k] != first_ctrl_col) {
                all_ctrl_cols_equal = false;
                break;
            }
        }

        // at least one control activated, compute the coefficient
        for (idx n1_ = 0; n1_ < DA; ++n1_) {
            internal::n2multiidx(n1_, subsyssize, CdimsA, CmidxArow);
            for (idx k = 0; k < subsyssize; ++k) {
                Cmidxrow[subsys[k]] = CmidxArow[k];
            }
            idx idxrowtmp = internal::multiidx2n(Cmidxrow, N, Cdims);

            if (all_ctrl_rows_equal) {
                lhs = Ai[first_ctrl_row](m1_, n1_);
            } else {
                lhs = (m1_ == n1_) ? 1 : 0; // identity matrix
            }

            for (idx n2_ = 0; n2_ < DA; ++n2_) {
                internal::n2multiidx(n2_, subsyssize, CdimsA, CmidxAcol);
                for (idx k = 0; k < subsyssize; ++k) {
                    Cmidxcol[subsys[k]] = CmidxAcol[k];
                }

                if (all_ctrl_cols_equal) {
                    rhs = Aidagger[first_ctrl_col](n2_, m2_);
                } else {
                    rhs = (n2_ == m2_) ? 1 : 0; // identity matrix
                }

                idx idxcoltmp = internal::multiidx2n(Cmidxcol, N, Cdims);

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
            throw exception::DimsMismatchCvector("qpp::applyCTRL()");
        if (D == 1)
            return rstate;

        dyn_mat<typename Derived1::Scalar> result = rstate;

#ifdef WITH_OPENMP_
#pragma omp parallel for collapse(2)
#endif // WITH_OPENMP_
        for (idx m = 0; m < DA; ++m)
            for (idx r = 0; r < DCTRLA_bar; ++r) {
                if (ctrlsize == 0) // no control
                {
                    result(coeff_idx_ket(1, m, r).second) =
                        coeff_idx_ket(1, m, r).first;
                } else
                    for (idx i = 0; i < d; ++i) {
                        result(coeff_idx_ket(i, m, r).second) =
                            coeff_idx_ket(i, m, r).first;
                    }
            }

        return result;
    }
    //************ density matrix ************//
    else if (internal::check_square_mat(rstate)) // we have a density operator
    {
        // check that dims match state matrix
        if (!internal::check_dims_match_mat(dims, rstate))
            throw exception::DimsMismatchMatrix("qpp::applyCTRL()");

        if (D == 1)
            return rstate;

        dyn_mat<typename Derived1::Scalar> result = rstate;

#ifdef WITH_OPENMP_
#pragma omp parallel for collapse(4)
#endif // WITH_OPENMP_
        for (idx m1 = 0; m1 < DA; ++m1)
            for (idx r1 = 0; r1 < DCTRLA_bar; ++r1)
                for (idx m2 = 0; m2 < DA; ++m2)
                    for (idx r2 = 0; r2 < DCTRLA_bar; ++r2)
                        if (ctrlsize == 0) // no control
                        {
                            auto coeff_idxes =
                                coeff_idx_rho(1, m1, r1, 1, m2, r2);
                            result(std::get<1>(coeff_idxes),
                                   std::get<2>(coeff_idxes)) =
                                std::get<0>(coeff_idxes);
                        } else {
                            for (idx i1 = 0; i1 < Dctrl; ++i1)
                                for (idx i2 = 0; i2 < Dctrl; ++i2) {
                                    auto coeff_idxes =
                                        coeff_idx_rho(i1, m1, r1, i2, m2, r2);
                                    result(std::get<1>(coeff_idxes),
                                           std::get<2>(coeff_idxes)) =
                                        std::get<0>(coeff_idxes);
                                }
                        }

        return result;
    }
    //************ Exception: not ket nor density matrix ************//
    else
        throw exception::MatrixNotSquareNorCvector("qpp::applyCTRL()");
}

/**
* \brief Applies the controlled-gate \a A to the part \a subsys
* of the multi-partite state vector or density matrix \a state
* \see qpp::Gates::CTRL()
*
* \note The dimension of the gate \a A must match
* the dimension of \a subsys
*
* \param state Eigen expression
* \param A Eigen expression
* \param ctrl Control subsystem indexes
* \param subsys Subsystem indexes where the gate \a A is applied
* \param d Subsystem dimensions
* \return CTRL-A gate applied to the part \a subsys of \a state
*/
template <typename Derived1, typename Derived2>
dyn_mat<typename Derived1::Scalar>
applyCTRL(const Eigen::MatrixBase<Derived1>& state,
          const Eigen::MatrixBase<Derived2>& A, const std::vector<idx>& ctrl,
          const std::vector<idx>& subsys, idx d = 2) {
    const typename Eigen::MatrixBase<Derived1>::EvalReturnType& rstate =
        state.derived();
    const dyn_mat<typename Derived1::Scalar>& rA = A.derived();

    // EXCEPTION CHECKS

    // check zero size
    if (!internal::check_nonzero_size(rstate))
        throw exception::ZeroSize("qpp::applyCTRL()");

    // check valid dims
    if (d < 2)
        throw exception::DimsInvalid("qpp::applyCTRL()");
    // END EXCEPTION CHECKS

    idx N = internal::get_num_subsys(static_cast<idx>(rstate.rows()), d);
    std::vector<idx> dims(N, d); // local dimensions vector

    return applyCTRL(rstate, rA, ctrl, subsys, dims);
}

/**
* \brief Applies the gate \a A to the part \a subsys
* of the multi-partite state vector or density matrix \a state
*
* \note The dimension of the gate \a A must match
* the dimension of \a subsys
*
* \param state Eigen expression
* \param A Eigen expression
* \param subsys Subsystem indexes where the gate \a A is applied
* \param dims Dimensions of the multi-partite system
* \return Gate \a A applied to the part \a subsys of \a state
*/
template <typename Derived1, typename Derived2>
dyn_mat<typename Derived1::Scalar>
apply(const Eigen::MatrixBase<Derived1>& state,
      const Eigen::MatrixBase<Derived2>& A, const std::vector<idx>& subsys,
      const std::vector<idx>& dims) {
    const typename Eigen::MatrixBase<Derived1>::EvalReturnType& rstate =
        state.derived();
    const dyn_mat<typename Derived2::Scalar>& rA = A.derived();

    // EXCEPTION CHECKS

    // check types
    if (!std::is_same<typename Derived1::Scalar,
                      typename Derived2::Scalar>::value)
        throw exception::TypeMismatch("qpp::apply()");

    // check zero sizes
    if (!internal::check_nonzero_size(rA))
        throw exception::ZeroSize("qpp::apply()");

    // check zero sizes
    if (!internal::check_nonzero_size(rstate))
        throw exception::ZeroSize("qpp::apply()");

    // check square matrix for the gate
    if (!internal::check_square_mat(rA))
        throw exception::MatrixNotSquare("qpp::apply()");

    // check that dimension is valid
    if (!internal::check_dims(dims))
        throw exception::DimsInvalid("qpp::apply()");

    // check subsys is valid w.r.t. dims
    if (!internal::check_subsys_match_dims(subsys, dims))
        throw exception::SubsysMismatchDims("qpp::apply()");

    // check that gate matches the dimensions of the subsys
    std::vector<idx> subsys_dims(subsys.size());
    for (idx i = 0; i < subsys.size(); ++i)
        subsys_dims[i] = dims[subsys[i]];
    if (!internal::check_dims_match_mat(subsys_dims, rA))
        throw exception::MatrixMismatchSubsys("qpp::apply()");
    // END EXCEPTION CHECKS

    //************ ket ************//
    if (internal::check_cvector(rstate)) // we have a ket
    {
        // check that dims match state vector
        if (!internal::check_dims_match_cvect(dims, rstate))
            throw exception::DimsMismatchCvector("qpp::apply()");

        return applyCTRL(rstate, rA, {}, subsys, dims);
    }
    //************ density matrix ************//
    else if (internal::check_square_mat(rstate)) // we have a density operator
    {

        // check that dims match state matrix
        if (!internal::check_dims_match_mat(dims, rstate))
            throw exception::DimsMismatchMatrix("qpp::apply()");

        return applyCTRL(rstate, rA, {}, subsys, dims);
    }
    //************ Exception: not ket nor density matrix ************//
    else
        throw exception::MatrixNotSquareNorCvector("qpp::apply()");
}

/**
* \brief Applies the gate \a A to the part \a subsys
* of the multi-partite state vector or density matrix \a state
*
* \note The dimension of the gate \a A must match
* the dimension of \a subsys
*
* \param state Eigen expression
* \param A Eigen expression
* \param subsys Subsystem indexes where the gate \a A is applied
* \param d Subsystem dimensions
* \return Gate \a A applied to the part \a subsys of \a state
*/
template <typename Derived1, typename Derived2>
dyn_mat<typename Derived1::Scalar>
apply(const Eigen::MatrixBase<Derived1>& state,
      const Eigen::MatrixBase<Derived2>& A, const std::vector<idx>& subsys,
      idx d = 2) {
    const typename Eigen::MatrixBase<Derived1>::EvalReturnType& rstate =
        state.derived();
    const dyn_mat<typename Derived1::Scalar>& rA = A.derived();

    // EXCEPTION CHECKS

    // check zero size
    if (!internal::check_nonzero_size(rstate))
        throw exception::ZeroSize("qpp::apply()");

    // check valid dims
    if (d < 2)
        throw exception::DimsInvalid("qpp::apply()");
    // END EXCEPTION CHECKS

    idx N = internal::get_num_subsys(static_cast<idx>(rstate.rows()), d);
    std::vector<idx> dims(N, d); // local dimensions vector

    return apply(rstate, rA, subsys, dims);
}

/**
* \brief Applies the channel specified by the set of Kraus operators \a Ks
* to the density matrix \a A
*
* \param A Eigen expression
* \param Ks Set of Kraus operators
* \return Output density matrix after the action of the channel
*/
template <typename Derived>
cmat apply(const Eigen::MatrixBase<Derived>& A, const std::vector<cmat>& Ks) {
    const cmat& rA = A.derived();

    // EXCEPTION CHECKS

    if (!internal::check_nonzero_size(rA))
        throw exception::ZeroSize("qpp::apply()");
    if (!internal::check_square_mat(rA))
        throw exception::MatrixNotSquare("qpp::apply()");
    if (Ks.size() == 0)
        throw exception::ZeroSize("qpp::apply()");
    if (!internal::check_square_mat(Ks[0]))
        throw exception::MatrixNotSquare("qpp::apply()");
    if (Ks[0].rows() != rA.rows())
        throw exception::DimsMismatchMatrix("qpp::apply()");
    for (auto&& it : Ks)
        if (it.rows() != Ks[0].rows() || it.cols() != Ks[0].rows())
            throw exception::DimsNotEqual("qpp::apply()");
    // END EXCEPTION CHECKS

    cmat result = cmat::Zero(rA.rows(), rA.rows());

#ifdef WITH_OPENMP_
#pragma omp parallel for
#endif // WITH_OPENMP_
    for (idx i = 0; i < Ks.size(); ++i) {
#ifdef WITH_OPENMP_
#pragma omp critical
#endif // WITH_OPENMP_
        { result += Ks[i] * rA * adjoint(Ks[i]); }
    }

    return result;
}

/**
* \brief Applies the channel specified by the set of Kraus operators \a Ks to
* the part \a subsys of the multi-partite density matrix \a A
*
* \param A Eigen expression
* \param Ks Set of Kraus operators
* \param subsys Subsystem indexes where the Kraus operators \a Ks are applied
* \param dims Dimensions of the multi-partite system
* \return Output density matrix after the action of the channel
*/
template <typename Derived>
cmat apply(const Eigen::MatrixBase<Derived>& A, const std::vector<cmat>& Ks,
           const std::vector<idx>& subsys, const std::vector<idx>& dims) {
    const cmat& rA = A.derived();

    // EXCEPTION CHECKS

    // check zero sizes
    if (!internal::check_nonzero_size(rA))
        throw exception::ZeroSize("qpp::apply()");

    // check square matrix for the A
    if (!internal::check_square_mat(rA))
        throw exception::MatrixNotSquare("qpp::apply()");

    // check that dimension is valid
    if (!internal::check_dims(dims))
        throw exception::DimsInvalid("qpp::apply()");

    // check that dims match A matrix
    if (!internal::check_dims_match_mat(dims, rA))
        throw exception::DimsMismatchMatrix("qpp::apply()");

    // check subsys is valid w.r.t. dims
    if (!internal::check_subsys_match_dims(subsys, dims))
        throw exception::SubsysMismatchDims("qpp::apply()");

    std::vector<idx> subsys_dims(subsys.size());
    for (idx i = 0; i < subsys.size(); ++i)
        subsys_dims[i] = dims[subsys[i]];

    // check the Kraus operators
    if (Ks.size() == 0)
        throw exception::ZeroSize("qpp::apply()");
    if (!internal::check_square_mat(Ks[0]))
        throw exception::MatrixNotSquare("qpp::apply()");
    if (!internal::check_dims_match_mat(subsys_dims, Ks[0]))
        throw exception::MatrixMismatchSubsys("qpp::apply()");
    for (auto&& it : Ks)
        if (it.rows() != Ks[0].rows() || it.cols() != Ks[0].rows())
            throw exception::DimsNotEqual("qpp::apply()");
    // END EXCEPTION CHECKS

    cmat result = cmat::Zero(rA.rows(), rA.rows());

    for (idx i = 0; i < Ks.size(); ++i)
        result += apply(rA, Ks[i], subsys, dims);

    return result;
}

/**
* \brief Applies the channel specified by the set of Kraus operators \a Ks to
* the part \a subsys of the multi-partite density matrix \a A
*
* \param A Eigen expression
* \param Ks Set of Kraus operators
* \param subsys Subsystem indexes where the Kraus operators \a Ks are applied
* \param d Subsystem dimensions
* \return Output density matrix after the action of the channel
*/
template <typename Derived>
cmat apply(const Eigen::MatrixBase<Derived>& A, const std::vector<cmat>& Ks,
           const std::vector<idx>& subsys, idx d = 2) {
    const cmat& rA = A.derived();

    // EXCEPTION CHECKS

    // check zero sizes
    if (!internal::check_nonzero_size(rA))
        throw exception::ZeroSize("qpp::apply()");

    // check valid dims
    if (d < 2)
        throw exception::DimsInvalid("qpp::apply()");
    // END EXCEPTION CHECKS

    idx N = internal::get_num_subsys(static_cast<idx>(rA.rows()), d);
    std::vector<idx> dims(N, d); // local dimensions vector

    return apply(rA, Ks, subsys, dims);
}

/**
* \brief Superoperator matrix
*
* Constructs the superoperator matrix of the channel specified by the set of
* Kraus operators \a Ks in the standard operator basis
* \f$\{|i\rangle\langle j|\}\f$ ordered in lexicographical order, i.e.
* \f$|0\rangle\langle 0|\f$, \f$|0\rangle\langle 1|\f$ etc.
*
* \param Ks Set of Kraus operators
* \return Superoperator matrix
*/
inline cmat kraus2super(const std::vector<cmat>& Ks) {
    // EXCEPTION CHECKS

    if (Ks.size() == 0)
        throw exception::ZeroSize("qpp::kraus2super()");
    if (!internal::check_nonzero_size(Ks[0]))
        throw exception::ZeroSize("qpp::kraus2super()");
    if (!internal::check_square_mat(Ks[0]))
        throw exception::MatrixNotSquare("qpp::kraus2super()");
    for (auto&& it : Ks)
        if (it.rows() != Ks[0].rows() || it.cols() != Ks[0].rows())
            throw exception::DimsNotEqual("qpp::kraus2super()");
    // END EXCEPTION CHECKS

    idx D = static_cast<idx>(Ks[0].rows());

    cmat result(D * D, D * D);
    cmat MN = cmat::Zero(D, D);
    bra A = bra::Zero(D);
    ket B = ket::Zero(D);
    cmat EMN = cmat::Zero(D, D);

#ifdef WITH_OPENMP_
#pragma omp parallel for collapse(2)
#endif // WITH_OPENMP_
    for (idx m = 0; m < D; ++m) {
        for (idx n = 0; n < D; ++n) {
#ifdef WITH_OPENMP_
#pragma omp critical
#endif // WITH_OPENMP_
            {
                // compute E(|m><n|)
                MN(m, n) = 1;
                for (idx i = 0; i < Ks.size(); ++i)
                    EMN += Ks[i] * MN * adjoint(Ks[i]);
                MN(m, n) = 0;

                for (idx a = 0; a < D; ++a) {
                    A(a) = 1;
                    for (idx b = 0; b < D; ++b) {
                        // compute result(ab,mn)=<a|E(|m><n)|b>
                        B(b) = 1;
                        result(a * D + b, m * D + n) =
                            static_cast<cmat>(A * EMN * B).value();
                        B(b) = 0;
                    }
                    A(a) = 0;
                }
                EMN = cmat::Zero(D, D);
            }
        }
    }

    return result;
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
* \note The superoperator matrix \f$S\f$ and the Choi matrix \f$ C\f$
* are related by \f$ S_{ab,mn} = C_{ma,nb}\f$
*
* \param Ks Set of Kraus operators
* \return Choi matrix
*/
inline cmat kraus2choi(const std::vector<cmat>& Ks) {
    // EXCEPTION CHECKS

    if (Ks.size() == 0)
        throw exception::ZeroSize("qpp::kraus2choi()");
    if (!internal::check_nonzero_size(Ks[0]))
        throw exception::ZeroSize("qpp::kraus2choi()");
    if (!internal::check_square_mat(Ks[0]))
        throw exception::MatrixNotSquare("qpp::kraus2choi()");
    for (auto&& it : Ks)
        if (it.rows() != Ks[0].rows() || it.cols() != Ks[0].rows())
            throw exception::DimsNotEqual("qpp::kraus2choi()");
    // END EXCEPTION CHECKS

    idx D = static_cast<idx>(Ks[0].rows());

    // construct the D x D \sum |jj> vector
    // (un-normalized maximally entangled state)
    cmat MES = cmat::Zero(D * D, 1);
    for (idx a = 0; a < D; ++a)
        MES(a * D + a) = 1;

    cmat Omega = MES * adjoint(MES);

    cmat result = cmat::Zero(D * D, D * D);

#ifdef WITH_OPENMP_
#pragma omp parallel for
#endif // WITH_OPENMP_
    for (idx i = 0; i < Ks.size(); ++i) {
#ifdef WITH_OPENMP_
#pragma omp critical
#endif // WITH_OPENMP_
        {
            result += kron(cmat::Identity(D, D), Ks[i]) * Omega *
                      adjoint(kron(cmat::Identity(D, D), Ks[i]));
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
* \note The Kraus operators satisfy \f$Tr(K_i^\dagger K_j)=\delta_{ij}\f$
* for all \f$i\neq j\f$
*
* \param A Choi matrix
* \return Set of orthogonal Kraus operators
*/
inline std::vector<cmat> choi2kraus(const cmat& A) {
    // EXCEPTION CHECKS

    if (!internal::check_nonzero_size(A))
        throw exception::ZeroSize("qpp::choi2kraus()");
    if (!internal::check_square_mat(A))
        throw exception::MatrixNotSquare("qpp::choi2kraus()");
    idx D = internal::get_dim_subsys(A.rows(), 2);
    // check equal dimensions
    if (D * D != static_cast<idx>(A.rows()))
        throw exception::DimsInvalid("qpp::choi2kraus()");
    // END EXCEPTION CHECKS

    dmat ev = hevals(A);
    cmat evec = hevects(A);
    std::vector<cmat> result;

    for (idx i = 0; i < D * D; ++i) {
        if (std::abs(ev(i)) > eps)
            result.push_back(std::sqrt(std::abs(ev(i))) *
                             reshape(evec.col(i), D, D));
    }

    return result;
}

/**
* \brief Converts Choi matrix to superoperator matrix
* \see qpp::super2choi()
*
* \param A Choi matrix
* \return Superoperator matrix
*/
inline cmat choi2super(const cmat& A) {
    // EXCEPTION CHECKS

    if (!internal::check_nonzero_size(A))
        throw exception::ZeroSize("qpp::choi2super()");
    if (!internal::check_square_mat(A))
        throw exception::MatrixNotSquare("qpp::choi2super()");
    idx D = internal::get_dim_subsys(static_cast<idx>(A.rows()), 2);
    // check equal dimensions
    if (D * D != static_cast<idx>(A.rows()))
        throw exception::DimsInvalid("qpp::choi2super()");
    // END EXCEPTION CHECKS

    cmat result(D * D, D * D);

#ifdef WITH_OPENMP_
#pragma omp parallel for collapse(4)
#endif // WITH_OPENMP_
    for (idx a = 0; a < D; ++a)
        for (idx b = 0; b < D; ++b)
            for (idx m = 0; m < D; ++m)
                for (idx n = 0; n < D; ++n)
                    result(a * D + b, m * D + n) = A(m * D + a, n * D + b);

    return result;
}

/**
* \brief Converts superoperator matrix to Choi matrix
* \see qpp::choi2super()
*
* \param A Superoperator matrix
* \return Choi matrix
*/
inline cmat super2choi(const cmat& A) {
    // EXCEPTION CHECKS

    if (!internal::check_nonzero_size(A))
        throw exception::ZeroSize("qpp::super2choi()");
    if (!internal::check_square_mat(A))
        throw exception::MatrixNotSquare("qpp::super2choi()");
    idx D = internal::get_dim_subsys(static_cast<idx>(A.rows()), 2);
    // check equal dimensions
    if (D * D != static_cast<idx>(A.rows()))
        throw exception::DimsInvalid("qpp::super2choi()");
    // END EXCEPTION CHECKS

    cmat result(D * D, D * D);

#ifdef WITH_OPENMP_
#pragma omp parallel for collapse(4)
#endif // WITH_OPENMP_
    for (idx a = 0; a < D; ++a)
        for (idx b = 0; b < D; ++b)
            for (idx m = 0; m < D; ++m)
                for (idx n = 0; n < D; ++n)
                    result(m * D + a, n * D + b) = A(a * D + b, m * D + n);

    return result;
}

/**
* \brief Partial trace
* \see qpp::ptrace2()
*
* Partial trace over the first subsystem
* of bi-partite state vector or density matrix
*
* \param A Eigen expression
* \param dims Dimensions of the bi-partite system
* \return Partial trace \f$Tr_{A}(\cdot)\f$ over the first subsytem \f$A\f$
* in a bi-partite system \f$A\otimes B\f$, as a dynamic matrix
* over the same scalar field as \a A
*/
template <typename Derived>
dyn_mat<typename Derived::Scalar> ptrace1(const Eigen::MatrixBase<Derived>& A,
                                          const std::vector<idx>& dims) {
    const typename Eigen::MatrixBase<Derived>::EvalReturnType& rA = A.derived();

    // EXCEPTION CHECKS

    // check zero-size
    if (!internal::check_nonzero_size(rA))
        throw exception::ZeroSize("qpp::ptrace1()");

    // check that dims is a valid dimension vector
    if (!internal::check_dims(dims))
        throw exception::DimsInvalid("qpp::ptrace1()");

    // check dims has only 2 elements
    if (dims.size() != 2)
        throw exception::NotBipartite("qpp::ptrace1()");
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
            throw exception::DimsMismatchCvector("qpp::ptrace1()");

        auto worker = [&](idx i, idx j) noexcept->typename Derived::Scalar {
            typename Derived::Scalar sum = 0;
            for (idx m = 0; m < DA; ++m)
                sum += rA(m * DB + i) * std::conj(rA(m * DB + j));

            return sum;
        }; /* end worker */

#ifdef WITH_OPENMP_
#pragma omp parallel for collapse(2)
#endif // WITH_OPENMP_
        // column major order for speed
        for (idx j = 0; j < DB; ++j)
            for (idx i = 0; i < DB; ++i)
                result(i, j) = worker(i, j);

        return result;
    }
    //************ density matrix ************//
    else if (internal::check_square_mat(rA)) // we have a density operator
    {
        // check that dims match the dimension of A
        if (!internal::check_dims_match_mat(dims, rA))
            throw exception::DimsMismatchMatrix("qpp::ptrace1()");

        auto worker = [&](idx i, idx j) noexcept->typename Derived::Scalar {
            typename Derived::Scalar sum = 0;
            for (idx m = 0; m < DA; ++m)
                sum += rA(m * DB + i, m * DB + j);

            return sum;
        }; /* end worker */

#ifdef WITH_OPENMP_
#pragma omp parallel for collapse(2)
#endif // WITH_OPENMP_
        // column major order for speed
        for (idx j = 0; j < DB; ++j)
            for (idx i = 0; i < DB; ++i)
                result(i, j) = worker(i, j);

        return result;
    }
    //************ Exception: not ket nor density matrix ************//
    else
        throw exception::MatrixNotSquareNorCvector("qpp::ptrace1()");
}

/**
* \brief Partial trace
* \see qpp::ptrace2()
*
* Partial trace over the first subsystem
* of bi-partite state vector or density matrix
*
* \param A Eigen expression
* \param d Subsystem dimensions
* \return Partial trace \f$Tr_{A}(\cdot)\f$ over the first subsytem \f$A\f$
* in a bi-partite system \f$A\otimes B\f$, as a dynamic matrix
* over the same scalar field as \a A
*/
template <typename Derived>
dyn_mat<typename Derived::Scalar> ptrace1(const Eigen::MatrixBase<Derived>& A,
                                          idx d = 2) {
    const dyn_mat<typename Derived::Scalar>& rA = A.derived();

    // EXCEPTION CHECKS

    // check zero size
    if (!internal::check_nonzero_size(A))
        throw exception::ZeroSize("qpp::ptrace1()");

    // check valid dims
    if (d == 0)
        throw exception::DimsInvalid("qpp::ptrace1()");
    // END EXCEPTION CHECKS

    std::vector<idx> dims(2, d); // local dimensions vector

    return ptrace1(A, dims);
}

/**
* \brief Partial trace
* \see qpp::ptrace1()
*
* Partial trace over the second subsystem
* of bi-partite state vector or density matrix
*
* \param A Eigen expression
* \param dims Dimensions of the bi-partite system
* \return Partial trace \f$Tr_{B}(\cdot)\f$ over the second subsytem \f$B\f$
* in a bi-partite system \f$A\otimes B\f$, as a dynamic matrix
* over the same scalar field as \a A
*/
template <typename Derived>
dyn_mat<typename Derived::Scalar> ptrace2(const Eigen::MatrixBase<Derived>& A,
                                          const std::vector<idx>& dims) {
    const typename Eigen::MatrixBase<Derived>::EvalReturnType& rA = A.derived();

    // EXCEPTION CHECKS

    // check zero-size
    if (!internal::check_nonzero_size(rA))
        throw exception::ZeroSize("qpp::ptrace2()");

    // check that dims is a valid dimension vector
    if (!internal::check_dims(dims))
        throw exception::DimsInvalid("qpp::ptrace2()");

    // check dims has only 2 elements
    if (dims.size() != 2)
        throw exception::NotBipartite("qpp::ptrace2()");
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
            throw exception::DimsMismatchCvector("qpp::ptrace2()");

        auto worker = [&](idx i, idx j) noexcept->typename Derived::Scalar {
            typename Derived::Scalar sum = 0;
            for (idx m = 0; m < DB; ++m)
                sum += rA(i * DB + m) * std::conj(rA(j * DB + m));

            return sum;
        }; /* end worker */

#ifdef WITH_OPENMP_
#pragma omp parallel for collapse(2)
#endif // WITH_OPENMP_
        // column major order for speed
        for (idx j = 0; j < DA; ++j)
            for (idx i = 0; i < DA; ++i)
                result(i, j) = worker(i, j);

        return result;
    }
    //************ density matrix ************//
    else if (internal::check_square_mat(rA)) // we have a density operator
    {
        // check that dims match the dimension of A
        if (!internal::check_dims_match_mat(dims, rA))
            throw exception::DimsMismatchMatrix("qpp::ptrace2()");

#ifdef WITH_OPENMP_
#pragma omp parallel for collapse(2)
#endif // WITH_OPENMP_
        // column major order for speed
        for (idx j = 0; j < DA; ++j)
            for (idx i = 0; i < DA; ++i)
                result(i, j) = trace(rA.block(i * DB, j * DB, DB, DB));

        return result;
    }
    //************ Exception: not ket nor density matrix ************//
    else
        throw exception::MatrixNotSquareNorCvector("qpp::ptrace1()");
}

/**
* \brief Partial trace
* \see qpp::ptrace1()
*
*  Partial trace over the second subsystem
*  of bi-partite state vector or density matrix
*
* \param A Eigen expression
* \param d Subsystem dimensions
* \return Partial trace \f$Tr_{B}(\cdot)\f$ over the second subsytem \f$B\f$
* in a bi-partite system \f$A\otimes B\f$, as a dynamic matrix
* over the same scalar field as \a A
*/
template <typename Derived>
dyn_mat<typename Derived::Scalar> ptrace2(const Eigen::MatrixBase<Derived>& A,
                                          idx d = 2) {
    const dyn_mat<typename Derived::Scalar>& rA = A.derived();

    // EXCEPTION CHECKS

    // check zero size
    if (!internal::check_nonzero_size(A))
        throw exception::ZeroSize("qpp::ptrace2()");

    // check valid dims
    if (d == 0)
        throw exception::DimsInvalid("qpp::ptrace2()");
    // END EXCEPTION CHECKS

    std::vector<idx> dims(2, d); // local dimensions vector

    return ptrace2(A, dims);
}

/**
* \brief Partial trace
* \see qpp::ptrace1(), qpp::ptrace2()
*
*  Partial trace of the multi-partite state vector or density matrix
*  over a list of subsystems
*
* \param A Eigen expression
* \param subsys Subsystem indexes
* \param dims Dimensions of the multi-partite system
* \return Partial trace \f$Tr_{subsys}(\cdot)\f$ over the subsytems \a subsys
* in a multi-partite system, as a dynamic matrix
* over the same scalar field as \a A
*/
template <typename Derived>
dyn_mat<typename Derived::Scalar> ptrace(const Eigen::MatrixBase<Derived>& A,
                                         const std::vector<idx>& subsys,
                                         const std::vector<idx>& dims) {
    const typename Eigen::MatrixBase<Derived>::EvalReturnType& rA = A.derived();

    // EXCEPTION CHECKS

    // check zero-size
    if (!internal::check_nonzero_size(rA))
        throw exception::ZeroSize("qpp::ptrace()");

    // check that dims is a valid dimension vector
    if (!internal::check_dims(dims))
        throw exception::DimsInvalid("qpp::ptrace()");

    // check that subsys are valid
    if (!internal::check_subsys_match_dims(subsys, dims))
        throw exception::SubsysMismatchDims("qpp::ptrace()");
    // END EXCEPTION CHECKS

    idx D = static_cast<idx>(rA.rows());
    idx N = dims.size();
    idx Nsubsys = subsys.size();
    idx Nsubsys_bar = N - Nsubsys;
    idx Dsubsys = 1;
    for (idx i = 0; i < Nsubsys; ++i)
        Dsubsys *= dims[subsys[i]];
    idx Dsubsys_bar = D / Dsubsys;

    idx Cdims[maxn];
    idx Csubsys[maxn];
    idx Cdimssubsys[maxn];
    idx Csubsys_bar[maxn];
    idx Cdimssubsys_bar[maxn];

    idx Cmidxcolsubsys_bar[maxn];

    std::vector<idx> subsys_bar = complement(subsys, N);
    std::copy(std::begin(subsys_bar), std::end(subsys_bar),
              std::begin(Csubsys_bar));

    for (idx i = 0; i < N; ++i) {
        Cdims[i] = dims[i];
    }
    for (idx i = 0; i < Nsubsys; ++i) {
        Csubsys[i] = subsys[i];
        Cdimssubsys[i] = dims[subsys[i]];
    }
    for (idx i = 0; i < Nsubsys_bar; ++i) {
        Cdimssubsys_bar[i] = dims[subsys_bar[i]];
    }

    dyn_mat<typename Derived::Scalar> result =
        dyn_mat<typename Derived::Scalar>(Dsubsys_bar, Dsubsys_bar);

    //************ ket ************//
    if (internal::check_cvector(rA)) // we have a ket
    {
        // check that dims match the dimension of A
        if (!internal::check_dims_match_cvect(dims, rA))
            throw exception::DimsMismatchCvector("qpp::ptrace()");

        if (subsys.size() == dims.size()) {
            result(0, 0) = (adjoint(rA) * rA).value();
            return result;
        }

        if (subsys.size() == 0)
            return rA * adjoint(rA);

        auto worker = [&](idx i) noexcept->typename Derived::Scalar {
            // use static allocation for speed!

            idx Cmidxrow[maxn];
            idx Cmidxcol[maxn];
            idx Cmidxrowsubsys_bar[maxn];
            idx Cmidxsubsys[maxn];

            /* get the row multi-indexes of the complement */
            internal::n2multiidx(i, Nsubsys_bar, Cdimssubsys_bar,
                                 Cmidxrowsubsys_bar);
            /* write them in the global row/col multi-indexes */
            for (idx k = 0; k < Nsubsys_bar; ++k) {
                Cmidxrow[Csubsys_bar[k]] = Cmidxrowsubsys_bar[k];
                Cmidxcol[Csubsys_bar[k]] = Cmidxcolsubsys_bar[k];
            }
            typename Derived::Scalar sm = 0;
            for (idx a = 0; a < Dsubsys; ++a) {
                // get the multi-index over which we do the summation
                internal::n2multiidx(a, Nsubsys, Cdimssubsys, Cmidxsubsys);
                // write it into the global row/col multi-indexes
                for (idx k = 0; k < Nsubsys; ++k)
                    Cmidxrow[Csubsys[k]] = Cmidxcol[Csubsys[k]] =
                        Cmidxsubsys[k];

                // now do the sum
                sm += rA(internal::multiidx2n(Cmidxrow, N, Cdims)) *
                      std::conj(rA(internal::multiidx2n(Cmidxcol, N, Cdims)));
            }

            return sm;
        }; /* end worker */

        for (idx j = 0; j < Dsubsys_bar; ++j) // column major order for speed
        {
            // compute the column multi-indexes of the complement
            internal::n2multiidx(j, Nsubsys_bar, Cdimssubsys_bar,
                                 Cmidxcolsubsys_bar);
#ifdef WITH_OPENMP_
#pragma omp parallel for
#endif // WITH_OPENMP_
            for (idx i = 0; i < Dsubsys_bar; ++i) {
                result(i, j) = worker(i);
            }
        }

        return result;
    }
    //************ density matrix ************//
    else if (internal::check_square_mat(rA)) // we have a density operator
    {
        // check that dims match the dimension of A
        if (!internal::check_dims_match_mat(dims, rA))
            throw exception::DimsMismatchMatrix("qpp::ptrace()");

        if (subsys.size() == dims.size()) {
            result(0, 0) = rA.trace();
            return result;
        }

        if (subsys.size() == 0)
            return rA;

        auto worker = [&](idx i) noexcept->typename Derived::Scalar {
            // use static allocation for speed!

            idx Cmidxrow[maxn];
            idx Cmidxcol[maxn];
            idx Cmidxrowsubsys_bar[maxn];
            idx Cmidxsubsys[maxn];

            /* get the row/col multi-indexes of the complement */
            internal::n2multiidx(i, Nsubsys_bar, Cdimssubsys_bar,
                                 Cmidxrowsubsys_bar);
            /* write them in the global row/col multi-indexes */
            for (idx k = 0; k < Nsubsys_bar; ++k) {
                Cmidxrow[Csubsys_bar[k]] = Cmidxrowsubsys_bar[k];
                Cmidxcol[Csubsys_bar[k]] = Cmidxcolsubsys_bar[k];
            }
            typename Derived::Scalar sm = 0;
            for (idx a = 0; a < Dsubsys; ++a) {
                // get the multi-index over which we do the summation
                internal::n2multiidx(a, Nsubsys, Cdimssubsys, Cmidxsubsys);
                // write it into the global row/col multi-indexes
                for (idx k = 0; k < Nsubsys; ++k)
                    Cmidxrow[Csubsys[k]] = Cmidxcol[Csubsys[k]] =
                        Cmidxsubsys[k];

                // now do the sum
                sm += rA(internal::multiidx2n(Cmidxrow, N, Cdims),
                         internal::multiidx2n(Cmidxcol, N, Cdims));
            }

            return sm;
        }; /* end worker */

        for (idx j = 0; j < Dsubsys_bar; ++j) // column major order for speed
        {
            // compute the column multi-indexes of the complement
            internal::n2multiidx(j, Nsubsys_bar, Cdimssubsys_bar,
                                 Cmidxcolsubsys_bar);
#ifdef WITH_OPENMP_
#pragma omp parallel for
#endif // WITH_OPENMP_
            for (idx i = 0; i < Dsubsys_bar; ++i) {
                result(i, j) = worker(i);
            }
        }

        return result;
    }
    //************ Exception: not ket nor density matrix ************//
    else
        throw exception::MatrixNotSquareNorCvector("qpp::ptrace()");
}

/**
* \brief Partial trace
* \see qpp::ptrace1(), qpp::ptrace2()
*
*  Partial trace of the multi-partite state vector or density matrix
*  over a list of subsystems
*
* \param A Eigen expression
* \param subsys Subsystem indexes
* \param d Subsystem dimensions
* \return Partial trace \f$Tr_{subsys}(\cdot)\f$ over the subsytems \a subsys
* in a multi-partite system, as a dynamic matrix
* over the same scalar field as \a A
*/
template <typename Derived>
dyn_mat<typename Derived::Scalar> ptrace(const Eigen::MatrixBase<Derived>& A,
                                         const std::vector<idx>& subsys,
                                         idx d = 2) {
    const typename Eigen::MatrixBase<Derived>::EvalReturnType& rA = A.derived();

    // EXCEPTION CHECKS

    // check zero size
    if (!internal::check_nonzero_size(rA))
        throw exception::ZeroSize("qpp::ptrace()");

    // check valid dims
    if (d < 2)
        throw exception::DimsInvalid("qpp::ptrace()");
    // END EXCEPTION CHECKS

    idx N = internal::get_num_subsys(static_cast<idx>(rA.rows()), d);
    std::vector<idx> dims(N, d); // local dimensions vector

    return ptrace(rA, subsys, dims);
}

/**
* \brief Partial transpose
*
*  Partial transpose of the multi-partite state vector or density matrix
*  over a list of subsystems
*
* \param A Eigen expression
* \param subsys Subsystem indexes
* \param dims Dimensions of the multi-partite system
* \return Partial transpose \f$(\cdot)^{T_{subsys}}\f$
* over the subsytems \a subsys in a multi-partite system, as a dynamic matrix
* over the same scalar field as \a A
*/
template <typename Derived>
dyn_mat<typename Derived::Scalar>
ptranspose(const Eigen::MatrixBase<Derived>& A, const std::vector<idx>& subsys,
           const std::vector<idx>& dims) {
    const typename Eigen::MatrixBase<Derived>::EvalReturnType& rA = A.derived();

    // EXCEPTION CHECKS

    // check zero-size
    if (!internal::check_nonzero_size(rA))
        throw exception::ZeroSize("qpp::ptranspose()");

    // check that dims is a valid dimension vector
    if (!internal::check_dims(dims))
        throw exception::DimsInvalid("qpp::ptranspose()");

    // check that subsys are valid
    if (!internal::check_subsys_match_dims(subsys, dims))
        throw exception::SubsysMismatchDims("qpp::ptranspose()");
    // END EXCEPTION CHECKS

    idx D = static_cast<idx>(rA.rows());
    idx N = dims.size();
    idx Nsubsys = subsys.size();
    idx Cdims[maxn];
    idx Cmidxcol[maxn];
    idx Csubsys[maxn];

    // copy dims in Cdims and subsys in Csubsys
    for (idx i = 0; i < N; ++i)
        Cdims[i] = dims[i];
    for (idx i = 0; i < Nsubsys; ++i)
        Csubsys[i] = subsys[i];

    dyn_mat<typename Derived::Scalar> result(D, D);

    //************ ket ************//
    if (internal::check_cvector(rA)) // we have a ket
    {
        // check that dims match the dimension of A
        if (!internal::check_dims_match_cvect(dims, rA))
            throw exception::DimsMismatchCvector("qpp::ptranspose()");

        if (subsys.size() == dims.size())
            return (rA * adjoint(rA)).transpose();

        if (subsys.size() == 0)
            return rA * adjoint(rA);

        auto worker = [&](idx i) noexcept->typename Derived::Scalar {
            // use static allocation for speed!
            idx midxcoltmp[maxn];
            idx midxrow[maxn];

            for (idx k = 0; k < N; ++k)
                midxcoltmp[k] = Cmidxcol[k];

            /* compute the row multi-index */
            internal::n2multiidx(i, N, Cdims, midxrow);

            for (idx k = 0; k < Nsubsys; ++k)
                std::swap(midxcoltmp[Csubsys[k]], midxrow[Csubsys[k]]);

            /* writes the result */
            return rA(internal::multiidx2n(midxrow, N, Cdims)) *
                   std::conj(rA(internal::multiidx2n(midxcoltmp, N, Cdims)));
        }; /* end worker */

        for (idx j = 0; j < D; ++j) {
            // compute the column multi-index
            internal::n2multiidx(j, N, Cdims, Cmidxcol);

#ifdef WITH_OPENMP_
#pragma omp parallel for
#endif // WITH_OPENMP_
            for (idx i = 0; i < D; ++i)
                result(i, j) = worker(i);
        }

        return result;
    }
    //************ density matrix ************//
    else if (internal::check_square_mat(rA)) // we have a density operator
    {
        // check that dims match the dimension of A
        if (!internal::check_dims_match_mat(dims, rA))
            throw exception::DimsMismatchMatrix("qpp::ptranspose()");

        if (subsys.size() == dims.size())
            return rA.transpose();

        if (subsys.size() == 0)
            return rA;

        auto worker = [&](idx i) noexcept->typename Derived::Scalar {
            // use static allocation for speed!
            idx midxcoltmp[maxn];
            idx midxrow[maxn];

            for (idx k = 0; k < N; ++k)
                midxcoltmp[k] = Cmidxcol[k];

            /* compute the row multi-index */
            internal::n2multiidx(i, N, Cdims, midxrow);

            for (idx k = 0; k < Nsubsys; ++k)
                std::swap(midxcoltmp[Csubsys[k]], midxrow[Csubsys[k]]);

            /* writes the result */
            return rA(internal::multiidx2n(midxrow, N, Cdims),
                      internal::multiidx2n(midxcoltmp, N, Cdims));
        }; /* end worker */

        for (idx j = 0; j < D; ++j) {
            // compute the column multi-index
            internal::n2multiidx(j, N, Cdims, Cmidxcol);

#ifdef WITH_OPENMP_
#pragma omp parallel for
#endif // WITH_OPENMP_
            for (idx i = 0; i < D; ++i)
                result(i, j) = worker(i);
        }

        return result;
    }
    //************ Exception: not ket nor density matrix ************//
    else
        throw exception::MatrixNotSquareNorCvector("qpp::ptranspose()");
}

/**
* \brief Partial transpose
*
*  Partial transpose of the multi-partite state vector or density matrix
*  over a list of subsystems
*
* \param A Eigen expression
* \param subsys Subsystem indexes
* \param d Subsystem dimensions
* \return Partial transpose \f$(\cdot)^{T_{subsys}}\f$
* over the subsytems \a subsys in a multi-partite system, as a dynamic matrix
* over the same scalar field as \a A
*/
template <typename Derived>
dyn_mat<typename Derived::Scalar>
ptranspose(const Eigen::MatrixBase<Derived>& A, const std::vector<idx>& subsys,
           idx d = 2) {
    const typename Eigen::MatrixBase<Derived>::EvalReturnType& rA = A.derived();

    // EXCEPTION CHECKS

    // check zero size
    if (!internal::check_nonzero_size(rA))
        throw exception::ZeroSize("qpp::ptranspose()");

    // check valid dims
    if (d < 2)
        throw exception::DimsInvalid("qpp::ptranspose()");
    // END EXCEPTION CHECKS

    idx N = internal::get_num_subsys(static_cast<idx>(rA.rows()), d);
    std::vector<idx> dims(N, d); // local dimensions vector

    return ptranspose(rA, subsys, dims);
}

/**
* \brief Subsystem permutation
*
* Permutes the subsystems of a state vector or density matrix.
* The qubit \a perm[\a i] is permuted to the location \a i.
*
* \param A Eigen expression
* \param perm Permutation
* \param dims Dimensions of the multi-partite system
* \return Permuted system, as a dynamic matrix
* over the same scalar field as \a A
*/
template <typename Derived>
dyn_mat<typename Derived::Scalar>
syspermute(const Eigen::MatrixBase<Derived>& A, const std::vector<idx>& perm,
           const std::vector<idx>& dims) {
    const typename Eigen::MatrixBase<Derived>::EvalReturnType& rA = A.derived();

    // EXCEPTION CHECKS

    // check zero-size
    if (!internal::check_nonzero_size(rA))
        throw exception::ZeroSize("qpp::syspermute()");

    // check that dims is a valid dimension vector
    if (!internal::check_dims(dims))
        throw exception::DimsInvalid("qpp::syspermute()");

    // check that we have a valid permutation
    if (!internal::check_perm(perm))
        throw exception::PermInvalid("qpp::syspermute()");

    // check that permutation match dimensions
    if (perm.size() != dims.size())
        throw exception::PermMismatchDims("qpp::syspermute()");
    // END EXCEPTION CHECKS

    idx D = static_cast<idx>(rA.rows());
    idx N = dims.size();

    dyn_mat<typename Derived::Scalar> result;

    //************ ket ************//
    if (internal::check_cvector(rA)) // we have a column vector
    {
        idx Cdims[maxn];
        idx Cperm[maxn];

        // check that dims match the dimension of rA
        if (!internal::check_dims_match_cvect(dims, rA))
            throw exception::DimsMismatchCvector("qpp::syspermute()");

        // copy dims in Cdims and perm in Cperm
        for (idx i = 0; i < N; ++i) {
            Cdims[i] = dims[i];
            Cperm[i] = perm[i];
        }
        result.resize(D, 1);

        auto worker = [&Cdims, &Cperm, N ](idx i) noexcept->idx {
            // use static allocation for speed,
            // double the size for matrices reshaped as vectors
            idx midx[maxn];
            idx midxtmp[maxn];
            idx permdims[maxn];

            /* compute the multi-index */
            internal::n2multiidx(i, N, Cdims, midx);

            for (idx k = 0; k < N; ++k) {
                permdims[k] = Cdims[Cperm[k]]; // permuted dimensions
                midxtmp[k] = midx[Cperm[k]];   // permuted multi-indexes
            }

            return internal::multiidx2n(midxtmp, N, permdims);
        }; /* end worker */

#ifdef WITH_OPENMP_
#pragma omp parallel for
#endif // WITH_OPENMP_
        for (idx i = 0; i < D; ++i)
            result(worker(i)) = rA(i);

        return result;
    }
    //************ density matrix ************//
    else if (internal::check_square_mat(rA)) // we have a density operator
    {
        idx Cdims[2 * maxn];
        idx Cperm[2 * maxn];

        // check that dims match the dimension of rA
        if (!internal::check_dims_match_mat(dims, rA))
            throw exception::DimsMismatchMatrix("qpp::syspermute()");

        // copy dims in Cdims and perm in Cperm
        for (idx i = 0; i < N; ++i) {
            Cdims[i] = dims[i];
            Cdims[i + N] = dims[i];
            Cperm[i] = perm[i];
            Cperm[i + N] = perm[i] + N;
        }
        result.resize(D * D, 1);
        // map A to a column vector
        dyn_mat<typename Derived::Scalar> vectA =
            Eigen::Map<dyn_mat<typename Derived::Scalar>>(
                const_cast<typename Derived::Scalar*>(rA.data()), D * D, 1);

        auto worker = [&Cdims, &Cperm, N ](idx i) noexcept->idx {
            // use static allocation for speed,
            // double the size for matrices reshaped as vectors
            idx midx[2 * maxn];
            idx midxtmp[2 * maxn];
            idx permdims[2 * maxn];

            /* compute the multi-index */
            internal::n2multiidx(i, 2 * N, Cdims, midx);

            for (idx k = 0; k < 2 * N; ++k) {
                permdims[k] = Cdims[Cperm[k]]; // permuted dimensions
                midxtmp[k] = midx[Cperm[k]];   // permuted multi-indexes
            }

            return internal::multiidx2n(midxtmp, 2 * N, permdims);
        }; /* end worker */

#ifdef WITH_OPENMP_
#pragma omp parallel for
#endif // WITH_OPENMP_
        for (idx i = 0; i < D * D; ++i)
            result(worker(i)) = rA(i);

        return reshape(result, D, D);
    }
    //************ Exception: not ket nor density matrix ************//
    else
        throw exception::MatrixNotSquareNorCvector("qpp::syspermute()");
}

/**
* \brief Subsystem permutation
*
* Permutes the subsystems of a state vector or density matrix.
* The qubit \a perm[\a i] is permuted to the location \a i.
*
* \param A Eigen expression
* \param perm Permutation
* \param d Subsystem dimensions
* \return Permuted system, as a dynamic matrix
* over the same scalar field as \a A
*/
template <typename Derived>
dyn_mat<typename Derived::Scalar>
syspermute(const Eigen::MatrixBase<Derived>& A, const std::vector<idx>& perm,
           idx d = 2) {
    const typename Eigen::MatrixBase<Derived>::EvalReturnType& rA = A.derived();

    // EXCEPTION CHECKS

    // check zero size
    if (!internal::check_nonzero_size(rA))
        throw exception::ZeroSize("qpp::syspermute()");

    // check valid dims
    if (d < 2)
        throw exception::DimsInvalid("qpp::syspermute()");
    // END EXCEPTION CHECKS

    idx N = internal::get_num_subsys(static_cast<idx>(rA.rows()), d);
    std::vector<idx> dims(N, d); // local dimensions vector

    return syspermute(rA, perm, dims);
}

} /* namespace qpp */

#endif /* OPERATIONS_H_ */
