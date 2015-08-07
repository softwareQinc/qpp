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
* \file operations.h
* \brief Quantum operation functions
*/

#ifndef OPERATIONS_H_
#define OPERATIONS_H_

//// silence g++4.8 bogus warning -Wunused-but-set-variable in lambda functions
#if (__GNUC__ && !__clang__)
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#endif

namespace qpp
{
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
template<typename Derived1, typename Derived2>
dyn_mat<typename Derived1::Scalar> applyCTRL(
        const Eigen::MatrixBase<Derived1>& state,
        const Eigen::MatrixBase<Derived2>& A,
        const std::vector<idx>& ctrl,
        const std::vector<idx>& subsys,
        const std::vector<idx>& dims)
{
    const dyn_mat<typename Derived1::Scalar>& rstate = state;
    const dyn_mat<typename Derived2::Scalar>& rA = A;

    // EXCEPTION CHECKS
    // check types
    if (!std::is_same<typename Derived1::Scalar,
            typename Derived2::Scalar>::value)
        throw Exception("qpp::applyCTRL()", Exception::Type::TYPE_MISMATCH);

    // check zero sizes
    if (!internal::_check_nonzero_size(rA))
        throw Exception("qpp::applyCTRL()", Exception::Type::ZERO_SIZE);

    // check zero sizes
    if (!internal::_check_nonzero_size(rstate))
        throw Exception("qpp::applyCTRL()", Exception::Type::ZERO_SIZE);

    // check square matrix for the gate
    if (!internal::_check_square_mat(rA))
        throw Exception("qpp::applyCTRL()",
                        Exception::Type::MATRIX_NOT_SQUARE);

    // check that all control subsystems have the same dimension
    idx d = ctrl.size() > 0 ? dims[ctrl[0]] : 1;
    for (idx i = 1; i < ctrl.size(); ++i)
        if (dims[ctrl[i]] != d)
            throw Exception("qpp::applyCTRL()",
                            Exception::Type::DIMS_NOT_EQUAL);

    // check that dimension is valid
    if (!internal::_check_dims(dims))
        throw Exception("qpp::applyCTRL()", Exception::Type::DIMS_INVALID);

    // check subsys is valid w.r.t. dims
    if (!internal::_check_subsys_match_dims(subsys, dims))
        throw Exception("qpp::applyCTRL()",
                        Exception::Type::SUBSYS_MISMATCH_DIMS);

    // check that gate matches the dimensions of the subsys
    std::vector<idx> subsys_dims(subsys.size());
    for (idx i = 0; i < subsys.size(); ++i)
        subsys_dims[i] = dims[subsys[i]];
    if (!internal::_check_dims_match_mat(subsys_dims, rA))
        throw Exception("qpp::applyCTRL()",
                        Exception::Type::MATRIX_MISMATCH_SUBSYS);

    std::vector<idx> ctrlgate = ctrl; // ctrl + gate subsystem vector
    ctrlgate.insert(std::end(ctrlgate), std::begin(subsys), std::end(subsys));
    std::sort(std::begin(ctrlgate), std::end(ctrlgate));

    // check that ctrl + gate subsystem is valid
    // with respect to local dimensions
    if (!internal::_check_subsys_match_dims(ctrlgate, dims))
        throw Exception("qpp::applyCTRL()",
                        Exception::Type::SUBSYS_MISMATCH_DIMS);

    // END EXCEPTION CHECKS

    // construct the table of A^i and (A^dagger)^i
    std::vector<dyn_mat<typename Derived1::Scalar>> Ai;
    std::vector<dyn_mat<typename Derived1::Scalar>> Aidagger;
    for (idx i = 0; i < std::max(d, static_cast<idx>(2)); ++i)
    {
        Ai.push_back(powm(rA, i));
        Aidagger.push_back(powm(adjoint(rA), i));
    }

    idx D = rstate.rows(); // total dimension
    idx n = dims.size();   // total number of subsystems
    idx ctrlsize = ctrl.size(); // dimension of ctrl subsystem
    idx DA = rA.rows(); // dimension of gate subsystem

    idx Cdims[maxn]; // local dimensions
    idx CdimsA[maxn]; // local dimensions
    idx CdimsCTRLAbar[maxn]; // local dimensions

    // compute the complementary subsystem of ctrlgate w.r.t. dims
    std::vector<idx> ctrlgatebar = complement(ctrlgate, n);

    idx DCTRLAbar = 1; // dimension of the rest
    for (idx i = 0; i < ctrlgatebar.size(); ++i)
        DCTRLAbar *= dims[ctrlgatebar[i]];

    for (idx k = 0; k < n; ++k)
        Cdims[k] = dims[k];
    for (idx k = 0; k < subsys.size(); ++k)
        CdimsA[k] = dims[subsys[k]];
    for (idx k = 0; k < ctrlgatebar.size(); ++k)
        CdimsCTRLAbar[k] = dims[ctrlgatebar[k]];


    // worker, computes the coefficient and the index for the ket case
    // used in #pragma omp parallel for collapse
    auto coeff_idx_ket = [=](idx _i, idx _m, idx _r) noexcept
            -> std::pair<typename Derived1::Scalar, idx>
    {
        idx indx = 0;
        typename Derived1::Scalar coeff = 0;

        idx Cmidx[maxn]; // the total multi-index
        idx CmidxA[maxn];// the gate part multi-index
        idx CmidxCTRLAbar[maxn];// the rest multi-index

        // compute the index

        // set the CTRL part
        for (idx k = 0; k < ctrl.size(); ++k)
        {
            Cmidx[ctrl[k]] = _i;
        }

        // set the rest
        internal::_n2multiidx(_r, n - ctrlgate.size(),
                              CdimsCTRLAbar, CmidxCTRLAbar);
        for (idx k = 0; k < n - ctrlgate.size(); ++k)
        {
            Cmidx[ctrlgatebar[k]] = CmidxCTRLAbar[k];
        }

        // set the A part
        internal::_n2multiidx(_m, subsys.size(), CdimsA, CmidxA);
        for (idx k = 0; k < subsys.size(); ++k)
        {
            Cmidx[subsys[k]] = CmidxA[k];
        }

        // we now got the total index
        indx = internal::_multiidx2n(Cmidx, n, Cdims);

        // compute the coefficient
        for (idx _n = 0; _n < DA; ++_n)
        {
            internal::_n2multiidx(_n, subsys.size(), CdimsA, CmidxA);
            for (idx k = 0; k < subsys.size(); ++k)
            {
                Cmidx[subsys[k]] = CmidxA[k];
            }
            coeff += Ai[_i](_m, _n) *
                     rstate(internal::_multiidx2n(Cmidx, n, Cdims));
        }

        return std::make_pair(coeff, indx);
    };

    // worker, computes the coefficient and the index
    // for the density matrix case
    // used in #pragma omp parallel for collapse
    auto coeff_idx_rho = [=](idx _i1, idx _m1,
                             idx _r1, idx _i2, idx _m2, idx _r2) noexcept
            -> std::tuple<typename Derived1::Scalar, idx, idx>
    {
        idx idxrow = 0;
        idx idxcol = 0;
        typename Derived1::Scalar coeff = 0;

        idx Cmidxrow[maxn]; // the total row multi-index
        idx Cmidxcol[maxn];// the total col multi-index
        idx CmidxArow[maxn];// the gate part row multi-index
        idx CmidxAcol[maxn];// the gate part col multi-index
        idx CmidxCTRLAbarrow[maxn];// the rest row multi-index
        idx CmidxCTRLAbarcol[maxn];// the rest col multi-index

        // compute the ket/bra indexes

        // set the CTRL part
        for (idx k = 0; k < ctrl.size(); ++k)
        {
            Cmidxrow[ctrl[k]] = _i1;
            Cmidxcol[ctrl[k]] = _i2;
        }

        // set the rest
        internal::_n2multiidx(_r1, n - ctrlgate.size(),
                              CdimsCTRLAbar, CmidxCTRLAbarrow);
        internal::_n2multiidx(_r2, n - ctrlgate.size(),
                              CdimsCTRLAbar, CmidxCTRLAbarcol);
        for (idx k = 0; k < n - ctrlgate.size(); ++k)
        {
            Cmidxrow[ctrlgatebar[k]] = CmidxCTRLAbarrow[k];
            Cmidxcol[ctrlgatebar[k]] = CmidxCTRLAbarcol[k];
        }

        // set the A part
        internal::_n2multiidx(_m1, subsys.size(), CdimsA, CmidxArow);
        internal::_n2multiidx(_m2, subsys.size(), CdimsA, CmidxAcol);
        for (idx k = 0; k < subsys.size(); ++k)
        {
            Cmidxrow[subsys[k]] = CmidxArow[k];
            Cmidxcol[subsys[k]] = CmidxAcol[k];
        }

        // we now got the total row/col indexes
        idxrow = internal::_multiidx2n(Cmidxrow, n, Cdims);
        idxcol = internal::_multiidx2n(Cmidxcol, n, Cdims);

        // compute the coefficient
        for (idx _n1 = 0; _n1 < DA; ++_n1)
        {
            internal::_n2multiidx(_n1, subsys.size(), CdimsA, CmidxArow);
            for (idx k = 0; k < subsys.size(); ++k)
            {
                Cmidxrow[subsys[k]] = CmidxArow[k];
            }
            for (idx _n2 = 0; _n2 < DA; ++_n2)
            {
                internal::_n2multiidx(_n2, subsys.size(), CdimsA, CmidxAcol);
                for (idx k = 0; k < subsys.size(); ++k)
                {
                    Cmidxcol[subsys[k]] = CmidxAcol[k];
                }
                coeff += Ai[_i1](_m1, _n1) *
                         rstate(internal::_multiidx2n(Cmidxrow, n, Cdims),
                                internal::_multiidx2n(Cmidxcol, n, Cdims)) *
                         Aidagger[_i2](_n2, _m2);
            }
        }

        return std::make_tuple(coeff, idxrow, idxcol);
    };

    //************ ket ************//
    if (internal::_check_cvector(rstate)) // we have a ket
    {
        // check that dims match state vector
        if (!internal::_check_dims_match_cvect(dims, rstate))
            throw Exception("qpp::applyCTRL()",
                            Exception::Type::DIMS_MISMATCH_CVECTOR);
        if (D == 1)
            return rstate;

        dyn_mat<typename Derived1::Scalar> result = rstate;

#pragma omp parallel for collapse(2)
        for (idx m = 0; m < DA; ++m)
            for (idx r = 0; r < DCTRLAbar; ++r)
            {
                if (ctrlsize == 0) // no control
                {
                    result(coeff_idx_ket(1, m, r).second) =
                            coeff_idx_ket(1, m, r).first;
                }
                else
                    for (idx i = 0; i < d; ++i)
                    {
                        result(coeff_idx_ket(i, m, r).second) =
                                coeff_idx_ket(i, m, r).first;
                    }
            }

        return result;
    }
        //************ density matrix ************//
    else if (internal::_check_square_mat(rstate)) // we have a density operator
    {
        // check that dims match state matrix
        if (!internal::_check_dims_match_mat(dims, rstate))
            throw Exception("qpp::applyCTRL()",
                            Exception::Type::DIMS_MISMATCH_MATRIX);

        if (D == 1)
            return rstate;

        dyn_mat<typename Derived1::Scalar> result = rstate;

#pragma omp parallel for collapse(4)
        for (idx m1 = 0; m1 < DA; ++m1)
            for (idx r1 = 0; r1 < DCTRLAbar; ++r1)
                for (idx m2 = 0; m2 < DA; ++m2)
                    for (idx r2 = 0; r2 < DCTRLAbar; ++r2)
                        if (ctrlsize == 0) // no control
                        {
                            auto coeff_idxes = coeff_idx_rho(1, m1, r1,
                                                             1, m2, r2);
                            result(std::get<1>(coeff_idxes),
                                   std::get<2>(coeff_idxes)) =
                                    std::get<0>(coeff_idxes);
                        }
                        else
                        {
                            for (idx i1 = 0; i1 < d; ++i1)
                                for (idx i2 = 0; i2 < d; ++i2)
                                {
                                    auto coeff_idxes = coeff_idx_rho(
                                            i1, m1, r1,
                                            i2, m2, r2);
                                    result(std::get<1>(coeff_idxes),
                                           std::get<2>(coeff_idxes)) =
                                            std::get<0>(coeff_idxes);
                                }
                        }

        return result;
    }
        //************ Exception: not ket nor density matrix ************//
    else
        throw Exception("qpp::applyCTRL()",
                        Exception::Type::MATRIX_NOT_SQUARE_OR_CVECTOR);
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
template<typename Derived1, typename Derived2>
dyn_mat<typename Derived1::Scalar> applyCTRL(
        const Eigen::MatrixBase<Derived1>& state,
        const Eigen::MatrixBase<Derived2>& A,
        const std::vector<idx>& ctrl,
        const std::vector<idx>& subsys,
        idx d = 2)
{
    const dyn_mat<typename Derived1::Scalar>& rstate = state;
    const dyn_mat<typename Derived1::Scalar>& rA = A;

    // check zero size
    if (!internal::_check_nonzero_size(rstate))
        throw Exception("qpp::applyCTRL()", Exception::Type::ZERO_SIZE);

    idx n =
            static_cast<idx>(std::llround(std::log2(rstate.rows()) /
                                          std::log2(d)));
    std::vector<idx> dims(n, d); // local dimensions vector

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
template<typename Derived1, typename Derived2>
dyn_mat<typename Derived1::Scalar> apply(
        const Eigen::MatrixBase<Derived1>& state,
        const Eigen::MatrixBase<Derived2>& A,
        const std::vector<idx>& subsys,
        const std::vector<idx>& dims)
{
    const dyn_mat<typename Derived1::Scalar>& rstate = state;
    const dyn_mat<typename Derived2::Scalar>& rA = A;

    // EXCEPTION CHECKS

    // check types
    if (!std::is_same<typename Derived1::Scalar,
            typename Derived2::Scalar>::value)
        throw Exception("qpp::apply()", Exception::Type::TYPE_MISMATCH);

    // check zero sizes
    if (!internal::_check_nonzero_size(rA))
        throw Exception("qpp::apply()", Exception::Type::ZERO_SIZE);

    // check zero sizes
    if (!internal::_check_nonzero_size(rstate))
        throw Exception("qpp::apply()", Exception::Type::ZERO_SIZE);

    // check square matrix for the gate
    if (!internal::_check_square_mat(rA))
        throw Exception("qpp::apply()", Exception::Type::MATRIX_NOT_SQUARE);

    // check that dimension is valid
    if (!internal::_check_dims(dims))
        throw Exception("qpp::apply()", Exception::Type::DIMS_INVALID);

    // check subsys is valid w.r.t. dims
    if (!internal::_check_subsys_match_dims(subsys, dims))
        throw Exception("qpp::apply()", Exception::Type::SUBSYS_MISMATCH_DIMS);

    // check that gate matches the dimensions of the subsys
    std::vector<idx> subsys_dims(subsys.size());
    for (idx i = 0; i < subsys.size(); ++i)
        subsys_dims[i] = dims[subsys[i]];
    if (!internal::_check_dims_match_mat(subsys_dims, rA))
        throw Exception("qpp::apply()",
                        Exception::Type::MATRIX_MISMATCH_SUBSYS);

    //************ ket ************//
    if (internal::_check_cvector(rstate)) // we have a ket
    {
        // check that dims match state vector
        if (!internal::_check_dims_match_cvect(dims, rstate))
            throw Exception("qpp::apply()",
                            Exception::Type::DIMS_MISMATCH_CVECTOR);

        return applyCTRL(rstate, rA, {}, subsys, dims);
    }
        //************ density matrix ************//
    else if (internal::_check_square_mat(rstate)) // we have a density operator
    {

        // check that dims match state matrix
        if (!internal::_check_dims_match_mat(dims, rstate))
            throw Exception("qpp::apply()",
                            Exception::Type::DIMS_MISMATCH_MATRIX);

        return applyCTRL(rstate, rA, {}, subsys, dims);
    }
        //************ Exception: not ket nor density matrix ************//
    else
        throw Exception("qpp::apply()",
                        Exception::Type::MATRIX_NOT_SQUARE_OR_CVECTOR);
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
template<typename Derived1, typename Derived2>
dyn_mat<typename Derived1::Scalar> apply(
        const Eigen::MatrixBase<Derived1>& state,
        const Eigen::MatrixBase<Derived2>& A,
        const std::vector<idx>& subsys,
        idx d = 2)
{
    const dyn_mat<typename Derived1::Scalar>& rstate = state;
    const dyn_mat<typename Derived1::Scalar>& rA = A;

    // check zero size
    if (!internal::_check_nonzero_size(rstate))
        throw Exception("qpp::apply()", Exception::Type::ZERO_SIZE);

    idx n =
            static_cast<idx>(std::llround(std::log2(rstate.rows()) /
                                          std::log2(d)));
    std::vector<idx> dims(n, d); // local dimensions vector

    return apply(rstate, rA, subsys, dims);
}

/**
* \brief Applies the channel specified by the set of Kraus operators \a Ks
* to the density matrix \a rho
*
* \param rho Eigen expression
* \param Ks Set of Kraus operators
* \return Output density matrix after the action of the channel
*/
template<typename Derived>
cmat apply(const Eigen::MatrixBase<Derived>& rho,
           const std::vector<cmat>& Ks)
{
    const cmat& rrho = rho;

    // EXCEPTION CHECKS
    if (!internal::_check_nonzero_size(rrho))
        throw Exception("qpp::apply()", Exception::Type::ZERO_SIZE);
    if (!internal::_check_square_mat(rrho))
        throw Exception("qpp::apply()", Exception::Type::MATRIX_NOT_SQUARE);
    if (Ks.size() == 0)
        throw Exception("qpp::apply()", Exception::Type::ZERO_SIZE);
    if (!internal::_check_square_mat(Ks[0]))
        throw Exception("qpp::apply()", Exception::Type::MATRIX_NOT_SQUARE);
    if (Ks[0].rows() != rrho.rows())
        throw Exception("qpp::apply()",
                        Exception::Type::DIMS_MISMATCH_MATRIX);
    for (auto&& it : Ks)
        if (it.rows() != Ks[0].rows() || it.cols() != Ks[0].rows())
            throw Exception("qpp::apply()", Exception::Type::DIMS_NOT_EQUAL);

    cmat result = cmat::Zero(rrho.rows(), rrho.rows());

#pragma omp parallel for
    for (idx i = 0; i < Ks.size(); ++i)
    {
#pragma omp critical
        {
            result += Ks[i] * rrho * adjoint(Ks[i]);
        }
    }

    return result;
}

/**
* \brief Applies the channel specified by the set of Kraus operators \a Ks to
* the part \a subsys of the multi-partite density matrix \a rho
*
* \param rho Eigen expression
* \param Ks Set of Kraus operators
* \param subsys Subsystem indexes where the Kraus operators \a Ks are applied
* \param dims Dimensions of the multi-partite system
* \return Output density matrix after the action of the channel
*/
template<typename Derived>
cmat apply(const Eigen::MatrixBase<Derived>& rho,
           const std::vector<cmat>& Ks,
           const std::vector<idx>& subsys,
           const std::vector<idx>& dims)
{
    const cmat& rrho = rho;

    // EXCEPTION CHECKS
    // check zero sizes
    if (!internal::_check_nonzero_size(rrho))
        throw Exception("qpp::apply()", Exception::Type::ZERO_SIZE);

    // check square matrix for the rho
    if (!internal::_check_square_mat(rrho))
        throw Exception("qpp::apply()", Exception::Type::MATRIX_NOT_SQUARE);

    // check that dimension is valid
    if (!internal::_check_dims(dims))
        throw Exception("qpp::apply()", Exception::Type::DIMS_INVALID);

    // check that dims match rho matrix
    if (!internal::_check_dims_match_mat(dims, rrho))
        throw Exception("qpp::apply()",
                        Exception::Type::DIMS_MISMATCH_MATRIX);

    // check subsys is valid w.r.t. dims
    if (!internal::_check_subsys_match_dims(subsys, dims))
        throw Exception("qpp::apply()",
                        Exception::Type::SUBSYS_MISMATCH_DIMS);

    std::vector<idx> subsys_dims(subsys.size());
    for (idx i = 0; i < subsys.size(); ++i)
        subsys_dims[i] = dims[subsys[i]];

    // check the Kraus operators
    if (Ks.size() == 0)
        throw Exception("qpp::apply()", Exception::Type::ZERO_SIZE);
    if (!internal::_check_square_mat(Ks[0]))
        throw Exception("qpp::apply()", Exception::Type::MATRIX_NOT_SQUARE);
    if (!internal::_check_dims_match_mat(subsys_dims, Ks[0]))
        throw Exception("qpp::apply()",
                        Exception::Type::MATRIX_MISMATCH_SUBSYS);
    for (auto&& it : Ks)
        if (it.rows() != Ks[0].rows() || it.cols() != Ks[0].rows())
            throw Exception("qpp::apply()", Exception::Type::DIMS_NOT_EQUAL);

    cmat result = cmat::Zero(rrho.rows(), rrho.rows());

    for (idx i = 0; i < Ks.size(); ++i)
        result += apply(rrho, Ks[i], subsys, dims);

    return result;
}

/**
* \brief Applies the channel specified by the set of Kraus operators \a Ks to
* the part \a subsys of the multi-partite density matrix \a rho
*
* \param rho Eigen expression
* \param Ks Set of Kraus operators
* \param subsys Subsystem indexes where the Kraus operators \a Ks are applied
* \param d Subsystem dimensions
* \return Output density matrix after the action of the channel
*/
template<typename Derived>
cmat apply(const Eigen::MatrixBase<Derived>& rho,
           const std::vector<cmat>& Ks,
           const std::vector<idx>& subsys,
           idx d = 2)
{
    const cmat& rrho = rho;

    // check zero sizes
    if (!internal::_check_nonzero_size(rrho))
        throw Exception("qpp::apply()", Exception::Type::ZERO_SIZE);

    idx n =
            static_cast<idx>(std::llround(std::log2(rrho.rows()) /
                                          std::log2(d)));
    std::vector<idx> dims(n, d); // local dimensions vector

    return apply(rrho, Ks, subsys, dims);
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
inline cmat kraus2super(const std::vector<cmat>& Ks)
{
    // EXCEPTION CHECKS
    if (Ks.size() == 0)
        throw Exception("qpp::kraus2super()", Exception::Type::ZERO_SIZE);
    if (!internal::_check_nonzero_size(Ks[0]))
        throw Exception("qpp::kraus2super()", Exception::Type::ZERO_SIZE);
    if (!internal::_check_square_mat(Ks[0]))
        throw Exception("qpp::kraus2super()",
                        Exception::Type::MATRIX_NOT_SQUARE);
    for (auto&& it : Ks)
        if (it.rows() != Ks[0].rows() || it.cols() != Ks[0].rows())
            throw Exception("qpp::kraus2super()",
                            Exception::Type::DIMS_NOT_EQUAL);
    idx D = static_cast<idx>(Ks[0].rows());

    cmat result(D * D, D * D);
    cmat MN = cmat::Zero(D, D);
    bra A = bra::Zero(D);
    ket B = ket::Zero(D);
    cmat EMN = cmat::Zero(D, D);

#pragma omp parallel for collapse(2)
    for (idx m = 0; m < D; ++m)
    {
        for (idx n = 0; n < D; ++n)
        {
#pragma omp critical
            {
                // compute E(|m><n|)
                MN(m, n) = 1;
                for (idx i = 0; i < Ks.size(); ++i)
                    EMN += Ks[i] * MN * adjoint(Ks[i]);
                MN(m, n) = 0;

                for (idx a = 0; a < D; ++a)
                {
                    A(a) = 1;
                    for (idx b = 0; b < D; ++b)
                    {
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
inline cmat kraus2choi(const std::vector<cmat>& Ks)
{
    // EXCEPTION CHECKS
    if (Ks.size() == 0)
        throw Exception("qpp::kraus2choi()", Exception::Type::ZERO_SIZE);
    if (!internal::_check_nonzero_size(Ks[0]))
        throw Exception("qpp::kraus2choi()", Exception::Type::ZERO_SIZE);
    if (!internal::_check_square_mat(Ks[0]))
        throw Exception("qpp::kraus2choi()",
                        Exception::Type::MATRIX_NOT_SQUARE);
    for (auto&& it : Ks)
        if (it.rows() != Ks[0].rows() || it.cols() != Ks[0].rows())
            throw Exception("qpp::kraus2choi()",
                            Exception::Type::DIMS_NOT_EQUAL);
    idx D = static_cast<idx>(Ks[0].rows());

    // construct the D x D \sum |jj> vector
    // (un-normalized maximally entangled state)
    cmat MES = cmat::Zero(D * D, 1);
    for (idx a = 0; a < D; ++a)
        MES(a * D + a) = 1;

    cmat Omega = MES * adjoint(MES);

    cmat result = cmat::Zero(D * D, D * D);

#pragma omp parallel for
    for (idx i = 0; i < Ks.size(); ++i)
    {
#pragma omp critical
        {
            result += kron(cmat::Identity(D, D), Ks[i]) * Omega
                      * adjoint(kron(cmat::Identity(D, D), Ks[i]));
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
inline std::vector<cmat> choi2kraus(const cmat& A)
{
    // EXCEPTION CHECKS
    if (!internal::_check_nonzero_size(A))
        throw Exception("qpp::choi2kraus()", Exception::Type::ZERO_SIZE);
    if (!internal::_check_square_mat(A))
        throw Exception("qpp::choi2kraus()",
                        Exception::Type::MATRIX_NOT_SQUARE);
    idx D = static_cast<idx>(std::llround(
            std::sqrt(static_cast<double>(A.rows()))));
    if (D * D != static_cast<idx>(A.rows()))
        throw Exception("qpp::choi2kraus()", Exception::Type::DIMS_INVALID);

    dmat ev = hevals(A);
    cmat evec = hevects(A);
    std::vector<cmat> result;

    for (idx i = 0; i < D * D; ++i)
    {
        if (std::abs(ev(i)) > eps)
            result.push_back(
                    std::sqrt(std::abs(ev(i))) * reshape(evec.col(i), D, D));
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
inline cmat choi2super(const cmat& A)
{
    // EXCEPTION CHECKS
    if (!internal::_check_nonzero_size(A))
        throw Exception("qpp::choi2super()", Exception::Type::ZERO_SIZE);
    if (!internal::_check_square_mat(A))
        throw Exception("qpp::choi2super()",
                        Exception::Type::MATRIX_NOT_SQUARE);
    idx D = static_cast<idx>(std::llround(
            std::sqrt(static_cast<double>(A.rows()))));
    if (D * D != static_cast<idx>(A.rows()))
        throw Exception("qpp::choi2super()", Exception::Type::DIMS_INVALID);

    cmat result(D * D, D * D);

#pragma omp parallel for collapse(4)
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
inline cmat super2choi(const cmat& A)
{
    // EXCEPTION CHECKS
    if (!internal::_check_nonzero_size(A))
        throw Exception("qpp::super2choi()", Exception::Type::ZERO_SIZE);
    if (!internal::_check_square_mat(A))
        throw Exception("qpp::super2choi()",
                        Exception::Type::MATRIX_NOT_SQUARE);
    idx D = static_cast<idx>(std::llround(
            std::sqrt(static_cast<double>(A.rows()))));
    if (D * D != static_cast<idx>(A.rows()))
        throw Exception("qpp::super2choi()", Exception::Type::DIMS_INVALID);

    cmat result(D * D, D * D);

#pragma omp parallel for collapse(4)
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
*  Partial trace over the first subsystem
*  of bi-partite state vector or density matrix
*
* \param A Eigen expression
* \param dims Dimensions of the bi-partite system
* (must be a std::vector with 2 elements)
* \return Partial trace \f$Tr_{A}(\cdot)\f$ over the first subsytem \f$A\f$
* in a bi-partite system \f$A\otimes B\f$, as a dynamic matrix
* over the same scalar field as \a A
*/
template<typename Derived>
dyn_mat<typename Derived::Scalar> ptrace1(const Eigen::MatrixBase<Derived>& A,
                                          const std::vector<idx>& dims)
{
    const dyn_mat<typename Derived::Scalar>& rA = A;

    // Error checks

    // check zero-size
    if (!internal::_check_nonzero_size(rA))
        throw Exception("qpp::ptrace1()", Exception::Type::ZERO_SIZE);

    // check that dims is a valid dimension vector
    if (!internal::_check_dims(dims))
        throw Exception("qpp::ptrace1()", Exception::Type::DIMS_INVALID);

    // check dims has only 2 elements
    if (dims.size() != 2)
        throw Exception("qpp::ptrace1()", Exception::Type::NOT_BIPARTITE);

    idx DA = dims[0];
    idx DB = dims[1];

    dyn_mat<typename Derived::Scalar> result =
            dyn_mat < typename Derived::Scalar > ::Zero(DB, DB);

    //************ ket ************//
    if (internal::_check_cvector(rA)) // we have a ket
    {
        // check that dims match the dimension of A
        if (!internal::_check_dims_match_cvect(dims, rA))
            throw Exception("qpp::ptrace1()",
                            Exception::Type::DIMS_MISMATCH_CVECTOR);

        auto worker = [=](idx i, idx j) noexcept
                -> typename Derived::Scalar
        {
            typename Derived::Scalar sum = 0;
            for (idx m = 0; m < DA; ++m)
                sum += rA(m * DB + i) * std::conj(rA(m * DB + j));

            return sum;
        };

#pragma omp parallel for collapse(2)
        for (idx j = 0; j < DB; ++j) // column major order for speed
            for (idx i = 0; i < DB; ++i)
                result(i, j) = worker(i, j);

        return result;
    }
        //************ density matrix ************//
    else if (internal::_check_square_mat(rA)) // we have a density operator
    {
        // check that dims match the dimension of A
        if (!internal::_check_dims_match_mat(dims, rA))
            throw Exception("qpp::ptrace1()",
                            Exception::Type::DIMS_MISMATCH_MATRIX);

        auto worker = [=](idx i, idx j) noexcept
                -> typename Derived::Scalar
        {
            typename Derived::Scalar sum = 0;
            for (idx m = 0; m < DA; ++m)
                sum += rA(m * DB + i, m * DB + j);

            return sum;
        };

#pragma omp parallel for collapse(2)
        for (idx j = 0; j < DB; ++j) // column major order for speed
            for (idx i = 0; i < DB; ++i)
                result(i, j) = worker(i, j);

        return result;
    }
        //************ Exception: not ket nor density matrix ************//
    else
        throw Exception("qpp::ptrace1()",
                        Exception::Type::MATRIX_NOT_SQUARE_OR_CVECTOR);
}

/**
* \brief Partial trace
* \see qpp::ptrace1()
*
*  Partial trace over the second subsystem
*  of bi-partite state vector or density matrix
*
* \param A Eigen expression
* \param dims Dimensions of the bi-partite system
* (must be a std::vector with 2 elements)
* \return Partial trace \f$Tr_{B}(\cdot)\f$ over the second subsytem \f$B\f$
* in a bi-partite system \f$A\otimes B\f$, as a dynamic matrix
* over the same scalar field as \a A
*/
template<typename Derived>
dyn_mat<typename Derived::Scalar> ptrace2(const Eigen::MatrixBase<Derived>& A,
                                          const std::vector<idx>& dims)
{
    const dyn_mat<typename Derived::Scalar>& rA = A;

    // Error checks

    // check zero-size
    if (!internal::_check_nonzero_size(rA))
        throw Exception("qpp::ptrace2()", Exception::Type::ZERO_SIZE);

    // check that dims is a valid dimension vector
    if (!internal::_check_dims(dims))
        throw Exception("qpp::ptrace2()", Exception::Type::DIMS_INVALID);

    // check dims has only 2 elements
    if (dims.size() != 2)
        throw Exception("qpp::ptrace2()", Exception::Type::NOT_BIPARTITE);

    idx DA = dims[0];
    idx DB = dims[1];

    dyn_mat<typename Derived::Scalar> result =
            dyn_mat < typename Derived::Scalar > ::Zero(DA, DA);

    //************ ket ************//
    if (internal::_check_cvector(rA)) // we have a ket
    {
        // check that dims match the dimension of A
        if (!internal::_check_dims_match_cvect(dims, rA))
            throw Exception("qpp::ptrace2()",
                            Exception::Type::DIMS_MISMATCH_CVECTOR);

        auto worker = [=](idx i, idx j) noexcept
                -> typename Derived::Scalar
        {
            typename Derived::Scalar sum = 0;
            for (idx m = 0; m < DB; ++m)
                sum += rA(i * DB + m) * std::conj(rA(j * DB + m));

            return sum;
        };

#pragma omp parallel for collapse(2)
        for (idx j = 0; j < DA; ++j) // column major order for speed
            for (idx i = 0; i < DA; ++i)
                result(i, j) = worker(i, j);

        return result;
    }
        //************ density matrix ************//
    else if (internal::_check_square_mat(rA)) // we have a density operator
    {
        // check that dims match the dimension of A
        if (!internal::_check_dims_match_mat(dims, rA))
            throw Exception("qpp::ptrace2()",
                            Exception::Type::DIMS_MISMATCH_MATRIX);

#pragma omp parallel for collapse(2)
        for (idx j = 0; j < DA; ++j) // column major order for speed
            for (idx i = 0; i < DA; ++i)
                result(i, j) = trace(rA.block(i * DB, j * DB, DB, DB));

        return result;
    }
        //************ Exception: not ket nor density matrix ************//
    else
        throw Exception("qpp::ptrace1()",
                        Exception::Type::MATRIX_NOT_SQUARE_OR_CVECTOR);
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
template<typename Derived>
dyn_mat<typename Derived::Scalar> ptrace(const Eigen::MatrixBase<Derived>& A,
                                         const std::vector<idx>& subsys,
                                         const std::vector<idx>& dims)
{
    const dyn_mat<typename Derived::Scalar>& rA = A;

    // error checks

    // check zero-size
    if (!internal::_check_nonzero_size(rA))
        throw Exception("qpp::ptrace()", Exception::Type::ZERO_SIZE);

    // check that dims is a valid dimension vector
    if (!internal::_check_dims(dims))
        throw Exception("qpp::ptrace()", Exception::Type::DIMS_INVALID);

    // check that subsys are valid
    if (!internal::_check_subsys_match_dims(subsys, dims))
        throw Exception("qpp::ptrace()",
                        Exception::Type::SUBSYS_MISMATCH_DIMS);

    idx D = static_cast<idx>(rA.rows());
    idx n = dims.size();
    idx nsubsys = subsys.size();
    idx nsubsysbar = n - nsubsys;
    idx dimsubsys = 1;
    for (idx i = 0; i < nsubsys; ++i)
        dimsubsys *= dims[subsys[i]];
    idx dimsubsysbar = D / dimsubsys;

    idx Cdims[maxn];
    idx Csubsys[maxn];
    idx Cdimssubsys[maxn];
    idx Csubsysbar[maxn];
    idx Cdimssubsysbar[maxn];

    idx Cmidxcolsubsysbar[maxn];

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

    dyn_mat<typename Derived::Scalar> result =
            dyn_mat < typename Derived::Scalar > (dimsubsysbar, dimsubsysbar);

    //************ ket ************//
    if (internal::_check_cvector(rA)) // we have a ket
    {
        // check that dims match the dimension of A
        if (!internal::_check_dims_match_cvect(dims, rA))
            throw Exception("qpp::ptrace()",
                            Exception::Type::DIMS_MISMATCH_CVECTOR);

        if (subsys.size() == dims.size())
        {
            result(0, 0) = (adjoint(rA) * rA).value();
            return result;
        }

        if (subsys.size() == 0)
            return rA * adjoint(rA);

        auto worker = [=, &Cmidxcolsubsysbar](idx i) noexcept
                -> typename Derived::Scalar
        {
            // use static allocation for speed!

            idx Cmidxrow[maxn];
            idx Cmidxcol[maxn];
            idx Cmidxrowsubsysbar[maxn];
            idx Cmidxsubsys[maxn];

            /* get the row multi-indexes of the complement */
            internal::_n2multiidx(i, nsubsysbar,
                                  Cdimssubsysbar, Cmidxrowsubsysbar);
            /* write them in the global row/col multi-indexes */
            for (idx k = 0; k < nsubsysbar; ++k)
            {
                Cmidxrow[Csubsysbar[k]] = Cmidxrowsubsysbar[k];
                Cmidxcol[Csubsysbar[k]] = Cmidxcolsubsysbar[k];
            }
            typename Derived::Scalar sm = 0;
            for (idx a = 0; a < dimsubsys; ++a)
            {
                // get the multi-index over which we do the summation
                internal::_n2multiidx(a, nsubsys, Cdimssubsys, Cmidxsubsys);
                // write it into the global row/col multi-indexes
                for (idx k = 0; k < nsubsys; ++k)
                    Cmidxrow[Csubsys[k]] = Cmidxcol[Csubsys[k]]
                            = Cmidxsubsys[k];

                // now do the sum
                sm += rA(internal::_multiidx2n(Cmidxrow, n, Cdims)) *
                      std::conj(rA(internal::_multiidx2n(Cmidxcol, n,
                                                         Cdims)));
            }

            return sm;
        };

        for (idx j = 0; j < dimsubsysbar; ++j) // column major order for speed
        {
            // compute the column multi-indexes of the complement
            internal::_n2multiidx(j, nsubsysbar,
                                  Cdimssubsysbar, Cmidxcolsubsysbar);
#pragma omp parallel for
            for (idx i = 0; i < dimsubsysbar; ++i)
            {
                result(i, j) = worker(i);
            }
        }

        return result;
    }
        //************ density matrix ************//
    else if (internal::_check_square_mat(rA)) // we have a density operator
    {
        // check that dims match the dimension of A
        if (!internal::_check_dims_match_mat(dims, rA))
            throw Exception("qpp::ptrace()",
                            Exception::Type::DIMS_MISMATCH_MATRIX);

        if (subsys.size() == dims.size())
        {
            result(0, 0) = rA.trace();
            return result;
        }

        if (subsys.size() == 0)
            return rA;

        auto worker = [=, &Cmidxcolsubsysbar](idx i) noexcept
                -> typename Derived::Scalar
        {
            // use static allocation for speed!

            idx Cmidxrow[maxn];
            idx Cmidxcol[maxn];
            idx Cmidxrowsubsysbar[maxn];
            idx Cmidxsubsys[maxn];

            /* get the row/col multi-indexes of the complement */
            internal::_n2multiidx(i, nsubsysbar,
                                  Cdimssubsysbar, Cmidxrowsubsysbar);
            /* write them in the global row/col multi-indexes */
            for (idx k = 0; k < nsubsysbar; ++k)
            {
                Cmidxrow[Csubsysbar[k]] = Cmidxrowsubsysbar[k];
                Cmidxcol[Csubsysbar[k]] = Cmidxcolsubsysbar[k];
            }
            typename Derived::Scalar sm = 0;
            for (idx a = 0; a < dimsubsys; ++a)
            {
                // get the multi-index over which we do the summation
                internal::_n2multiidx(a, nsubsys, Cdimssubsys, Cmidxsubsys);
                // write it into the global row/col multi-indexes
                for (idx k = 0; k < nsubsys; ++k)
                    Cmidxrow[Csubsys[k]] = Cmidxcol[Csubsys[k]]
                            = Cmidxsubsys[k];

                // now do the sum
                sm += rA(internal::_multiidx2n(Cmidxrow, n, Cdims),
                         internal::_multiidx2n(Cmidxcol, n, Cdims));
            }

            return sm;
        };

        for (idx j = 0; j < dimsubsysbar; ++j) // column major order for speed
        {
            // compute the column multi-indexes of the complement
            internal::_n2multiidx(j, nsubsysbar,
                                  Cdimssubsysbar, Cmidxcolsubsysbar);
#pragma omp parallel for
            for (idx i = 0; i < dimsubsysbar; ++i)
            {
                result(i, j) = worker(i);
            }
        }

        return result;
    }
        //************ Exception: not ket nor density matrix ************//
    else
        throw Exception("qpp::ptrace()",
                        Exception::Type::MATRIX_NOT_SQUARE_OR_CVECTOR);
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
template<typename Derived>
dyn_mat<typename Derived::Scalar> ptrace(const Eigen::MatrixBase<Derived>& A,
                                         const std::vector<idx>& subsys,
                                         idx d = 2)
{
    const dyn_mat<typename Derived::Scalar>& rA = A;

    // check zero size
    if (!internal::_check_nonzero_size(rA))
        throw Exception("qpp::ptrace()", Exception::Type::ZERO_SIZE);

    idx n =
            static_cast<idx>(std::llround(std::log2(rA.rows()) /
                                          std::log2(d)));
    std::vector<idx> dims(n, d); // local dimensions vector

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
template<typename Derived>
dyn_mat<typename Derived::Scalar> ptranspose(
        const Eigen::MatrixBase<Derived>& A,
        const std::vector<idx>& subsys,
        const std::vector<idx>& dims)
{
    const dyn_mat<typename Derived::Scalar>& rA = A;

    // error checks

    // check zero-size
    if (!internal::_check_nonzero_size(rA))
        throw Exception("qpp::ptranspose()", Exception::Type::ZERO_SIZE);

    // check that dims is a valid dimension vector
    if (!internal::_check_dims(dims))
        throw Exception("qpp::ptranspose()", Exception::Type::DIMS_INVALID);

    // check that subsys are valid
    if (!internal::_check_subsys_match_dims(subsys, dims))
        throw Exception("qpp::ptranspose()",
                        Exception::Type::SUBSYS_MISMATCH_DIMS);

    idx D = static_cast<idx>(rA.rows());
    idx numdims = dims.size();
    idx numsubsys = subsys.size();
    idx Cdims[maxn];
    idx Cmidxcol[maxn];
    idx Csubsys[maxn];

    // copy dims in Cdims and subsys in Csubsys
    for (idx i = 0; i < numdims; ++i)
        Cdims[i] = dims[i];
    for (idx i = 0; i < numsubsys; ++i)
        Csubsys[i] = subsys[i];

    dyn_mat<typename Derived::Scalar> result(D, D);

    //************ ket ************//
    if (internal::_check_cvector(rA)) // we have a ket
    {
        // check that dims match the dimension of A
        if (!internal::_check_dims_match_cvect(dims, rA))
            throw Exception("qpp::ptranspose()",
                            Exception::Type::DIMS_MISMATCH_CVECTOR);

        if (subsys.size() == dims.size())
            return (rA * adjoint(rA)).transpose();

        if (subsys.size() == 0)
            return rA * adjoint(rA);

        auto worker = [=, &Cmidxcol](idx i) noexcept
                -> typename Derived::Scalar
        {
            // use static allocation for speed!
            idx midxcoltmp[maxn];
            idx midxrow[maxn];

            for (idx k = 0; k < numdims; ++k)
                midxcoltmp[k] = Cmidxcol[k];

            /* compute the row multi-index */
            internal::_n2multiidx(i, numdims, Cdims, midxrow);

            for (idx k = 0; k < numsubsys; ++k)
                std::swap(midxcoltmp[Csubsys[k]], midxrow[Csubsys[k]]);

            /* writes the result */
            return rA(internal::_multiidx2n(midxrow, numdims, Cdims)) *
                   std::conj(rA(internal::_multiidx2n(midxcoltmp, numdims,
                                                      Cdims)));
        };

        for (idx j = 0; j < D; ++j)
        {
            // compute the column multi-index
            internal::_n2multiidx(j, numdims, Cdims, Cmidxcol);
#pragma omp parallel for
            for (idx i = 0; i < D; ++i)
                result(i, j) = worker(i);
        }

        return result;
    }
        //************ density matrix ************//
    else if (internal::_check_square_mat(rA)) // we have a density operator
    {
        // check that dims match the dimension of A
        if (!internal::_check_dims_match_mat(dims, rA))
            throw Exception("qpp::ptranspose()",
                            Exception::Type::DIMS_MISMATCH_MATRIX);

        if (subsys.size() == dims.size())
            return rA.transpose();

        if (subsys.size() == 0)
            return rA;

        auto worker = [=, &Cmidxcol](idx i) noexcept
                -> typename Derived::Scalar
        {
            // use static allocation for speed!
            idx midxcoltmp[maxn];
            idx midxrow[maxn];

            for (idx k = 0; k < numdims; ++k)
                midxcoltmp[k] = Cmidxcol[k];

            /* compute the row multi-index */
            internal::_n2multiidx(i, numdims, Cdims, midxrow);

            for (idx k = 0; k < numsubsys; ++k)
                std::swap(midxcoltmp[Csubsys[k]], midxrow[Csubsys[k]]);

            /* writes the result */
            return rA(internal::_multiidx2n(midxrow, numdims, Cdims),
                      internal::_multiidx2n(midxcoltmp, numdims, Cdims));
        };

        for (idx j = 0; j < D; ++j)
        {
            // compute the column multi-index
            internal::_n2multiidx(j, numdims, Cdims, Cmidxcol);
#pragma omp parallel for
            for (idx i = 0; i < D; ++i)
                result(i, j) = worker(i);
        }

        return result;
    }
        //************ Exception: not ket nor density matrix ************//
    else
        throw Exception("qpp::ptranspose()",
                        Exception::Type::MATRIX_NOT_SQUARE_OR_CVECTOR);
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
template<typename Derived>
dyn_mat<typename Derived::Scalar> ptranspose(
        const Eigen::MatrixBase<Derived>& A,
        const std::vector<idx>& subsys,
        idx d = 2)
{
    const dyn_mat<typename Derived::Scalar>& rA = A;

    // check zero size
    if (!internal::_check_nonzero_size(rA))
        throw Exception("qpp::ptranspose()", Exception::Type::ZERO_SIZE);

    idx n =
            static_cast<idx>(std::llround(std::log2(rA.rows()) /
                                          std::log2(d)));
    std::vector<idx> dims(n, d); // local dimensions vector

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
template<typename Derived>
dyn_mat<typename Derived::Scalar> syspermute(
        const Eigen::MatrixBase<Derived>& A,
        const std::vector<idx>& perm,
        const std::vector<idx>& dims)
{
    const dyn_mat<typename Derived::Scalar>& rA = A;

    // Error checks

    // check zero-size
    if (!internal::_check_nonzero_size(rA))
        throw Exception("qpp::syspermute()", Exception::Type::ZERO_SIZE);

    // check that dims is a valid dimension vector
    if (!internal::_check_dims(dims))
        throw Exception("qpp::syspermute()", Exception::Type::DIMS_INVALID);

    // check that we have a valid permutation
    if (!internal::_check_perm(perm))
        throw Exception("qpp::syspermute()", Exception::Type::PERM_INVALID);

    // check that permutation match dimensions
    if (perm.size() != dims.size())
        throw Exception("qpp::syspermute()",
                        Exception::Type::PERM_MISMATCH_DIMS);

    idx D = static_cast<idx>(rA.rows());
    idx numdims = dims.size();

    dyn_mat<typename Derived::Scalar> result;

    //************ ket ************//
    if (internal::_check_cvector(rA)) // we have a column vector
    {
        idx Cdims[maxn];
        idx Cperm[maxn];

        // check that dims match the dimension of rA
        if (!internal::_check_dims_match_cvect(dims, rA))
            throw Exception("qpp::syspermute()",
                            Exception::Type::DIMS_MISMATCH_CVECTOR);

        // copy dims in Cdims and perm in Cperm
        for (idx i = 0; i < numdims; ++i)
        {
            Cdims[i] = dims[i];
            Cperm[i] = perm[i];
        }
        result.resize(D, 1);

        auto worker = [&Cdims, &Cperm, numdims](idx i) noexcept
                -> idx
        {
            // use static allocation for speed,
            // double the size for matrices reshaped as vectors
            idx midx[maxn];
            idx midxtmp[maxn];
            idx permdims[maxn];

            /* compute the multi-index */
            internal::_n2multiidx(i, numdims, Cdims, midx);

            for (idx k = 0; k < numdims; ++k)
            {
                permdims[k] = Cdims[Cperm[k]]; // permuted dimensions
                midxtmp[k] = midx[Cperm[k]];// permuted multi-indexes
            }

            return internal::_multiidx2n(midxtmp, numdims, permdims);
        };

#pragma omp parallel for
        for (idx i = 0; i < D; ++i)
            result(worker(i)) = rA(i);

        return result;
    }
        //************ density matrix ************//
    else if (internal::_check_square_mat(rA)) // we have a density operator
    {
        idx Cdims[2 * maxn];
        idx Cperm[2 * maxn];

        // check that dims match the dimension of rA
        if (!internal::_check_dims_match_mat(dims, rA))
            throw Exception("qpp::syspermute()",
                            Exception::Type::DIMS_MISMATCH_MATRIX);

        // copy dims in Cdims and perm in Cperm
        for (idx i = 0; i < numdims; ++i)
        {
            Cdims[i] = dims[i];
            Cdims[i + numdims] = dims[i];
            Cperm[i] = perm[i];
            Cperm[i + numdims] = perm[i] + numdims;
        }
        result.resize(D * D, 1);
        // map A to a column vector
        dyn_mat<typename Derived::Scalar> vectA =
                Eigen::Map<dyn_mat<typename Derived::Scalar>>(
                        const_cast<typename Derived::Scalar*>(rA.data()), D * D,
                        1);

        auto worker = [&Cdims, &Cperm, numdims](idx i) noexcept
                -> idx
        {
            // use static allocation for speed,
            // double the size for matrices reshaped as vectors
            idx midx[2 * maxn];
            idx midxtmp[2 * maxn];
            idx permdims[2 * maxn];

            /* compute the multi-index */
            internal::_n2multiidx(i, 2 * numdims, Cdims, midx);

            for (idx k = 0; k < 2 * numdims; ++k)
            {
                permdims[k] = Cdims[Cperm[k]]; // permuted dimensions
                midxtmp[k] = midx[Cperm[k]];// permuted multi-indexes
            }

            return internal::_multiidx2n(midxtmp, 2 * numdims, permdims);
        };

#pragma omp parallel for
        for (idx i = 0; i < D * D; ++i)
            result(worker(i)) = rA(i);

        return reshape(result, D, D);
    }
        //************ Exception: not ket nor density matrix ************//
    else
        throw Exception("qpp::syspermute()",
                        Exception::Type::MATRIX_NOT_SQUARE_OR_CVECTOR);
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
template<typename Derived>
dyn_mat<typename Derived::Scalar> syspermute(
        const Eigen::MatrixBase<Derived>& A,
        const std::vector<idx>& perm,
        idx d = 2)
{
    const dyn_mat<typename Derived::Scalar>& rA = A;

    // check zero size
    if (!internal::_check_nonzero_size(rA))
        throw Exception("qpp::syspermute()", Exception::Type::ZERO_SIZE);

    idx n =
            static_cast<idx>(std::llround(std::log2(rA.rows()) /
                                          std::log2(d)));
    std::vector<idx> dims(n, d); // local dimensions vector

    return syspermute(rA, perm, dims);
}

} /* namespace qpp */

#endif /* OPERATIONS_H_ */
