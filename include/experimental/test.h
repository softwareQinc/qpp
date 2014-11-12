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

#ifndef INCLUDE_EXPERIMENTAL_TEST_H_
#define INCLUDE_EXPERIMENTAL_TEST_H_

// testing functions, do not use/modify
/**
* \namespace qpp::experimental Experimental/test functions,
* do not use/modify these functions/classes
*/
namespace qpp
{
    namespace experimental
    {

        /**
        * \brief Applies the gate \a A to the part \a subsys
        * of a multi-partite state vector or density matrix
        *
        * \note The dimension of the gate \a A must match
        * the dimension of \a subsys
        *
        * \param state Eigen expression
        * \param A Eigen expression
        * \param subsys Subsystem indexes where the gate \a A is applied
        * \param dims Local dimensions of all local Hilbert spaces (can be different)
        * \return Gate \a A applied to the part \a subsys of \a state
        */
        template<typename Derived1, typename Derived2>
        DynMat<typename Derived1::Scalar> apply_old(
                const Eigen::MatrixBase<Derived1> &state,
                const Eigen::MatrixBase<Derived2> &A,
                const std::vector<std::size_t> &subsys,
                const std::vector<std::size_t> &dims)
        {
            const DynMat<typename Derived1::Scalar> &rstate = state;
            const DynMat<typename Derived2::Scalar> &rA = A;

            // EXCEPTION CHECKS

            // check types
            if (!std::is_same<typename Derived1::Scalar, typename Derived2::Scalar>::value)
                throw Exception("qpp::experimental::apply()", Exception::Type::TYPE_MISMATCH);

            // check zero sizes
            if (!internal::_check_nonzero_size(rA))
                throw Exception("qpp::experimental::apply()", Exception::Type::ZERO_SIZE);

            // check zero sizes
            if (!internal::_check_nonzero_size(rstate))
                throw Exception("qpp::experimental::apply()", Exception::Type::ZERO_SIZE);

            // check square matrix for the gate
            if (!internal::_check_square_mat(rA))
                throw Exception("qpp::experimental::apply()", Exception::Type::MATRIX_NOT_SQUARE);

            // check that dims is a valid dimension vector
            if (!internal::_check_dims(dims))
                throw Exception("qpp::experimental::apply()", Exception::Type::DIMS_INVALID);

            // check subsys is valid w.r.t. dims
            if (!internal::_check_subsys_match_dims(subsys, dims))
                throw Exception("qpp::experimental::apply()", Exception::Type::SUBSYS_MISMATCH_DIMS);

            // Use static allocation for speed!
            std::size_t Cdims[maxn];
            std::size_t midx_row[maxn];
            std::size_t midx_rho_row[maxn];

            std::size_t CdimsA[maxn];
            std::size_t CsubsysA[maxn];
            std::size_t midxA_row[maxn];
            std::size_t midxA_rho_row[maxn];

            std::size_t CdimsA_bar[maxn];
            std::size_t CsubsysA_bar[maxn];
            std::size_t midxA_bar_row[maxn];

            std::size_t n = dims.size();
            std::size_t nA = subsys.size();
            std::size_t nA_bar = n - nA;

            std::size_t D = 1;
            std::size_t DA_bar = 1;

            for (std::size_t k = 0, cnt = 0; k < n; ++k)
            {
                midx_row[k] = midx_rho_row[k] = 0;
                Cdims[k] = dims[k];
                D *= dims[k];

                // compute the complement of subsys w.r.t. dims
                if (std::find(std::begin(subsys), std::end(subsys), k)
                        == std::end(subsys))
                {
                    CsubsysA_bar[cnt] = k;
                    CdimsA_bar[cnt] = dims[k];
                    midxA_bar_row[cnt] = 0;
                    DA_bar *= dims[k];
                    cnt++;
                }
            }

            std::size_t DA = 1;
            for (std::size_t k = 0; k < nA; ++k)
            {
                midxA_row[k] = midxA_rho_row[k] = 0;
                CdimsA[k] = dims[subsys[k]];
                CsubsysA[k] = subsys[k];
                DA *= dims[subsys[k]];
            }

            // check that gate mathches the dimensions of the subsys
            if (static_cast<std::size_t>(rA.rows()) != DA)
                throw Exception("qpp::experimental::apply()", Exception::Type::DIMS_MISMATCH_MATRIX);

            if (internal::_check_col_vector(rstate)) // we have a ket
            {
                // check that dims match state vector
                if (!internal::_check_dims_match_cvect(dims, rstate))
                    throw Exception("qpp::experimental::apply()", Exception::Type::DIMS_MISMATCH_CVECTOR);

                DynMat<typename Derived1::Scalar> result(D, 1);

                // run over the subsys's row multi-index
                for (std::size_t a = 0; a < DA; ++a)
                {
                    // get the subsys's row multi-index
                    internal::_n2multiidx(a, nA, CdimsA, midxA_row);
                    // compute subsys part of the result's row multi-index
                    for (std::size_t k = 0; k < nA; ++k)
                        midx_row[CsubsysA[k]] = midxA_row[k];

                    // run over the complement's row multi-index
                    for (std::size_t i = 0; i < DA_bar; ++i)
                    {
                        // get the complement's row multi-index
                        internal::_n2multiidx(i, nA_bar, CdimsA_bar, midxA_bar_row);
                        // now compute the complement part of the
                        // result's row multi-index
                        // and the complement part
                        // of the state's total row multi-index
                        for (std::size_t k = 0; k < nA_bar; ++k)
                            midx_row[CsubsysA_bar[k]] = midx_rho_row[CsubsysA_bar[k]] =
                                    midxA_bar_row[k];
                        // compute the results's row index
                        std::size_t result_row_idx = internal::_multiidx2n(midx_row, n,
                                Cdims);

                        // compute the coefficient
                        typename Derived1::Scalar coeff = 0;
                        for (std::size_t c = 0; c < DA; ++c)
                        {
                            // compute the subsys part state's row multi-index
                            internal::_n2multiidx(c, nA, CdimsA, midxA_rho_row);
                            // now we have the total state's row multi-index
                            for (std::size_t k = 0; k < nA; ++k)
                                midx_rho_row[CsubsysA[k]] = midxA_rho_row[k];

                            coeff += rA(a, c)
                                    * rstate(
                                    internal::_multiidx2n(midx_rho_row, n,
                                            Cdims));
                        }
                        // write down the result
                        result(result_row_idx) = coeff;
                    }
                }
                return result;
            }
            else if (internal::_check_square_mat(rstate)) // we have a density matrix
            {

                // check that dims match state matrix
                if (!internal::_check_dims_match_mat(dims, rstate))
                    throw Exception("qpp::experimental::apply()", Exception::Type::DIMS_MISMATCH_MATRIX);

                DynMat<typename Derived1::Scalar> result(D, D);

                // run over the subsys's row multi-index
                for (std::size_t a = 0; a < DA; ++a)
                {
                    // get the subsys's row multi-index
                    internal::_n2multiidx(a, nA, CdimsA, midxA_row);
                    // compute subsys part of the result's row multi-index
                    for (std::size_t k = 0; k < nA; ++k)
                        midx_row[CsubsysA[k]] = midxA_row[k];

                    // run over the complement's row multi-index
                    for (std::size_t i = 0; i < DA_bar; ++i)
                    {
                        // get the complement's row multi-index
                        internal::_n2multiidx(i, nA_bar, CdimsA_bar, midxA_bar_row);
                        // now compute the complement part
                        // of the result's row multi-index
                        // and the complement part of the
                        // state's total row multi-index
                        for (std::size_t k = 0; k < nA_bar; ++k)
                            midx_row[CsubsysA_bar[k]] = midx_rho_row[CsubsysA_bar[k]] =
                                    midxA_bar_row[k];
                        // compute the results's row index
                        std::size_t result_row_idx = internal::_multiidx2n(midx_row, n,
                                Cdims);

                        // run over the col index
                        for (std::size_t j = 0; j < D; ++j)
                        {
                            // compute the coefficient
                            typename Derived1::Scalar coeff = 0;
                            for (std::size_t c = 0; c < DA; ++c)
                            {
                                // compute the subsys part state's row multi-index
                                internal::_n2multiidx(c, nA, CdimsA, midxA_rho_row);
                                // now we have the total state's row multi-index
                                for (std::size_t k = 0; k < nA; ++k)
                                    midx_rho_row[CsubsysA[k]] = midxA_rho_row[k];

                                coeff += rA(a, c)
                                        * rstate(
                                        internal::_multiidx2n(midx_rho_row, n,
                                                Cdims), j);

                            }
                            // write down the result
                            result(result_row_idx, j) = coeff;
                        }
                    }
                }
                return result;
            }
            else
                throw Exception("qpp::experimental::apply()", Exception::Type::MATRIX_NOT_SQUARE_OR_CVECTOR);
        }

        /**
        * \brief Applies the channel specified by the set of Kraus operators \a Ks to
        * the part of the density matrix \a rho specified by \a subsys
        *
        * \param rho Eigen expression
        * \param Ks Set of Kraus operators
        * \param subsys Subsystems' indexes
        * \param dims Local dimensions of all local Hilbert spaces (can be different)
        * \return Output density matrix after the action of the channel
        */
        template<typename Derived>
        cmat channel_old(const Eigen::MatrixBase<Derived> &rho, const std::vector<cmat> &Ks,
                const std::vector<std::size_t> &subsys,
                const std::vector<std::size_t> &dims)
        {
            const cmat &rrho = rho;

            // EXCEPTION CHECKS
            // check zero sizes
            if (!internal::_check_nonzero_size(rrho))
                throw Exception("qpp::experimental::channel()", Exception::Type::ZERO_SIZE);

            // check square matrix for the rho
            if (!internal::_check_square_mat(rrho))
                throw Exception("qpp::experimental::channel()", Exception::Type::MATRIX_NOT_SQUARE);

            // check that dims is a valid dimension vector
            if (!internal::_check_dims(dims))
                throw Exception("qpp::experimental::channel()", Exception::Type::DIMS_INVALID);

            // check that dims match rho matrix
            if (!internal::_check_dims_match_mat(dims, rrho))
                throw Exception("qpp::experimental::channel()", Exception::Type::DIMS_MISMATCH_MATRIX);

            // check subsys is valid w.r.t. dims
            if (!internal::_check_subsys_match_dims(subsys, dims))
                throw Exception("qpp::experimental::channel()", Exception::Type::SUBSYS_MISMATCH_DIMS);

            // check the Kraus operators
            if (!internal::_check_nonzero_size(Ks))
                throw Exception("qpp::experimental::channel()", Exception::Type::ZERO_SIZE);
            if (!internal::_check_square_mat(Ks[0]))
                throw Exception("qpp::experimental::channel()", Exception::Type::MATRIX_NOT_SQUARE);
            for (auto &&it : Ks)
                if (it.rows() != Ks[0].rows() || it.cols() != Ks[0].rows())
                    throw Exception("qpp::experimental::channel()", Exception::Type::DIMS_NOT_EQUAL);

            // Use static allocation for speed!
            std::size_t Cdims[maxn];
            std::size_t midx_row[maxn];
            std::size_t midx_col[maxn];
            std::size_t midx_rho_row[maxn];
            std::size_t midx_rho_col[maxn];

            std::size_t CdimsA[maxn];
            std::size_t CsubsysA[maxn];
            std::size_t midxA_row[maxn];
            std::size_t midxA_col[maxn];
            std::size_t midxA_rho_row[maxn];
            std::size_t midxA_rho_col[maxn];

            std::size_t CsubsysA_bar[maxn];
            std::size_t midxA_bar_row[maxn];
            std::size_t midxA_bar_col[maxn];

            std::size_t n = dims.size();
            std::size_t nA = subsys.size();
            std::size_t nA_bar = n - nA;

            std::size_t D = 1;
            std::size_t DA_bar = 1;

            for (std::size_t k = 0, cnt = 0; k < n; ++k)
            {
                midx_row[k] = midx_col[k] = midx_rho_row[k] = midx_rho_col[k] = 0;
                Cdims[k] = dims[k];
                D *= dims[k];

                // compute the complementary subsystem of subsys w.r.t. dims
                if (std::find(std::begin(subsys), std::end(subsys), k)
                        == std::end(subsys))
                {
                    CsubsysA_bar[cnt] = k;
                    midxA_bar_row[cnt] = midxA_bar_col[cnt] = 0;
                    DA_bar *= dims[k];
                    cnt++;
                }
            }

            std::size_t DA = 1;
            for (std::size_t k = 0; k < nA; ++k)
            {
                midxA_row[k] = midxA_col[k] = midxA_rho_row[k] = midxA_rho_col[k] = 0;
                CdimsA[k] = dims[subsys[k]];
                CsubsysA[k] = subsys[k];
                DA *= dims[subsys[k]];
            }

            // check that dimension of Kraus matches the dimension of the subsys
            if (static_cast<std::size_t>(Ks[0].rows()) != DA)
                throw Exception("qpp::experimental::channel()", Exception::Type::DIMS_MISMATCH_MATRIX);

            // get the superoperator matrix of the channel
            cmat sop = super(Ks);

            cmat result(D, D);

            // run over rows
            for (std::size_t i = 0; i < D; ++i)
            {
                // get the result's row multi-index
                internal::_n2multiidx(i, n, Cdims, midx_row);
                // get the subsys' complement row multi-index
                for (std::size_t k = 0; k < nA_bar; ++k)
                    midxA_bar_row[k] = midx_row[CsubsysA_bar[k]];
                // get the subsys' row multi-index
                for (std::size_t k = 0; k < nA; ++k)
                    midxA_row[k] = midx_row[CsubsysA[k]];

                // run over cols
                for (std::size_t j = 0; j < D; ++j)
                {
                    // get the result's col multi-index
                    internal::_n2multiidx(j, n, Cdims, midx_col);
                    // get the subsys' complement col multi-index
                    for (std::size_t k = 0; k < nA_bar; ++k)
                        midxA_bar_col[k] = midx_col[CsubsysA_bar[k]];
                    // get the subsys' col multi-index
                    for (std::size_t k = 0; k < nA; ++k)
                        midxA_col[k] = midx_col[CsubsysA[k]];

                    // now compute the coefficient
                    cplx coeff = 0;
                    for (std::size_t a = 0; a < DA; ++a)
                    {
                        // get the subsys part of row multi-index for rho
                        internal::_n2multiidx(a, nA, CdimsA, midxA_rho_row);
                        for (std::size_t b = 0; b < DA; ++b)
                        {
                            // get the subsys part of col multi-index for rho
                            internal::_n2multiidx(b, nA, CdimsA, midxA_rho_col);

                            // get the total row/col multi-index for rho
                            for (std::size_t k = 0; k < nA; ++k)
                            {
                                midx_rho_row[CsubsysA[k]] = midxA_rho_row[k];
                                midx_rho_col[CsubsysA[k]] = midxA_rho_col[k];
                            }
                            for (std::size_t k = 0; k < nA_bar; ++k)
                            {
                                midx_rho_row[CsubsysA_bar[k]] = midxA_bar_row[k];
                                midx_rho_col[CsubsysA_bar[k]] = midxA_bar_col[k];
                            }

                            std::size_t midx_sop_col[2]; // index the superop using 2 indices
                            std::size_t midx_sop_row[2];
                            std::size_t sop_dims[2];
                            sop_dims[0] = sop_dims[1] = DA;
                            midx_sop_row[0] = internal::_multiidx2n(midxA_row, nA,
                                    CdimsA);
                            midx_sop_row[1] = internal::_multiidx2n(midxA_col, nA,
                                    CdimsA);
                            midx_sop_col[0] = a;
                            midx_sop_col[1] = b;

                            coeff += sop(
                                    internal::_multiidx2n(midx_sop_row, 2, sop_dims),
                                    internal::_multiidx2n(midx_sop_col, 2, sop_dims))
                                    * rrho(
                                    internal::_multiidx2n(midx_rho_row, n,
                                            Cdims),
                                    internal::_multiidx2n(midx_rho_col, n,
                                            Cdims));
                        }
                    }
                    result(i, j) = coeff;
                }
            }
            return result;
        }

        /**
        * \brief Superoperator matrix representation
        *
        * Constructs the superoperator matrix of the channel specified by the set of
        * Kraus operators \a Ks in the standard operator basis
        * \f$\{|i\rangle\langle j|\}\f$ ordered in lexicographical order, i.e.
        * \f$|0\rangle\langle 0|\f$, \f$|0\rangle\langle 1|\f$ etc.
        *
        * \param Ks Set of Kraus operators
        * \return Superoperator matrix representation
        */
        cmat super(const std::vector<cmat> &Ks)
        {
            // EXCEPTION CHECKS
            if (!internal::_check_nonzero_size(Ks))
                throw Exception("qpp::experimental::super()", Exception::Type::ZERO_SIZE);
            if (!internal::_check_nonzero_size(Ks[0]))
                throw Exception("qpp::experimental::super()", Exception::Type::ZERO_SIZE);
            if (!internal::_check_square_mat(Ks[0]))
                throw Exception("qpp::experimental::super()", Exception::Type::MATRIX_NOT_SQUARE);
            for (auto &&it : Ks)
                if (it.rows() != Ks[0].rows() || it.cols() != Ks[0].rows())
                    throw Exception("qpp::experimental::super()", Exception::Type::DIMS_NOT_EQUAL);
            std::size_t D = static_cast<std::size_t>(Ks[0].rows());

            std::size_t midx_row[2] = {0, 0};
            std::size_t midx_col[2] = {0, 0};
            std::size_t dims[2];
            dims[0] = dims[1] = D;

            cmat result(D * D, D * D);
            cmat MN = cmat::Zero(D, D);
            cmat BA = cmat::Zero(D, D);
            cmat EMN = cmat::Zero(D, D);

            for (std::size_t m = 0; m < D; ++m)
            {
                midx_col[0] = m;
                for (std::size_t n = 0; n < D; ++n)
                {
                    midx_col[1] = n;
                    MN(m, n) = 1;
                    // compute E(|m><n|)
                    for (std::size_t i = 0; i < Ks.size(); ++i)
                        EMN += Ks[i] * MN * adjoint(Ks[i]);
                    MN(m, n) = 0;
                    for (std::size_t a = 0; a < D; ++a)
                    {
                        midx_row[0] = a;
                        for (std::size_t b = 0; b < D; ++b)
                        {
                            midx_row[1] = b;
                            BA(b, a) = 1;

                            // compute result(ab,mn)=<a|E(|m><n)|b>

                            result(internal::_multiidx2n(midx_row, 2, dims),
                                    internal::_multiidx2n(midx_col, 2, dims)) = (EMN
                                    * BA).trace();
                            BA(b, a) = 0;
                        }
                    }
                    EMN = cmat::Zero(D, D);
                }
            }
            return result;
        }

        /**
        * \brief Generates the multi-partite multiple-controlled-\a A gate
        * in matrix form
        *
        * \note The dimension of the gate \a A must match
        * the dimension of \a subsys
        *
        * \param A Eigen expression
        * \param ctrl Control subsystem indexes
        * \param subsys Subsystem indexes where the gate \a A is applied
        * \param n Total number of subsystes
        * \param d Local dimensions of all local Hilbert spaces (must all be equal)
        * \return CTRL-A gate, as a matrix over the same scalar field as \a A
        */
        // Parallel version, seems slower than the sequential qpp::Gates::CTRL
        template<typename Derived>
        DynMat<typename Derived::Scalar> CTRL(const Eigen::MatrixBase<Derived> &A,
                const std::vector<std::size_t> &ctrl,
                const std::vector<std::size_t> &subsys, std::size_t n,
                std::size_t d = 2)
        {
            const DynMat<typename Derived::Scalar> &rA = A;

            // EXCEPTION CHECKS
            // check matrix zero size
            if (!internal::_check_nonzero_size(rA))
                throw Exception("qpp::experimental::CTRL()", Exception::Type::ZERO_SIZE);

            // check square matrix
            if (!internal::_check_square_mat(rA))
                throw Exception("qpp::experimental::CTRL()", Exception::Type::MATRIX_NOT_SQUARE);

            // check lists zero size
            if (subsys.size() == 0)
                throw Exception("qpp::experimental::CTRL()", Exception::Type::ZERO_SIZE);

            // check out of range
            if (n == 0)
                throw Exception("qpp::experimental::CTRL()", Exception::Type::OUT_OF_RANGE);

            // check valid local dimension
            if (d == 0)
                throw Exception("qpp::experimental::CTRL()", Exception::Type::DIMS_INVALID);

            std::vector<std::size_t> ctrlgate = ctrl;    // ctrl + gate subsystem vector
            ctrlgate.insert(std::end(ctrlgate), std::begin(subsys), std::end(subsys));
            std::sort(std::begin(ctrlgate), std::end(ctrlgate));

            std::vector<std::size_t> dims(n, d); // local dimensions vector

            // check that ctrl + gate subsystem is valid
            // with respect to local dimensions
            if (!internal::_check_subsys_match_dims(ctrlgate, dims))
                throw Exception("qpp::experimental::CTRL()", Exception::Type::SUBSYS_MISMATCH_DIMS);

            // check that subsys list match the dimension of the matrix
            if (rA.cols() != std::pow(d, subsys.size()))
                throw Exception("qpp::experimental::CTRL()", Exception::Type::DIMS_MISMATCH_MATRIX);
            // END EXCEPTION CHECKS

            if (d == 1)
                return rA;

            std::size_t ctrlsize = ctrl.size();

            // construct the table of A^i
            std::vector<DynMat<typename Derived::Scalar>> Ai;
            std::vector<DynMat<typename Derived::Scalar>> Aidagger;
            for (std::size_t i = 0; i < d; ++i)
            {
                Ai.push_back(powm(rA, i));
            }

            std::size_t D = std::pow(d, n);
            std::size_t DA = rA.rows();
            std::size_t DCTRLAbar = static_cast<std::size_t>(std::pow(d,
                    n - ctrlgate.size()));

            std::size_t Cdims[maxn]; // local dimensions
            std::size_t CdimsA[maxn]; // local dimensions
            std::size_t CdimsCTRLAbar[maxn]; // local dimensions

            std::vector<std::size_t> ctrlgatebar(n - ctrlgate.size()); // rest
            std::vector<std::size_t> allsubsys(n);
            std::iota(std::begin(allsubsys), std::end(allsubsys), 0);
            std::set_difference(std::begin(allsubsys), std::end(allsubsys),
                    std::begin(ctrlgate), std::end(ctrlgate), std::begin(ctrlgatebar));

            for (std::size_t k = 0; k < n; ++k)
                Cdims[k] = d;
            for (std::size_t k = 0; k < subsys.size(); ++k)
                CdimsA[k] = d;
            for (std::size_t k = 0; k < n - ctrlgate.size(); ++k)
                CdimsCTRLAbar[k] = d;

            auto coeff =
                    [=](std::size_t _i, std::size_t _m, std::size_t _n, std::size_t _r)
                            -> std::pair<std::size_t, std::size_t>
                    {
                        std::size_t idxrow = 0;
                        std::size_t idxcol = 0;
                        std::size_t Cmidxrow[maxn]; // the total row multi-index
                        std::size_t Cmidxcol[maxn];// the total col multi-index
                        std::size_t CmidxArow[maxn];// the gate part row multi-index
                        std::size_t CmidxAcol[maxn];// the gate part col multi-index
                        std::size_t CmidxCTRLAbar[maxn];// the rest multi-index

                        // compute the index

                        // set the CTRL part
                        for (std::size_t k = 0; k < ctrl.size(); ++k)
                        {
                            Cmidxrow[ctrl[k]] = Cmidxcol[ctrl[k]] = _i;
                        }

                        // set the rest
                        internal::_n2multiidx(_r, n - ctrlgate.size(),
                                CdimsCTRLAbar, CmidxCTRLAbar);
                        for (std::size_t k = 0; k < n - ctrlgate.size(); ++k)
                        {
                            Cmidxrow[ctrlgatebar[k]] = Cmidxcol[ctrlgatebar[k]] = CmidxCTRLAbar[k];
                        }

                        // set the A part
                        internal::_n2multiidx(_m, subsys.size(), CdimsA, CmidxArow);
                        internal::_n2multiidx(_n, subsys.size(), CdimsA, CmidxAcol);
                        for (std::size_t k = 0; k < subsys.size(); ++k)
                        {
                            Cmidxrow[subsys[k]] = CmidxArow[k];
                            Cmidxcol[subsys[k]] = CmidxAcol[k];
                        }

                        // we now got the total row/col indexes
                        idxrow = internal::_multiidx2n(Cmidxrow, n, Cdims);
                        idxcol = internal::_multiidx2n(Cmidxcol, n, Cdims);

                        return std::make_pair(idxrow, idxcol);
                    };

            DynMat<typename Derived::Scalar> result =
                    DynMat<typename Derived::Scalar>::Identity(D, D);

#pragma omp parallel for collapse(4)
            for (std::size_t m = 0; m < DA; ++m)
                for (std::size_t n = 0; n < DA; ++n)
                    for (std::size_t r = 0; r < DCTRLAbar; ++r)
                        for (std::size_t i = 0; i < d; ++i)
                            if (ctrlsize == 0) // no control
                            {
                                result(coeff(i, m, n, r).first,
                                        coeff(i, m, n, r).second) = Ai[1](m, n);
                            }
                            else
                            {
                                result(coeff(i, m, n, r).first,
                                        coeff(i, m, n, r).second) = Ai[i](m, n);
                            }

            return result;
        }

        /**
        * \brief Choi matrix representation
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
        * \return Choi matrix representation
        */
        cmat choi(const std::vector<cmat> &Ks)
        {
            // EXCEPTION CHECKS
            if (!internal::_check_nonzero_size(Ks))
                throw Exception("qpp::experimental::choi()", Exception::Type::ZERO_SIZE);
            if (!internal::_check_nonzero_size(Ks[0]))
                throw Exception("qpp::experimental::choi()", Exception::Type::ZERO_SIZE);
            if (!internal::_check_square_mat(Ks[0]))
                throw Exception("qpp::experimental::choi()", Exception::Type::MATRIX_NOT_SQUARE);
            for (auto &&it : Ks)
                if (it.rows() != Ks[0].rows() || it.cols() != Ks[0].rows())
                    throw Exception("qpp::experimental::choi()", Exception::Type::DIMS_NOT_EQUAL);
            std::size_t D = static_cast<std::size_t>(Ks[0].rows());

            // construct the D x D \sum |jj> vector
            // (un-normalized maximally entangled state)
            cmat MES = cmat::Zero(D * D, 1);
            for (std::size_t a = 0; a < D; ++a)
                MES(a * D + a) = 1;

            cmat Omega = static_cast<cmat>(MES * adjoint(MES));

            cmat result = cmat::Zero(D * D, D * D);

            for (std::size_t i = 0; i < Ks.size(); ++i)
            {
                result += kron(cmat::Identity(D, D), Ks[i]) * Omega
                        * adjoint(kron(cmat::Identity(D, D), Ks[i]));
            }

            return result;
        }

        /**
        * \brief Generates a set of random Kraus operators
        *
        * \note The set of Kraus operators satisfy the closure condition
        * \f$ \sum_i K_i^\dagger K_i = I\f$
        *
        * \param n Number of Kraus operators
        * \param D Dimension of the Hilbert space
        * \return Set of \a n Kraus operators satisfying the closure condition
        */
        std::vector<cmat> randkraus(std::size_t n, std::size_t D)
        {
            if (n == 0)
                throw Exception("qpp::experimental::randkraus()", Exception::Type::OUT_OF_RANGE);
            if (D == 0)
                throw Exception("qpp::experimental::randkraus()", Exception::Type::DIMS_INVALID);

            std::vector<cmat> result;
            cmat Fk(D, D);
            cmat U = randU(n * D);
            std::size_t dims[2];
            dims[0] = D;
            dims[1] = n;
            std::size_t midx_row[2] = {0, 0};
            std::size_t midx_col[2] = {0, 0};

            for (std::size_t k = 0; k < n; ++k)
            {
                midx_row[1] = k;
                for (std::size_t a = 0; a < D; ++a)
                {
                    midx_row[0] = a;
                    for (std::size_t b = 0; b < D; ++b)
                    {
                        midx_col[0] = b;
                        Fk(a, b) = U(internal::_multiidx2n(midx_row, 2, dims),
                                internal::_multiidx2n(midx_col, 2, dims));
                    }
                }
                result.push_back(Fk);
            }

            return result;
        }

        /**
        * \brief Renyi-\f$\infty\f$ entropy (min entropy) of the
        * probability distribution/density matrix \a A
        *
        * \param A Eigen expression, representing a probability distribution
        * (real dynamic column vector) or a density matrix (complex dynamic matrix)
        * \return Renyi-\f$\infty\f$ entropy (min entropy),
        * with the logarithm in base 2
        */
        template<typename Derived>
        double renyi_inf(const Eigen::MatrixBase<Derived> &A)
        {
            const DynMat<typename Derived::Scalar> &rA = A;

            // check zero-size
            if (!internal::_check_nonzero_size(rA))
                throw Exception("qpp::experimental::renyi_inf()", Exception::Type::ZERO_SIZE);

            // input is a vector
            if (internal::_check_vector(rA))
            {
                double max = 0;
                for (std::size_t i = 0; i < static_cast<std::size_t>(rA.size()); ++i)
                    if (std::abs(rA(i)) > max)
                        max = std::abs(rA(i));

                return -std::log2(max);
            }

            // input is a matrix

            // check square matrix
            if (!internal::_check_square_mat(rA))
                throw Exception("qpp::experimental::renyi_inf()", Exception::Type::MATRIX_NOT_SQUARE);

            // get the eigenvalues
            dmat ev = hevals(rA);
            double max = 0;
            // take the absolut value to get rid of tiny negatives
            for (std::size_t i = 0; i < static_cast<std::size_t>(ev.size()); ++i)
                if (std::abs((cplx) ev(i)) > max)
                    max = std::abs((cplx) ev(i));

            return -std::log2(max);
        }

        /**
        * \brief Displays a range. Does not add a newline.
        *
        * \see qpp::experimental::displn()
        *
        * \param first Iterator to the first element of the range
        * \param last  Iterator to the last element of the range
        * \param separator Separator
        * \param start Left marking
        * \param end Right marking
        * \param os Output stream
        * \return Output stream
        */
        template<typename InputIterator>
        void disp(const InputIterator &first, const InputIterator &last,
                const std::string &separator, const std::string &start = "[",
                const std::string &end = "]", std::ostream &os = std::cout)
        {
            os << start;

            auto it = first;
            auto it_end = last;

            if (it != it_end)
            {
                // the iterator just before the end, need this for containers
                // that do not have backwards iterators
                decltype(it_end) it_before_end = it;
                while (it_before_end = it, ++it != it_end);

                it = first;
                for (; it != it_before_end; ++it)
                    os << *it << separator;
                os << *it;
            }

            os << end;
        }

        /**
        * \brief Displays a range. Adds a newline.
        *
        * \see qpp::experimental::disp()
        *
        * \param first Iterator to the first element of the range
        * \param last  Iterator to the last element of the range
        * \param separator Separator
        * \param start Left marking
        * \param end Right marking
        * \param os Output stream
        * \return Output stream
        */
        template<typename InputIterator>
        std::ostream &displn(const InputIterator &first, const InputIterator &last,
                const std::string &separator, const std::string &start = "[",
                const std::string &end = "]", std::ostream &os = std::cout)
        {
            disp(first, last, separator, start, end, os);
            os << std::endl;
            return os;
        }

        /**
        * \brief Displays a standard container that supports std::begin, std::end
        * and forward iteration. Does not add a newline.
        *
        * \see qpp::experimental::displn()
        *
        * \param x Container
        * \param separator Separator
        * \param start Left marking
        * \param end Right marking
        * \param os Output stream
        * \return Output stream
        */
        template<typename T>
        std::ostream &disp(const T &x, const std::string &separator,
                const std::string &start = "[", const std::string &end = "]",
                std::ostream &os = std::cout)
        {
            disp(std::begin(x), std::end(x), separator, start, end, os);
            return os;
        }

        /**
        * \brief Displays a standard container that supports std::begin, std::end
        * and forward iteration. Adds a newline.
        *
        * \see qpp::experimental::disp()
        *
        * \param x Container
        * \param separator Separator
        * \param start Left marking
        * \param end Right marking
        * \param os Output stream
        * \return Output stream
        */
        template<typename T>
        std::ostream &displn(const T &x, const std::string &separator,
                const std::string &start = "[", const std::string &end = "]",
                std::ostream &os = std::cout)
        {
            disp(x, separator, start, end, os);
            os << std::endl;
            return os;
        }

        /**
        * \brief Displays a C-style array. Does not add a newline.
        *
        * \see qpp::experimental::displn()
        *
        * \param x Pointer to the first element
        * \param n Number of elements to be displayed
        * \param separator Separator
        * \param start Left marking
        * \param end Right marking
        * \param os Output stream
        * \return Output stream
        */
        template<typename T>
        std::ostream &disp(const T *x, const std::size_t n,
                const std::string &separator, const std::string &start = "[",
                const std::string &end = "]", std::ostream &os = std::cout)
        {
            os << start;

            for (std::size_t i = 0; i < n - 1; ++i)
                os << x[i] << separator;
            if (n > 0)
                os << x[n - 1];

            os << end;
            return os;
        }

        /**
        * \brief Displays a C-style array. Adds a newline.
        *
        * \see qpp::experimental::disp()
        *
        * \param x Pointer to the first element
        * \param n Number of elements to be displayed
        * \param separator Separator
        * \param start Left marking
        * \param end Right marking
        * \param os Output stream
        * \return Output stream
        */
        template<typename T>
        std::ostream &displn(const T *x, const std::size_t n,
                const std::string &separator, const std::string &start = "[",
                const std::string &end = "]", std::ostream &os = std::cout)
        {
            disp(x, n, separator, start, end, os);
            os << std::endl;
            return os;
        }

        /**
        * \brief Displays an Eigen expression in matrix friendly form. Does not add a
        * new line.
        *
        * \see qpp::experimental::displn()
        *
        * \param A Eigen expression
        * \param chop Set to zero the elements smaller in absolute value
        * than \a chop
        * \param os Output stream
        * \return Output stream
        */
        template<typename Derived>
        std::ostream &disp(const Eigen::MatrixBase<Derived> &A, double chop = qpp::chop,
                std::ostream &os = std::cout)
        {
            const DynMat<typename Derived::Scalar> &rA = A;

            if (rA.size() == 0)
            {
                os << "Empty [" << rA.rows() << " x " << rA.cols() << "] matrix";
                return os;
            };

            std::ostringstream ostr;
            ostr.copyfmt(os); // copy os' state

            std::vector<std::string> vstr;
            std::string strA;

            for (std::size_t i = 0; i < static_cast<std::size_t>(rA.rows()); ++i)
            {
                for (std::size_t j = 0; j < static_cast<std::size_t>(rA.cols()); ++j)
                {
                    strA.clear(); // clear the temporary string
                    ostr.clear();
                    ostr.str(std::string {}); // clear the ostringstream

                    // convert to complex
                    double re = static_cast<cplx>(rA(i, j)).real();
                    double im = static_cast<cplx>(rA(i, j)).imag();

                    if (std::abs(re) < chop && std::abs(im) < chop)
                    {
                        ostr << "0 "; // otherwise segfault on destruction
                        // if using only vstr.push_back("0 ");
                        // bug in MATLAB's libmx
                        vstr.push_back(ostr.str());
                    }
                    else if (std::abs(re) < chop)
                    {
                        ostr << im;
                        vstr.push_back(ostr.str() + "i");
                    }
                    else if (std::abs(im) < chop)
                    {
                        ostr << re;
                        vstr.push_back(ostr.str() + " ");
                    }
                    else
                    {
                        ostr << re;
                        strA = ostr.str();

                        strA += (im > 0 ? " + " : " - ");
                        ostr.clear();
                        ostr.str(std::string()); // clear
                        ostr << std::abs(im);
                        strA += ostr.str();
                        strA += "i";
                        vstr.push_back(strA);
                    }
                }
            }

            // determine the maximum lenght of the entries in each column
            std::vector<std::size_t> maxlengthcols(rA.cols(), 0);

            for (std::size_t i = 0; i < static_cast<std::size_t>(rA.rows()); ++i)
                for (std::size_t j = 0; j < static_cast<std::size_t>(rA.cols()); ++j)
                    if (vstr[i * rA.cols() + j].size() > maxlengthcols[j])
                        maxlengthcols[j] = vstr[i * rA.cols() + j].size();

            // finally display it!
            for (std::size_t i = 0; i < static_cast<std::size_t>(rA.rows()); ++i)
            {
                os << std::setw(static_cast<int>(maxlengthcols[0])) << std::right
                        << vstr[i * rA.cols()]; // display first column
                // then the rest
                for (std::size_t j = 1; j < static_cast<std::size_t>(rA.cols()); ++j)
                    os << std::setw(static_cast<int>(maxlengthcols[j] + 2))
                            << std::right << vstr[i * rA.cols() + j];

                if (i < static_cast<std::size_t>(rA.rows()) - 1)
                    os << std::endl;
            }
            return os;
        }

        /**
        * \brief Displays an Eigen expression in matrix friendly form. Adds a newline.
        *
        * \see qpp::experimental::disp()
        *
        * \param A Eigen expression
        * \param chop Set to zero the elements smaller in absolute value
        * than \a chop
        * \param os Output stream
        * \return Output stream
        */
        template<typename Derived>
        std::ostream &displn(const Eigen::MatrixBase<Derived> &A, double chop = qpp::chop,
                std::ostream &os = std::cout)
        {
            disp(A, chop, os);
            os << std::endl;
            return os;
        }

        /**
        * \brief Displays a number (implicitly converted to std::complex<double>)
        * in friendly form. Does not add a new line.
        *
        * \see qpp::experimental::displn()
        *
        * \param z Real/complex number
        * \param chop Set to zero the elements smaller in absolute value
        * than \a chop
        * \param os Output stream
        * \return Output stream
        */
        std::ostream &disp(const cplx z, double chop = qpp::chop, std::ostream &os = std::cout)
        {
            // put the complex number inside an Eigen matrix
            cmat A(1, 1);
            A(0, 0) = z;
            disp(A, chop, os);
            return os;
        }

        /**
        * \brief Displays a number (implicitly converted to std::complex<double>)
        * in friendly form. Adds a new line.
        *
        * \see qpp::experimental::disp()
        *
        * \param z Real/complex number
        * \param chop Set to zero the elements smaller in absolute value
        * than \a chop
        * \param os Output stream
        * \return Output stream
        */
        std::ostream &displn(const cplx z, double chop = qpp::chop, std::ostream &os = std::cout)
        {
            disp(z, chop, os);
            os << std::endl;
            return os;
        }

        /**
        * \brief Applies the controlled-gate \a A to the part \a subsys
        * of a multi-partite state vector or density matrix
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
        DynMat<typename Derived1::Scalar> applyCTRL(
                const Eigen::MatrixBase<Derived1> &state,
                const Eigen::MatrixBase<Derived2> &A,
                const std::vector<std::size_t> &ctrl,
                const std::vector<std::size_t> &subsys,
                const std::vector<std::size_t> &dims)
        {
            const DynMat<typename Derived1::Scalar> &rstate = state;
            const DynMat<typename Derived2::Scalar> &rA = A;

            // EXCEPTION CHECKS
            // check types
            if (!std::is_same<typename Derived1::Scalar, typename Derived2::Scalar>::value)
                throw Exception("qpp::experimental::applyCTRL()", Exception::Type::TYPE_MISMATCH);

            // check zero sizes
            if (!internal::_check_nonzero_size(rA))
                throw Exception("qpp::experimental::applyCTRL()", Exception::Type::ZERO_SIZE);

            // check zero sizes
            if (!internal::_check_nonzero_size(rstate))
                throw Exception("qpp::experimental::applyCTRL()", Exception::Type::ZERO_SIZE);

            // check square matrix for the gate
            if (!internal::_check_square_mat(rA))
                throw Exception("qpp::experimental::applyCTRL()", Exception::Type::MATRIX_NOT_SQUARE);

            // check that all control subsystems have the same dimension
            std::size_t d = ctrl.size() > 0 ? ctrl[0] : 1;
            for (std::size_t i = 1; i < ctrl.size(); ++i)
                if (ctrl[i] != d)
                    throw Exception("qpp::experimental::applyCTRL()", Exception::Type::DIMS_NOT_EQUAL);

            // check that dimension is valid
            if (!internal::_check_dims(dims))
                throw Exception("qpp::experimental::applyCTRL()", Exception::Type::DIMS_INVALID);

            // check subsys is valid w.r.t. dims
            if (!internal::_check_subsys_match_dims(subsys, dims))
                throw Exception("qpp::experimental::applyCTRL()", Exception::Type::SUBSYS_MISMATCH_DIMS);

            // check that gate matches the dimensions of the subsys
            std::vector<std::size_t> subsys_dims(subsys.size());
            for (std::size_t i = 0; i < subsys.size(); ++i)
                subsys_dims[i] = dims[subsys[i]];
            if (!internal::_check_dims_match_mat(subsys_dims, rA))
                throw Exception("qpp::experimental::applyCTRL()", Exception::Type::MATRIX_MISMATCH_SUBSYS);

            std::vector<std::size_t> ctrlgate = ctrl; // ctrl + gate subsystem vector
            ctrlgate.insert(std::end(ctrlgate), std::begin(subsys), std::end(subsys));
            std::sort(std::begin(ctrlgate), std::end(ctrlgate));

            // check that ctrl + gate subsystem is valid
            // with respect to local dimensions
            if (!internal::_check_subsys_match_dims(ctrlgate, dims))
                throw Exception("qpp::experimental::applyCTRL()", Exception::Type::SUBSYS_MISMATCH_DIMS);

            // END EXCEPTION CHECKS

            // construct the table of A^i and (A^dagger)^i
            std::vector<DynMat<typename Derived1::Scalar>> Ai;
            std::vector<DynMat<typename Derived1::Scalar>> Aidagger;
            for (std::size_t i = 0; i < std::max(d, (std::size_t) 2); ++i)
            {
                Ai.push_back(powm(rA, i));
                Aidagger.push_back(powm(adjoint(rA), i));
            }

            std::size_t D = rstate.rows(); // total dimension
            std::size_t n = dims.size();   // total number of subsystems

            std::size_t ctrlsize = ctrl.size();

            std::size_t DA = rA.rows();
            std::size_t DCTRLAbar = static_cast<std::size_t>(D / (DA * std::pow(d, ctrl.size())));

            std::size_t Cdims[maxn]; // local dimensions
            std::size_t CdimsA[maxn]; // local dimensions
            std::size_t CdimsCTRLAbar[maxn]; // local dimensions

            std::vector<std::size_t> ctrlgatebar(n - ctrlgate.size()); // rest
            std::vector<std::size_t> allsubsys(n); // all subsystems
            std::iota(std::begin(allsubsys), std::end(allsubsys), 0);
            // compute the complementary subsystem of ctrlgate w.r.t. dims
            std::set_difference(std::begin(allsubsys), std::end(allsubsys),
                    std::begin(ctrlgate), std::end(ctrlgate), std::begin(ctrlgatebar));

            for (std::size_t k = 0; k < n; ++k)
                Cdims[k] = dims[k];
            for (std::size_t k = 0; k < subsys.size(); ++k)
                CdimsA[k] = dims[subsys[k]];
            for (std::size_t k = 0; k < n - ctrlgate.size(); ++k)
                CdimsCTRLAbar[k] = dims[ctrlgatebar[k]];

            // worker, computes the coefficient and the index for the ket case
            // used in #pragma omp parallel for collapse
            auto coeff_idx_ket = [=](std::size_t _i, std::size_t _m, std::size_t _r)
                    -> std::pair<typename Derived1::Scalar, std::size_t>
            {
                std::size_t idx = 0;
                typename Derived1::Scalar coeff = 0;

                std::size_t Cmidx[maxn]; // the total multi-index
                std::size_t CmidxA[maxn];// the gate part multi-index
                std::size_t CmidxCTRLAbar[maxn];// the rest multi-index

                // compute the index

                // set the CTRL part
                for (std::size_t k = 0; k < ctrl.size(); ++k)
                {
                    Cmidx[ctrl[k]] = _i;
                }

                // set the rest
                internal::_n2multiidx(_r, n - ctrlgate.size(),
                        CdimsCTRLAbar, CmidxCTRLAbar);
                for (std::size_t k = 0; k < n - ctrlgate.size(); ++k)
                {
                    Cmidx[ctrlgatebar[k]] = CmidxCTRLAbar[k];
                }

                // set the A part
                internal::_n2multiidx(_m, subsys.size(), CdimsA, CmidxA);
                for (std::size_t k = 0; k < subsys.size(); ++k)
                {
                    Cmidx[subsys[k]] = CmidxA[k];
                }

                // we now got the total index
                idx = internal::_multiidx2n(Cmidx, n, Cdims);

                // compute the coefficient
                for (std::size_t _n = 0; _n < DA; ++_n)
                {
                    internal::_n2multiidx(_n, subsys.size(), CdimsA, CmidxA);
                    for (std::size_t k = 0; k < subsys.size(); ++k)
                    {
                        Cmidx[subsys[k]] = CmidxA[k];
                    }
                    coeff += Ai[_i](_m, _n) *
                            rstate(internal::_multiidx2n(Cmidx, n, Cdims));
                }

                return std::make_pair(coeff, idx);
            };

            // worker, computes the coefficient and the index
            // for the density matrix case
            // used in #pragma omp parallel for collapse
            auto coeff_idx_rho = [=](std::size_t _i1, std::size_t _m1,
                    std::size_t _r1, std::size_t _i2, std::size_t _m2,
                    std::size_t _r2)
                    -> std::tuple<typename Derived1::Scalar, std::size_t, std::size_t>
            {
                std::size_t idxrow = 0;
                std::size_t idxcol = 0;
                typename Derived1::Scalar coeff = 0;

                std::size_t Cmidxrow[maxn]; // the total row multi-index
                std::size_t Cmidxcol[maxn];// the total col multi-index
                std::size_t CmidxArow[maxn];// the gate part row multi-index
                std::size_t CmidxAcol[maxn];// the gate part col multi-index
                std::size_t CmidxCTRLAbarrow[maxn];// the rest row multi-index
                std::size_t CmidxCTRLAbarcol[maxn];// the rest col multi-index

                // compute the ket/bra indexes

                // set the CTRL part
                for (std::size_t k = 0; k < ctrl.size(); ++k)
                {
                    Cmidxrow[ctrl[k]] = _i1;
                    Cmidxcol[ctrl[k]] = _i2;
                }

                // set the rest
                internal::_n2multiidx(_r1, n - ctrlgate.size(),
                        CdimsCTRLAbar, CmidxCTRLAbarrow);
                internal::_n2multiidx(_r2, n - ctrlgate.size(),
                        CdimsCTRLAbar, CmidxCTRLAbarcol);
                for (std::size_t k = 0; k < n - ctrlgate.size(); ++k)
                {
                    Cmidxrow[ctrlgatebar[k]] = CmidxCTRLAbarrow[k];
                    Cmidxcol[ctrlgatebar[k]] = CmidxCTRLAbarcol[k];
                }

                // set the A part
                internal::_n2multiidx(_m1, subsys.size(), CdimsA, CmidxArow);
                internal::_n2multiidx(_m2, subsys.size(), CdimsA, CmidxAcol);
                for (std::size_t k = 0; k < subsys.size(); ++k)
                {
                    Cmidxrow[subsys[k]] = CmidxArow[k];
                    Cmidxcol[subsys[k]] = CmidxAcol[k];
                }

                // we now got the total row/col indexes
                idxrow = internal::_multiidx2n(Cmidxrow, n, Cdims);
                idxcol = internal::_multiidx2n(Cmidxcol, n, Cdims);

                // compute the coefficient
                for (std::size_t _n1 = 0; _n1 < DA; ++_n1)
                {
                    internal::_n2multiidx(_n1, subsys.size(), CdimsA, CmidxArow);
                    for (std::size_t k = 0; k < subsys.size(); ++k)
                    {
                        Cmidxrow[subsys[k]] = CmidxArow[k];
                    }
                    for (std::size_t _n2 = 0; _n2 < DA; ++_n2)
                    {
                        internal::_n2multiidx(_n2, subsys.size(), CdimsA, CmidxAcol);
                        for (std::size_t k = 0; k < subsys.size(); ++k)
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
            if (internal::_check_col_vector(rstate)) // we have a ket
            {
                // check that dims match state vector
                if (!internal::_check_dims_match_cvect(dims, rstate))
                    throw Exception("qpp::experimental::applyCTRL()", Exception::Type::DIMS_MISMATCH_CVECTOR);

                if (D == 1)
                    return rstate;

                DynMat<typename Derived1::Scalar> result = rstate;

#pragma omp parallel for collapse(2)
                for (std::size_t m = 0; m < DA; ++m)
                    for (std::size_t r = 0; r < DCTRLAbar; ++r)
                        if (ctrlsize == 0) // no control
                        {
                            result(coeff_idx_ket(1, m, r).second) = coeff_idx_ket(1, m,
                                    r).first;
                        }
                        else
                            for (std::size_t i = 0; i < d; ++i)
                            {
                                result(coeff_idx_ket(i, m, r).second) = coeff_idx_ket(i,
                                        m, r).first;
                            }

                return result;
            }
                //************ density matrix ************//
            else if (internal::_check_square_mat(rstate)) // we have a density matrix
            {
                // check that dims match state matrix
                if (!internal::_check_dims_match_mat(dims, rstate))
                    throw Exception("qpp::experimental::applyCTRL()", Exception::Type::DIMS_MISMATCH_MATRIX);

                if (D == 1)
                    return rstate;

                DynMat<typename Derived1::Scalar> result = rstate;

#pragma omp parallel for collapse(4)
                for (std::size_t m1 = 0; m1 < DA; ++m1)
                    for (std::size_t r1 = 0; r1 < DCTRLAbar; ++r1)
                        for (std::size_t m2 = 0; m2 < DA; ++m2)
                            for (std::size_t r2 = 0; r2 < DCTRLAbar; ++r2)
                                if (ctrlsize == 0) // no control
                                {
                                    auto coeff_idxes = coeff_idx_rho(1, m1, r1, 1, m2,
                                            r2);
                                    result(std::get<1>(coeff_idxes),
                                            std::get<2>(coeff_idxes)) = std::get<0>(
                                            coeff_idxes);
                                }
                                else
                                {
                                    for (std::size_t i1 = 0; i1 < d; ++i1)
                                        for (std::size_t i2 = 0; i2 < d; ++i2)
                                        {
                                            auto coeff_idxes = coeff_idx_rho(i1, m1, r1,
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
                throw Exception("qpp::experimental::applyCTRL()",
                        Exception::Type::MATRIX_NOT_SQUARE_OR_CVECTOR);
        }

        /**
        * \brief Applies the gate \a A to the part \a subsys
        * of a multi-partite state vector or density matrix
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
        DynMat<typename Derived1::Scalar> apply(
                const Eigen::MatrixBase<Derived1> &state,
                const Eigen::MatrixBase<Derived2> &A,
                const std::vector<std::size_t> &subsys,
                const std::vector<std::size_t> &dims)
        {
            const DynMat<typename Derived1::Scalar> &rstate = state;
            const DynMat<typename Derived2::Scalar> &rA = A;

            // EXCEPTION CHECKS

            // check types
            if (!std::is_same<typename Derived1::Scalar, typename Derived2::Scalar>::value)
                throw Exception("qpp::experimental::apply()", Exception::Type::TYPE_MISMATCH);

            // check zero sizes
            if (!internal::_check_nonzero_size(rA))
                throw Exception("qpp::experimental::apply()", Exception::Type::ZERO_SIZE);

            // check zero sizes
            if (!internal::_check_nonzero_size(rstate))
                throw Exception("qpp::experimental::apply()", Exception::Type::ZERO_SIZE);

            // check square matrix for the gate
            if (!internal::_check_square_mat(rA))
                throw Exception("qpp::experimental::apply()", Exception::Type::MATRIX_NOT_SQUARE);

            // check that dimension is valid
            if (!internal::_check_dims(dims))
                throw Exception("qpp::experimental::apply()", Exception::Type::DIMS_INVALID);

            // check subsys is valid w.r.t. dims
            if (!internal::_check_subsys_match_dims(subsys, dims))
                throw Exception("qpp::experimental::apply()", Exception::Type::SUBSYS_MISMATCH_DIMS);

            // check that gate matches the dimensions of the subsys
            std::vector<std::size_t> subsys_dims(subsys.size());
            for (std::size_t i = 0; i < subsys.size(); ++i)
                subsys_dims[i] = dims[subsys[i]];
            if (!internal::_check_dims_match_mat(subsys_dims, rA))
                throw Exception("qpp::experimental::apply()", Exception::Type::MATRIX_MISMATCH_SUBSYS);

            if (internal::_check_col_vector(rstate)) // we have a ket
            {
                // check that dims match state vector
                if (!internal::_check_dims_match_cvect(dims, rstate))
                    throw Exception("qpp::experimental::apply()", Exception::Type::DIMS_MISMATCH_CVECTOR);

                return applyCTRL(rstate, rA, {}, subsys, dims);
            }
            else if (internal::_check_square_mat(rstate)) // we have a density matrix
            {

                // check that dims match state matrix
                if (!internal::_check_dims_match_mat(dims, rstate))
                    throw Exception("qpp::experimental::apply()", Exception::Type::DIMS_MISMATCH_MATRIX);

                return applyCTRL(rstate, rA, {}, subsys, dims);;
            }
            else
                throw Exception("qpp::experimental::apply()", Exception::Type::MATRIX_NOT_SQUARE_OR_CVECTOR);
        }

        /**
        * \brief Applies the channel specified by the set of Kraus operators \a Ks to
        * the part of the density matrix \a rho specified by \a subsys
        *
        * \param rho Eigen expression
        * \param Ks Set of Kraus operators
        * \param subsys Subsystems' indexes where the Kraus operators \a Ks are applied
        * \param dims Dimensions of the multi-partite system
        * \return Output density matrix after the action of the channel
        */
        template<typename Derived>
        cmat channel(const Eigen::MatrixBase<Derived> &rho, const std::vector<cmat> &Ks,
                const std::vector<std::size_t> &subsys, const std::vector<std::size_t> &dims)
        {
            const cmat &rrho = rho;

            // EXCEPTION CHECKS
            // check zero sizes
            if (!internal::_check_nonzero_size(rrho))
                throw Exception("qpp::experimental::channel()", Exception::Type::ZERO_SIZE);

            // check square matrix for the rho
            if (!internal::_check_square_mat(rrho))
                throw Exception("qpp::experimental::channel()", Exception::Type::MATRIX_NOT_SQUARE);

            // check that dimension is valid
            if (!internal::_check_dims(dims))
                throw Exception("qpp::experimental::channel()", Exception::Type::DIMS_INVALID);

            // check that dims match rho matrix
            if (!internal::_check_dims_match_mat(dims, rrho))
                throw Exception("qpp::experimental::channel()", Exception::Type::DIMS_MISMATCH_MATRIX);

            // check subsys is valid w.r.t. dims
            if (!internal::_check_subsys_match_dims(subsys, dims))
                throw Exception("qpp::experimental::channel()", Exception::Type::SUBSYS_MISMATCH_DIMS);

            std::vector<std::size_t> subsys_dims(subsys.size());
            for (std::size_t i = 0; i < subsys.size(); ++i)
                subsys_dims[i] = dims[subsys[i]];

            // check the Kraus operators
            if (!internal::_check_nonzero_size(Ks))
                throw Exception("qpp::experimental::channel()", Exception::Type::ZERO_SIZE);
            if (!internal::_check_square_mat(Ks[0]))
                throw Exception("qpp::experimental::channel()", Exception::Type::MATRIX_NOT_SQUARE);
            if (!internal::_check_dims_match_mat(subsys_dims, Ks[0]))
                throw Exception("qpp::experimental::channel()", Exception::Type::MATRIX_MISMATCH_SUBSYS);
            for (auto &&it : Ks)
                if (it.rows() != Ks[0].rows() || it.cols() != Ks[0].rows())
                    throw Exception("qpp::experimental::channel()", Exception::Type::DIMS_NOT_EQUAL);

            cmat result = cmat::Zero(rrho.rows(), rrho.rows());

            for (std::size_t i = 0; i < Ks.size(); ++i)
            {
                result += apply(rrho, Ks[i], subsys, dims);
            }
            return result;
        }

    } /* namespace experimental */
} /* namespace qpp */

#endif /* INCLUDE_EXPERIMENTAL_TEST_H_ */
