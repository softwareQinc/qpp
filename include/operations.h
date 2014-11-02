/*
 * operations.h
 *
 *  Created on: Apr 6, 2014
 *      Author: vlad
 */

#ifndef OPERATIONS_H_
#define OPERATIONS_H_

namespace qpp
{

/**
 * \brief Applies the controlled-gate \a A to the part \a subsys
 * of a multipartite state vector or density matrix
 *
 * \note The dimension of the gate \a A must match
 * the dimension of \a subsys
 *
 * \param state Eigen expression
 * \param A Eigen expression
 * \param ctrl Control subsystem indexes
 * \param subsys Subsystem indexes where the gate \a A is applied
 * \param n Total number of subsystems
 * \param d Local dimensions of all local Hilbert spaces (must all be equal)
 * \return CTRL-A gate applied to the part \a subsys of \a state
 */
template<typename Derived1, typename Derived2>
DynMat<typename Derived1::Scalar> applyCTRL(
		const Eigen::MatrixBase<Derived1>& state,
		const Eigen::MatrixBase<Derived2>& A,
		const std::vector<std::size_t>& ctrl,
		const std::vector<std::size_t>& subsys, std::size_t n,
		std::size_t d = 2)
{
	const DynMat<typename Derived1::Scalar> & rstate = state;
	const DynMat<typename Derived2::Scalar> & rA = A;

	// EXCEPTION CHECKS
	// check types
	if (!std::is_same<typename Derived1::Scalar, typename Derived2::Scalar>::value)
		throw Exception("applyCTRL", Exception::Type::TYPE_MISMATCH);

	// check zero sizes
	if (!internal::_check_nonzero_size(rA))
		throw Exception("applyCTRL", Exception::Type::ZERO_SIZE);

	// check zero sizes
	if (!internal::_check_nonzero_size(rstate))
		throw Exception("applyCTRL", Exception::Type::ZERO_SIZE);

	// check square matrix for the gate
	if (!internal::_check_square_mat(rA))
		throw Exception("applyCTRL", Exception::Type::MATRIX_NOT_SQUARE);

	// check that gate matches the dimensions of the subsys
	if (static_cast<std::size_t>(rA.rows()) != std::pow(d, subsys.size()))
		throw Exception("applyCTRL", Exception::Type::DIMS_MISMATCH_MATRIX);

	// check out of range
	if (n == 0)
		throw Exception("applyCTRL", Exception::Type::OUT_OF_RANGE);

	// check that dimension is valid
	if (d == 0)
		throw Exception("applyCTRL", Exception::Type::DIMS_INVALID);

	std::vector<std::size_t> dims(n, d); // local dimensions vector

	// check subsys is valid w.r.t. dims
	if (!internal::_check_subsys_match_dims(subsys, dims))
		throw Exception("applyCTRL", Exception::Type::SUBSYS_MISMATCH_DIMS);

	std::vector<std::size_t> ctrlgate = ctrl; // ctrl + gate subsystem vector
	ctrlgate.insert(std::end(ctrlgate), std::begin(subsys), std::end(subsys));
	std::sort(std::begin(ctrlgate), std::end(ctrlgate));

	// check that ctrl + gate subsystem is valid
	// with respect to local dimensions
	if (!internal::_check_subsys_match_dims(ctrlgate, dims))
		throw Exception("applyCTRL", Exception::Type::SUBSYS_MISMATCH_DIMS);

	// END EXCEPTION CHECKS

	std::size_t ctrlsize = ctrl.size();

	// construct the table of A^i and (A^dagger)^i
	std::vector<DynMat<typename Derived1::Scalar>> Ai;
	std::vector<DynMat<typename Derived1::Scalar>> Aidagger;
	for (std::size_t i = 0; i < d; i++)
	{
		Ai.push_back(powm(rA, i));
		Aidagger.push_back(powm(adjoint(rA), i));
	}

	std::size_t D = static_cast<std::size_t>(std::pow(d, n));
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

	for (std::size_t k = 0; k < n; k++)
		Cdims[k] = d;
	for (std::size_t k = 0; k < subsys.size(); k++)
		CdimsA[k] = d;
	for (std::size_t k = 0; k < n - ctrlgate.size(); k++)
		CdimsCTRLAbar[k] = d;

	// worker, computes the coefficient and the index for the ket case
	// used in #pragma omp parallel for collapse
	auto coeff_idx_ket = [&](std::size_t _i, std::size_t _m, std::size_t _r)
	-> std::pair<typename Derived1::Scalar, std::size_t>
	{
		std::size_t idx = 0;
		typename Derived1::Scalar coeff = 0;

		std::size_t Cmidx[maxn]; // the total multi-index
			std::size_t CmidxA[maxn];// the gate multi-index
			std::size_t CmidxCTRLAbar[maxn];// the rest multi-index

			// compute the index

			// set the CTRL part
			for(std::size_t k = 0; k < ctrl.size(); k++)
			{
				Cmidx[ctrl[k]] = _i;
			}

			// set the rest
			internal::_n2multiidx(_r, n - ctrlgate.size(),
					CdimsCTRLAbar, CmidxCTRLAbar);
			for(std::size_t k = 0; k < n - ctrlgate.size(); k++)
			{
				Cmidx[ctrlgatebar[k]] = CmidxCTRLAbar[k];
			}

			// set the A part
			internal::_n2multiidx(_m, subsys.size(), CdimsA, CmidxA);
			for(std::size_t k = 0; k < subsys.size(); k++)
			{
				Cmidx[subsys[k]] = CmidxA[k];
			}

			idx = internal::_multiidx2n(Cmidx, n, Cdims); // we got the index

			// compute the coefficient
			for(std::size_t _n = 0; _n < DA; _n++)
			{
				internal::_n2multiidx(_n, subsys.size(), CdimsA, CmidxA);
				for(std::size_t k = 0; k < subsys.size(); k++)
				{
					Cmidx[subsys[k]] = CmidxA[k];
				}
				coeff += Ai[_i](_m,_n) *
				rstate(internal::_multiidx2n(Cmidx, n, Cdims));
			}

			return std::make_pair(coeff, idx);
		};

	// worker, computes the coefficient and the index
	// for the density matrix case
	// used in #pragma omp parallel for collapse
	auto coeff_idx_rho = [&](std::size_t _i1, std::size_t _m1,
			std::size_t _r1, std::size_t _i2, std::size_t _m2,
			std::size_t _r2 )
	-> std::tuple<typename Derived1::Scalar, std::size_t, std::size_t>
	{
		std::size_t idxrow = 0;
		std::size_t idxcol = 0;
		typename Derived1::Scalar coeff = 0;

		std::size_t Cmidxrow[maxn]; // the total row multi-index
			std::size_t Cmidxcol[maxn];// the total col multi-index
			std::size_t CmidxArow[maxn];// the gate row multi-index
			std::size_t CmidxAcol[maxn];// the gate col multi-index
			std::size_t CmidxCTRLAbarrow[maxn];// the rest row multi-index
			std::size_t CmidxCTRLAbarcol[maxn];// the rest col multi-index

			// compute the ket/bra indexes

			// set the CTRL part
			for(std::size_t k = 0; k < ctrl.size(); k++)
			{
				Cmidxrow[ctrl[k]] = _i1;
				Cmidxcol[ctrl[k]] = _i2;
			}

			// set the rest
			internal::_n2multiidx(_r1, n - ctrlgate.size(),
					CdimsCTRLAbar, CmidxCTRLAbarrow);
			internal::_n2multiidx(_r2, n - ctrlgate.size(),
					CdimsCTRLAbar, CmidxCTRLAbarcol);
			for(std::size_t k = 0; k < n - ctrlgate.size(); k++)
			{
				Cmidxrow[ctrlgatebar[k]] = CmidxCTRLAbarrow[k];
				Cmidxcol[ctrlgatebar[k]] = CmidxCTRLAbarcol[k];
			}

			// set the A part
			internal::_n2multiidx(_m1, subsys.size(), CdimsA, CmidxArow);
			internal::_n2multiidx(_m2, subsys.size(), CdimsA, CmidxAcol);
			for(std::size_t k = 0; k < subsys.size(); k++)
			{
				Cmidxrow[subsys[k]] = CmidxArow[k];
				Cmidxcol[subsys[k]] = CmidxAcol[k];
			}

			// we got the ket/bra indexes
			idxrow = internal::_multiidx2n(Cmidxrow, n, Cdims);
			idxcol = internal::_multiidx2n(Cmidxcol, n, Cdims);

			// compute the coefficient
			for(std::size_t _n1 = 0; _n1 < DA; _n1++)
			{
				internal::_n2multiidx(_n1, subsys.size(), CdimsA, CmidxArow);
				for(std::size_t k = 0; k < subsys.size(); k++)
				{
					Cmidxrow[subsys[k]] = CmidxArow[k];
				}
				for(std::size_t _n2 = 0; _n2 < DA; _n2++)
				{
					internal::_n2multiidx(_n2, subsys.size(), CdimsA, CmidxAcol);
					for(std::size_t k = 0; k < subsys.size(); k++)
					{
						Cmidxcol[subsys[k]] = CmidxAcol[k];
					}
					coeff += Ai[_i1](_m1,_n1)*
					rstate(internal::_multiidx2n(Cmidxrow, n, Cdims),
							internal::_multiidx2n(Cmidxcol, n, Cdims))*
					Aidagger[_i2](_n2,_m2);
				}
			}

			return std::make_tuple(coeff, idxrow, idxcol);
		};

	//************ ket ************//
	if (internal::_check_col_vector(rstate)) // we have a ket
	{
		// check that dims match state vector
		if (!internal::_check_dims_match_cvect(dims, rstate))
			throw Exception("applyCTRL",
					Exception::Type::DIMS_MISMATCH_CVECTOR);

		DynMat<typename Derived1::Scalar> result = DynMat<
				typename Derived1::Scalar>::Zero(D, 1);

#pragma omp parallel for collapse(2)
		for (std::size_t m = 0; m < DA; m++)
			for (std::size_t r = 0; r < DCTRLAbar; r++)
				if (ctrlsize == 0) // no control
				{
					result(coeff_idx_ket(1, m, r).second) = coeff_idx_ket(1, m,
							r).first;
				}
				else
					for (std::size_t i = 0; i < d; i++)
					{
						result(coeff_idx_ket(i, m, r).second) = coeff_idx_ket(i,
								m, r).first;
					}

		return result;
	}

	//************ density matrix ************//
	else if (internal::_check_square_mat(rstate)) // we have a matrix
	{
		// check that dims match state matrix
		if (!internal::_check_dims_match_mat(dims, rstate))
			throw Exception("applyCTRL", Exception::Type::DIMS_MISMATCH_MATRIX);

		DynMat<typename Derived1::Scalar> result = DynMat<
				typename Derived1::Scalar>::Zero(D, D);

#pragma omp parallel for collapse(4)
		for (std::size_t m1 = 0; m1 < DA; m1++)
			for (std::size_t r1 = 0; r1 < DCTRLAbar; r1++)
				for (std::size_t m2 = 0; m2 < DA; m2++)
					for (std::size_t r2 = 0; r2 < DCTRLAbar; r2++)
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
							for (std::size_t i1 = 0; i1 < d; i1++)
								for (std::size_t i2 = 0; i2 < d; i2++)

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
		throw Exception("applyCTRL",
				Exception::Type::MATRIX_NOT_SQUARE_OR_CVECTOR);
}

/**
 * \brief Applies the gate \a A to the part \a subsys
 * of a multipartite state vector or density matrix
 *
 * \note The dimension of the gate \a A must match
 * the dimension of \a subsys
 *
 * \param state Eigen expression
 * \param A Eigen expression
 * \param subsys Subsystem indexes where the gate \a A is applied
 * \param n Total number of subsystems
 * \param d Local dimensions of all local Hilbert spaces (must all be equal)
 * \return Gate \a A applied to the part \a subsys of \a state
 */
template<typename Derived1, typename Derived2>
DynMat<typename Derived1::Scalar> apply(
		const Eigen::MatrixBase<Derived1>& state,
		const Eigen::MatrixBase<Derived2>& A,
		const std::vector<std::size_t>& subsys, std::size_t n,
		std::size_t d = 2)
{
	const DynMat<typename Derived1::Scalar> & rstate = state;
	const DynMat<typename Derived2::Scalar> & rA = A;

	// EXCEPTION CHECKS

	// check types
	if (!std::is_same<typename Derived1::Scalar, typename Derived2::Scalar>::value)
		throw Exception("apply", Exception::Type::TYPE_MISMATCH);

	// check zero sizes
	if (!internal::_check_nonzero_size(rA))
		throw Exception("apply", Exception::Type::ZERO_SIZE);

	// check zero sizes
	if (!internal::_check_nonzero_size(rstate))
		throw Exception("apply", Exception::Type::ZERO_SIZE);

	// check square matrix for the gate
	if (!internal::_check_square_mat(rA))
		throw Exception("apply", Exception::Type::MATRIX_NOT_SQUARE);

	// check out of range
	if (n == 0)
		throw Exception("applyCTRL", Exception::Type::OUT_OF_RANGE);

	// check that dimension is valid
	if (d == 0)
		throw Exception("applyCTRL", Exception::Type::DIMS_INVALID);

	std::vector<std::size_t> dims(n, d); // local dimensions vector

	// check subsys is valid w.r.t. dims
	if (!internal::_check_subsys_match_dims(subsys, dims))
		throw Exception("apply", Exception::Type::SUBSYS_MISMATCH_DIMS);

	// check that gate matches the dimensions of the subsys
	if (static_cast<std::size_t>(rA.rows()) != std::pow(d, subsys.size()))
		throw Exception("apply", Exception::Type::DIMS_MISMATCH_MATRIX);

	if (internal::_check_col_vector(rstate)) // we have a ket
	{
		// check that dims match state vector
		if (!internal::_check_dims_match_cvect(dims, rstate))
			throw Exception("apply", Exception::Type::DIMS_MISMATCH_CVECTOR);

		return applyCTRL(rstate, rA, { }, subsys, n, d);
	}
	else if (internal::_check_square_mat(rstate)) // we have a matrix
	{

		// check that dims match state matrix
		if (!internal::_check_dims_match_mat(dims, rstate))
			throw Exception("apply", Exception::Type::DIMS_MISMATCH_MATRIX);

		return applyCTRL(rstate, rA, { }, subsys, n, d);;
	}
	else
		throw Exception("apply", Exception::Type::MATRIX_NOT_SQUARE_OR_CVECTOR);
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
cmat super(const std::vector<cmat>& Ks)
{
	// EXCEPTION CHECKS
	if (!internal::_check_nonzero_size(Ks))
		throw Exception("super", Exception::Type::ZERO_SIZE);
	if (!internal::_check_nonzero_size(Ks[0]))
		throw Exception("super", Exception::Type::ZERO_SIZE);
	if (!internal::_check_square_mat(Ks[0]))
		throw Exception("super", Exception::Type::MATRIX_NOT_SQUARE);
	for (auto&& it : Ks)
		if (it.rows() != Ks[0].rows() || it.cols() != Ks[0].rows())
			throw Exception("super", Exception::Type::DIMS_NOT_EQUAL);
	std::size_t D = static_cast<std::size_t>(Ks[0].rows());

	std::size_t midx_row[2] = { 0, 0 };
	std::size_t midx_col[2] = { 0, 0 };
	std::size_t dims[2];
	dims[0] = dims[1] = D;

	cmat result(D * D, D * D);
	cmat MN = cmat::Zero(D, D);
	cmat BA = cmat::Zero(D, D);
	cmat EMN = cmat::Zero(D, D);

	for (std::size_t m = 0; m < D; m++)
	{
		midx_col[0] = m;
		for (std::size_t n = 0; n < D; n++)
		{
			midx_col[1] = n;
			MN(m, n) = 1;
			// compute E(|m><n|)
			for (std::size_t i = 0; i < Ks.size(); i++)
				EMN += Ks[i] * MN * adjoint(Ks[i]);
			MN(m, n) = 0;
			for (std::size_t a = 0; a < D; a++)
			{
				midx_row[0] = a;
				for (std::size_t b = 0; b < D; b++)
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
 * \brief Applies the channel specified by the set of Kraus operators \a Ks
 * to the density matrix \a rho
 *
 * \param rho Eigen expression
 * \param Ks Set of Kraus operators
 * \return Output density matrix after the action of the channel
 */
template<typename Derived>
cmat channel(const Eigen::MatrixBase<Derived>& rho, const std::vector<cmat>& Ks)
{
	const cmat & rrho = rho;

	// EXCEPTION CHECKS
	if (!internal::_check_nonzero_size(rrho))
		throw Exception("channel", Exception::Type::ZERO_SIZE);
	if (!internal::_check_square_mat(rrho))
		throw Exception("channel", Exception::Type::MATRIX_NOT_SQUARE);
	if (!internal::_check_nonzero_size(Ks))
		throw Exception("channel", Exception::Type::ZERO_SIZE);
	if (!internal::_check_square_mat(Ks[0]))
		throw Exception("channel", Exception::Type::MATRIX_NOT_SQUARE);
	if (Ks[0].rows() != rrho.rows())
		throw Exception("channel", Exception::Type::DIMS_MISMATCH_MATRIX);
	for (auto&& it : Ks)
		if (it.rows() != Ks[0].rows() || it.cols() != Ks[0].rows())
			throw Exception("channel", Exception::Type::DIMS_NOT_EQUAL);

	cmat result = cmat::Zero(rrho.rows(), rrho.cols());

#pragma omp parallel for
	for (std::size_t i = 0; i < Ks.size(); i++)
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
 * the part of the density matrix \a rho specified by \a subsys
 *
 * \param rho Eigen expression
 * \param Ks Set of Kraus operators
 * \param subsys Subsystems' indexes
 * \param n Total number of subsystems
 * \param d Local dimensions of all local Hilbert spaces (must all be equal)
 * \return Output density matrix after the action of the channel
 */
template<typename Derived>
cmat channel(const Eigen::MatrixBase<Derived>& rho, const std::vector<cmat>& Ks,
		const std::vector<std::size_t>& subsys, std::size_t n,
		std::size_t d = 2)
{
	const cmat & rrho = rho;

	// EXCEPTION CHECKS
	// check zero sizes
	if (!internal::_check_nonzero_size(rrho))
		throw Exception("channel", Exception::Type::ZERO_SIZE);

	// check square matrix for the rho
	if (!internal::_check_square_mat(rrho))
		throw Exception("channel", Exception::Type::MATRIX_NOT_SQUARE);

	// check out of range
	if (n == 0)
		throw Exception("applyCTRL", Exception::Type::OUT_OF_RANGE);

	// check that dimension is valid
	if (d == 0)
		throw Exception("applyCTRL", Exception::Type::DIMS_INVALID);

	std::vector<std::size_t> dims(n, d); // local dimensions vector

	// check that dims match rho matrix
	if (!internal::_check_dims_match_mat(dims, rrho))
		throw Exception("channel", Exception::Type::DIMS_MISMATCH_MATRIX);

	// check subsys is valid w.r.t. dims
	if (!internal::_check_subsys_match_dims(subsys, dims))
		throw Exception("channel", Exception::Type::SUBSYS_MISMATCH_DIMS);

	// check the Kraus operators
	if (!internal::_check_nonzero_size(Ks))
		throw Exception("channel", Exception::Type::ZERO_SIZE);
	if (!internal::_check_square_mat(Ks[0]))
		throw Exception("channel", Exception::Type::MATRIX_NOT_SQUARE);
	for (auto&& it : Ks)
		if (it.rows() != Ks[0].rows() || it.cols() != Ks[0].rows())
			throw Exception("channel", Exception::Type::DIMS_NOT_EQUAL);

	cmat result = cmat::Zero(std::pow(d, n), std::pow(d, n));

#pragma omp parallel for
	for (size_t i = 0; i < Ks.size(); i++)
	{
#pragma omp critical
		{
			result += apply(rrho, Ks[i], subsys, n, d);
		}
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
cmat choi(const std::vector<cmat>& Ks)
{
	// EXCEPTION CHECKS
	if (!internal::_check_nonzero_size(Ks))
		throw Exception("choi", Exception::Type::ZERO_SIZE);
	if (!internal::_check_nonzero_size(Ks[0]))
		throw Exception("choi", Exception::Type::ZERO_SIZE);
	if (!internal::_check_square_mat(Ks[0]))
		throw Exception("choi", Exception::Type::MATRIX_NOT_SQUARE);
	for (auto&& it : Ks)
		if (it.rows() != Ks[0].rows() || it.cols() != Ks[0].rows())
			throw Exception("choi", Exception::Type::DIMS_NOT_EQUAL);
	std::size_t D = static_cast<std::size_t>(Ks[0].rows());

	// construct the D x D \sum |jj> vector
	// (un-normalized maximally entangled state)
	cmat MES = cmat::Zero(D * D, 1);
	for (std::size_t a = 0; a < D; a++)
		MES(a * D + a) = 1;

	cmat Omega = static_cast<cmat>(MES * adjoint(MES));

	cmat result = cmat::Zero(D * D, D * D);

#pragma omp parallel for
	for (std::size_t i = 0; i < Ks.size(); i++)
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
 * \brief Extracts orthogonal Kraus operators from Choi matrix
 *
 * Extracts a set of orthogonal (under Hilbert-Schmidt operator norm) Kraus
 * operators from the Choi representation \a A of the channel
 *
 * \note The Kraus operators satisfy \f$Tr(K_i^\dagger K_j)=\delta_{ij}\f$
 * for all \f$i\neq j\f$
 *
 * \param A Choi matrix
 * \return Set of Kraus operators
 */
std::vector<cmat> choi2kraus(const cmat& A)
{
	// EXCEPTION CHECKS
	if (!internal::_check_nonzero_size(A))
		throw Exception("choi2kraus", Exception::Type::ZERO_SIZE);
	if (!internal::_check_square_mat(A))
		throw Exception("choi2kraus", Exception::Type::MATRIX_NOT_SQUARE);
	std::size_t D = static_cast<std::size_t>(std::sqrt((double) A.rows()));
	if (D * D != static_cast<std::size_t>(A.rows()))
		throw Exception("choi2kraus", Exception::Type::DIMS_INVALID);

	dmat ev = hevals(A);
	cmat evec = hevects(A);
	std::vector<cmat> result;

	for (std::size_t i = 0; i < D * D; i++)
	{
		// take the absolute value to get rid of tiny negatives
		if (std::abs((double) ev(i)) > eps)
			result.push_back(
					(cmat) (std::sqrt((double) ev(i))
							* reshape((cmat) evec.col(i), D, D)));
	}
	return result;
}

/**
 * \brief Partial trace
 *
 *  Partial trace of density matrix
 *  over the first subsystem in a bi-partite system
 *
 * \param A Eigen expression
 * \param dims Dimensions of bi-partite system
 * (must be a std::vector with 2 elements)
 * \return Partial trace \f$Tr_{A}(\cdot)\f$ over the first subsytem \f$A\f$
 * in a bi-partite system \f$A\otimes B\f$, as a dynamic matrix
 * over the same scalar field as \a A
 */
template<typename Derived>
DynMat<typename Derived::Scalar> ptrace1(const Eigen::MatrixBase<Derived>& A,
		const std::vector<std::size_t>& dims)
{
	const DynMat<typename Derived::Scalar> & rA = A;

// Error checks

// check zero-size
	if (!internal::_check_nonzero_size(rA))
		throw Exception("ptrace1", Exception::Type::ZERO_SIZE);

// check that dims is a valid dimension vector
	if (!internal::_check_dims(dims))
		throw Exception("ptrace1", Exception::Type::DIMS_INVALID);

// check square matrix
	if (!internal::_check_square_mat(rA))
		throw Exception("ptrace1", Exception::Type::MATRIX_NOT_SQUARE);

// check dims has only 2 elements
	if (dims.size() != 2)
		throw Exception("ptrace1", Exception::Type::NOT_BIPARTITE);

// check that dims match the dimension of A
	if (!internal::_check_dims_match_mat(dims, rA))
		throw Exception("ptrace1", Exception::Type::DIMS_MISMATCH_MATRIX);

	std::size_t DA = dims[0];
	std::size_t DB = dims[1];

	DynMat<typename Derived::Scalar> result =
			DynMat<typename Derived::Scalar>::Zero(DB, DB);

	auto worker = [&](std::size_t i, std::size_t j)
	{
		typename Derived::Scalar sum = 0;
		for (std::size_t m = 0; m < DA; m++)
		sum += rA(m * DB + i, m * DB + j);
		return sum;
	};

#pragma omp parallel for collapse(2)
	for (std::size_t j = 0; j < DB; j++) // column major order for speed
		for (std::size_t i = 0; i < DB; i++)
			result(i, j) = worker(i, j);

	return result;
}

/**
 * \brief Partial trace
 *
 * \param A Eigen expression
 * \param dims Dimensions of bi-partite system
 * (must be a std::vector with 2 elements)
 * \return Partial trace \f$Tr_{B}(\cdot)\f$ over the second subsytem \f$B\f$
 * in a bi-partite system \f$A\otimes B\f$, as a dynamic matrix
 * over the same scalar field as \a A
 */
template<typename Derived>
DynMat<typename Derived::Scalar> ptrace2(const Eigen::MatrixBase<Derived>& A,
		const std::vector<std::size_t>& dims)
{
	const DynMat<typename Derived::Scalar> & rA = A;

// Error checks

// check zero-size
	if (!internal::_check_nonzero_size(rA))
		throw Exception("ptrace2", Exception::Type::ZERO_SIZE);

// check that dims is a valid dimension vector
	if (!internal::_check_dims(dims))
		throw Exception("ptrace2", Exception::Type::DIMS_INVALID);

// check square matrix
	if (!internal::_check_square_mat(rA))
		throw Exception("ptrace2", Exception::Type::MATRIX_NOT_SQUARE);

// check dims has only 2 elements
	if (dims.size() != 2)
		throw Exception("ptrace2", Exception::Type::NOT_BIPARTITE);

// check that dims match the dimension of A
	if (!internal::_check_dims_match_mat(dims, rA))
		throw Exception("ptrace2", Exception::Type::DIMS_MISMATCH_MATRIX);

	std::size_t DA = dims[0];
	std::size_t DB = dims[1];

	DynMat<typename Derived::Scalar> result =
			DynMat<typename Derived::Scalar>::Zero(DA, DA);

#pragma omp parallel for collapse(2)
	for (std::size_t j = 0; j < DA; j++) // column major order for speed
		for (std::size_t i = 0; i < DA; i++)
			result(i, j) = trace(rA.block(i * DB, j * DB, DB, DB));

	return result;
}

/**
 * \brief Partial trace
 *
 *  Partial trace of the multi-partite density matrix
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
DynMat<typename Derived::Scalar> ptrace(const Eigen::MatrixBase<Derived>& A,
		const std::vector<std::size_t>& subsys,
		const std::vector<std::size_t>& dims)

{
	const DynMat<typename Derived::Scalar> & rA = A;

// error checks

// check zero-size
	if (!internal::_check_nonzero_size(rA))
		throw Exception("ptrace", Exception::Type::ZERO_SIZE);

// check that dims is a valid dimension vector
	if (!internal::_check_dims(dims))
		throw Exception("ptrace", Exception::Type::DIMS_INVALID);

// check square matrix
	if (!internal::_check_square_mat(rA))
		throw Exception("ptrace", Exception::Type::MATRIX_NOT_SQUARE);

// check that dims match the dimension of A
	if (!internal::_check_dims_match_mat(dims, rA))
		throw Exception("ptrace", Exception::Type::DIMS_MISMATCH_MATRIX);

	if (subsys.size() == dims.size())
	{
		DynMat<typename Derived::Scalar> result = DynMat<
				typename Derived::Scalar>(1, 1);
		result(0, 0) = rA.trace();
		return result;
	}
	if (subsys.size() == 0)
		return rA;
	// check that subsys are valid
	if (!internal::_check_subsys_match_dims(subsys, dims))
		throw Exception("ptrace", Exception::Type::SUBSYS_MISMATCH_DIMS);

	std::size_t D = static_cast<std::size_t>(rA.rows());
	std::size_t n = dims.size();
	std::size_t nsubsys = subsys.size();
	std::size_t nsubsysbar = n - nsubsys;
	std::size_t dimsubsys = 1;
	for (std::size_t i = 0; i < nsubsys; i++)
		dimsubsys *= dims[subsys[i]];
	std::size_t dimsubsysbar = D / dimsubsys;

	std::size_t Cdims[maxn];
	std::size_t Csubsys[maxn];
	std::size_t Cdimssubsys[maxn];
	std::size_t Csubsysbar[maxn];
	std::size_t Cdimssubsysbar[maxn];

	for (std::size_t i = 0; i < n; i++)
		Cdims[i] = dims[i];
	for (std::size_t i = 0; i < nsubsys; i++)
	{
		Csubsys[i] = subsys[i];
		Cdimssubsys[i] = dims[subsys[i]];
	}
	// construct the complement of subsys
	std::size_t cnt = 0;
	for (std::size_t i = 0; i < n; i++)
	{
		bool found = false;
		for (std::size_t m = 0; m < nsubsys; m++)
			if (subsys[m] == i)
			{
				found = true;
				break;
			}
		if (!found)
		{
			Csubsysbar[cnt] = i;
			Cdimssubsysbar[cnt] = dims[i];
			cnt++;
		}
	}

	DynMat<typename Derived::Scalar> result = DynMat<typename Derived::Scalar>(
			dimsubsysbar, dimsubsysbar);

	auto worker = [&](std::size_t i, std::size_t j)
	{
		// use static allocation for speed!

			std::size_t Cmidxrow[maxn];
			std::size_t Cmidxcol[maxn];
			std::size_t Cmidxrowsubsysbar[maxn];
			std::size_t Cmidxcolsubsysbar[maxn];
			std::size_t Cmidxsubsys[maxn];

			/* get the row/col multi-indexes of the complement */
			internal::_n2multiidx(i, nsubsysbar, Cdimssubsysbar, Cmidxrowsubsysbar);
			internal::_n2multiidx(j, nsubsysbar, Cdimssubsysbar, Cmidxcolsubsysbar);
			/* write them in the global row/col multi-indexes */
			for(std::size_t k=0;k<nsubsysbar;k++)
			{
				Cmidxrow[Csubsysbar[k]]=Cmidxrowsubsysbar[k];
				Cmidxcol[Csubsysbar[k]]=Cmidxcolsubsysbar[k];
			}
			typename Derived::Scalar sm = 0;
			for(std::size_t a=0; a<dimsubsys; a++)
			{
				// get the multi-index over which we do the summation
				internal::_n2multiidx(a, nsubsys, Cdimssubsys, Cmidxsubsys);
				// write it into the global row/col multi-indexes
				for(std::size_t k=0;k<nsubsys;k++)
				Cmidxrow[Csubsys[k]]=Cmidxcol[Csubsys[k]]=Cmidxsubsys[k];

				// now do the sum
				sm+= rA(internal::_multiidx2n(Cmidxrow,n,Cdims),
						internal::_multiidx2n(Cmidxcol,n,Cdims));
			}

			return sm;
		};

#pragma omp parallel for collapse(2)
	for (std::size_t i = 0; i < dimsubsysbar; i++)
		for (std::size_t j = 0; j < dimsubsysbar; j++)
			result(i, j) = worker(i, j);

	return result;
}

/**
 * \brief Partial transpose
 *
 *  Partial transpose of the multi-partite density matrix
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
DynMat<typename Derived::Scalar> ptranspose(const Eigen::MatrixBase<Derived>& A,
		const std::vector<std::size_t>& subsys,
		const std::vector<std::size_t>& dims)

{
	const DynMat<typename Derived::Scalar> & rA = A;

// error checks

// check zero-size
	if (!internal::_check_nonzero_size(rA))
		throw Exception("ptranspose", Exception::Type::ZERO_SIZE);

// check that dims is a valid dimension vector
	if (!internal::_check_dims(dims))
		throw Exception("ptranspose", Exception::Type::DIMS_INVALID);

// check square matrix
	if (!internal::_check_square_mat(rA))
		throw Exception("ptranspose", Exception::Type::MATRIX_NOT_SQUARE);

// check that dims match the dimension of A
	if (!internal::_check_dims_match_mat(dims, rA))
		throw Exception("ptranspose", Exception::Type::DIMS_MISMATCH_MATRIX);

	if (subsys.size() == dims.size())
		return rA.transpose();
	if (subsys.size() == 0)
		return rA;
// check that subsys are valid
	if (!internal::_check_subsys_match_dims(subsys, dims))
		throw Exception("ptranspose", Exception::Type::SUBSYS_MISMATCH_DIMS);

	std::size_t D = static_cast<std::size_t>(rA.rows());
	std::size_t numdims = dims.size();
	std::size_t numsubsys = subsys.size();
	std::size_t cdims[maxn];
	std::size_t midxcol[maxn];
	std::size_t csubsys[maxn];

// copy dims in cdims and subsys in csubsys
	for (std::size_t i = 0; i < numdims; i++)
		cdims[i] = dims[i];
	for (std::size_t i = 0; i < numsubsys; i++)
		csubsys[i] = subsys[i];

	DynMat<typename Derived::Scalar> result(D, D);

	auto worker = [&](std::size_t i, std::size_t j)
	{
		// use static allocation for speed!
			std::size_t midxcoltmp[maxn];
			std::size_t midxrow[maxn];

			for (std::size_t k = 0; k < numdims; k++)
			midxcoltmp[k] = midxcol[k];

			/* compute the row multi-index */
			internal::_n2multiidx(i, numdims, cdims, midxrow);

			for (std::size_t k = 0; k < numsubsys; k++)
			std::swap(midxcoltmp[csubsys[k]], midxrow[csubsys[k]]);

			/* writes the result */
			result(i, j)=rA(internal::_multiidx2n(midxrow, numdims, cdims),
					internal::_multiidx2n(midxcoltmp, numdims, cdims));

		};

	for (std::size_t j = 0; j < D; j++)
	{
		// compute the column multi-index
		internal::_n2multiidx(j, numdims, cdims, midxcol);
#pragma omp parallel for
		for (std::size_t i = 0; i < D; i++)
			worker(i, j);
	}

	return result;
}
/**
 * \brief System permutation
 *
 * Permutes the subsystems in a state vector or density matrix\n
 * The qubit \a perm[\a i] is permuted to the location \a i
 *
 * \param A Eigen expression
 * \param perm Permutation
 * \param dims Subsystems' dimensions
 * \return Permuted system, as a dynamic matrix
 * over the same scalar field as \a A
 */
template<typename Derived>
DynMat<typename Derived::Scalar> syspermute(const Eigen::MatrixBase<Derived>& A,
		const std::vector<std::size_t>& perm,
		const std::vector<std::size_t>& dims)

{
	const DynMat<typename Derived::Scalar> & rA = A;

// Error checks

	// check zero-size
	if (!internal::_check_nonzero_size(rA))
		throw Exception("syspermute", Exception::Type::ZERO_SIZE);

	// check that dims is a valid dimension vector
	if (!internal::_check_dims(dims))
		throw Exception("syspermute", Exception::Type::DIMS_INVALID);

	// check that we have a valid permutation
	if (!internal::_check_perm(perm))
		throw Exception("syspermute", Exception::Type::PERM_INVALID);

	// check permutation size
	if (perm.size() != dims.size())
		throw Exception("syspermute", Exception::Type::PERM_INVALID);

	std::size_t D = static_cast<std::size_t>(rA.rows());
	std::size_t numdims = dims.size();

	DynMat<typename Derived::Scalar> result;

	auto worker =
			[](std::size_t i, std::size_t numdims, const std::size_t* cdims,
					const std::size_t* cperm)
			{
				// use static allocation for speed,
				// double the size for matrices reshaped as vectors
				std::size_t midx[2 * maxn];
				std::size_t midxtmp[2 * maxn];
				std::size_t permdims[2 * maxn];

				/* compute the multi-index */
				internal::_n2multiidx(i, numdims, cdims, midx);

				for (std::size_t k = 0; k < numdims; k++)
				{
					permdims[k] = cdims[cperm[k]]; // permuted dimensions
					midxtmp[k] = midx[cperm[k]];// permuted multi-indexes
				}
				return internal::_multiidx2n(midxtmp, numdims, permdims);
			};

// check column vector
	if (internal::_check_col_vector(rA)) // we have a column vector
	{
		std::size_t cdims[maxn];
		std::size_t cperm[maxn];

		// check that dims match the dimension of rA
		if (!internal::_check_dims_match_cvect(dims, rA))
			throw Exception("syspermute",
					Exception::Type::DIMS_MISMATCH_CVECTOR);

		// copy dims in cdims and perm in cperm
		for (std::size_t i = 0; i < numdims; i++)
		{
			cdims[i] = dims[i];
			cperm[i] = perm[i];
		}
		result.resize(D, 1);

#pragma omp parallel for
		for (std::size_t i = 0; i < D; i++)
			result(worker(i, numdims, cdims, cperm)) = rA(i);

		return result;
	}

	else if (internal::_check_square_mat(rA)) // we have a square matrix
	{
		std::size_t cdims[2 * maxn];
		std::size_t cperm[2 * maxn];

		// check that dims match the dimension of rA
		if (!internal::_check_dims_match_mat(dims, rA))
			throw Exception("syspermute",
					Exception::Type::DIMS_MISMATCH_MATRIX);

		// copy dims in cdims and perm in cperm
		for (std::size_t i = 0; i < numdims; i++)
		{
			cdims[i] = dims[i];
			cdims[i + numdims] = dims[i];
			cperm[i] = perm[i];
			cperm[i + numdims] = perm[i] + numdims;
		}
		result.resize(D * D, 1);
		// map A to a column vector
		DynMat<typename Derived::Scalar> vectA = Eigen::Map<
				DynMat<typename Derived::Scalar>>(
				const_cast<typename Derived::Scalar*>(rA.data()), D * D, 1);

#pragma omp parallel for
		for (std::size_t i = 0; i < D * D; i++)
			result(worker(i, 2 * numdims, cdims, cperm)) = rA(i);

		return reshape(result, D, D);
	}

	else
		throw Exception("syspermute",
				Exception::Type::MATRIX_NOT_SQUARE_OR_CVECTOR);
}

} /* namespace qpp */

#endif /* OPERATIONS_H_ */
