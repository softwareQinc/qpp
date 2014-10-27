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
 * @param state Eigen expression
 * @param A Eigen expression
 * @param ctrl Control subsystem indexes
 * @param subsys Subsystem indexes where the gate \a A is applied
 * @param n Total number of subsystes
 * @param d Local dimensions of all local Hilbert spaces (must all be equal)
 * @return CTRL-A gate applied to the part \a subsys of \a state
 */
template<typename Derived1, typename Derived2>
DynMat<typename Derived1::Scalar> applyCTRL(
		const Eigen::MatrixBase<Derived1>& state,
		const Eigen::MatrixBase<Derived2>& A,
		const std::vector<std::size_t>& ctrl,
		const std::vector<std::size_t>& subsys, std::size_t n,
		std::size_t d = 2)
{
}

/**
 * \brief Applies the gate \a A to the part \a subsys
 * of a multipartite state vector or density matrix
 *
 * \note The dimension of the gate \a A must match
 * the dimension of \a subsys
 *
 * @param state Eigen expression
 * @param A Eigen expression
 * @param subsys Subsystem indexes where the gate \a A is applied
 * @param dims Local dimensions of all local Hilbert spaces (can be different)
 * @return Gate \a A applied to the part \a subsys of \a state
 */
template<typename Derived1, typename Derived2>
DynMat<typename Derived1::Scalar> apply(
		const Eigen::MatrixBase<Derived1>& state,
		const Eigen::MatrixBase<Derived2>& A,
		const std::vector<std::size_t>& subsys,
		const std::vector<std::size_t>& dims)
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

	// check square matrix for the subsys
	if (!internal::_check_square_mat(rA))
		throw Exception("apply", Exception::Type::MATRIX_NOT_SQUARE);

	// check that dims is a valid dimension vector
	if (!internal::_check_dims(dims))
		throw Exception("apply", Exception::Type::DIMS_INVALID);

	// check subsys is valid w.r.t. dims
	if (!internal::_check_subsys_match_dims(subsys, dims))
		throw Exception("apply", Exception::Type::SUBSYS_MISMATCH_DIMS);

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

	for (std::size_t k = 0, cnt = 0; k < n; k++)
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
	for (std::size_t k = 0; k < nA; k++)
	{
		midxA_row[k] = midxA_rho_row[k] = 0;
		CdimsA[k] = dims[subsys[k]];
		CsubsysA[k] = subsys[k];
		DA *= dims[subsys[k]];
	}

	if (internal::_check_col_vector(rstate)) // we have a ket
	{
		// check that dims match state vector
		if (!internal::_check_dims_match_cvect(dims, rstate))
			throw Exception("apply", Exception::Type::DIMS_MISMATCH_CVECTOR);

		// check that state vector matches the dimensions of the subsys
		if (static_cast<std::size_t>(rA.cols()) != DA)
			throw Exception("apply", Exception::Type::DIMS_MISMATCH_CVECTOR);

		DynMat<typename Derived1::Scalar> result(D, 1);

		// run over the subsys's row multi-index
		for (std::size_t a = 0; a < DA; a++)
		{
			// get the subsys's row multi-index
			internal::_n2multiidx(a, nA, CdimsA, midxA_row);
			// compute subsys part of the result's row multi-index
			for (std::size_t k = 0; k < nA; k++)
				midx_row[CsubsysA[k]] = midxA_row[k];

			// run over the complement's row multi-index
			for (std::size_t i = 0; i < DA_bar; i++)
			{
				// get the complement's row multi-index
				internal::_n2multiidx(i, nA_bar, CdimsA_bar, midxA_bar_row);
				// now compute the complement part of the
				// result's row multi-index
				// and the complement part
				// of the state's total row multi-index
				for (std::size_t k = 0; k < nA_bar; k++)
					midx_row[CsubsysA_bar[k]] = midx_rho_row[CsubsysA_bar[k]] =
							midxA_bar_row[k];
				// compute the results's row index
				std::size_t result_row_idx = internal::_multiidx2n(midx_row, n,
						Cdims);

				// compute the coefficient
				typename Derived1::Scalar coeff = 0;
				for (std::size_t c = 0; c < DA; c++)
				{
					// compute the subsys part state's row multi-index
					internal::_n2multiidx(c, nA, CdimsA, midxA_rho_row);
					// now we have the total state's row multi-index
					for (std::size_t k = 0; k < nA; k++)
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
	else if (internal::_check_square_mat(rstate)) // we have a matrix
	{

		// check that dims match state matrix
		if (!internal::_check_dims_match_mat(dims, rstate))
			throw Exception("apply", Exception::Type::DIMS_MISMATCH_MATRIX);

		// check that state matrix matches the dimensions of the subsys
		if (static_cast<std::size_t>(rA.cols()) != DA)
			throw Exception("apply", Exception::Type::DIMS_MISMATCH_MATRIX);

		DynMat<typename Derived1::Scalar> result(D, D);

		// run over the subsys's row multi-index
		for (std::size_t a = 0; a < DA; a++)
		{
			// get the subsys's row multi-index
			internal::_n2multiidx(a, nA, CdimsA, midxA_row);
			// compute subsys part of the result's row multi-index
			for (std::size_t k = 0; k < nA; k++)
				midx_row[CsubsysA[k]] = midxA_row[k];

			// run over the complement's row multi-index
			for (std::size_t i = 0; i < DA_bar; i++)
			{
				// get the complement's row multi-index
				internal::_n2multiidx(i, nA_bar, CdimsA_bar, midxA_bar_row);
				// now compute the complement part
				// of the result's row multi-index
				// and the complement part of the
				// state's total row multi-index
				for (std::size_t k = 0; k < nA_bar; k++)
					midx_row[CsubsysA_bar[k]] = midx_rho_row[CsubsysA_bar[k]] =
							midxA_bar_row[k];
				// compute the results's row index
				std::size_t result_row_idx = internal::_multiidx2n(midx_row, n,
						Cdims);

				// run over the col index
				for (std::size_t j = 0; j < D; j++)
				{
					// compute the coefficient
					typename Derived1::Scalar coeff = 0;
					for (std::size_t c = 0; c < DA; c++)
					{
						// compute the subsys part state's row multi-index
						internal::_n2multiidx(c, nA, CdimsA, midxA_rho_row);
						// now we have the total state's row multi-index
						for (std::size_t k = 0; k < nA; k++)
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
	for (auto it : Ks)
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
	for (auto it : Ks)
		if (it.rows() != Ks[0].rows() || it.cols() != Ks[0].rows())
			throw Exception("channel", Exception::Type::DIMS_NOT_EQUAL);

	cmat result = cmat::Zero(rrho.rows(), rrho.cols());
	for (auto it : Ks)
		result += it * rrho * adjoint(it);

	return result;
}

/**
 * \brief Applies the channel specified by the set of Kraus operators \a Ks to
 * the part of the density matrix \a rho specified by \a subsys
 *
 * \param rho Eigen expression
 * \param Ks Set of Kraus operators
 * \param subsys Subsystems' indexes
 * \param dims Dimensions of the multi-partite system
 * \return Output density matrix after the action of the channel
 */
template<typename Derived>
cmat channel(const Eigen::MatrixBase<Derived>& rho, const std::vector<cmat>& Ks,
		const std::vector<std::size_t>& subsys,
		const std::vector<std::size_t>& dims)
{
	const cmat & rrho = rho;

	// EXCEPTION CHECKS
	// check zero sizes
	if (!internal::_check_nonzero_size(rrho))
		throw Exception("channel", Exception::Type::ZERO_SIZE);

	// check square matrix for the rho
	if (!internal::_check_square_mat(rrho))
		throw Exception("channel", Exception::Type::MATRIX_NOT_SQUARE);

	// check that dims is a valid dimension vector
	if (!internal::_check_dims(dims))
		throw Exception("channel", Exception::Type::DIMS_INVALID);

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
	for (auto it : Ks)
		if (it.rows() != Ks[0].rows() || it.cols() != Ks[0].rows())
			throw Exception("channel", Exception::Type::DIMS_NOT_EQUAL);

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

	for (std::size_t k = 0, cnt = 0; k < n; k++)
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
	for (std::size_t k = 0; k < nA; k++)
	{
		midxA_row[k] = midxA_col[k] = midxA_rho_row[k] = midxA_rho_col[k] = 0;
		CdimsA[k] = dims[subsys[k]];
		CsubsysA[k] = subsys[k];
		DA *= dims[subsys[k]];
	}

	// check that dimension of Kraus matches the dimension of the subsys
	if (static_cast<std::size_t>(Ks[0].rows()) != DA)
		throw Exception("channel", Exception::Type::DIMS_MISMATCH_MATRIX);

	// get the superoperator matrix of the channel
	cmat sop = super(Ks);

	cmat result(D, D);

	// run over rows
	for (std::size_t i = 0; i < D; i++)
	{
		// get the result's row multi-index
		internal::_n2multiidx(i, n, Cdims, midx_row);
		// get the subsys' complement row multi-index
		for (std::size_t k = 0; k < nA_bar; k++)
			midxA_bar_row[k] = midx_row[CsubsysA_bar[k]];
		// get the subsys' row multi-index
		for (std::size_t k = 0; k < nA; k++)
			midxA_row[k] = midx_row[CsubsysA[k]];

		// run over cols
		for (std::size_t j = 0; j < D; j++)
		{
			// get the result's col multi-index
			internal::_n2multiidx(j, n, Cdims, midx_col);
			// get the subsys' complement col multi-index
			for (std::size_t k = 0; k < nA_bar; k++)
				midxA_bar_col[k] = midx_col[CsubsysA_bar[k]];
			// get the subsys' col multi-index
			for (std::size_t k = 0; k < nA; k++)
				midxA_col[k] = midx_col[CsubsysA[k]];

			// now compute the coefficient
			cplx coeff = 0;
			for (std::size_t a = 0; a < DA; a++)
			{
				// get the subsys part of row multi-index for rho
				internal::_n2multiidx(a, nA, CdimsA, midxA_rho_row);
				for (std::size_t b = 0; b < DA; b++)
				{
					// get the subsys part of col multi-index for rho
					internal::_n2multiidx(b, nA, CdimsA, midxA_rho_col);

					// get the total row/col multi-index for rho
					for (std::size_t k = 0; k < nA; k++)
					{
						midx_rho_row[CsubsysA[k]] = midxA_rho_row[k];
						midx_rho_col[CsubsysA[k]] = midxA_rho_col[k];
					}
					for (std::size_t k = 0; k < nA_bar; k++)
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
	for (auto it : Ks)
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
	for (auto it : Ks)
		result += kron(cmat::Identity(D, D), it) * Omega
				* adjoint(kron(cmat::Identity(D, D), it));

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

} /* namespace qpp */

#endif /* OPERATIONS_H_ */
