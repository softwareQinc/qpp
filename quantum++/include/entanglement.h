/*
 * entanglement.h
 *
 *  Created on: Apr 13, 2014
 *      Author: vlad
 */

#ifndef INCLUDE_ENTANGLEMENT_H_
#define INCLUDE_ENTANGLEMENT_H_

// entanglement

namespace qpp
{

/**
 * \brief Schmidt coefficients of the bi-partite pure state \a A
 *
 * \note The sum of the squares of the Schmidt coefficients equals 1
 * \see \a qpp::schmidtprob()
 *
 * \param A Eigen expression
 * \param dims Subsystems' dimensions
 * \return Schmidt coefficients of \a A, as a complex dynamic matrix,
 *  with the Schmidt coefficients on the diagonal
 */
template<typename Derived>
cmat schmidtcoeff(const Eigen::MatrixBase<Derived>& A,
		const std::vector<std::size_t>& dims)
{
	const DynMat<typename Derived::Scalar> & rA = A;
	// check zero-size
	if (!internal::_check_nonzero_size(rA))
		throw Exception("schmidtcoeff", Exception::Type::ZERO_SIZE);
	// check bi-partite
	if (dims.size() != 2)
		throw Exception("schmidtcoeff", Exception::Type::NOT_BIPARTITE);
	// check column vector
	if (!internal::_check_col_vector(rA))
		throw Exception("schmidtcoeff", Exception::Type::MATRIX_NOT_CVECTOR);
	// check matching dimensions
	if (!internal::_check_dims_match_mat(dims, A))
		throw Exception("schmidtcoeff", Exception::Type::DIMS_MISMATCH_MATRIX);

	Eigen::JacobiSVD<DynMat<typename Derived::Scalar>> svd(
			transpose(reshape(rA, dims[1], dims[0])));
	return svd.singularValues().template cast<cplx>();
}

/**
 * \brief Schmidt basis on Alice's side
 *
 * \param A Eigen expression
 * \param dims Subsystems' dimensions
 * \return Unitary matrix \f$ U \f$ representing the Schmidt basis
 * on Alice's side, as a complex dynamic matrix, acting on
 * the computational basis as \f$ U|j\rangle = |\bar j\rangle\f$
 * (Schmidt vector)
 */
template<typename Derived>
cmat schmidtU(const Eigen::MatrixBase<Derived>& A,
		const std::vector<std::size_t>& dims)
{
	const DynMat<typename Derived::Scalar> & rA = A;
	// check zero-size
	if (!internal::_check_nonzero_size(rA))
		throw Exception("schmidtU", Exception::Type::ZERO_SIZE);
	// check bi-partite
	if (dims.size() != 2)
		throw Exception("schmidtU", Exception::Type::NOT_BIPARTITE);
	// check column vector
	if (!internal::_check_col_vector(rA))
		throw Exception("schmidtU", Exception::Type::MATRIX_NOT_CVECTOR);
	// check matching dimensions
	if (!internal::_check_dims_match_mat(dims, A))
		throw Exception("schmidtU", Exception::Type::DIMS_MISMATCH_MATRIX);

	Eigen::JacobiSVD<DynMat<typename Derived::Scalar>> svd(
			transpose(reshape(rA, dims[1], dims[0])),
			Eigen::DecompositionOptions::ComputeFullU);
	return svd.matrixU();
}

/**
 * \brief Schmidt basis on Bob's side
 *
 * \param A Eigen expression
 * \param dims Subsystems' dimensions
 * \return Unitary matrix \f$ V \f$ representing the Schmidt basis
 * on Bob's side, as a complex dynamic matrix, acting on
 * the computational basis as \f$ V|j\rangle = |\bar j\rangle\f$
 * (Schmidt vector)
 */
template<typename Derived>
cmat schmidtV(const Eigen::MatrixBase<Derived>& A,
		const std::vector<std::size_t>& dims)
{
	const DynMat<typename Derived::Scalar> & rA = A;
	// check zero-size
	if (!internal::_check_nonzero_size(rA))
		throw Exception("schmidtV", Exception::Type::ZERO_SIZE);
	// check bi-partite
	if (dims.size() != 2)
		throw Exception("schmidtV", Exception::Type::NOT_BIPARTITE);
	// check column vector
	if (!internal::_check_col_vector(rA))
		throw Exception("schmidtV", Exception::Type::MATRIX_NOT_CVECTOR);
	// check matching dimensions
	if (!internal::_check_dims_match_mat(dims, A))
		throw Exception("schmidtV", Exception::Type::DIMS_MISMATCH_MATRIX);

	Eigen::JacobiSVD<DynMat<typename Derived::Scalar>> svd(
			transpose(reshape(rA, dims[1], dims[0])),
			Eigen::DecompositionOptions::ComputeFullV);
	// by default returns V^*, we need V, i.e. the complex conjugate,
	// i.e. adjoint(transpose(V))
	return adjoint(transpose(svd.matrixV()));
}

/**
 * \brief Schmidt probabilities of the bi-partite pure state \a A
 *
 * Defined as the squares of the Schmidt coefficients\n
 * The sum of the Schmidt probabilities equals 1
 * \see \a qpp::schmidtcoeff()
 *
 * \param A Eigen expression
 * \param dims Subsystems' dimensions
 * \return Schmidt probabilites of \a A, as a complex dynamic matrix,
 * with the Schmidt probabilities on the diagonal
 */
template<typename Derived>
cmat schmidtprob(const Eigen::MatrixBase<Derived>& A,
		const std::vector<std::size_t>& dims)
{
	const DynMat<typename Derived::Scalar> & rA = A;
	// check zero-size
	if (!internal::_check_nonzero_size(rA))
		throw Exception("schmidtprob", Exception::Type::ZERO_SIZE);
	// check bi-partite
	if (dims.size() != 2)
		throw Exception("schmidtprob", Exception::Type::NOT_BIPARTITE);
	// check column vector
	if (!internal::_check_col_vector(rA))
		throw Exception("schmidtprob", Exception::Type::MATRIX_NOT_CVECTOR);
	// check matching dimensions
	if (!internal::_check_dims_match_mat(dims, A))
		throw Exception("schmidtprob", Exception::Type::DIMS_MISMATCH_MATRIX);

	Eigen::JacobiSVD<DynMat<typename Derived::Scalar>> svd(
			transpose(reshape(rA, dims[1], dims[0])));

	return powm(
			static_cast<cmat>((svd.singularValues().template cast<cplx>()).asDiagonal()),
			2).diagonal();
}

/**
 * \brief Entanglement of the bi-partite pure state \a A
 *
 * Defined as the von-Neumann entropy of the reduced density matrix
 * of one of the subsystems
 * \see qpp::shannon()
 *
 * \param A Eigen expression
 * \param dims Subsystems' dimensions
 * \return Entanglement, with the logarithm in base 2
 */
template<typename Derived>
double entanglement(const Eigen::MatrixBase<Derived>& A,
		const std::vector<std::size_t>& dims)
{
	const DynMat<typename Derived::Scalar> & rA = A;
	// check zero-size
	if (!internal::_check_nonzero_size(rA))
		throw Exception("entanglement", Exception::Type::ZERO_SIZE);
	// check bi-partite
	if (dims.size() != 2)
		throw Exception("entanglement", Exception::Type::NOT_BIPARTITE);
	// check column vector
	if (!internal::_check_col_vector(rA))
		throw Exception("entanglement", Exception::Type::MATRIX_NOT_CVECTOR);
	// check matching dimensions
	if (!internal::_check_dims_match_mat(dims, A))
		throw Exception("entanglement", Exception::Type::DIMS_MISMATCH_MATRIX);

	return shannon(schmidtprob(rA, dims));
}

/**
 * \brief G-concurrence of the bi-partite pure state \a A
 *
 * Uses \a qpp::logdet() to avoid overflows
 * \see qpp::logdet()
 *
 * \param A Eigen expression
 * \param dims Subsystems' dimensions
 * \return G-concurrence
 */
template<typename Derived> // the G-concurrence
double gconcurrence(const Eigen::MatrixBase<Derived>& A)
{
	const DynMat<typename Derived::Scalar> & rA = A;

	// check zero-size
	if (!internal::_check_nonzero_size(rA))
		throw Exception("gconcurrence", Exception::Type::ZERO_SIZE);
	// check column vector
	if (!internal::_check_col_vector(rA))
		throw Exception("gconcurrence", Exception::Type::MATRIX_NOT_CVECTOR);

	std::size_t D = static_cast<std::size_t>(std::sqrt((double) A.rows()));
	if (D * D != static_cast<std::size_t>(A.rows()))
		throw Exception("gconcurrence", Exception::Type::DIMS_NOT_EQUAL);

	// we compute exp(logdet()) to avoid underflow
	return D * std::abs(std::exp(2. / D * logdet(reshape(rA, D, D))));
}

} /* namespace qpp */

#endif /* INCLUDE_ENTANGLEMENT_H_ */
