/*
 * entanglement.h
 *
 *  Created on: Apr 13, 2014
 *      Author: vlad
 */

#ifndef ENTANGLEMENT_H_
#define ENTANGLEMENT_H_

// entanglement

namespace qpp
{

// schmidt coefficients
template<typename Derived>
types::cmat schmidtcoeff(const Eigen::MatrixBase<Derived>& A,
		const std::vector<std::size_t>& dims)
{
	const types::DynMat<typename Derived::Scalar> & rA = A;
	// check zero-size
	if (!internal::_check_nonzero_size(rA))
		throw Exception("schmidtcoeff", Exception::Type::ZERO_SIZE);
	// check bipartite
	if (dims.size() != 2)
		throw Exception("schmidtcoeff", Exception::Type::NOT_BIPARTITE);
	// check column vector
	if (!internal::_check_col_vector(rA))
		throw Exception("schmidtcoeff", Exception::Type::MATRIX_NOT_CVECTOR);
	// check matching dimensions
	if (!internal::_check_dims_match_mat(dims, A))
		throw Exception("schmidtcoeff", Exception::Type::DIMS_MISMATCH_MATRIX);

	Eigen::JacobiSVD<types::DynMat<typename Derived::Scalar>> svd(
			transpose(reshape(rA, dims[1], dims[0])));
	return svd.singularValues().template cast<types::cplx>();
}

// schmidt U (basis on Alice's side, i.e. U|j> = |\bar j> (schmidt vector))
template<typename Derived>
types::cmat schmidtU(const Eigen::MatrixBase<Derived>& A,
		const std::vector<std::size_t>& dims)
{
	const types::DynMat<typename Derived::Scalar> & rA = A;
	// check zero-size
	if (!internal::_check_nonzero_size(rA))
		throw Exception("schmidtU", Exception::Type::ZERO_SIZE);
	// check bipartite
	if (dims.size() != 2)
		throw Exception("schmidtU", Exception::Type::NOT_BIPARTITE);
	// check column vector
	if (!internal::_check_col_vector(rA))
		throw Exception("schmidtU", Exception::Type::MATRIX_NOT_CVECTOR);
	// check matching dimensions
	if (!internal::_check_dims_match_mat(dims, A))
		throw Exception("schmidtU", Exception::Type::DIMS_MISMATCH_MATRIX);

	Eigen::JacobiSVD<types::DynMat<typename Derived::Scalar>> svd(
			transpose(reshape(rA, dims[1], dims[0])),
			Eigen::DecompositionOptions::ComputeFullU);
	return svd.matrixU();
}

// schmidt V (basis on Bob's side, i.e. V|j> = |\bar j> (schmidt vector))
template<typename Derived>
types::cmat schmidtV(const Eigen::MatrixBase<Derived>& A,
		const std::vector<std::size_t>& dims)
{
	const types::DynMat<typename Derived::Scalar> & rA = A;
	// check zero-size
	if (!internal::_check_nonzero_size(rA))
		throw Exception("schmidtV", Exception::Type::ZERO_SIZE);
	// check bipartite
	if (dims.size() != 2)
		throw Exception("schmidtV", Exception::Type::NOT_BIPARTITE);
	// check column vector
	if (!internal::_check_col_vector(rA))
		throw Exception("schmidtV", Exception::Type::MATRIX_NOT_CVECTOR);
	// check matching dimensions
	if (!internal::_check_dims_match_mat(dims, A))
		throw Exception("schmidtV", Exception::Type::DIMS_MISMATCH_MATRIX);

	Eigen::JacobiSVD<types::DynMat<typename Derived::Scalar>> svd(
			transpose(reshape(rA, dims[1], dims[0])),
			Eigen::DecompositionOptions::ComputeFullV);
	// by default returns V^*, we need V, i.e. the complex conjugate,
	// i.e. adjoint(transpose(V))
	return adjoint(transpose(svd.matrixV()));
}

// schmidt probabilities (sum up to one)
template<typename Derived>
types::cmat schmidtprob(const Eigen::MatrixBase<Derived>& A,
		const std::vector<std::size_t>& dims)
{
	const types::DynMat<typename Derived::Scalar> & rA = A;
	// check zero-size
	if (!internal::_check_nonzero_size(rA))
		throw Exception("schmidtprob", Exception::Type::ZERO_SIZE);
	// check bipartite
	if (dims.size() != 2)
		throw Exception("schmidtprob", Exception::Type::NOT_BIPARTITE);
	// check column vector
	if (!internal::_check_col_vector(rA))
		throw Exception("schmidtprob", Exception::Type::MATRIX_NOT_CVECTOR);
	// check matching dimensions
	if (!internal::_check_dims_match_mat(dims, A))
		throw Exception("schmidtprob", Exception::Type::DIMS_MISMATCH_MATRIX);

	Eigen::JacobiSVD<types::DynMat<typename Derived::Scalar>> svd(
			transpose(reshape(rA, dims[1], dims[0])));

	return powm(
			static_cast<types::cmat>((svd.singularValues().template cast<
					types::cplx>()).asDiagonal()), 2).diagonal();
}

template<typename Derived>
double entanglement(const Eigen::MatrixBase<Derived>& A,
		const std::vector<std::size_t>& dims)
{
	const types::DynMat<typename Derived::Scalar> & rA = A;
	// check zero-size
	if (!internal::_check_nonzero_size(rA))
		throw Exception("entanglement", Exception::Type::ZERO_SIZE);
	// check bipartite
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

template<typename Derived> // the G-concurrence
double gconcurrence(const Eigen::MatrixBase<Derived>& A)
{
	const types::DynMat<typename Derived::Scalar> & rA = A;

	// check zero-size
	if (!internal::_check_nonzero_size(rA))
		throw Exception("gconcurrence", Exception::Type::ZERO_SIZE);
	// check column vector
	if (!internal::_check_col_vector(rA))
		throw Exception("gconcurrence", Exception::Type::MATRIX_NOT_CVECTOR);

	std::size_t D = static_cast<std::size_t>((double)std::sqrt(A.rows()));
	if (D * D != static_cast<std::size_t>(A.rows()))
		throw Exception("gconcurrence", Exception::Type::DIMS_NOT_EQUAL);

	// we compute exp(logdet()) to avoid underflow
	return D * std::abs(std::exp(2. / D * logdet(reshape(rA, D, D))));
}

} /* namespace qpp */
#endif /* ENTANGLEMENT_H_ */
