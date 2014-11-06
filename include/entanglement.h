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
 * \return Schmidt coefficients of \a A, as a complex dynamic column vector
 */
template<typename Derived>
DynColVect<cplx> schmidtcoeff(const Eigen::MatrixBase<Derived>& A,
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
	if (!internal::_check_dims_match_mat(dims, rA))
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
 * \return Unitary matrix \f$ U \f$ whose columns represent
 * the Schmidt basis vectors on Alice's side. Acts on the computational basis
 * as \f$ U|j\rangle = |\bar j\rangle\f$, where \f$|\bar j\rangle\f$ denotes
 * the Schmidt vector.
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
	if (!internal::_check_dims_match_mat(dims, rA))
		throw Exception("schmidtU", Exception::Type::DIMS_MISMATCH_MATRIX);

	return svdU(transpose(reshape(rA, dims[1], dims[0])));
}

/**
 * \brief Schmidt basis on Bob's side
 *
 * \param A Eigen expression
 * \param dims Subsystems' dimensions
 * \return Unitary matrix \f$ V \f$ whose columns represent
 * the Schmidt basis vectors on Bob's side. Acts on the computational basis
 * as \f$ V|j\rangle = |\bar j\rangle\f$, where \f$|\bar j\rangle\f$ denotes
 * the Schmidt vector.
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
	if (!internal::_check_dims_match_mat(dims, rA))
		throw Exception("schmidtV", Exception::Type::DIMS_MISMATCH_MATRIX);

	// by default returns V^*, we need V, i.e. the complex conjugate,
	// i.e. adjoint(transpose(V))
	return svdV(transpose(reshape(conjugate(rA), dims[1], dims[0])));
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
 * \return Schmidt probabilites of \a A, as a real dynamic column vector
 */
template<typename Derived>
DynColVect<double> schmidtprob(const Eigen::MatrixBase<Derived>& A,
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
	if (!internal::_check_dims_match_mat(dims, rA))
		throw Exception("schmidtprob", Exception::Type::DIMS_MISMATCH_MATRIX);

	Eigen::JacobiSVD<DynMat<typename Derived::Scalar>> svd(
			transpose(reshape(rA, dims[1], dims[0])));

	return powm((dmat) svd.singularValues().asDiagonal(), 2).diagonal();
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
	if (!internal::_check_dims_match_mat(dims, rA))
		throw Exception("entanglement", Exception::Type::DIMS_MISMATCH_MATRIX);

	return shannon(schmidtprob(rA, dims));
}

/**
 * \brief G-concurrence of the bi-partite pure state \a A
 *
 * \note Both local dimensions must be equal
 *
 * Uses \a qpp::logdet() to avoid overflows
 * \see qpp::logdet()
 *
 * \param A Eigen expression
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

// check equal local dimensions
	std::size_t D = static_cast<std::size_t>(std::sqrt((double) rA.rows()));
	if (D * D != static_cast<std::size_t>(rA.rows()))
		throw Exception("gconcurrence", Exception::Type::DIMS_NOT_EQUAL);

// we compute exp(logdet()) to avoid underflow
	return D * std::abs(std::exp(2. / D * logdet(reshape(rA, D, D))));
}

/**
 * \brief Negativity of the bi-partite mixed state \a A
 *
 * \param A Eigen expression
 * \param dims Subsystems' dimensions
 * \return Negativity
 */
template<typename Derived>
double negativity(const Eigen::MatrixBase<Derived>& A,
		const std::vector<std::size_t>& dims)
{
	const DynMat<typename Derived::Scalar> & rA = A;
// check zero-size
	if (!internal::_check_nonzero_size(rA))
		throw Exception("negativity", Exception::Type::ZERO_SIZE);
// check bi-partite
	if (dims.size() != 2)
		throw Exception("negativity", Exception::Type::NOT_BIPARTITE);
// check square matrix vector
	if (!internal::_check_square_mat(rA))
		throw Exception("negativity", Exception::Type::MATRIX_NOT_SQUARE);
// check matching dimensions
	if (!internal::_check_dims_match_mat(dims, rA))
		throw Exception("negativity", Exception::Type::DIMS_MISMATCH_MATRIX);

	return (schatten(ptranspose(rA, { 0 }, dims), 1) - 1.) / 2.;
}

/**
 * \brief Logarithmic negativity of the bi-partite mixed state \a A
 *
 * \param A Eigen expression
 * \param dims Subsystems' dimensions
 * \return Logarithmic negativity, with the logarithm in base 2
 */
template<typename Derived>
double lognegativity(const Eigen::MatrixBase<Derived>& A,
		const std::vector<std::size_t>& dims)
{
	const DynMat<typename Derived::Scalar> & rA = A;
// check zero-size
	if (!internal::_check_nonzero_size(rA))
		throw Exception("lognegativity", Exception::Type::ZERO_SIZE);
// check bi-partite
	if (dims.size() != 2)
		throw Exception("lognegativity", Exception::Type::NOT_BIPARTITE);
// check square matrix vector
	if (!internal::_check_square_mat(rA))
		throw Exception("lognegativity", Exception::Type::MATRIX_NOT_SQUARE);
// check matching dimensions
	if (!internal::_check_dims_match_mat(dims, rA))
		throw Exception("lognegativity", Exception::Type::DIMS_MISMATCH_MATRIX);

	return std::log2(2 * negativity(rA, dims) + 1);
}

/**
 * \brief Wootters concurrence of the bi-partite qubit mixed state \a A
 *
 * \param A Eigen expression
 * \return Wootters concurrence
 */
template<typename Derived>
double concurrence(const Eigen::MatrixBase<Derived>& A)
{
	const DynMat<typename Derived::Scalar> & rA = A;
// check zero-size
	if (!internal::_check_nonzero_size(rA))
		throw Exception("concurrence", Exception::Type::ZERO_SIZE);
// check square matrix vector
	if (!internal::_check_square_mat(rA))
		throw Exception("concurrence", Exception::Type::MATRIX_NOT_SQUARE);
// check that the state is a 2-qubit state
	if (rA.rows() != 4)
		throw Exception("concurrence", Exception::Type::NOT_QUBIT_SUBSYS);

	cmat sigmaY = Gates::get_instance().Y;
	DynColVect<double> lambdas =
			evals(
					rA * kron(sigmaY, sigmaY) * conjugate(rA)
							* kron(sigmaY, sigmaY)).real();

	std::vector<double> lambdas_sorted(lambdas.data(),
			lambdas.data() + lambdas.size());

	std::sort(std::begin(lambdas_sorted), std::end(lambdas_sorted),
			std::greater<double>());
	std::transform(std::begin(lambdas_sorted), std::end(lambdas_sorted),
			std::begin(lambdas_sorted), [](double elem)
			{	return std::sqrt(std::abs(elem));}); // chop tiny negatives

	return std::max(0.,
			lambdas_sorted[0] - lambdas_sorted[1] - lambdas_sorted[2]
					- lambdas_sorted[3]);
}

} /* namespace qpp */

#endif /* INCLUDE_ENTANGLEMENT_H_ */
