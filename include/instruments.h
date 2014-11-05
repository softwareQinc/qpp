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

#ifndef INCLUDE_INSTRUMENTS_H_
#define INCLUDE_INSTRUMENTS_H_

// measurements
namespace qpp
{

/**
 * \brief Measures the state \a A using the set of Kraus operators \a Ks
 *
 * \param A Eigen expression
 * \param Ks Set of Kraus operators
 * \return \return Pair of vector of probabilities and vector of
 * resulting normalized states
 */
template<typename Derived>
std::pair<std::vector<double>, std::vector<cmat>> measure(
		const Eigen::MatrixBase<Derived>& A, const std::vector<cmat>& Ks)
{
	const DynMat<typename Derived::Scalar> &rA = A;

	// EXCEPTION CHECKS
	// check zero-size
	if (!internal::_check_nonzero_size(rA))
		throw Exception("measure", Exception::Type::ZERO_SIZE);

	// check the Kraus operators
	if (!internal::_check_nonzero_size(Ks))
		throw Exception("measure", Exception::Type::ZERO_SIZE);
	if (!internal::_check_square_mat(Ks[0]))
		throw Exception("measure", Exception::Type::MATRIX_NOT_SQUARE);
	if (Ks[0].rows() != rA.rows())
		throw Exception("measure", Exception::Type::DIMS_MISMATCH_MATRIX);
	for (auto && it : Ks)
		if (it.rows() != Ks[0].rows() || it.cols() != Ks[0].rows())
			throw Exception("measure", Exception::Type::DIMS_NOT_EQUAL);
	// END EXCEPTION CHECKS

	// probabilities
	std::vector<double> prob(Ks.size());
	// resulting states
	std::vector<cmat> outstates;

	if (internal::_check_square_mat(rA)) // square matrix
	{
		for (std::size_t i = 0; i < Ks.size(); i++)
		{
			cmat tmp;
			tmp = Ks[i] * rA * adjoint(Ks[i]); // un-normalized
			prob[i] = std::abs(trace(tmp)); // probability
			if (prob[i] > eps)
				outstates.push_back(tmp / prob[i]); // normalized
		}
	}
	else if (internal::_check_col_vector(rA)) // column vector
	{
		for (std::size_t i = 0; i < Ks.size(); i++)
		{
			cmat tmp;
			tmp = Ks[i] * rA; // un-normalized
			// probability
			prob[i] = std::abs((adjoint(tmp) * tmp).value());
			if (prob[i] > eps)
				outstates.push_back(tmp / std::sqrt(prob[i])); // normalized
		}
	}
	else
		throw Exception("measure",
				Exception::Type::MATRIX_NOT_SQUARE_OR_CVECTOR);

	return std::make_pair(prob, outstates);
}

/**
 * \brief Measures the state \a A in the basis specified by the unitary
 * matrix \a U
 *
 * \param A Eigen expression
 * \param U Unitary matrix representing the measurement basis
 * \return Pair of vector of probabilities and vector of
 * resulting normalized states
 */
template<typename Derived>
std::pair<std::vector<double>, std::vector<cmat>> measure(
		const Eigen::MatrixBase<Derived>& A, const cmat& U)
{
	const DynMat<typename Derived::Scalar> &rA = A;

	// EXCEPTION CHECKS
	// check zero-size
	if (!internal::_check_nonzero_size(rA))
		throw Exception("measure", Exception::Type::ZERO_SIZE);

	// check the gate U
	if (!internal::_check_nonzero_size(U))
		throw Exception("measure", Exception::Type::ZERO_SIZE);
	if (!internal::_check_square_mat(U))
		throw Exception("measure", Exception::Type::MATRIX_NOT_SQUARE);
	if (U.rows() != rA.rows())
		throw Exception("measure", Exception::Type::DIMS_MISMATCH_MATRIX);
	// END EXCEPTION CHECKS

	// probabilities
	std::vector<double> prob(U.rows());
	// resulting states
	std::vector<cmat> outstates;

	if (internal::_check_square_mat(rA)) // square matrix
	{
		for (std::size_t i = 0; i < static_cast<std::size_t>(U.rows()); i++)
		{
			cmat tmp;
			// un-normalized
			tmp = evects(U).col(i) * rA * adjoint((ket) evects(U).col(i));
			prob[i] = std::abs(trace(tmp)); // probability
			if (prob[i] > eps)
				outstates.push_back(tmp / prob[i]); // normalized
		}
	}
	else if (internal::_check_col_vector(rA)) // column vector
	{
		for (std::size_t i = 0; i < static_cast<std::size_t>(U.rows()); i++)
		{
			cmat tmp;
			tmp = prj((ket) evects(U).col(i)) * rA;
			// probability
			prob[i] = std::abs((adjoint(tmp) * tmp).value());
			if (prob[i] > eps)
				outstates.push_back(tmp / std::sqrt(prob[i])); // normalized
		}
	}
	else
		throw Exception("measure",
				Exception::Type::MATRIX_NOT_SQUARE_OR_CVECTOR);

	return std::make_pair(prob, outstates);
}

}
/* namespace qpp */

#endif /* INCLUDE_INSTRUMENTS_H_ */
