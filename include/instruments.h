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
 * \brief
 *
 * \param A
 * \param Ks
 * \return
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
	std::vector<double> prob(rA.size());
	// resulting states
	std::vector<cmat> rhos(rA.size());

	if (internal::_check_square_mat(rA)) // square matrix
	{

	}
	else if (internal::_check_col_vector(rA)) // column vector
	{

	}
	else
		throw Exception("measure",
				Exception::Type::MATRIX_NOT_SQUARE_OR_CVECTOR);

	//return rA;
}

/**
 * \brief
 *
 * \param A
 * \param U
 * \return
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
	std::vector<double> prob(rA.size());
	// resulting states
	std::vector<cmat> rhos(rA.size());

	if (internal::_check_square_mat(rA)) // square matrix
	{

	}
	else if (internal::_check_col_vector(rA)) // column vector
	{

	}
	else
		throw Exception("measure",
				Exception::Type::MATRIX_NOT_SQUARE_OR_CVECTOR);

	//return rA;
}

}
/* namespace qpp */

#endif /* INCLUDE_INSTRUMENTS_H_ */
