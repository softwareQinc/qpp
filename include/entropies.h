/*
 * entropy.h
 *
 *  Created on: Mar 27, 2014
 *      Author: vlad
 */

#ifndef ENTROPY_H_
#define ENTROPY_H_

#include <cmath>

#include "functions.h"
#include "internal.h"
#include "types.h"
#include "classes/exception.h"
#include "io.h"

// various entropies, assume as input either
// a normalized hermitian matrix or a probability vector

namespace qpp
{

// Shannon/von-Neumann entropy
template<typename Derived>
double shannon(const Eigen::MatrixBase<Derived>& A)
{
	const types::DynMat<typename Derived::Scalar> & rA = A;

	// check zero-size
	if (!internal::_check_nonzero_size(rA))
		throw Exception("shannon", Exception::Type::ZERO_SIZE);

	// input is a vector
	if (internal::_check_vector(rA))
	{
		double result = 0;
		// take the absolut value to get rid of tiny negatives
		for (size_t i = 0; i < static_cast<size_t>(rA.size()); i++)
			if (std::abs(rA(i)) != 0) // not identically zero
				result -= std::abs(rA(i)) * std::log2(std::abs(rA(i)));

		return result;
	}

	// input is a matrix

	// check square matrix
	if (!internal::_check_square_mat(rA))
		throw Exception("shannon", Exception::Type::MATRIX_NOT_SQUARE);

	// get the eigenvalues
	types::dmat ev = hevals(rA);
	double result = 0;
	// take the absolut value to get rid of tiny negatives
	for (size_t i = 0; i < static_cast<size_t>(ev.rows()); i++)
		if (std::abs((types::cplx) ev(i)) != 0) // not identically zero
			result -= std::abs((types::cplx) ev(i))
					* std::log2(std::abs((types::cplx) ev(i)));

	return result;
}

// Renyi-alpha entropy (alpha>=0)
template<typename Derived>
double renyi(const double alpha, const Eigen::MatrixBase<Derived>& A)
{
	const types::DynMat<typename Derived::Scalar> & rA = A;

	if (alpha < 0)
		throw Exception("renyi", Exception::Type::OUT_OF_RANGE);

	if (alpha == 1) // Shannon/von Neumann
		return shannon(rA);

	// check zero-size
	if (!internal::_check_nonzero_size(rA))
		throw Exception("renyi", Exception::Type::ZERO_SIZE);

	// input is a vector
	if (internal::_check_vector(rA))
	{
		if (alpha == 0)
			return std::log2(rA.size());

		double result = 0;
		// take the absolut value to get rid of tiny negatives
		for (size_t i = 0; i < static_cast<size_t>(rA.size()); i++)
			if (std::abs((types::cplx) rA(i)) != 0) // not identically zero
				result += std::pow(std::abs(rA(i)), alpha);

		return std::log2(result) / (1 - alpha);
	}

	// input is a matrix

	// check square matrix
	if (!internal::_check_square_mat(rA))
		throw Exception("renyi", Exception::Type::MATRIX_NOT_SQUARE);

	if (alpha == 0)
		return std::log2(rA.rows());

	// get the eigenvalues
	types::dmat ev = hevals(rA);
	double result = 0;
	// take the absolut value to get rid of tiny negatives
	for (size_t i = 0; i < static_cast<size_t>(ev.rows()); i++)
		if (std::abs((types::cplx) ev(i)) != 0) // not identically zero
			result += std::pow(std::abs((types::cplx) ev(i)), alpha);

	return std::log2(result) / (1 - alpha);
}

// Renyi-infinity entropy (min entropy)
template<typename Derived>
double renyi_inf(const Eigen::MatrixBase<Derived>& A)
{
	const types::DynMat<typename Derived::Scalar> & rA = A;

	// check zero-size
	if (!internal::_check_nonzero_size(rA))
		throw Exception("renyi_inf", Exception::Type::ZERO_SIZE);

	// input is a vector
	if (internal::_check_vector(rA))
	{
		double max = 0;
		for (size_t i = 0; i < static_cast<size_t>(rA.size()); i++)
			if (std::abs(rA(i)) > max)
				max = std::abs(rA(i));

		return -std::log2(max);
	}

	// input is a matrix

	// check square matrix
	if (!internal::_check_square_mat(rA))
		throw Exception("renyi_inf", Exception::Type::MATRIX_NOT_SQUARE);

	// get the eigenvalues
	types::dmat ev = hevals(rA);
	double max = 0;
	// take the absolut value to get rid of tiny negatives
	for (size_t i = 0; i < static_cast<size_t>(ev.size()); i++)
		if (std::abs((types::cplx) ev(i)) > max)
			max = std::abs((types::cplx) ev(i));

	return -std::log2(max);
}

// Tsallis-alpha entropy (alpha >=0)
// when alpha->1 converges to Shannon/von Neumann with base e logarithm
template<typename Derived>
double tsallis(const double alpha, const Eigen::MatrixBase<Derived>& A)
{
	const types::DynMat<typename Derived::Scalar> & rA = A;

	if (alpha < 0)
		throw Exception("tsallis", Exception::Type::OUT_OF_RANGE);

	if (alpha == 1) // Shannon/von Neumann with base e logarithm
		return shannon(rA) * std::log(2);

	// check zero-size
	if (!internal::_check_nonzero_size(rA))
		throw Exception("tsallis", Exception::Type::ZERO_SIZE);

	// input is a vector
	if (internal::_check_vector(rA))
	{
		double result = 0;
		// take the absolut value to get rid of tiny negatives
		for (size_t i = 0; i < static_cast<size_t>(rA.size()); i++)
			if (std::abs((types::cplx) rA(i)) != 0) // not identically zero
				result += std::pow(std::abs(rA(i)), alpha);

		return (result - 1) / (1 - alpha);
	}

	// input is a matrix

	// check square matrix
	if (!internal::_check_square_mat(rA))
		throw Exception("tsallis", Exception::Type::MATRIX_NOT_SQUARE);

	// get the eigenvalues
	types::dmat ev = hevals(rA);
	double result = 0;
	// take the absolut values of the entries
	//of tiny negativesginary parts
	for (size_t i = 0; i < static_cast<size_t>(ev.rows()); i++)
		if (std::abs((types::cplx) ev(i)) != 0) // not identically zero
			result += std::pow(std::abs((types::cplx) ev(i)), alpha);

	return (result - 1) / (1 - alpha);
}

// quantum mutual information between 2 subsystems
template<typename Derived>
double qmutualinfo(const Eigen::MatrixBase<Derived>& A,
		const std::vector<size_t>& subsys, const std::vector<size_t>& dims)

{
	const types::DynMat<typename Derived::Scalar> & rA = A;

// error checks

	// check that there are only 2 subsystems
	if (subsys.size() != 2)
		throw Exception("mutualinfo", Exception::Type::NOT_BIPARTITE);

	// check zero-size
	if (!internal::_check_nonzero_size(rA))
		throw Exception("mutualinfo", Exception::Type::ZERO_SIZE);

	// check that dims is a valid dimension vector
	if (!internal::_check_dims(dims))
		throw Exception("mutualinfo", Exception::Type::DIMS_INVALID);

	// check square matrix
	if (!internal::_check_square_mat(rA))
		throw Exception("mutualinfo", Exception::Type::MATRIX_NOT_SQUARE);

// check that dims match the dimension of A
	if (!internal::_check_dims_match_mat(dims, rA))
		throw Exception("mutualinfo", Exception::Type::DIMS_MISMATCH_MATRIX);

// check that subsys are valid
	if (!internal::_check_subsys_match_dims(subsys, dims))
		throw Exception("mutualinfo", Exception::Type::SUBSYS_MISMATCH_DIMS);

// construct the complement of subsys
	std::vector<size_t> subsysbarA;
	std::vector<size_t> subsysbarB;
	std::vector<size_t> subsysbarAB;
	for (size_t i = 0; i < dims.size(); i++)
	{
		if (subsys[0] != i)
			subsysbarA.push_back(i);
		if (subsys[1] != i)
			subsysbarB.push_back(i);
		if (subsys[0] != i && subsys[1] != i)
			subsysbarAB.push_back(i);

	};

	types::cmat rhoA;
	types::cmat rhoB;
	types::cmat rhoAB;
	if (dims.size() == 2) // bipartite state
	{
		rhoA = ptrace2(rA, dims);
		rhoB = ptrace1(rA, dims);
		rhoAB = rA;
	}
	else
	{
		rhoA = ptrace(rA, subsysbarA, dims);
		rhoB = ptrace(rA, subsysbarB, dims);
		rhoAB = ptrace(rA, subsysbarAB, dims);
	}

	return shannon(rhoA) + shannon(rhoB) - shannon(rhoAB);
}

} /* namespace qpp */
#endif /* ENTROPY_H_ */
