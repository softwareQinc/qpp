/*
 * qubit.h
 *
 *  Created on: Apr 9, 2014
 *      Author: vlad
 */

#ifndef QUBIT_H_
#define QUBIT_H_

#include "types.h"
#include "gates.h"
#include "internal.h"
#include "exception.h"
#include "functions.h"
#include "stat.h"

namespace qpp
{

class Qubit
{
	types::ket _psi;
public:
	Qubit(const types::ket& psi = Gates::getInstance().z0)
	{
		_psi = psi;
	}
	Qubit(std::initializer_list<types::cplx> ampl)
	{
		_psi.resize(2);
		size_t cnt = 0;
		for (auto i : ampl)
		{
			_psi(cnt) = i;
			cnt++;
		}
	}

	size_t measure(const types::cmat& U = types::cmat::Identity(2, 2),
			bool destructive = false)
	{
		if (!internal::_check_square_mat(U))
			throw Exception("Qubit::measure",
					Exception::Type::MATRIX_NOT_SQUARE);
		if (U.rows() != 2)
			throw Exception("Qubit::measure", Exception::Type::DIMS_INVALID);
		std::vector<double> p(2);
		for (size_t i = 0; i < 2; i++)
			p[i] = std::abs(trace(proj(evects(U).col(i)) * proj(_psi)));

		stat::DiscreteDistribution dd(p);
		size_t result = dd.sample();
		if (destructive) // von Neumann
			_psi = proj(evects(U).col(result)) * _psi
					 / std::sqrt(p[result]);

		return result;

	}
	virtual ~Qubit() = default;
};

} /* namespace qpp */

#endif /* QUBIT_H_ */
