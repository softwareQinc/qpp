/*
 * qudit.h
 *
 *  Created on: Apr 9, 2014
 *      Author: vlad
 */

#ifndef QUDIT_H_
#define QUDIT_H_

#include "exception.h"
#include "functions.h"
#include "gates.h"
#include "internal.h"
#include "types.h"
#include "classes/stat.h"

namespace qpp
{

class Qudit
{
	types::cmat _rho;
	size_t _D;
public:
	// by default we have a standard qubit in state |0>
	Qudit(const types::cmat& rho = Gates::getInstance().pz0) :
			_rho(types::cmat::Zero(2, 2)), _D(2) // qubit by default
	{
		if (!internal::_check_nonzero_size(rho))
			throw Exception("Qudit::Qudit", Exception::Type::ZERO_SIZE);
		if (!internal::_check_square_mat(rho))
			throw Exception("Qudit::Qudit", Exception::Type::MATRIX_NOT_SQUARE);
		_D = rho.rows();

		_rho = rho;
	}

	size_t measure(const types::cmat& U, bool destructive = false)
	{
		if (!internal::_check_square_mat(U))
			throw Exception("Qudit::measure",
					Exception::Type::MATRIX_NOT_SQUARE);
		if (static_cast<size_t>(U.rows()) != _D)
			throw Exception("Qudit::measure", Exception::Type::DIMS_INVALID);

		std::vector<double> p(_D);
		for (size_t i = 0; i < _D; i++)
			p[i] = std::abs(trace(proj(evects(U).col(i)) * _rho));

		stat::DiscreteDistribution dd(p);
		size_t result = dd.sample();
		if (destructive) // von Neumann
			_rho = proj(evects(U).col(result)) * _rho
					* proj(evects(U).col(result)) / p[result];

		return result;
	}

	// measure in the standard basis {|j>}_{j=0}^{_D-1}
	size_t measure(bool destructive = false)
	{
		return measure(types::cmat::Identity(_D, _D), destructive);
	}

	types::cmat getRho() const
	{
		return _rho;
	}
	size_t getD() const
	{
		return _D;
	}

	virtual ~Qudit() = default;
};

} /* namespace qpp */

#endif /* QUDIT_H_ */
