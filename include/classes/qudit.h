/*
 * qudit.h
 *
 *  Created on: Apr 9, 2014
 *      Author: vlad
 */

#ifndef QUDIT_H_
#define QUDIT_H_

namespace qpp
{

class Qudit
{
	cmat _rho;
	std::size_t _D;
public:
	// by default we have a standard qubit in state |0>
	Qudit(const cmat& rho = States::get_instance().pz0) :
			_rho(cmat::Zero(2, 2)), _D(2) // qubit by default
	{
		if (!internal::_check_nonzero_size(rho))
			throw Exception("Qudit::Qudit", Exception::Type::ZERO_SIZE);
		if (!internal::_check_square_mat(rho))
			throw Exception("Qudit::Qudit", Exception::Type::MATRIX_NOT_SQUARE);
		_D = rho.rows();

		_rho = rho;
	}

	std::size_t measure(const cmat& U, bool destructive = false)
	{
		if (!internal::_check_square_mat(U))
			throw Exception("Qudit::measure",
					Exception::Type::MATRIX_NOT_SQUARE);
		if (static_cast<std::size_t>(U.rows()) != _D)
			throw Exception("Qudit::measure", Exception::Type::DIMS_INVALID);

		std::vector<double> p(_D);
		for (std::size_t i = 0; i < _D; i++)
			p[i] = std::abs(
					(cplx) trace(
							prj((cmat) evects(U).col(i)) * _rho));

		DiscreteDistribution dd(p);
		std::size_t result = dd.sample();
		if (destructive) // von Neumann
			_rho = prj((cmat) evects(U).col(result)) * _rho
					* prj((cmat) evects(U).col(result)) / p[result];

		return result;
	}

	// measure in the standard basis {|j>}_{j=0}^{_D-1}
	std::size_t measure(bool destructive = false)
	{
		return measure(cmat::Identity(_D, _D), destructive);
	}

	cmat getRho() const
	{
		return _rho;
	}
	std::size_t getD() const
	{
		return _D;
	}
};

} /* namespace qpp */

#endif /* QUDIT_H_ */
