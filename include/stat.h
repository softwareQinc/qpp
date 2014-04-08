/*
 * stat.h
 *
 *  Created on: Dec 17, 2013
 *      Author: vlad
 */

#ifndef STAT_H_
#define STAT_H_

#include <algorithm>
#include <functional>
#include "types.h"
#include "internal.h"
#include "exception.h"
#include "randevs.h"

// statistical distributions etc

namespace qpp
{
namespace stat
{

// light wrappers around C++11 statistical distributions

class NormalDistribution
{
protected:
	std::normal_distribution<> _d;

public:
	NormalDistribution(double mean = 0, double sigma = 1) :
			_d(std::normal_distribution<>(mean, sigma))
	{
	}

	double sample()
	{
		return _d(RandomDevices::getInstance()->_rng);
	}
};

class UniformRealDistribution
{
protected:
	std::uniform_real_distribution<> _d;

public:
	UniformRealDistribution(double a = 0, double b = 1) :
			_d(std::uniform_real_distribution<>(a, b))
	{
	}

	double sample()
	{
		return _d(RandomDevices::getInstance()->_rng);
	}
};

class DiscreteDistribution
{
protected:
	std::discrete_distribution<size_t> _d;

public:
	template<typename InputIterator>
	DiscreteDistribution(InputIterator first, InputIterator last) :
			_d(first, last)
	{
	}

	DiscreteDistribution(std::initializer_list<double> weights) :
			_d(weights)
	{
	}

	DiscreteDistribution(std::vector<double> weights) :
			_d(std::begin(weights), std::end(weights))
	{
	}

	size_t sample()
	{
		return _d(RandomDevices::getInstance()->_rng);
	}

	std::vector<double> probabilities()
	{
		return _d.probabilities();
	}
};

class DiscreteDistributionFromComplex
{
protected:
	std::discrete_distribution<size_t> _d;

	template<typename InputIterator>
	std::vector<double> cplx2amplitudes(InputIterator first, InputIterator last)
	{
		std::vector<double> weights(last - first);
		std::transform(first, last, std::begin(weights),
				[](const types::cplx &z)->double
				{	return std::pow(std::abs(z),2);});
		return weights;
	}

public:
	template<typename InputIterator>
	DiscreteDistributionFromComplex(InputIterator first, InputIterator last) :
			_d()
	{
		std::vector<double> weights = cplx2amplitudes(first, last);
		std::discrete_distribution<size_t> tmp(std::begin(weights),
				std::end(weights));
		_d = tmp;
	}

	DiscreteDistributionFromComplex(
			std::initializer_list<types::cplx> amplitudes) :
			_d()
	{
		std::vector<double> weights = cplx2amplitudes(std::begin(amplitudes),
				std::end(amplitudes));
		std::discrete_distribution<size_t> tmp(std::begin(weights),
				std::end(weights));
		_d = tmp;
	}

	DiscreteDistributionFromComplex(std::vector<types::cplx> amplitudes) :
			_d()
	{
		std::vector<double> weights = cplx2amplitudes(std::begin(amplitudes),
				std::end(amplitudes));
		std::discrete_distribution<size_t> tmp(std::begin(weights),
				std::end(weights));
		_d = tmp;
	}

	DiscreteDistributionFromComplex(const types::cmat& V) :
			_d()
	{
		// check zero-size
		if (!internal::_check_nonzero_size(V))
			throw Exception("DiscreteDistributionFromComplex::"
					"DiscreteDistributionFromComplex",
					Exception::Type::ZERO_SIZE);

		// check vector
		if (!internal::_check_vector(V))
			throw Exception("DiscreteDistributionFromComplex::"
					"DiscreteDistributionFromComplex",
					Exception::Type::MATRIX_NOT_VECTOR);

		std::vector<double> weights = cplx2amplitudes(V.data(),
				V.data() + V.size());
		std::discrete_distribution<size_t> tmp(std::begin(weights),
				std::end(weights));
		_d = tmp;

	}

	size_t sample()
	{
		return _d(RandomDevices::getInstance()->_rng);
	}

	std::vector<double> probabilities()
	{
		return _d.probabilities();
	}
};

} /* namespace stat */
} /* namespace qpp */

#endif /* STAT_H_ */
