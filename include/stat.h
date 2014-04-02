/*
 * stat.h
 *
 *  Created on: Dec 17, 2013
 *      Author: vlad
 */

#ifndef STAT_H_
#define STAT_H_

#include <random>
#include <algorithm>
#include <functional>
#include <stdexcept>
#include "types.h"
#include "internal.h"

// statistical distributions etc

namespace qpp
{
namespace stat
{

extern std::random_device _rd; // use for seeding
extern std::mt19937 _rng; // our random number generator

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
		return _d(_rng);
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
		return _d(_rng);
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
			_d(weights.begin(), weights.end())
	{
	}

	size_t sample()
	{
		return _d(_rng);
	}

	std::vector<double> probabilities()
	{
		return _d.probabilities();
	}
};

//TODO: Add constructor for cmat
class DiscreteDistributionFromComplex
{
protected:
	std::discrete_distribution<size_t> _d;

	template<typename InputIterator>
	std::vector<double> cplx2double(InputIterator first, InputIterator last)
	{
		std::vector<double> weights(last - first);
		std::transform(first, last, weights.begin(),
				[](const types::cplx &z)->double
				{	return std::abs(z);});
		return weights;
	}

public:
	template<typename InputIterator>
	DiscreteDistributionFromComplex(InputIterator first, InputIterator last) :
			_d()
	{
		std::vector<double> weights = cplx2double(first, last);
		std::discrete_distribution<size_t> tmp(weights.begin(), weights.end());
		_d = tmp;
	}

	DiscreteDistributionFromComplex(
			std::initializer_list<types::cplx> amplitudes) :
			_d()
	{
		std::vector<double> weights = cplx2double(amplitudes.begin(),
				amplitudes.end());
		std::discrete_distribution<size_t> tmp(weights.begin(), weights.end());
		_d = tmp;
	}

	DiscreteDistributionFromComplex(std::vector<types::cplx> amplitudes) :
			_d()
	{
		std::vector<double> weights = cplx2double(amplitudes.begin(),
				amplitudes.end());
		std::discrete_distribution<size_t> tmp(weights.begin(), weights.end());
		_d = tmp;
	}

	DiscreteDistributionFromComplex(const types::cmat &v) :
			_d()
	{
		if (!internal::_check_nonzero_size(v))
			throw std::runtime_error(
					"DiscreteDistributionFromComplex::DiscreteDistributionFromComplex: "
							"zero-sized matrix!");
		if (!internal::_check_vector(v))
			throw std::runtime_error(
					"DiscreteDistributionFromComplex::DiscreteDistributionFromComplex: "
							"input must be a row/column vector!");
		std::vector<double> weights = cplx2double(v.data(),
				v.data() + v.size());
		std::discrete_distribution<size_t> tmp(weights.begin(), weights.end());
		_d = tmp;

	}

	size_t sample()
	{
		return _d(_rng);
	}

	std::vector<double> probabilities()
	{
		return _d.probabilities();
	}
};

}
}

#endif /* STAT_H_ */
