/*
 * stat.h
 *
 *  Created on: Dec 17, 2013
 *      Author: vlad
 */

#ifndef STAT_H_
#define STAT_H_

#include <random>

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
public:
	std::normal_distribution<> _d;
	NormalDistribution(double mean = 0, double sigma = 1) :
			_d(std::normal_distribution<>(mean, sigma))
	{
	}
	;
	double sample()
	{
		return _d(_rng);
	}
};

class UniformRealDistribution
{
public:
	std::uniform_real_distribution<> _d;
	UniformRealDistribution(double a = 0, double b = 1) :
			_d(std::uniform_real_distribution<>(a, b))
	{
	}
	;
	double sample()
	{
		return _d(_rng);
	}
};

class DiscreteDistribution
{
public:
	std::discrete_distribution<size_t> _d;
	DiscreteDistribution(std::initializer_list<double> weights) :
			_d(weights)
	{
	}
	;
	DiscreteDistribution(std::vector<double> weights) :
			_d(weights.begin(), weights.end())
	{
	}
	;
	template<typename InputIterator>
	DiscreteDistribution(InputIterator first, InputIterator last) :
			_d(first, last)
	{
	}
	;
	size_t sample()
	{
		return _d(_rng);
	}
};

/*
 class DiscreteDistributionFromComplex
 {
 std::discrete_distribution<size_t> d;
 public:
 DiscreteDistribution(std::initializer_list<types::cplx> amplitudes)
 {

 };
 size_t sample()
 {
 return d(_rng);
 }
 };*/

}
}

#endif /* STAT_H_ */
