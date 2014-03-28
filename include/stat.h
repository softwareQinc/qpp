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

class NormalDistribution
{
	std::normal_distribution<> d;
public:
	NormalDistribution(double mean = 0, double sigma = 1) :
			d(std::normal_distribution<>(mean, sigma))
	{
	}
	;
	double sample()
	{
		return d(_rng);
	}
};

class UniformRealDistribution
{
	std::uniform_real_distribution<> d;
public:
	UniformRealDistribution(double a = 0, double b = 1) :
			d(std::uniform_real_distribution<>(a, b))
	{
	}
	;
	double sample()
	{
		return d(_rng);
	}
};

}
}

#endif /* STAT_H_ */
