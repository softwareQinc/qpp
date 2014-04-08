/*
 * randevs.h
 *
 *  Created on: Apr 7, 2014
 *      Author: vlad
 */

#ifndef RANDEVS_H_
#define RANDEVS_H_

#include <cstdlib>
#include <random>

// Random devices singleton class
// public members:
// std::random_device _rd
// std::mt19937 _rng

namespace qpp
{

class RandomDevices // make it a singleton
{
private:
	RandomDevices() :
			_rd(), _rng(_rd())
	{
		// seed the standard C number generator, used by Eigen
		std::srand(_rd());
	}
public:
	RandomDevices(const RandomDevices&) = delete;
	RandomDevices& operator=(const RandomDevices&) = delete;

	std::random_device _rd;
	std::mt19937 _rng;

	static RandomDevices& getInstance()
	{
		static RandomDevices instance; // Guaranteed to be destroyed.
		// Instantiated on first use.
		return instance;
	}

	virtual ~RandomDevices() = default;
};

} /* namespace qpp */

#endif /* RANDEVS_H_ */
