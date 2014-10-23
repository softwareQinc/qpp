/*
 * randevs.h
 *
 *  Created on: Apr 7, 2014
 *      Author: vlad
 */

#ifndef RANDEVS_H_
#define RANDEVS_H_

// Random devices singleton class
// public members:
// std::random_device _rd
// std::mt19937 _rng

namespace qpp
{

class RandomDevices: public Singleton<RandomDevices> // Singleton
{
	friend class Singleton<RandomDevices> ;
	std::random_device _rd; // used to seed std::mt19937 _rng
public:
	std::mt19937 _rng; // Mersenne twister random number generator engine
private:
	RandomDevices() :
			_rd(), _rng(_rd())
	{
		// seeds the standard C number generator, used by Eigen
		std::srand(_rng());
	}
};
/* class RandomDevices */

} /* namespace qpp */

#endif /* RANDEVS_H_ */
