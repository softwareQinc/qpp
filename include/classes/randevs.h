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

// Singleton, cannot be const
class RandomDevices: public Singleton<RandomDevices>
{
	friend class Singleton<RandomDevices> ;
	~RandomDevices()
	{
	}
	;
public:
	std::random_device _rd;
	std::mt19937 _rng;
private:
	RandomDevices() :
			_rd(), _rng(_rd())
	{
		// seed the standard C number generator, used by Eigen
		std::srand(_rd());
	}
};

} /* namespace qpp */

#endif /* RANDEVS_H_ */
