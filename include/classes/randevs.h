/*
 * randevs.h
 *
 *  Created on: Apr 7, 2014
 *      Author: vlad
 */

#ifndef INCLUDE_CLASSES_RANDEVS_H_
#define INCLUDE_CLASSES_RANDEVS_H_

namespace qpp
{

/**
 * \class qpp::RandomDevices
 * \brief Singeleton class that manages the source of randomness in the library
 *
 * It consists of a wrapper around an std::mt19937 Mersenne twister
 * random number generator engine and an std::random_device engine. The latter
 * is used to seed the Mersenne twister. The class also seeds the standard
 * std::srand C number generator, as it is used by Eigen.
 */
class RandomDevices: public internal::Singleton<RandomDevices> // Singleton
{
	friend class internal::Singleton<RandomDevices>;
	std::random_device _rd; ///< used to seed std::mt19937 _rng
public:
	std::mt19937 _rng; ///< Mersenne twister random number generator engine
private:
	/**
	 * \brief Initializes and seeds the random number generators.
	 */
	RandomDevices() :
			_rd(), _rng(_rd())
	{
		std::srand(_rng());
	}
};
/* class RandomDevices */

} /* namespace qpp */

#endif /* INCLUDE_CLASSES_RANDEVS_H_ */
