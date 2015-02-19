/*
 * Quantum++
 *
 * Copyright (c) 2013 - 2015 Vlad Gheorghiu (vgheorgh@gmail.com)
 *
 * This file is part of Quantum++.
 *
 * Quantum++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Quantum++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Quantum++.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
* \file classes/random_devices.h
* \brief Random devices
*/

#ifndef CLASSES_RANDOM_DEVICES_H_
#define CLASSES_RANDOM_DEVICES_H_

namespace qpp
{

/**
* \class qpp::RandomDevices
* \brief Singeleton class that manages the source of randomness in the library
*
* Consists of a wrapper around an std::mt19937 Mersenne twister
* random number generator engine and an std::random_device engine. The latter
* is used to seed the Mersenne twister. This class also seeds the standard
* C number generator, as it is used by Eigen.
*/
class RandomDevices : public internal::Singleton<RandomDevices> // Singleton
{
    friend class internal::Singleton<RandomDevices>;

    std::random_device _rd; ///< used to seed std::mt19937 _rng
public:
    std::mt19937 _rng;      ///< Mersenne twister random number generator
private:
    /**
    * \brief Initializes and seeds the random number generators
    */
    RandomDevices() :
            _rd{}, _rng{_rd()}
    {
        std::srand(_rng());
    }
}; /* class RandomDevices */

} /* namespace qpp */

#endif /* CLASSES_RANDOM_DEVICES_H_ */
