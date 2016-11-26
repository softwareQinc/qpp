/*
 * Quantum++
 *
 * Copyright (c) 2013 - 2017 Vlad Gheorghiu (vgheorgh@gmail.com)
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
* is used to seed the Mersenne twister.
*
* \warning This class DOES NOT seed the standard C number generator used by
* Eigen::Matrix::Random(), since it is not thread safe. Do not use
* Eigen::Matrix::Random() or functions that depend on the C style random
* number engine, but use qpp::rand() instead!
*/
class RandomDevices final : public internal::Singleton<RandomDevices> //
// Singleton
{
    friend class internal::Singleton<RandomDevices>;

    std::random_device rd_; ///< used to seed std::mt19937 rng_
    std::mt19937 prng_;     ///< Mersenne twister random number generator
public:
    /**
    * \brief Returns a reference to the internal PRNG object
    * \return Reference to the internal PRNG object
    */
    std::mt19937& get_prng()
    {
        return prng_;
    }

    /**
    * \brief Loads the state of the PRNG from an input stream
    * \param is Input stream
    * \return The input stream
    */
    std::istream& load(std::istream& is)
    {
        return is >> prng_;
    }

    /**
    * \brief Saves the state of the PRNG to an output stream
    * \param os Output stream
    * \return The output stream
    */
    std::ostream& save(std::ostream& os) const
    {
        return os << prng_;
    }

private:
    /**
    * \brief Initializes and seeds the random number generators
    */
    RandomDevices() : rd_{}, prng_{rd_()}
    {
    }

    /**
    * \brief Default destructor
    */
    ~RandomDevices() = default;
}; /* class RandomDevices */

} /* namespace qpp */

#endif /* CLASSES_RANDOM_DEVICES_H_ */
