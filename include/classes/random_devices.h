/*
 * This file is part of Quantum++.
 *
 * MIT License
 *
 * Copyright (c) 2013 - 2018 Vlad Gheorghiu (vgheorgh@gmail.com)
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

/**
* \file classes/random_devices.h
* \brief Random devices
*/

#ifndef CLASSES_RANDOM_DEVICES_H_
#define CLASSES_RANDOM_DEVICES_H_

namespace qpp {
/**
* \class qpp::RandomDevices
* \brief Singleton class that manages the source of randomness in the library
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

    std::random_device rd_; ///< used to seed std::mt19937 prng_
    std::mt19937 prng_;     ///< Mersenne twister random number generator
  public:
    /**
    * \brief Returns a reference to the internal PRNG object
    * \return Reference to the internal PRNG object
    */
    std::mt19937& get_prng() { return prng_; }

    /**
    * \brief Loads the state of the PRNG from an input stream
    * \param is Input stream
    * \return The input stream
    */
    std::istream& load(std::istream& is) { return is >> prng_; }

    /**
    * \brief Saves the state of the PRNG to an output stream
    * \param os Output stream
    * \return The output stream
    */
    std::ostream& save(std::ostream& os) const { return os << prng_; }

  private:
    /**
    * \brief Initializes and seeds the random number generators
    */
    RandomDevices() : rd_{}, prng_{rd_()} {}

    /**
    * \brief Default destructor
    */
    ~RandomDevices() = default;
}; /* class RandomDevices */

} /* namespace qpp */

#endif /* CLASSES_RANDOM_DEVICES_H_ */
