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
* \file classes/codes.h
* \brief Quantum error correcting codes
*/

#ifndef CLASSES_CODES_H_
#define CLASSES_CODES_H_

namespace qpp {
/**
* \class qpp::Codes
* \brief const Singleton class that defines quantum error correcting codes
*/
class Codes final : public internal::Singleton<const Codes> // const Singleton
{
    friend class internal::Singleton<const Codes>;

  public:
    /**
    * \brief Code types, add more codes here if needed
    * \see qpp::Codes::codeword()
    */
    enum class Type       // exception types
    { FIVE_QUBIT = 1,     ///< [[5,1,3]] qubit code
      SEVEN_QUBIT_STEANE, ///< [[7,1,3]] Steane qubit code
      NINE_QUBIT_SHOR     ///< [[9,1,3]] Shor qubit code
    };

  private:
    /**
    * \brief Default constructor
    */
    Codes() {} // = default; // clang++ spits the error below if defaulted:
    //    error:
    //    default initialization of an object of const type 'const qpp::Codes'
    //    requires a user-provided default constructor

    /**
    * \brief Default destructor
    */
    ~Codes() = default;

  public:
    /**
    * \brief Returns the codeword of the specified code type
    * \see qpp::Codes::Type
    *
    * \param type Code type
    * \param i Codeword index
    * \return \a i-th codeword  of the code \a type
    */
    ket codeword(Type type, idx i) const {
        ket result;
        switch (type) {
        // [[5,1,3]] Five qubit code (as in Nielsen and Chuang)
        case Type::FIVE_QUBIT:
            switch (i) {
            case 0:
                result = (mket({0, 0, 0, 0, 0}) + mket({1, 0, 0, 1, 0}) +
                          mket({0, 1, 0, 0, 1}) + mket({1, 0, 1, 0, 0}) +
                          mket({0, 1, 0, 1, 0}) - mket({1, 1, 0, 1, 1}) -
                          mket({0, 0, 1, 1, 0}) - mket({1, 1, 0, 0, 0}) -
                          mket({1, 1, 1, 0, 1}) - mket({0, 0, 0, 1, 1}) -
                          mket({1, 1, 1, 1, 0}) - mket({0, 1, 1, 1, 1}) -
                          mket({1, 0, 0, 0, 1}) - mket({0, 1, 1, 0, 0}) -
                          mket({1, 0, 1, 1, 1}) + mket({0, 0, 1, 0, 1})) /
                         4.;
                break;
            case 1:
                result = (mket({1, 1, 1, 1, 1}) + mket({0, 1, 1, 0, 1}) +
                          mket({1, 0, 1, 1, 0}) + mket({0, 1, 0, 1, 1}) +
                          mket({1, 0, 1, 0, 1}) - mket({0, 0, 1, 0, 0}) -
                          mket({1, 1, 0, 0, 1}) - mket({0, 0, 1, 1, 1}) -
                          mket({0, 0, 0, 1, 0}) - mket({1, 1, 1, 0, 0}) -
                          mket({0, 0, 0, 0, 1}) - mket({1, 0, 0, 0, 0}) -
                          mket({0, 1, 1, 1, 0}) - mket({1, 0, 0, 1, 1}) -
                          mket({0, 1, 0, 0, 0}) + mket({1, 1, 0, 1, 0})) /
                         4.;
                break;
            default:
                throw exception::NoCodeword("qpp::Codes::codeword()");
            }
            break;
        // [[7,1,3]] Steane code according to Nielsen and Chuang)
        case Type::SEVEN_QUBIT_STEANE:
            switch (i) {
            case 0:
                result =
                    (mket({0, 0, 0, 0, 0, 0, 0}) + mket({1, 0, 1, 0, 1, 0, 1}) +
                     mket({0, 1, 1, 0, 0, 1, 1}) + mket({1, 1, 0, 0, 1, 1, 0}) +
                     mket({0, 0, 0, 1, 1, 1, 1}) + mket({1, 0, 1, 1, 0, 1, 0}) +
                     mket({0, 1, 1, 1, 1, 0, 0}) +
                     mket({1, 1, 0, 1, 0, 0, 1})) /
                    std::sqrt(8.);

                break;
            case 1:
                result =
                    (mket({1, 1, 1, 1, 1, 1, 1}) + mket({0, 1, 0, 1, 0, 1, 0}) +
                     mket({1, 0, 0, 1, 1, 0, 0}) + mket({0, 0, 1, 1, 0, 0, 1}) +
                     mket({1, 1, 1, 0, 0, 0, 0}) + mket({0, 1, 0, 0, 1, 0, 1}) +
                     mket({1, 0, 0, 0, 0, 1, 1}) +
                     mket({0, 0, 1, 0, 1, 1, 0})) /
                    std::sqrt(8.);
                break;
            default:
                throw exception::NoCodeword("qpp::Codes::codeword()");
            }
            break;
        // [[9,1,3]] Shor code
        case Type::NINE_QUBIT_SHOR:
            ket shora, shorb;
            shora = mket({0, 0, 0}) + mket({
                                          1, 1, 1,
                                      });
            shorb = mket({0, 0, 0}) - mket({
                                          1, 1, 1,
                                      });
            switch (i) {
            case 0:
                result = kron(shora, kron(shora, shora)) / std::sqrt(8.);
                break;
            case 1:
                result = kron(shorb, kron(shorb, shorb)) / std::sqrt(8.);
                break;
            default:
                throw exception::NoCodeword("qpp::Codes::codeword()");
            }
        }

        return result;
    }
}; /* class Codes */

} /* namespace qpp */

#endif /* CLASSES_CODES_H_ */
