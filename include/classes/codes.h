/*
 * Quantum++
 *
 * Copyright (c) 2013 - 2014 Vlad Gheorghiu (vgheorgh@gmail.com)
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

#ifndef INCLUDE_CLASSES_CODES_H_
#define INCLUDE_CLASSES_CODES_H_

namespace qpp
{
/**
 * \class qpp::Codes
 * \brief const Singleton class that defines quantum error correcting codes
 */
class Codes: public internal::Singleton<const Codes> // const Singleton
{
	friend class internal::Singleton<const Codes>;
public:
	/**
	 * \brief Code types, add more codes here if needed
	 *
	 * \see \a qpp::Codes::codeword()
	 */
	enum class Type // exception types
	{
		FIVE_QUBIT = 1, 	//!< [[5,1,3]] qubit code
		SEVEN_QUBIT_STEANE, //!< [[7,1,3]] Steane qubit code
		NINE_QUBIT_SHOR //!< [[9,1,3]] Shor qubit code
	};
private:
	/**
	 * \brief Default constructor
	 */
	Codes(){};
public:
	/**
	 * \brief Returns the codeword of the specified code
	 *
	 * \param type Code type, defined in the enum \a qpp::Codes::Types
	 * \param i Codeword index
	 * \return \a i-th codeword  of the code \a type
	 */
	ket codeword(Type type, std::size_t i) const
	{
		ket result;
		switch (type)
		{
		case Type::FIVE_QUBIT:  // [[5,1,3]] Five qubit code according to Nielsen and Chuang)
			switch (i)
			{
			case 0:
				result = (  mket({0,0,0,0,0}) + mket({1,0,0,1,0}) + mket({0,1,0,0,1}) + mket({1,0,1,0,0})
						  + mket({0,1,0,1,0}) - mket({1,1,0,1,1}) - mket({0,0,1,1,0}) - mket({1,1,0,0,0})
						  - mket({1,1,1,0,1}) - mket({0,0,0,1,1}) - mket({1,1,1,1,0}) - mket({0,1,1,1,1})
						  - mket({1,0,0,0,1}) - mket({0,1,1,0,0}) - mket({1,0,1,1,1}) + mket({0,0,1,0,1}))/4.;
				break;
			case 1:
				result = (  mket({1,1,1,1,1}) + mket({0,1,1,0,1}) + mket({1,0,1,1,0}) + mket({0,1,0,1,1})
				          + mket({1,0,1,0,1}) - mket({0,0,1,0,0}) - mket({1,1,0,0,1}) - mket({0,0,1,1,1})
				          - mket({0,0,0,1,0}) - mket({1,1,1,0,0}) - mket({0,0,0,0,1}) - mket({1,0,0,0,0})
				          - mket({0,1,1,1,0}) - mket({1,0,0,1,1}) - mket({0,1,0,0,0}) + mket({1,1,0,1,0}))/4.;
				break;
			default:
				throw Exception("codes", Exception::Type::NO_CODEWORD);
			}
			break;
		case Type::SEVEN_QUBIT_STEANE: // [[7,1,3]] Steane code according to Nielsen and Chuang)
			switch (i)
			{
			case 0:
				result =  (  mket({0,0,0,0,0,0,0}) + mket({1,0,1,0,1,0,1}) + mket({0,1,1,0,0,1,1}) + mket({1,1,0,0,1,1,0})
				           + mket({0,0,0,1,1,1,1}) + mket({1,0,1,1,0,1,0}) + mket({0,1,1,1,1,0,0}) + mket({1,1,0,1,0,0,1}))/std::sqrt(8.);

				break;
			case 1:
				result =  (  mket({1,1,1,1,1,1,1}) + mket({0,1,0,1,0,1,0}) + mket({1,0,0,1,1,0,0}) + mket({0,0,1,1,0,0,1})
				           + mket({1,1,1,0,0,0,0}) + mket({0,1,0,0,1,0,1}) + mket({1,0,0,0,0,1,1}) + mket({0,0,1,0,1,1,0}))/std::sqrt(8.);
				break;
			default:
				throw Exception("codes", Exception::Type::NO_CODEWORD);
			}
			break;
		case Type::NINE_QUBIT_SHOR: // [[9,1,3]] Shor code
			ket shora, shorb;
			shora = mket({0,0,0}) + mket({1,1,1,});
			shorb = mket({0,0,0}) - mket({1,1,1,});
			switch (i)
			{
			case 0:
				result = kron(shora, kron(shora, shora))/std::sqrt(8.);
				break;
			case 1:
				result = kron(shorb, kron(shorb, shorb))/std::sqrt(8.);
				break;
			default:
				throw Exception("codes", Exception::Type::NO_CODEWORD);
			}
		}
		return result;
	}
};
/* class Codes */

}

#endif /* INCLUDE_CLASSES_CODES_H_ */
