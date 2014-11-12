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

#ifndef INCLUDE_CLASSES_EXCEPTION_H_
#define INCLUDE_CLASSES_EXCEPTION_H_

namespace qpp
{

    /**
    * \class qpp::Exception
    * \brief Generates custom exceptions, used when validating function parameters
    *
    * Customize this class if more exceptions are needed
    */
    class Exception : public std::exception
    {
    public:
        /**
        * \brief Exception types, add more exceptions here if needed
        *
        * \see qpp:Exception::_construct_exception_msg()
        */
        enum class Type // exception types
        {
            UNKNOWN_EXCEPTION = 1,       //!< UNKNOWN_EXCEPTION
            /*!< Unknown exception */
                    ZERO_SIZE,                   //!< ZERO_SIZE
            /*!< Zero sized object, e.g. empty Eigen::Matrix
             * or std::vector with no elements */
                    MATRIX_NOT_SQUARE,           //!< MATRIX_NOT_SQUARE
            /*!< Eigen::Matrix is not square */
                    MATRIX_NOT_CVECTOR,          //!< MATRIX_NOT_CVECTOR
            /*!< Eigen::Matrix is not a column vector */
                    MATRIX_NOT_RVECTOR,          //!< MATRIX_NOT_RVECTOR
            /*!< Eigen::Matrix is not a row vector */
                    MATRIX_NOT_VECTOR,           //!< MATRIX_NOT_VECTOR
            /*!< Eigen::Matrix is not a row/column vector */
                    MATRIX_NOT_SQUARE_OR_CVECTOR,          //!< MATRIX_NOT_SQUARE_OR_CVECTOR
            /*!< Eigen::Matrix is not square nor a column vector */
                    MATRIX_NOT_SQUARE_OR_RVECTOR,          //!< MATRIX_NOT_SQUARE_OR_RVECTOR
            /*!< Eigen::Matrix is not square nor a row vector */
                    MATRIX_NOT_SQUARE_OR_VECTOR, //!< MATRIX_NOT_SQUARE_OR_VECTOR
            /*!< Eigen::Matrix is not square nor a row/column vector */
                    MATRIX_MISMATCH_SUBSYS,        //!< MATRIX_MISMATCH_SUBSYS
            /*! Matrix size mismatch subsystem sizes (e.g. in
             * qpp::apply(), or qpp::channel() */
                    DIMS_INVALID,                //!< DIMS_INVALID
            /*!< std::vector<std::size_t> representing the dimensions
             * has zero size or contains zeros */
                    DIMS_NOT_EQUAL,              //!< DIMS_NOT_EQUAL
            /*!< std::vector<std::size_t> representing the dimensions
             * contains non-equal elements */
                    DIMS_MISMATCH_MATRIX,        //!< DIMS_MISMATCH_MATRIX
            /*!< Product of the dimenisons' std::vector<std::size_t>
             * is not equal to the number of rows of Eigen::Matrix
             * (assumed to be square) */
                    DIMS_MISMATCH_CVECTOR,       //!< DIMS_MISMATCH_CVECTOR
            /*!< Product of the dimenisons' std::vector<std::size_t>
             * is not equal to the number of cols of Eigen::Matrix
             * (assumed to be a column vector) */
                    DIMS_MISMATCH_RVECTOR,       //!< DIMS_MISMATCH_RVECTOR
            /*!< Product of the dimenisons' std::vector<std::size_t>
             * is not equal to the number of cols of Eigen::Matrix
             * (assumed to be a row vector) */
                    DIMS_MISMATCH_VECTOR,        //!< DIMS_MISMATCH_VECTOR
            /*!< Product of the dimenisons' std::vector<std::size_t>
             * is not equal to the number of cols of Eigen::Matrix
             * (assumed to be a row/column vector) */
                    SUBSYS_MISMATCH_DIMS,        //!< SUBSYS_MISMATCH_DIMS
            /*!< std::vector<std::size_t> representing the subsystem
             * labels has duplicatates, or has entries that are larger than
             * the size of the std::vector<std::size_t> representing the
             * dimensions */
                    NOT_QUBIT_GATE,              //!< NOT_QUBIT_GATE
            /*!<  Eigen::Matrix is not 2 x 2 */
                    NOT_QUBIT_SUBSYS,            //!< NOT_QUBIT_SUBSYS
            /*!< Subsystems are not 2-dimensional */
                    NOT_BIPARTITE,               //!< NOT_BIPARTITE
            /*!< std::vector<std::size_t> representing the dimensions
             * has size different from 2 */
                    NO_CODEWORD,                 //!< NO_CODEWORD
            /*!< Codeword does not exist, thrown when calling
             * qpp::Codes::codeword() with invalid \a i */
                    PERM_INVALID,                //!< PERM_INVALID
            /*!< Invalid std::vector<std::size_t> permutation */
                    OUT_OF_RANGE,                //!< OUT_OF_RANGE
            /*!< Parameter out of range */
                    TYPE_MISMATCH,               //!< TYPE_MISMATCH
            /*!< Types do not match (i.e. Matrix<double> vs Matrix<cplx>) */
                    UNDEFINED_TYPE,              //!< UNDEFINED_TYPE
            /*!< Templated function not defined for this type */
                    CUSTOM_EXCEPTION             //!< CUSTOM_EXCEPTION
            /*!< Custom exception, user must provide a custom message */
        };

        /**
        * \brief Constructs an exception
        *
        * \param where Text representing where the exception occured
        * \param type Exception's type, see the strong enumeration
        * qpp::Exception::Type
        */
        Exception(const std::string &where, const Type &type) :
                _where(where), _msg(), _type(type), _custom()
        {
            _construct_exception_msg();
        }

        /**
        * \brief Constructs an exception
        *
        * \overload
        *
        * \param where Text representing where the exception occured
        * \param custom Exception's description
        */
        Exception(const std::string &where, const std::string &custom) :
                _where(where), _msg(), _type(Type::CUSTOM_EXCEPTION), _custom(
                custom)
        {
            _construct_exception_msg();
            _msg += custom; // add the custom message at the end
        }

        /**
        * \brief Overrides std::exception::what()
        *
        * \return Exception's description
        */
        virtual const char *what() const noexcept override
        {
            return _msg.c_str();
        }

    private:
        std::string _where, _msg;
        Type _type;
        std::string _custom;

        /**
        * \brief Constructs the exception's description from its type
        *
        * Must modify the code of this function if more exceptions are added
        *
        * \return Exception's description
        */
        std::string _construct_exception_msg()
        {
            _msg += "IN ";
            _msg += _where;
            _msg += ": ";

            switch (_type)
            {
                case Type::UNKNOWN_EXCEPTION:
                    _msg += "UNKNOWN EXCEPTION!";
                    break;
                case Type::ZERO_SIZE:
                    _msg += "Object has zero size!";
                    break;
                case Type::MATRIX_NOT_SQUARE:
                    _msg += "Matrix is not square!";
                    break;
                case Type::MATRIX_NOT_CVECTOR:
                    _msg += "Matrix is not column vector!";
                    break;
                case Type::MATRIX_NOT_RVECTOR:
                    _msg += "Matrix is not row vector!";
                    break;
                case Type::MATRIX_NOT_VECTOR:
                    _msg += "Matrix is not vector!";
                    break;
                case Type::MATRIX_NOT_SQUARE_OR_CVECTOR:
                    _msg += "Matrix is not square nor column vector!";
                    break;
                case Type::MATRIX_NOT_SQUARE_OR_RVECTOR:
                    _msg += "Matrix is not square nor row vector!";
                    break;
                case Type::MATRIX_NOT_SQUARE_OR_VECTOR:
                    _msg += "Matrix is not square nor vector!";
                    break;
                case Type::MATRIX_MISMATCH_SUBSYS:
                    _msg += "Matrix mismatch subsystems!";
                    break;
                case Type::DIMS_INVALID:
                    _msg += "Invalid dimension(s)!";
                    break;
                case Type::DIMS_NOT_EQUAL:
                    _msg += "Dimensions not equal!";
                    break;
                case Type::DIMS_MISMATCH_MATRIX:
                    _msg += "Dimension(s) mismatch matrix size!";
                    break;
                case Type::DIMS_MISMATCH_CVECTOR:
                    _msg += "Dimension(s) mismatch column vector!";
                    break;
                case Type::DIMS_MISMATCH_RVECTOR:
                    _msg += "Dimension(s) mismatch row vector!";
                    break;
                case Type::DIMS_MISMATCH_VECTOR:
                    _msg += "Dimension(s) mismatch vector!";
                    break;
                case Type::SUBSYS_MISMATCH_DIMS:
                    _msg += "Subsystems mismatch dimensions!";
                    break;
                case Type::NOT_QUBIT_GATE:
                    _msg += "Matrix is not qubit gate!";
                    break;
                case Type::NOT_QUBIT_SUBSYS:
                    _msg += "Subsystems are not qubits!";
                    break;
                case Type::NOT_BIPARTITE:
                    _msg += "Not bi-partite!";
                    break;
                case Type::NO_CODEWORD:
                    _msg += "Codeword does not exist!";
                    break;
                case Type::PERM_INVALID:
                    _msg += "Invalid permutation!";
                    break;
                case Type::OUT_OF_RANGE:
                    _msg += "Parameter out of range!";
                    break;
                case Type::TYPE_MISMATCH:
                    _msg += "Type mismatch!";
                    break;
                case Type::UNDEFINED_TYPE:
                    _msg += "Not defined for this type!";
                    break;
                case Type::CUSTOM_EXCEPTION:
                    _msg += "CUSTOM EXCEPTION ";
                    break;
            }
            return _msg;
        }
    }; /* class Exception */


} /* namespace qpp */

#endif /* INCLUDE_CLASSES_EXCEPTION_H_ */
