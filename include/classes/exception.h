/*
 * Quantum++
 *
 * Copyright (c) 2013 - 2016 Vlad Gheorghiu (vgheorgh@gmail.com)
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
* \file classes/exception.h
* \brief Exceptions
*/

#ifndef CLASSES_EXCEPTION_H_
#define CLASSES_EXCEPTION_H_

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
    * \brief Exception types, add more here if needed
    * \see qpp::Exception::construct_exception_msg_()
    */
    enum class Type // exception types
    {
        UNKNOWN_EXCEPTION = 1,
        /*!< Unknown exception */
                ZERO_SIZE,
        /*!< Zero sized object, e.g. empty Eigen::Matrix
         * or std::vector<> with no elements */
                MATRIX_NOT_SQUARE,
        /*!< Eigen::Matrix is not square */
                MATRIX_NOT_CVECTOR,
        /*!< Eigen::Matrix is not a column vector */
                MATRIX_NOT_RVECTOR,
        /*!< Eigen::Matrix is not a row vector */
                MATRIX_NOT_VECTOR,
        /*!< Eigen::Matrix is not a row/column vector */
                MATRIX_NOT_SQUARE_OR_CVECTOR,
        /*!< Eigen::Matrix is not square nor a column vector */
                MATRIX_NOT_SQUARE_OR_RVECTOR,
        /*!< Eigen::Matrix is not square nor a row vector */
                MATRIX_NOT_SQUARE_OR_VECTOR,
        /*!< Eigen::Matrix is not square nor a row/column vector */
                MATRIX_MISMATCH_SUBSYS,
        /*!< Matrix size mismatch subsystem sizes (e.g. in qpp::apply()) */
                DIMS_INVALID,
        /*!< std::vector<idx> of dimensions has zero size or contains zeros */
                DIMS_NOT_EQUAL,
        /*!< Local/global dimensions are not equal */
                DIMS_MISMATCH_MATRIX,
        /*!< Product of the elements of std::vector<idx> of dimensions
         * is not equal to the number of rows of Eigen::Matrix
         * (assumed to be a square matrix) */
                DIMS_MISMATCH_CVECTOR,
        /*!< Product of the elements of std::vector<idx> of dimensions
         * is not equal to the number of elements of Eigen::Matrix
         * (assumed to be a column vector) */
                DIMS_MISMATCH_RVECTOR,
        /*!< Product of the elements of std::vector<idx> of dimensions
         * is not equal to the number of elements of Eigen::Matrix
         * (assumed to be a row vector) */
                DIMS_MISMATCH_VECTOR,
        /*!< Product of the elements of std::vector<idx> of dimensions
         * is not equal to the number of elements of Eigen::Matrix
         * (assumed to be a row/column vector) */
                SUBSYS_MISMATCH_DIMS,
        /*!< std::vector<idx> of subsystem labels has duplicates,
         * or has entries that are larger than the size of
         * the std::vector<idx> of dimensions */
                PERM_INVALID,
        /*!< std::vector<idx> does note represent a valid permutation */
                PERM_MISMATCH_DIMS,
        /*!< Size of the std::vector<idx> representing the permutation
         * is different from the size of the std::vector<idx> of dimensions */
                NOT_QUBIT_MATRIX,
        /*!<  Eigen::Matrix is not 2 x 2*/
                NOT_QUBIT_CVECTOR,
        /*!<  Eigen::Matrix is not 2 x 1*/
                NOT_QUBIT_RVECTOR,
        /*!<  Eigen::Matrix is not 1 x 2*/
                NOT_QUBIT_VECTOR,
        /*!<  Eigen::Matrix is not 1 x 2 nor 2 x 1*/
                NOT_QUBIT_SUBSYS,
        /*!< Subsystems are not 2-dimensional */
                NOT_BIPARTITE,
        /*!< std::vector<idx> of dimensions has size different from 2 */
                NO_CODEWORD,
        /*!< Codeword does not exist, thrown when calling
         * qpp::Codes::codeword() with invalid index \a i */
                OUT_OF_RANGE,
        /*!< Parameter out of range */
                TYPE_MISMATCH,
        /*!< Scalar types do not match */
                SIZE_MISMATCH,
        /*!< Sizes do not match */
                UNDEFINED_TYPE,
        /*!< Templated specialization not defined for this type */
                CUSTOM_EXCEPTION
        /*!< Custom exception, user must provide a custom message */
    };

    /**
    * \brief Constructs an exception
    *
    * \param where Text representing where the exception occured
    * \param type Exception type, defined in qpp::Exception::Type
    */
    Exception(const std::string& where, const Type& type) :
            where_{where}, msg_{}, type_{type}, custom_{}
    {
        construct_exception_msg_();
    }

    /**
    * \brief Constructs an exception
    *
    * \overload
    *
    * \param where Text representing where the exception occured
    * \param custom Exception description
    */
    Exception(const std::string& where, const std::string& custom) :
            where_{where}, msg_{}, type_{Type::CUSTOM_EXCEPTION},
            custom_{custom}
    {
        construct_exception_msg_();
        msg_ += custom; // add the custom message at the end
    }

    /**
    * \brief Overrides std::exception::what()
    *
    * \return Exception description
    */
    virtual const char* what() const noexcept override
    {
        return msg_.c_str();
    }

private:
    std::string where_, msg_;
    Type type_;
    std::string custom_;

    /**
    * \brief Constructs the exception description from its type
    * \see qpp::Exception::Type
    *
    * Must modify the code of this function if more exceptions are added
    */
    void construct_exception_msg_()
    {
        msg_ += "IN ";
        msg_ += where_;
        msg_ += ": ";

        switch (type_)
        {
            case Type::UNKNOWN_EXCEPTION:
                msg_ += "UNKNOWN EXCEPTION!";
                break;
            case Type::ZERO_SIZE:
                msg_ += "Object has zero size!";
                break;
            case Type::MATRIX_NOT_SQUARE:
                msg_ += "Matrix is not square!";
                break;
            case Type::MATRIX_NOT_CVECTOR:
                msg_ += "Matrix is not column vector!";
                break;
            case Type::MATRIX_NOT_RVECTOR:
                msg_ += "Matrix is not row vector!";
                break;
            case Type::MATRIX_NOT_VECTOR:
                msg_ += "Matrix is not vector!";
                break;
            case Type::MATRIX_NOT_SQUARE_OR_CVECTOR:
                msg_ += "Matrix is not square nor column vector!";
                break;
            case Type::MATRIX_NOT_SQUARE_OR_RVECTOR:
                msg_ += "Matrix is not square nor row vector!";
                break;
            case Type::MATRIX_NOT_SQUARE_OR_VECTOR:
                msg_ += "Matrix is not square nor vector!";
                break;
            case Type::MATRIX_MISMATCH_SUBSYS:
                msg_ += "Matrix mismatch subsystems!";
                break;
            case Type::DIMS_INVALID:
                msg_ += "Invalid dimension(s)!";
                break;
            case Type::DIMS_NOT_EQUAL:
                msg_ += "Dimensions not equal!";
                break;
            case Type::DIMS_MISMATCH_MATRIX:
                msg_ += "Dimension(s) mismatch matrix size!";
                break;
            case Type::DIMS_MISMATCH_CVECTOR:
                msg_ += "Dimension(s) mismatch column vector!";
                break;
            case Type::DIMS_MISMATCH_RVECTOR:
                msg_ += "Dimension(s) mismatch row vector!";
                break;
            case Type::DIMS_MISMATCH_VECTOR:
                msg_ += "Dimension(s) mismatch vector!";
                break;
            case Type::SUBSYS_MISMATCH_DIMS:
                msg_ += "Subsystems mismatch dimensions!";
                break;
            case Type::PERM_INVALID:
                msg_ += "Invalid permutation!";
                break;
            case Type::PERM_MISMATCH_DIMS:
                msg_ += "Permutation mismatch dimensions!";
                break;
            case Type::NOT_QUBIT_MATRIX:
                msg_ += "Matrix is not 2 x 2!";
                break;
            case Type::NOT_QUBIT_CVECTOR:
                msg_ += "Column vector is not 2 x 1!";
                break;
            case Type::NOT_QUBIT_RVECTOR:
                msg_ += "Row vector is not 1 x 2!";
                break;
            case Type::NOT_QUBIT_VECTOR:
                msg_ += "Vector is not 2 x 1 nor 1 x 2!";
                break;
            case Type::NOT_QUBIT_SUBSYS:
                msg_ += "Subsystems are not qubits!";
                break;
            case Type::NOT_BIPARTITE:
                msg_ += "Not bi-partite!";
                break;
            case Type::NO_CODEWORD:
                msg_ += "Codeword does not exist!";
                break;
            case Type::OUT_OF_RANGE:
                msg_ += "Parameter out of range!";
                break;
            case Type::TYPE_MISMATCH:
                msg_ += "Type mismatch!";
                break;
            case Type::SIZE_MISMATCH:
                msg_ += "Size mismatch!";
                break;
            case Type::UNDEFINED_TYPE:
                msg_ += "Not defined for this type!";
                break;
            case Type::CUSTOM_EXCEPTION:
                msg_ += "CUSTOM EXCEPTION ";
                break;
        }
    }
}; /* class Exception */

} /* namespace qpp */

#endif /* CLASSES_EXCEPTION_H_ */
