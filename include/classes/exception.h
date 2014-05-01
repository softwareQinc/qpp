/*
 * exception.h
 *
 *  Created on: Apr 2, 2014
 *      Author: vlad
 */

#ifndef EXCEPTION_H_
#define EXCEPTION_H_

#include <stdexcept>
#include <string>

// exception class, function parameter checking

namespace qpp
{

class Exception: public std::exception
{
public:
	enum class Type // exception types
	{
		// unknown exception
		UNKNOWN_EXCEPTION = 1,
		// zero sized object, e.g. non-initialized Eigen::Matrix
		// or std::vector with no elements
		ZERO_SIZE,
		// Eigen::Matrix is not square
		MATRIX_NOT_SQUARE,
		// Eigen::Matrix is not a column vector
		MATRIX_NOT_CVECTOR,
		// Eigen::Matrix is not a row vector
		MATRIX_NOT_RVECTOR,
		// Eigen::Matrix is not a row/column vector
		MATRIX_NOT_VECTOR,
		// Eigen::Matrix is not square nor a column vector
		MATRIX_NOT_SQUARE_OR_CVECTOR,
		// Eigen::Matrix is not square nor a row vector
		MATRIX_NOT_SQUARE_OR_RVECTOR,
		// Eigen::Matrix is not square nor a row/column vector
		MATRIX_NOT_SQUARE_OR_VECTOR,
		// std::vector<size_t> of dimensions has zero size or contains zeros
		DIMS_INVALID,
		// std::vector<size_t> of dimensions contains un-equal elements
		DIMS_NOT_EQUAL,
		// product of dimenison's std::vector<size_t> not equal to
		// the number of rows of Eigen::Matrix (assumed to be square)
		DIMS_MISMATCH_MATRIX,
		// product of dimenison's std::vector<size_t> not equal to
		// the number of rows of Eigen::Matrix (assumed to be column vector)
		DIMS_MISMATCH_CVECTOR,
		// product of dimenison's std::vector<size_t> not equal to
		// the number of cols of Eigen::Matrix (assumed to be row vector)
		DIMS_MISMATCH_RVECTOR,
		// product of dimenison's std::vector<size_t> not equal to
		// the size of Eigen::Matrix (assumed to be row/column vector)
		DIMS_MISMATCH_VECTOR,
		// std::vector<size_t> subsystem vector has duplicatates, or
		// has entries that are larger than the size of std::vector<size_t>
		// of dimensions
		SUBSYS_MISMATCH_DIMS,
		// invalid std::vector<size_t> permutation
		PERM_INVALID,
		// Eigen::Matrix is not 2 x 2
		NOT_QUBIT_GATE,
		// not 2-dimensional subsystems
		NOT_QUBIT_SUBSYS,
		// std::vector<size_t> of dimensions has size different from 2
		NOT_BIPARTITE,
		// parameter out of range
		OUT_OF_RANGE,
		// template function not defined for this type
		UNDEFINED_TYPE,
		// custom exception, user must provide a custom message
		CUSTOM_EXCEPTION
	};

	Exception(const std::string & where, const Type& type) :
			_where(where), _msg(), _type(type), _custom()
	{
		_construct_exception_msg();
	}

	Exception(const std::string & where, const std::string & custom) :
			_where(where), _msg(), _type(Type::CUSTOM_EXCEPTION), _custom(
					custom)
	{
		_construct_exception_msg();
		_msg += custom; // add the custom message at the end
	}

	virtual const char* what() const noexcept override
	{
		return _msg.c_str();
	}

	virtual ~Exception() noexcept
	{
	}
private:
	std::string _where, _msg;
	Type _type;
	std::string _custom;

	// construct exception messages
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
		case Type::PERM_INVALID:
			_msg += "Invalid permutation!";
			break;
		case Type::NOT_QUBIT_GATE:
			_msg += "Matrix is not qubit gate!";
			break;
		case Type::NOT_QUBIT_SUBSYS:
			_msg += "Subsystems are not qubits!";
			break;
		case Type::NOT_BIPARTITE:
			_msg += "Not bipartite!";
			break;
		case Type::OUT_OF_RANGE:
			_msg += "Parameter out of range!";
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
};

} /* namespace qpp */

#endif /* EXCEPTION_H_ */
