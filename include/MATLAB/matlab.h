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
* \file MATLAB/matlab.h
* \brief Input/output interfacing with MATLAB
*/

#ifndef MATLAB_MATLAB_H_
#define MATLAB_MATLAB_H_

// MATLAB I/O interfacing
// add the path to $MATLAB_INSTALLATION_FOLDER/extern/include in include path

#include "mat.h"  // path to this file is defined in the Makefile
#include "mex.h"  // path to this file is defined in the Makefile

namespace qpp
{

/**
* \brief Loads an Eigen dynamic matrix from a MATLAB .mat file, generic version
* \see qpp::saveMATLABmatrix()
*
* This is the generic version that always throws
* qpp::Exception::Type::UNDEFINED_TYPE. It is specialized only for
* qpp::dmat and qpp::cmat (the only matrix types that can be loaded)
*/
template<typename Derived>
Derived loadMATLABmatrix(const std::string& mat_file,
                         const std::string& var_name)
{
    throw Exception("qpp::loadMATLABmatrix()",
                    Exception::Type::UNDEFINED_TYPE);
}

/**
* \brief Loads an Eigen dynamic matrix from a MATLAB .mat file,
* specialization for double matrices (qpp::dmat)
* \see qpp::saveMATLABmatrix()
*
* The template parameter cannot be automatically deduced and
* must be explicitly provided
*
* Example:
* \code
* // loads a previously saved Eigen dynamic double matrix
* // from the MATLAB file "input.mat"
* auto mat = loadMATLABmatrix<dmat>("input.mat");
* \endcode
*
* \note If \a var_name is a complex matrix, only the real part is loaded
*
* \param mat_file MATALB .mat file
* \param var_name Variable name in the .mat file representing
* the matrix to be loaded
*
* \return Eigen double dynamic matrix (qpp::dmat)
*/
template<>
inline dmat loadMATLABmatrix(const std::string& mat_file,
                             const std::string& var_name)
{
    MATFile* pmat = matOpen(mat_file.c_str(), "r");
    if (!pmat)
    {
        throw std::runtime_error(
            "qpp::loadMATLABmatrix(): Can not open MATLAB file "
            + mat_file + "!");
    }

    mxArray* pa = matGetVariable(pmat, var_name.c_str());
    if (!pa)
        throw std::runtime_error(
            "qpp::loadMATLABmatrix(): Can not load the variable "
            + var_name + " from MATLAB file " + mat_file + "!");

    if (mxGetNumberOfDimensions(pa) != 2) // not a matrix
        throw std::runtime_error(
            "qpp::loadMATLABmatrix(): Loaded variable " + var_name
            + " is not 2-dimensional!");

    if (!mxIsDouble(pa))
        throw std::runtime_error(
            "qpp::loadMATLABmatrix(): Loaded variable " + var_name
            + " is not in double-precision format!");

    idx rows = mxGetM(pa);
    idx cols = mxGetN(pa);

    dmat result(rows, cols);

    std::memcpy(result.data(), mxGetPr(pa),
                sizeof(double) * mxGetNumberOfElements(pa));

    mxDestroyArray(pa);
    matClose(pmat);

    return result;
}

/**
* \brief Loads an Eigen dynamic matrix from a MATLAB .mat file,
* specialization for complex matrices (qpp::cmat)
* \see qpp::saveMATLABmatrix()
*
* The template parameter cannot be automatically deduced and
* must be explicitly provided
*
* Example:
* \code
* // loads a previously saved Eigen dynamic complex matrix
* // from the MATLAB file "input.mat"
* auto mat = loadMATLABmatrix<cmat>("input.mat");
* \endcode
*
* \param mat_file MATALB .mat file
* \param var_name Variable name in the .mat file representing
* the matrix to be loaded
* \return Eigen complex dynamic matrix (qpp::cmat)
*/
template<>
inline cmat loadMATLABmatrix(const std::string& mat_file,
                             const std::string& var_name)
{
    MATFile* pmat = matOpen(mat_file.c_str(), "r");
    if (!pmat)
    {
        throw std::runtime_error(
            "qpp::loadMATLABmatrix(): Can not open MATLAB file "
            + mat_file + "!");
    }

    mxArray* pa = matGetVariable(pmat, var_name.c_str());
    if (!pa)
        throw std::runtime_error(
            "qpp::loadMATLABmatrix(): Can not load the variable "
            + var_name + " from MATLAB file " + mat_file + "!");

    if (mxGetNumberOfDimensions(pa) != 2) // not a matrix
        throw std::runtime_error(
            "qpp::loadMATLABmatrix(): Loaded variable " + var_name
            + " is not 2-dimensional!");

    if (!mxIsDouble(pa))
        throw std::runtime_error(
            "qpp::loadMATLABmatrix(): Loaded variable " + var_name
            + " is not in double-precision format!");

    idx rows = mxGetM(pa);
    idx cols = mxGetN(pa);

    dmat result_re(rows, cols);
    dmat result_im(rows, cols);

    // real part and imaginary part pointers
    double* pa_re = nullptr, * pa_im = nullptr;

    // Populate the real part of the created array.
    pa_re = static_cast<double*>(mxGetPr(pa));
    std::memcpy(result_re.data(), pa_re,
                sizeof(double) * mxGetNumberOfElements(pa));

    if (mxIsComplex(pa)) // populate the imaginary part if exists
    {
        pa_im = static_cast<double*>(mxGetPi(pa));
        std::memcpy(result_im.data(), pa_im,
                    sizeof(double) * mxGetNumberOfElements(pa));
    }
    else // set to zero the imaginary part
    {
        std::memset(result_im.data(), 0,
                    sizeof(double) * mxGetNumberOfElements(pa));
    }

    mxDestroyArray(pa);
    matClose(pmat);

    return (result_re.cast<cplx>()) + 1_i * (result_im.cast<cplx>());
}

/**
* \brief Saves an Eigen dynamic matrix to a MATLAB .mat file, generic version
* \see qpp::loadMATLABmatrix()
*
* This is the generic version that always throws
* qpp::Exception::Type::UNDEFINED_TYPE. It is specialized only for
* qpp::dmat and qpp::cmat (the only matrix types that can be saved)
*/
template<typename Derived>
void saveMATLABmatrix(const Eigen::MatrixBase <Derived>& A,
                      const std::string& mat_file, const std::string& var_name,
                      const std::string& mode)
{
    throw Exception("qpp::saveMATLABmatrix()",
                    Exception::Type::UNDEFINED_TYPE);
}

/**
* \brief Saves an Eigen dynamic matrix to a MATLAB .mat file,
* specialization for double matrices (qpp::dmat)
* \see qpp::loadMATLABmatrix()
*
* \param A Eigen expression over the complex field
* \param mat_file MATALB .mat file
* \param var_name Variable name in the .mat file representing
* the matrix to be saved
* \param mode Saving mode (append, overwrite etc.),
* see MATLAB \a matOpen() documentation for details
*/
template<>
// Eigen::MatrixXd specialization
inline void saveMATLABmatrix(const Eigen::MatrixBase <dmat>& A,
                             const std::string& mat_file,
                             const std::string& var_name,
                             const std::string& mode)
{
    const dmat& rA = A;

    // check zero-size
    if (!internal::_check_nonzero_size(rA))
        throw Exception("qpp::saveMATLABmatrix()", Exception::Type::ZERO_SIZE);

    MATFile* pmat = matOpen(mat_file.c_str(), mode.c_str());
    if (!pmat)
        throw std::runtime_error(
            "qpp::saveMATLABmatrix(): Can not open/create MATLAB file "
            + mat_file + "!");

    mxArray* pa = mxCreateDoubleMatrix(rA.rows(), rA.cols(), mxREAL);
    if (!pa)
        throw std::runtime_error(
            "qpp::saveMATLABmatrix(): mxCreateDoubleMatrix failed!");

    std::memcpy(mxGetPr(pa), rA.data(), sizeof(double) * rA.size());

    if (matPutVariable(pmat, var_name.c_str(), pa))
        throw std::runtime_error(
            "qpp::saveMATLABmatrix(): Can not write the variable "
            + var_name + " to MATLAB file " + mat_file + "!");

    mxDestroyArray(pa);
    matClose(pmat);
}

/**
* \brief Saves an Eigen dynamic matrix to a MATLAB .mat file,
* specialization for complex matrices (qpp::cmat)
* \see qpp::loadMATLABmatrix()
*
* \param A Eigen expression over the complex field
* \param mat_file MATALB .mat file
* \param var_name Variable name in the .mat file representing
* the matrix to be saved
* \param mode Saving mode (append, overwrite etc.),
* see MATLAB \a matOpen() documentation for details
*/
template<>
// Eigen::MatrixXcd specialization
inline void saveMATLABmatrix(const Eigen::MatrixBase <cmat>& A,
                             const std::string& mat_file,
                             const std::string& var_name,
                             const std::string& mode)
{
    const cmat& rA = A;

    // check zero-size
    if (!internal::_check_nonzero_size(rA))
        throw Exception("qpp::saveMATLABmatrix()", Exception::Type::ZERO_SIZE);

    // cast the input to a double (internal MATLAB format)
    dmat tmp_re = rA.real();
    dmat tmp_im = rA.imag();

    MATFile* pmat = matOpen(mat_file.c_str(), mode.c_str());
    if (!pmat)
        throw std::runtime_error(
            "qpp::saveMATLABmatrix(): Can not open/create MATLAB file "
            + mat_file + "!");

    mxArray* pa = mxCreateDoubleMatrix(
        tmp_re.rows(), tmp_re.cols(), mxCOMPLEX);
    if (!pa)
        throw std::runtime_error(
            "qpp::saveMATLABmatrix(): mxCreateDoubleMatrix failed!");

    double* pa_re, * pa_im;

    /* Populate the real part of the created array. */
    pa_re = static_cast<double*>(mxGetPr(pa));
    std::memcpy(pa_re, tmp_re.data(), sizeof(double) * tmp_re.size());

    /* Populate the imaginary part of the created array. */
    pa_im = static_cast<double*>(mxGetPi(pa));
    std::memcpy(pa_im, tmp_im.data(), sizeof(double) * tmp_im.size());

    if (matPutVariable(pmat, var_name.c_str(), pa))
        throw std::runtime_error(
            "qpp::saveMATLABmatrix(): Can not write the variable "
            + var_name + " to MATLAB file " + mat_file + "!");

    mxDestroyArray(pa);
    matClose(pmat);
}

} /* namespace qpp */
#endif /* MATLAB_MATLAB_H_ */
