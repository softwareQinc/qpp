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
* \file MATLAB/matlab.h
* \brief Input/output interfacing with MATLAB
*/

#ifndef MATLAB_MATLAB_H_
#define MATLAB_MATLAB_H_

// MATLAB I/O interfacing
// add the path to $MATLAB_INSTALLATION_FOLDER/extern/include in include path

#include "mat.h" // path to this file is defined in the Makefile
#include "mex.h" // path to this file is defined in the Makefile

namespace qpp {
/**
* \brief Loads a complex Eigen dynamic matrix from a MATLAB .mat file,
* \see qpp::saveMATLAB()
*
* The template parameter cannot be automatically deduced and
* must be explicitly provided
*
* Example:
* \code
* // loads a previously saved Eigen ket
* // from the MATLAB file "input.mat"
* ket psi = loadMATLAB<ket>("input.mat");
* \endcode
*
* \tparam Derived Complex Eigen type
* \param mat_file MATALB .mat file
* \param var_name Variable name in the .mat file representing the matrix
* to be loaded
*
* \return Eigen dynamic matrix
*/
template <typename Derived> // double
typename std::enable_if<std::is_same<typename Derived::Scalar, cplx>::value,
                        dyn_mat<cplx>>::type
loadMATLAB(const std::string& mat_file, const std::string& var_name) {
    MATFile* pmat = matOpen(mat_file.c_str(), "r");

    // EXCEPTION CHECKS

    if (!pmat) {
        throw std::runtime_error(
            "qpp::loadMATLAB(): Can not open MATLAB file " + mat_file + "!");
    }

    mxArray* pa = matGetVariable(pmat, var_name.c_str());
    if (!pa)
        throw std::runtime_error(
            "qpp::loadMATLAB(): Can not load the variable " + var_name +
            " from MATLAB file " + mat_file + "!");

    if (mxGetNumberOfDimensions(pa) != 2) // not a matrix
        throw std::runtime_error("qpp::loadMATLAB(): Loaded variable " +
                                 var_name + " is not 2-dimensional!");

    if (!mxIsDouble(pa))
        throw std::runtime_error("qpp::loadMATLAB(): Loaded variable " +
                                 var_name +
                                 " is not in double-precision format!");
    // END EXCEPTION CHECKS

    idx rows = mxGetM(pa);
    idx cols = mxGetN(pa);

    dyn_mat<double> result_re(rows, cols);
    dyn_mat<double> result_im(rows, cols);

    // real part and imaginary part pointers
    double* pa_re = nullptr;
    double* pa_im = nullptr;

    // Populate the real part of the created array.
    pa_re = mxGetPr(pa);
    std::memcpy(result_re.data(), pa_re,
                sizeof(double) * mxGetNumberOfElements(pa));

    if (mxIsComplex(pa)) // populate the imaginary part if exists
    {
        pa_im = mxGetPi(pa);
        std::memcpy(result_im.data(), pa_im,
                    sizeof(double) * mxGetNumberOfElements(pa));
    } else // set to zero the imaginary part
    {
        std::memset(result_im.data(), 0,
                    sizeof(double) * mxGetNumberOfElements(pa));
    }

    mxDestroyArray(pa);
    matClose(pmat);

    return (result_re.cast<cplx>()) + 1_i * (result_im.cast<cplx>());
}

/**
* \brief Loads a non-complex Eigen dynamic matrix from a MATLAB .mat file,
* \see qpp::saveMATLAB()
*
* The template parameter cannot be automatically deduced and
* must be explicitly provided
*
* Example:
* \code
* // loads a previously saved Eigen dynamic double matrix
* // from the MATLAB file "input.mat"
* dmat mat = loadMATLAB<dmat>("input.mat");
* \endcode
*
* \tparam Derived Non-complex Eigen type
* \param mat_file MATALB .mat file
* \param var_name Variable name in the .mat file representing
* the matrix to be loaded
*
* \return Eigen dynamic matrix
*/
template <typename Derived> // cplx
typename std::enable_if<!std::is_same<typename Derived::Scalar, cplx>::value,
                        dyn_mat<typename Derived::Scalar>>::type
loadMATLAB(const std::string& mat_file, const std::string& var_name) {
    MATFile* pmat = matOpen(mat_file.c_str(), "r");

    // EXCEPTION CHECKS

    if (!pmat) {
        throw std::runtime_error(
            "qpp::loadMATLAB(): Can not open MATLAB file " + mat_file + "!");
    }

    mxArray* pa = matGetVariable(pmat, var_name.c_str());
    if (!pa)
        throw std::runtime_error(
            "qpp::loadMATLAB(): Can not load the variable " + var_name +
            " from MATLAB file " + mat_file + "!");

    if (mxGetNumberOfDimensions(pa) != 2) // not a matrix
        throw std::runtime_error("qpp::loadMATLAB(): Loaded variable " +
                                 var_name + " is not 2-dimensional!");

    if (!mxIsDouble(pa))
        throw std::runtime_error("qpp::loadMATLAB(): Loaded variable " +
                                 var_name +
                                 " is not in double-precision format!");
    // END EXCEPTION CHECKS

    idx rows = mxGetM(pa);
    idx cols = mxGetN(pa);

    dyn_mat<double> result(rows, cols);

    std::memcpy(result.data(), mxGetPr(pa),
                sizeof(double) * mxGetNumberOfElements(pa));

    mxDestroyArray(pa);
    matClose(pmat);

    // cast back to the original type
    return result.cast<typename Derived::Scalar>();
}

/**
* \brief Saves a complex Eigen dynamic matrix to a MATLAB .mat file,
* \see qpp::loadMATLAB()
*
* \tparam Complex Eigen type
* \param A Eigen expression over the complex field
* \param mat_file MATALB .mat file
* \param var_name Variable name in the .mat file representing
* the matrix to be saved
* \param mode Saving mode (append, overwrite etc.),
* see MATLAB \a matOpen() documentation for details
*/
template <typename Derived>
// double
typename std::enable_if<
    std::is_same<typename Derived::Scalar, cplx>::value>::type
saveMATLAB(const Eigen::MatrixBase<Derived>& A, const std::string& mat_file,
           const std::string& var_name, const std::string& mode) {
    const dyn_mat<cplx>& rA = A.derived();

    // EXCEPTION CHECKS

    // check zero-size
    if (!internal::check_nonzero_size(rA))
        throw exception::ZeroSize("qpp::saveMATLAB()");

    // cast the input to a double (internal MATLAB format)
    dyn_mat<double> tmp_re = rA.real();
    dyn_mat<double> tmp_im = rA.imag();

    MATFile* pmat = matOpen(mat_file.c_str(), mode.c_str());
    if (!pmat)
        throw std::runtime_error(
            "qpp::saveMATLAB(): Can not open/create MATLAB file " + mat_file +
            "!");

    mxArray* pa = mxCreateDoubleMatrix(tmp_re.rows(), tmp_re.cols(), mxCOMPLEX);
    if (!pa)
        throw std::runtime_error(
            "qpp::saveMATLAB(): mxCreateDoubleMatrix failed!");
    // END EXCEPTION CHECKS

    // real part and imaginary part pointers
    double* pa_re = nullptr;
    double* pa_im = nullptr;

    /* Populate the real part of the created array. */
    pa_re = mxGetPr(pa);
    std::memcpy(pa_re, tmp_re.data(), sizeof(double) * tmp_re.size());

    /* Populate the imaginary part of the created array. */
    pa_im = mxGetPi(pa);
    std::memcpy(pa_im, tmp_im.data(), sizeof(double) * tmp_im.size());

    if (matPutVariable(pmat, var_name.c_str(), pa))
        throw std::runtime_error(
            "qpp::saveMATLAB(): Can not write the variable " + var_name +
            " to MATLAB file " + mat_file + "!");

    mxDestroyArray(pa);
    matClose(pmat);
}

/**
* \brief Saves a non-complex Eigen dynamic matrix to a MATLAB .mat file,
* \see qpp::loadMATLAB()
*
* \tparam Npn-complex Eigen type
* \param A Non-complex Eigen expression
* \param mat_file MATALB .mat file
* \param var_name Variable name in the .mat file representing
* the matrix to be saved
* \param mode Saving mode (append, overwrite etc.),
* see MATLAB \a matOpen() documentation for details
*/
template <typename Derived>
// cplx
typename std::enable_if<
    !std::is_same<typename Derived::Scalar, cplx>::value>::type
saveMATLAB(const Eigen::MatrixBase<Derived>& A, const std::string& mat_file,
           const std::string& var_name, const std::string& mode) {
    // cast to double, as MATLAB doesn't work with other types
    const dyn_mat<double>& rA = A.template cast<double>();

    // EXCEPTION CHECKS

    // check zero-size
    if (!internal::check_nonzero_size(rA))
        throw exception::ZeroSize("qpp::saveMATLAB()");

    MATFile* pmat = matOpen(mat_file.c_str(), mode.c_str());
    if (!pmat)
        throw std::runtime_error(
            "qpp::saveMATLAB(): Can not open/create MATLAB file " + mat_file +
            "!");

    mxArray* pa = mxCreateDoubleMatrix(rA.rows(), rA.cols(), mxREAL);
    if (!pa)
        throw std::runtime_error(
            "qpp::saveMATLAB(): mxCreateDoubleMatrix failed!");
    // END EXCEPTION CHECKS

    std::memcpy(mxGetPr(pa), rA.data(), sizeof(double) * rA.size());

    if (matPutVariable(pmat, var_name.c_str(), pa))
        throw std::runtime_error(
            "qpp::saveMATLAB(): Can not write the variable " + var_name +
            " to MATLAB file " + mat_file + "!");

    mxDestroyArray(pa);
    matClose(pmat);
}

} /* namespace qpp */
#endif /* MATLAB_MATLAB_H_ */
