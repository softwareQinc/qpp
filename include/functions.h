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
* \file functions.h
* \brief Generic quantum computing functions
*/

#ifndef FUNCTIONS_H_
#define FUNCTIONS_H_

namespace qpp {
// Eigen function wrappers
/**
* \brief Transpose
*
* \param A Eigen expression
* \return Transpose of \a A, as a dynamic matrix
* over the same scalar field as \a A
*/
template <typename Derived>
dyn_mat<typename Derived::Scalar>
transpose(const Eigen::MatrixBase<Derived>& A) {
    const dyn_mat<typename Derived::Scalar>& rA = A.derived();

    // check zero-size
    if (!internal::check_nonzero_size(rA))
        throw exception::ZeroSize("qpp::transpose()");

    return rA.transpose();
}

/**
* \brief Complex conjugate
*
* \param A Eigen expression
* \return Complex conjugate of \a A, as a dynamic matrix
* over the same scalar field as \a A
*/
template <typename Derived>
dyn_mat<typename Derived::Scalar>
conjugate(const Eigen::MatrixBase<Derived>& A) {
    const dyn_mat<typename Derived::Scalar>& rA = A.derived();

    // check zero-size
    if (!internal::check_nonzero_size(rA))
        throw exception::ZeroSize("qpp::conjugate()");

    return rA.conjugate();
}

/**
* \brief Adjoint
*
* \param A Eigen expression
* \return Adjoint (Hermitian conjugate) of \a A, as a dynamic matrix
* over the same scalar field as \a A
*/
template <typename Derived>
dyn_mat<typename Derived::Scalar> adjoint(const Eigen::MatrixBase<Derived>& A) {
    const dyn_mat<typename Derived::Scalar>& rA = A.derived();

    // EXCEPTION CHECKS

    // check zero-size
    if (!internal::check_nonzero_size(rA))
        throw exception::ZeroSize("qpp::adjoint()");
    // END EXCEPTION CHECKS

    return rA.adjoint();
}

/**
* \brief Inverse
*
* \param A Eigen expression
* \return Inverse of \a A, as a dynamic matrix
* over the same scalar field as \a A
*/
template <typename Derived>
dyn_mat<typename Derived::Scalar> inverse(const Eigen::MatrixBase<Derived>& A) {
    const dyn_mat<typename Derived::Scalar>& rA = A.derived();

    // EXCEPTION CHECKS

    // check zero-size
    if (!internal::check_nonzero_size(rA))
        throw exception::ZeroSize("qpp::inverse()");
    // END EXCEPTION CHECKS

    return rA.inverse();
}

/**
* \brief Trace
*
* \param A Eigen expression
* \return Trace of \a A, as a scalar over the same scalar field as \a A
*/
template <typename Derived>
typename Derived::Scalar trace(const Eigen::MatrixBase<Derived>& A) {
    const dyn_mat<typename Derived::Scalar>& rA = A.derived();

    // EXCEPTION CHECKS

    // check zero-size
    if (!internal::check_nonzero_size(rA))
        throw exception::ZeroSize("qpp::trace()");
    // END EXCEPTION CHECKS

    return rA.trace();
}

/**
* \brief Determinant
*
* \param A Eigen expression
* \return Determinant of \a A, as a scalar over the same scalar field as \a A.
* Returns \f$\pm \infty\f$ when the determinant overflows/underflows.
*/
template <typename Derived>
typename Derived::Scalar det(const Eigen::MatrixBase<Derived>& A) {
    const dyn_mat<typename Derived::Scalar>& rA = A.derived();

    // EXCEPTION CHECKS

    // check zero-size
    if (!internal::check_nonzero_size(rA))
        throw exception::ZeroSize("qpp::det()");
    // END EXCEPTION CHECKS

    return rA.determinant();
}

/**
* \brief Logarithm of the determinant
*
* Useful when the determinant overflows/underflows
*
* \param A Eigen expression
* \return Logarithm of the determinant of \a A, as a scalar
* over the same scalar field as \a A
*/
template <typename Derived>
typename Derived::Scalar logdet(const Eigen::MatrixBase<Derived>& A) {
    const dyn_mat<typename Derived::Scalar>& rA = A.derived();

    // EXCEPTION CHECKS

    // check zero-size
    if (!internal::check_nonzero_size(rA))
        throw exception::ZeroSize("qpp::logdet()");

    // check square matrix
    if (!internal::check_square_mat(rA))
        throw exception::MatrixNotSquare("qpp::logdet()");
    // END EXCEPTION CHECKS

    Eigen::PartialPivLU<dyn_mat<typename Derived::Scalar>> lu(rA);
    dyn_mat<typename Derived::Scalar> U =
        lu.matrixLU().template triangularView<Eigen::Upper>();
    typename Derived::Scalar result = std::log(U(0, 0));

    for (idx i = 1; i < static_cast<idx>(rA.rows()); ++i)
        result += std::log(U(i, i));

    return result;
}

/**
* \brief Element-wise sum of \a A
*
* \param A Eigen expression
* \return Element-wise sum of \a A, as a scalar
* over the same scalar field as \a A
*/
template <typename Derived>
typename Derived::Scalar sum(const Eigen::MatrixBase<Derived>& A) {
    const dyn_mat<typename Derived::Scalar>& rA = A.derived();

    // EXCEPTION CHECKS

    // check zero-size
    if (!internal::check_nonzero_size(rA))
        throw exception::ZeroSize("qpp::sum()");
    // END EXCEPTION CHECKS

    return rA.sum();
}

/**
* \brief Element-wise product of \a A
*
* \param A Eigen expression
* \return Element-wise product of \a A, as a scalar
* over the same scalar field as \a A
*/
template <typename Derived>
typename Derived::Scalar prod(const Eigen::MatrixBase<Derived>& A) {
    const dyn_mat<typename Derived::Scalar>& rA = A.derived();

    // EXCEPTION CHECKS

    // check zero-size
    if (!internal::check_nonzero_size(rA))
        throw exception::ZeroSize("qpp::prod()");
    // END EXCEPTION CHECKS

    return rA.prod();
}

/**
* \brief Frobenius norm
*
* \param A Eigen expression
* \return Frobenius norm of \a A
*/
template <typename Derived>
double norm(const Eigen::MatrixBase<Derived>& A) {
    const dyn_mat<typename Derived::Scalar>& rA = A.derived();

    // EXCEPTION CHECKS

    // check zero-size
    if (!internal::check_nonzero_size(rA))
        throw exception::ZeroSize("qpp::norm()");
    // END EXCEPTION CHECKS

    // convert matrix to complex then return its norm
    return (rA.template cast<cplx>()).norm();
}

/**
* \brief Full eigen decomposition
* \see qpp::heig()
*
* \param A Eigen expression
* \return Pair of:  1. Eigenvalues of \a A, as a complex dynamic column vector,
* and 2. Eigenvectors of \a A, as columns of a complex dynamic matrix
*/
template <typename Derived>
std::pair<dyn_col_vect<cplx>, cmat>

eig(const Eigen::MatrixBase<Derived>& A) {
    const dyn_mat<typename Derived::Scalar>& rA = A.derived();

    // EXCEPTION CHECKS

    // check zero-size
    if (!internal::check_nonzero_size(rA))
        throw exception::ZeroSize("qpp::eig()");

    // check square matrix
    if (!internal::check_square_mat(rA))
        throw exception::MatrixNotSquare("qpp::eig()");
    // END EXCEPTION CHECKS

    Eigen::ComplexEigenSolver<cmat> es(rA.template cast<cplx>());

    return std::make_pair(es.eigenvalues(), es.eigenvectors());
}

/**
* \brief Eigenvalues
* \see qpp::hevals()
*
* \param A Eigen expression
* \return Eigenvalues of \a A, as a complex dynamic column vector
*/
template <typename Derived>
dyn_col_vect<cplx> evals(const Eigen::MatrixBase<Derived>& A) {
    const dyn_mat<typename Derived::Scalar>& rA = A.derived();

    // EXCEPTION CHECKS

    // check zero-size
    if (!internal::check_nonzero_size(rA))
        throw exception::ZeroSize("qpp::evals()");

    // check square matrix
    if (!internal::check_square_mat(rA))
        throw exception::MatrixNotSquare("qpp::evals()");
    // END EXCEPTION CHECKS

    return eig(rA).first;
}

/**
* \brief Eigenvectors
* \see qpp::hevects()
*
* \param A Eigen expression
* \return Eigenvectors of \a A, as columns of a complex dynamic matrix
*/
template <typename Derived>
cmat evects(const Eigen::MatrixBase<Derived>& A) {
    const dyn_mat<typename Derived::Scalar>& rA = A.derived();

    // EXCEPTION CHECKS

    // check zero-size
    if (!internal::check_nonzero_size(rA))
        throw exception::ZeroSize("qpp::evects()");

    // check square matrix
    if (!internal::check_square_mat(rA))
        throw exception::MatrixNotSquare("qpp::evects()");
    // END EXCEPTION CHECKS

    Eigen::ComplexEigenSolver<cmat> es(rA.template cast<cplx>());

    return eig(rA).second;
}

/**
* \brief Full eigen decomposition of Hermitian expression
* \see qpp::eig()
*
* \param A Eigen expression
* \return Pair of:  1. Eigenvalues of \a A, as a real dynamic column vector,
* and 2. Eigenvectors of \a A, as columns of a complex dynamic matrix
*/
template <typename Derived>
std::pair<dyn_col_vect<double>, cmat>

heig(const Eigen::MatrixBase<Derived>& A) {
    const dyn_mat<typename Derived::Scalar>& rA = A.derived();

    // EXCEPTION CHECKS

    // check zero-size
    if (!internal::check_nonzero_size(rA))
        throw exception::ZeroSize("qpp::heig()");

    // check square matrix
    if (!internal::check_square_mat(rA))
        throw exception::MatrixNotSquare("qpp::heig()");
    // END EXCEPTION CHECKS

    Eigen::SelfAdjointEigenSolver<cmat> es(rA.template cast<cplx>());

    return std::make_pair(es.eigenvalues(), es.eigenvectors());
}

/**
* \brief Hermitian eigenvalues
* \see qpp::evals()
*
* \param A Eigen expression
* \return Eigenvalues of Hermitian \a A, as a real dynamic column vector
*/
template <typename Derived>
dyn_col_vect<double> hevals(const Eigen::MatrixBase<Derived>& A) {
    const dyn_mat<typename Derived::Scalar>& rA = A.derived();

    // EXCEPTION CHECKS

    // check zero-size
    if (!internal::check_nonzero_size(rA))
        throw exception::ZeroSize("qpp::hevals()");

    // check square matrix
    if (!internal::check_square_mat(rA))
        throw exception::MatrixNotSquare("qpp::hevals()");
    // END EXCEPTION CHECKS

    return heig(rA).first;
}

/**
* \brief Hermitian eigenvectors
* \see qpp::evects()
*
* \param A Eigen expression
* \return Eigenvectors of Hermitian \a A, as columns of a complex matrix
*/
template <typename Derived>
cmat hevects(const Eigen::MatrixBase<Derived>& A) {
    const dyn_mat<typename Derived::Scalar>& rA = A.derived();

    // EXCEPTION CHECKS

    // check zero-size
    if (!internal::check_nonzero_size(rA))
        throw exception::ZeroSize("qpp::hevects()");

    // check square matrix
    if (!internal::check_square_mat(rA))
        throw exception::MatrixNotSquare("qpp::hevects()");
    // END EXCEPTION CHECKS

    return heig(rA).second;
}

/**
* \brief Full singular value decomposition
*
* \param A Eigen expression
* \return Tuple of: 1. Left sigular vectors of \a A, as columns of a complex
* dynamic matrix, 2. Singular values of \a A, ordered in decreasing order,
* as a real dynamic column vector, and 3. Right singular vectors of \a A,
* as columns of a complex dynamic matrix
*/
template <typename Derived>
std::tuple<cmat, dyn_col_vect<double>, cmat>

svd(const Eigen::MatrixBase<Derived>& A) {
    const dyn_mat<typename Derived::Scalar>& rA = A.derived();

    // EXCEPTION CHECKS

    // check zero-size
    if (!internal::check_nonzero_size(rA))
        throw exception::ZeroSize("qpp::svd()");
    // END EXCEPTION CHECKS

    Eigen::JacobiSVD<dyn_mat<typename Derived::Scalar>> sv(
        rA, Eigen::DecompositionOptions::ComputeFullU |
                Eigen::DecompositionOptions::ComputeFullV);

    return std::make_tuple(sv.matrixU(), sv.singularValues(), sv.matrixV());
}

/**
* \brief Singular values
*
* \param A Eigen expression
* \return Singular values of \a A, ordered in decreasing order,
* as a real dynamic column vector
*/
template <typename Derived>
dyn_col_vect<double> svals(const Eigen::MatrixBase<Derived>& A) {
    const dyn_mat<typename Derived::Scalar>& rA = A.derived();

    // EXCEPTION CHECKS

    // check zero-size
    if (!internal::check_nonzero_size(rA))
        throw exception::ZeroSize("qpp::svals()");
    // END EXCEPTION CHECKS

    Eigen::JacobiSVD<dyn_mat<typename Derived::Scalar>> sv(rA);

    return sv.singularValues();
}

/**
* \brief Left singular vectors
*
* \param A Eigen expression
* \return Complex dynamic matrix, whose columns are the left singular
* vectors of \a A
*/
template <typename Derived>
cmat svdU(const Eigen::MatrixBase<Derived>& A) {
    const dyn_mat<typename Derived::Scalar>& rA = A.derived();

    // EXCEPTION CHECKS

    // check zero-size
    if (!internal::check_nonzero_size(rA))
        throw exception::ZeroSize("qpp::svdU()");
    // END EXCEPTION CHECKS

    Eigen::JacobiSVD<dyn_mat<typename Derived::Scalar>> sv(
        rA, Eigen::DecompositionOptions::ComputeFullU);

    return sv.matrixU();
}

/**
* \brief Right singular vectors
*
* \param A Eigen expression
* \return Complex dynamic matrix, whose columns are the right singular
* vectors of \a A
*/
template <typename Derived>
cmat svdV(const Eigen::MatrixBase<Derived>& A) {
    const dyn_mat<typename Derived::Scalar>& rA = A.derived();

    // EXCEPTION CHECKS

    // check zero-size
    if (!internal::check_nonzero_size(rA))
        throw exception::ZeroSize("qpp::svdV()");
    // END EXCEPTION CHECKS

    Eigen::JacobiSVD<dyn_mat<typename Derived::Scalar>> sv(
        rA, Eigen::DecompositionOptions::ComputeFullV);

    return sv.matrixV();
}

// Matrix functional calculus

/**
* \brief Functional calculus f(A)
*
* \param A Eigen expression
* \param f Pointer-to-function from complex to complex
* \return \a \f$f(A)\f$
*/
template <typename Derived>
cmat funm(const Eigen::MatrixBase<Derived>& A, cplx (*f)(const cplx&)) {
    const dyn_mat<typename Derived::Scalar>& rA = A.derived();

    // EXCEPTION CHECKS

    // check zero-size
    if (!internal::check_nonzero_size(rA))
        throw exception::ZeroSize("qpp::funm()");

    // check square matrix
    if (!internal::check_square_mat(rA))
        throw exception::MatrixNotSquare("qpp::funm()");
    // END EXCEPTION CHECKS

    Eigen::ComplexEigenSolver<cmat> es(rA.template cast<cplx>());
    cmat evects = es.eigenvectors();
    cmat evals = es.eigenvalues();
    for (idx i = 0; i < static_cast<idx>(evals.rows()); ++i)
        evals(i) = (*f)(evals(i)); // apply f(x) to each eigenvalue

    cmat evalsdiag = evals.asDiagonal();

    return evects * evalsdiag * evects.inverse();
}

/**
* \brief Matrix square root
*
* \param A Eigen expression
* \return Matrix square root of \a A
*/
template <typename Derived>
cmat sqrtm(const Eigen::MatrixBase<Derived>& A) {
    const dyn_mat<typename Derived::Scalar>& rA = A.derived();

    // EXCEPTION CHECKS

    // check zero-size
    if (!internal::check_nonzero_size(rA))
        throw exception::ZeroSize("qpp::sqrtm()");

    // check square matrix
    if (!internal::check_square_mat(rA))
        throw exception::MatrixNotSquare("qpp::sqrtm()");
    // END EXCEPTION CHECKS

    return funm(rA, &std::sqrt);
}

/**
* \brief Matrix absolute value
*
* \param A Eigen expression
* \return Matrix absolute value of \a A
*/
template <typename Derived>
cmat absm(const Eigen::MatrixBase<Derived>& A) {
    const dyn_mat<typename Derived::Scalar>& rA = A.derived();

    // EXCEPTION CHECKS

    // check zero-size
    if (!internal::check_nonzero_size(rA))
        throw exception::ZeroSize("qpp::absm()");

    // check square matrix
    if (!internal::check_square_mat(rA))
        throw exception::MatrixNotSquare("qpp::absm()");
    // END EXCEPTION CHECKS

    return sqrtm(adjoint(rA) * rA);
}

/**
* \brief Matrix exponential
*
* \param A Eigen expression
* \return Matrix exponential of \a A
*/
template <typename Derived>
cmat expm(const Eigen::MatrixBase<Derived>& A) {
    const dyn_mat<typename Derived::Scalar>& rA = A.derived();

    // EXCEPTION CHECKS

    // check zero-size
    if (!internal::check_nonzero_size(rA))
        throw exception::ZeroSize("qpp::expm()");

    // check square matrix
    if (!internal::check_square_mat(rA))
        throw exception::MatrixNotSquare("qpp::expm()");
    // END EXCEPTION CHECKS

    return funm(rA, &std::exp);
}

/**
* \brief Matrix logarithm
*
* \param A Eigen expression
* \return Matrix logarithm of \a A
*/
template <typename Derived>
cmat logm(const Eigen::MatrixBase<Derived>& A) {
    const dyn_mat<typename Derived::Scalar>& rA = A.derived();

    // EXCEPTION CHECKS

    // check zero-size
    if (!internal::check_nonzero_size(rA))
        throw exception::ZeroSize("qpp::logm()");

    // check square matrix
    if (!internal::check_square_mat(rA))
        throw exception::MatrixNotSquare("qpp::logm()");
    // END EXCEPTION CHECKS

    return funm(rA, &std::log);
}

/**
* \brief Matrix sin
*
* \param A Eigen expression
* \return Matrix sine of \a A
*/
template <typename Derived>
cmat sinm(const Eigen::MatrixBase<Derived>& A) {
    const dyn_mat<typename Derived::Scalar>& rA = A.derived();

    // EXCEPTION CHECKS

    // check zero-size
    if (!internal::check_nonzero_size(rA))
        throw exception::ZeroSize("qpp::sinm()");

    // check square matrix
    if (!internal::check_square_mat(rA))
        throw exception::MatrixNotSquare("qpp::sinm()");
    // END EXCEPTION CHECKS

    return funm(rA, &std::sin);
}

/**
* \brief Matrix cos
*
* \param A Eigen expression
* \return Matrix cosine of \a A
*/
template <typename Derived>
cmat cosm(const Eigen::MatrixBase<Derived>& A) {
    const dyn_mat<typename Derived::Scalar>& rA = A.derived();

    // EXCEPTION CHECKS

    // check zero-size
    if (!internal::check_nonzero_size(rA))
        throw exception::ZeroSize("qpp::cosm()");

    // check square matrix
    if (!internal::check_square_mat(rA))
        throw exception::MatrixNotSquare("qpp::cosm()");
    // END EXCEPTION CHECKS

    return funm(rA, &std::cos);
}

/**
* \brief Matrix power
* \see qpp::powm()
*
* Uses the spectral decomposition of \a A to compute the matrix power.
* By convention \f$A^0 = I\f$.
*
* \param A Eigen expression
* \param z Complex number
* \return Matrix power \f$A^z\f$
*/
template <typename Derived>
cmat spectralpowm(const Eigen::MatrixBase<Derived>& A, const cplx z) {
    const dyn_mat<typename Derived::Scalar>& rA = A.derived();

    // EXCEPTION CHECKS

    // check zero-size
    if (!internal::check_nonzero_size(rA))
        throw exception::ZeroSize("qpp::spectralpowm()");

    // check square matrix
    if (!internal::check_square_mat(rA))
        throw exception::MatrixNotSquare("qpp::spectralpowm()");
    // END EXCEPTION CHECKS

    // Define A^0 = Id, for z IDENTICALLY zero
    if (real(z) == 0 && imag(z) == 0)
        return cmat::Identity(rA.rows(), rA.rows());

    Eigen::ComplexEigenSolver<cmat> es(rA.template cast<cplx>());
    cmat evects = es.eigenvectors();
    cmat evals = es.eigenvalues();
    for (idx i = 0; i < static_cast<idx>(evals.rows()); ++i)
        evals(i) = std::pow(evals(i), z);

    cmat evalsdiag = evals.asDiagonal();

    return evects * evalsdiag * evects.inverse();
}

/**
* \brief Fast matrix power based on the SQUARE-AND-MULTIPLY algorithm
* \see qpp::spectralpowm()
*
* Explicitly multiplies the matrix \a A with itself \a n times.
* By convention \f$A^0 = I\f$.
*
* \param A Eigen expression
* \param n Non-negative integer
* \return Matrix power \f$A^n\f$, as a dynamic matrix
* over the same scalar field as \a A
*/
template <typename Derived>
dyn_mat<typename Derived::Scalar> powm(const Eigen::MatrixBase<Derived>& A,
                                       idx n) {
    // EXCEPTION CHECKS

    // check zero-size
    if (!internal::check_nonzero_size(A))
        throw exception::ZeroSize("qpp::powm()");

    // check square matrix
    if (!internal::check_square_mat(A))
        throw exception::MatrixNotSquare("qpp::powm()");
    // END EXCEPTION CHECKS

    // if n = 1, return the matrix unchanged
    if (n == 1)
        return A;

    dyn_mat<typename Derived::Scalar> result =
        dyn_mat<typename Derived::Scalar>::Identity(A.rows(), A.rows());

    // if n = 0, return the identity (as just prepared in result)
    if (n == 0)
        return result;

    dyn_mat<typename Derived::Scalar> cA = A.derived(); // copy

    // fast matrix power
    for (; n > 0; n /= 2) {
        if (n % 2)
            result = (result * cA).eval();
        cA = (cA * cA).eval();
    }

    return result;
}

/**
* \brief Schatten matrix norm
*
* \param A Eigen expression
* \param p Real number, greater or equal to 1,
* use qpp::infty for \f$p = \infty\f$
* \return Schatten-\a p matrix norm of \a A
*/
template <typename Derived>
double schatten(const Eigen::MatrixBase<Derived>& A, double p) {
    const dyn_mat<typename Derived::Scalar>& rA = A.derived();

    // EXCEPTION CHECKS

    // check zero-size
    if (!internal::check_nonzero_size(rA))
        throw exception::ZeroSize("qpp::schatten()");
    if (p < 1)
        throw exception::OutOfRange("qpp::schatten()");
    // END EXCEPTION CHECKS

    if (p == infty) // infinity norm (largest singular value)
        return svals(rA)(0);

    const dyn_col_vect<double> sv = svals(rA);
    double result = 0;
    for (idx i = 0; i < static_cast<idx>(sv.rows()); ++i)
        result += std::pow(sv[i], p);

    return std::pow(result, 1. / p);
}

// other functions

/**
* \brief Functor
*
* \param A Eigen expression
* \param f Pointer-to-function from scalars of \a A to \a OutputScalar
* \return Component-wise \f$f(A)\f$, as a dynamic matrix
* over the \a OutputScalar scalar field
*/
template <typename OutputScalar, typename Derived>
dyn_mat<OutputScalar>
cwise(const Eigen::MatrixBase<Derived>& A,
      OutputScalar (*f)(const typename Derived::Scalar&)) {
    const dyn_mat<typename Derived::Scalar>& rA = A.derived();

    // EXCEPTION CHECKS

    // check zero-size
    if (!internal::check_nonzero_size(rA))
        throw exception::ZeroSize("qpp::cwise()");
    // END EXCEPTION CHECKS

    dyn_mat<OutputScalar> result(rA.rows(), rA.cols());

#ifdef WITH_OPENMP_
#pragma omp parallel for collapse(2)
#endif // WITH_OPENMP_
    // column major order for speed
    for (idx j = 0; j < static_cast<idx>(rA.cols()); ++j)
        for (idx i = 0; i < static_cast<idx>(rA.rows()); ++i)
            result(i, j) = (*f)(rA(i, j));

    return result;
}

// Kronecker product of multiple matrices, preserve return type
// variadic template
/**
* \brief Kronecker product
* \see qpp::kronpow()
*
* Used to stop the recursion for the variadic template version of
* qpp::kron()
*
* \param head Eigen expression
* \return Its argument \a head
*/
template <typename T>
dyn_mat<typename T::Scalar> kron(const T& head) {
    return head;
}

/**
* \brief Kronecker product
* \see qpp::kronpow()
*
* \param head Eigen expression
* \param tail Variadic Eigen expression (zero or more parameters)
* \return Kronecker product of all input parameters,
* evaluated from left to right, as a dynamic matrix
* over the same scalar field as its arguments
*/
template <typename T, typename... Args>
dyn_mat<typename T::Scalar> kron(const T& head, const Args&... tail) {
    return internal::kron2(head, kron(tail...));
}

/**
* \brief Kronecker product
* \see qpp::kronpow()
*
* \param As std::vector of Eigen expressions
* \return Kronecker product of all elements in \a As,
* evaluated from left to right, as a dynamic matrix
* over the same scalar field as its arguments
*/
template <typename Derived>
dyn_mat<typename Derived::Scalar> kron(const std::vector<Derived>& As) {
    // EXCEPTION CHECKS

    if (As.size() == 0)
        throw exception::ZeroSize("qpp::kron()");

    for (auto&& it : As)
        if (!internal::check_nonzero_size(it))
            throw exception::ZeroSize("qpp::kron()");
    // END EXCEPTION CHECKS

    dyn_mat<typename Derived::Scalar> result = As[0].derived();
    for (idx i = 1; i < As.size(); ++i) {
        result = kron(result, As[i]);
    }

    return result;
}

// Kronecker product of a list of matrices, preserve return type
// deduce the template parameters from initializer_list
/**
* \brief Kronecker product
* \see qpp::kronpow()
*
* \param As std::initializer_list of Eigen expressions,
* such as \a {A1, A2, ... ,Ak}
* \return Kronecker product of all elements in \a As,
* evaluated from left to right, as a dynamic matrix
* over the same scalar field as its arguments
*/
template <typename Derived>
dyn_mat<typename Derived::Scalar>
kron(const std::initializer_list<Derived>& As) {
    return kron(std::vector<Derived>(As));
}

/**
* \brief Kronecker power
* \see qpp::kron()
*
* \param A Eigen expression
* \param n Non-negative integer
* \return Kronecker product of \a A with itself \a n times \f$A^{\otimes n}\f$,
* as a dynamic matrix over the same scalar field as \a A
*/
template <typename Derived>
dyn_mat<typename Derived::Scalar> kronpow(const Eigen::MatrixBase<Derived>& A,
                                          idx n) {
    const dyn_mat<typename Derived::Scalar>& rA = A.derived();

    // EXCEPTION CHECKS

    // check zero-size
    if (!internal::check_nonzero_size(rA))
        throw exception::ZeroSize("qpp::kronpow()");

    // check out of range
    if (n == 0)
        throw exception::OutOfRange("qpp::kronpow()");
    // END EXCEPTION CHECKS

    std::vector<dyn_mat<typename Derived::Scalar>> As(n, rA);

    return kron(As);
}

// Direct sum of multiple matrices, preserve return type
// variadic template
/**
* \brief Direct sum
* \see qpp::dirsumpow()
*
* Used to stop the recursion for the variadic template version of
* qpp::dirsum()
*
* \param head Eigen expression
* \return Its argument \a head
*/
template <typename T>
dyn_mat<typename T::Scalar> dirsum(const T& head) {
    return head;
}

/**
* \brief Direct sum
* \see qpp::dirsumpow()
*
* \param head Eigen expression
* \param tail Variadic Eigen expression (zero or more parameters)
* \return Direct sum of all input parameters,
* evaluated from left to right, as a dynamic matrix
* over the same scalar field as its arguments
*/
template <typename T, typename... Args>
dyn_mat<typename T::Scalar> dirsum(const T& head, const Args&... tail) {
    return internal::dirsum2(head, dirsum(tail...));
}

/**
* \brief Direct sum
* \see qpp::dirsumpow()
*
* \param As std::vector of Eigen expressions
* \return Direct sum of all elements in \a As,
* evaluated from left to right, as a dynamic matrix
* over the same scalar field as its arguments
*/
template <typename Derived>
dyn_mat<typename Derived::Scalar> dirsum(const std::vector<Derived>& As) {
    // EXCEPTION CHECKS

    if (As.size() == 0)
        throw exception::ZeroSize("qpp::dirsum()");

    for (auto&& it : As)
        if (!internal::check_nonzero_size(it))
            throw exception::ZeroSize("qpp::dirsum()");
    // END EXCEPTION CHECKS

    idx total_rows = 0, total_cols = 0;
    for (idx i = 0; i < As.size(); ++i) {
        total_rows += static_cast<idx>(As[i].rows());
        total_cols += static_cast<idx>(As[i].cols());
    }
    dyn_mat<typename Derived::Scalar> result =
        dyn_mat<typename Derived::Scalar>::Zero(total_rows, total_cols);

    idx cur_row = 0, cur_col = 0;
    for (idx i = 0; i < As.size(); ++i) {
        result.block(cur_row, cur_col, As[i].rows(), As[i].cols()) = As[i];
        cur_row += static_cast<idx>(As[i].rows());
        cur_col += static_cast<idx>(As[i].cols());
    }

    return result;
}

// Direct sum of a list of matrices, preserve return type
// deduce the template parameters from initializer_list
/**
* \brief Direct sum
* \see qpp::dirsumpow()
*
* \param As std::initializer_list of Eigen expressions,
* such as \a {A1, A2, ... ,Ak}
* \return Direct sum of all elements in \a As,
* evaluated from left to right, as a dynamic matrix
* over the same scalar field as its arguments
*/
template <typename Derived>
dyn_mat<typename Derived::Scalar>
dirsum(const std::initializer_list<Derived>& As) {
    return dirsum(std::vector<Derived>(As));
}

/**
* \brief Direct sum power
* \see qpp::dirsum()
*
* \param A Eigen expression
* \param n Non-negative integer
* \return Direct sum of \a A with itself \a n times \f$A^{\oplus n}\f$,
* as a dynamic matrix over the same scalar field as \a A
*/
template <typename Derived>
dyn_mat<typename Derived::Scalar> dirsumpow(const Eigen::MatrixBase<Derived>& A,
                                            idx n) {
    const dyn_mat<typename Derived::Scalar>& rA = A.derived();

    // EXCEPTION CHECKS

    // check zero-size
    if (!internal::check_nonzero_size(rA))
        throw exception::ZeroSize("qpp::dirsumpow()");

    // check out of range
    if (n == 0)
        throw exception::OutOfRange("qpp::dirsumpow()");
    // END EXCEPTION CHECKS

    std::vector<dyn_mat<typename Derived::Scalar>> As(n, rA);

    return dirsum(As);
}

/**
* \brief Reshape
*
*  Uses column-major order when reshaping (same as MATLAB)
*
* \param A Eigen expression
* \param rows Number of rows of the reshaped matrix
* \param cols Number of columns of the reshaped matrix
* \return Reshaped matrix with \a rows rows and \a cols columns,
* as a dynamic matrix over the same scalar field as \a A
*/
template <typename Derived>
dyn_mat<typename Derived::Scalar> reshape(const Eigen::MatrixBase<Derived>& A,
                                          idx rows, idx cols) {
    const dyn_mat<typename Derived::Scalar>& rA = A.derived();

    idx Arows = static_cast<idx>(rA.rows());
    idx Acols = static_cast<idx>(rA.cols());

    // EXCEPTION CHECKS

    // check zero-size
    if (!internal::check_nonzero_size(rA))
        throw exception::ZeroSize("qpp::reshape()");

    if (Arows * Acols != rows * cols)
        throw exception::DimsMismatchMatrix("qpp::reshape()");
    // END EXCEPTION CHECKS

    return Eigen::Map<dyn_mat<typename Derived::Scalar>>(
        const_cast<typename Derived::Scalar*>(rA.data()), rows, cols);
}

/**
* \brief Commutator
* \see qpp::anticomm()
*
*  Commutator \f$ [A,B] = AB - BA \f$.
*  Both \a A and \a B must be Eigen expressions over the same scalar field.
*
* \param A Eigen expression
* \param B Eigen expression
* \return Commutator \f$AB -BA\f$, as a dynamic matrix
* over the same scalar field as \a A
*/
template <typename Derived1, typename Derived2>
dyn_mat<typename Derived1::Scalar> comm(const Eigen::MatrixBase<Derived1>& A,
                                        const Eigen::MatrixBase<Derived2>& B) {
    const dyn_mat<typename Derived1::Scalar>& rA = A.derived();
    const dyn_mat<typename Derived2::Scalar>& rB = B.derived();

    // EXCEPTION CHECKS

    // check types
    if (!std::is_same<typename Derived1::Scalar,
                      typename Derived2::Scalar>::value)
        throw exception::TypeMismatch("qpp::comm()");

    // check zero-size
    if (!internal::check_nonzero_size(rA) || !internal::check_nonzero_size(rB))
        throw exception::ZeroSize("qpp::comm()");

    // check square matrices
    if (!internal::check_square_mat(rA) || !internal::check_square_mat(rB))
        throw exception::MatrixNotSquare("qpp::comm()");

    // check equal dimensions
    if (rA.rows() != rB.rows())
        throw exception::DimsNotEqual("qpp::comm()");
    // END EXCEPTION CHECKS

    return rA * rB - rB * rA;
}

/**
* \brief Anti-commutator
* \see qpp::comm()
*
*  Anti-commutator \f$ \{A,B\} = AB + BA \f$.
*  Both \a A and \a B must be Eigen expressions over the same scalar field.
*
* \param A Eigen expression
* \param B Eigen expression
* \return Anti-commutator \f$AB +BA\f$, as a dynamic matrix
* over the same scalar field as \a A
*/
template <typename Derived1, typename Derived2>
dyn_mat<typename Derived1::Scalar>
anticomm(const Eigen::MatrixBase<Derived1>& A,
         const Eigen::MatrixBase<Derived2>& B) {
    const dyn_mat<typename Derived1::Scalar>& rA = A.derived();
    const dyn_mat<typename Derived2::Scalar>& rB = B.derived();

    // EXCEPTION CHECKS

    // check types
    if (!std::is_same<typename Derived1::Scalar,
                      typename Derived2::Scalar>::value)
        throw exception::TypeMismatch("qpp::anticomm()");

    // check zero-size
    if (!internal::check_nonzero_size(rA) || !internal::check_nonzero_size(rB))
        throw exception::ZeroSize("qpp::anticomm()");

    // check square matrices
    if (!internal::check_square_mat(rA) || !internal::check_square_mat(rB))
        throw exception::MatrixNotSquare("qpp::anticomm()");

    // check equal dimensions
    if (rA.rows() != rB.rows())
        throw exception::DimsNotEqual("qpp::anticomm()");
    // END EXCEPTION CHECKS

    return rA * rB + rB * rA;
}

/**
* \brief Projector
*
*  Normalized projector onto state vector
*
* \param A Eigen expression
* \return Projector onto the state vector \a A, or the matrix \a Zero
* if \a A has norm zero (i.e. smaller than qpp::eps),
* as a dynamic matrix over the same scalar field as \a A
*/
template <typename Derived>
dyn_mat<typename Derived::Scalar> prj(const Eigen::MatrixBase<Derived>& A) {
    const dyn_mat<typename Derived::Scalar>& rA = A.derived();

    // EXCEPTION CHECKS

    // check zero-size
    if (!internal::check_nonzero_size(rA))
        throw exception::ZeroSize("qpp::prj()");

    // check column vector
    if (!internal::check_cvector(rA))
        throw exception::MatrixNotCvector("qpp::prj()");
    // END EXCEPTION CHECKS

    double normA = norm(rA);
    if (normA > eps)
        return rA * adjoint(rA) / (normA * normA);
    else
        return dyn_mat<typename Derived::Scalar>::Zero(rA.rows(), rA.rows());
}

/**
* \brief Gram-Schmidt orthogonalization
*
* \param As std::vector of Eigen expressions as column vectors
* \return Gram-Schmidt vectors of \a As as columns of a dynamic matrix
* over the same scalar field as its arguments
*/
template <typename Derived>
dyn_mat<typename Derived::Scalar> grams(const std::vector<Derived>& As) {
    // EXCEPTION CHECKS

    // check empty list
    if (!internal::check_nonzero_size(As))
        throw exception::ZeroSize("qpp::grams()");

    for (auto&& it : As)
        if (!internal::check_nonzero_size(it))
            throw exception::ZeroSize("qpp::grams()");

    // check that As[0] is a column vector
    if (!internal::check_cvector(As[0]))
        throw exception::MatrixNotCvector("qpp::grams()");

    // now check that all the rest match the size of the first vector
    for (auto&& it : As)
        if (it.rows() != As[0].rows() || it.cols() != 1)
            throw exception::DimsNotEqual("qpp::grams()");
    // END EXCEPTION CHECKS

    dyn_mat<typename Derived::Scalar> cut =
        dyn_mat<typename Derived::Scalar>::Identity(As[0].rows(), As[0].rows());

    dyn_mat<typename Derived::Scalar> vi =
        dyn_mat<typename Derived::Scalar>::Zero(As[0].rows(), 1);

    std::vector<dyn_mat<typename Derived::Scalar>> outvecs;
    // find the first non-zero vector in the list
    idx pos = 0;
    for (pos = 0; pos < As.size(); ++pos) {
        if (norm(As[pos]) > eps) // add it as the first element
        {
            outvecs.push_back(As[pos]);
            break;
        }
    }

    // start the process
    for (idx i = pos + 1; i < As.size(); ++i) {
        cut -= prj(outvecs[i - 1 - pos]);
        vi = cut * As[i];
        outvecs.push_back(vi);
    }

    dyn_mat<typename Derived::Scalar> result(As[0].rows(), outvecs.size());

    idx cnt = 0;
    for (auto&& it : outvecs) {
        double normA = norm(it);
        if (normA > eps) // we add only the non-zero vectors
        {
            result.col(cnt) = it / normA;
            cnt++;
        }
    }

    return result.block(0, 0, As[0].rows(), cnt);
}

// deduce the template parameters from initializer_list
/**
* \brief Gram-Schmidt orthogonalization
*
* \param As std::initializer_list of Eigen expressions as column vectors
* \return Gram-Schmidt vectors of \a As as columns of a dynamic matrix
* over the same scalar field as its arguments
*/
template <typename Derived>
dyn_mat<typename Derived::Scalar>
grams(const std::initializer_list<Derived>& As) {
    return grams(std::vector<Derived>(As));
}

/**
* \brief Gram-Schmidt orthogonalization
*
* \param A Eigen expression, the input vectors are the columns of \a A
* \return Gram-Schmidt vectors of the columns of \a A,
* as columns of a dynamic matrix over the same scalar field as \a A
*/
template <typename Derived>
dyn_mat<typename Derived::Scalar> grams(const Eigen::MatrixBase<Derived>& A) {
    const dyn_mat<typename Derived::Scalar>& rA = A.derived();

    // EXCEPTION CHECKS

    if (!internal::check_nonzero_size(rA))
        throw exception::ZeroSize("qpp::grams()");
    // END EXCEPTION CHECKS

    std::vector<dyn_mat<typename Derived::Scalar>> input;

    for (idx i = 0; i < static_cast<idx>(rA.cols()); ++i)
        input.push_back(rA.col(i));

    return grams<dyn_mat<typename Derived::Scalar>>(input);
}

/**
* \brief Non-negative integer index to multi-index
* \see qpp::multiidx2n()
*
* Uses standard lexicographical order, i.e. 00...0, 00...1 etc.
*
* \param n Non-negative integer index
* \param dims Dimensions of the multi-partite system
* \return Multi-index of the same size as \a dims
*/
inline std::vector<idx> n2multiidx(idx n, const std::vector<idx>& dims) {
    // EXCEPTION CHECKS

    if (!internal::check_dims(dims))
        throw exception::DimsInvalid("qpp::n2multiidx()");

    if (n >= std::accumulate(std::begin(dims), std::end(dims),
                             static_cast<idx>(1), std::multiplies<idx>()))
        throw exception::OutOfRange("qpp::n2multiidx()");
    // END EXCEPTION CHECKS

    // double the size for matrices reshaped as vectors
    idx result[2 * maxn];
    internal::n2multiidx(n, dims.size(), dims.data(), result);

    return std::vector<idx>(result, result + dims.size());
}

/**
* \brief Multi-index to non-negative integer index
* \see qpp::n2multiidx()
*
* Uses standard lexicographical order, i.e. 00...0, 00...1 etc.
*
* \param midx Multi-index
* \param dims Dimensions of the multi-partite system
* \return Non-negative integer index
*/
inline idx multiidx2n(const std::vector<idx>& midx,
                      const std::vector<idx>& dims) {
    // EXCEPTION CHECKS

    if (!internal::check_dims(dims))
        throw exception::DimsInvalid("qpp::multiidx2n()");

    for (idx i = 0; i < dims.size(); ++i)
        if (midx[i] >= dims[i])
            throw exception::OutOfRange("qpp::multiidx2n()");
    // END EXCEPTION CHECKS

    return internal::multiidx2n(midx.data(), dims.size(), dims.data());
}

/**
* \brief Multi-partite qudit ket
* \see ket template<char... Bits> qpp::operator "" _ket()
*
*
* Constructs the multi-partite qudit ket \f$|\mathrm{mask}\rangle\f$,
* where \a mask is a std::vector of non-negative integers.
* Each element in \a mask has to be smaller than the corresponding element
* in \a dims.
*
* \param mask std::vector of non-negative integers
* \param dims Dimensions of the multi-partite system
* \return Multi-partite qudit state vector, as a complex dynamic column vector
*/
inline ket mket(const std::vector<idx>& mask, const std::vector<idx>& dims) {
    idx N = mask.size();

    idx D = std::accumulate(std::begin(dims), std::end(dims),
                            static_cast<idx>(1), std::multiplies<idx>());

    // EXCEPTION CHECKS

    // check zero size
    if (N == 0)
        throw exception::ZeroSize("qpp::mket()");
    // check valid dims
    if (!internal::check_dims(dims))
        throw exception::DimsInvalid("qpp::mket()");
    // check mask and dims have the same size
    if (mask.size() != dims.size())
        throw exception::SubsysMismatchDims("qpp::mket()");
    // check mask is a valid vector
    for (idx i = 0; i < N; ++i)
        if (mask[i] >= dims[i])
            throw exception::SubsysMismatchDims("qpp::mket()");
    // END EXCEPTION CHECKS

    ket result = ket::Zero(D);
    idx pos = multiidx2n(mask, dims);
    result(pos) = 1;

    return result;
}

/**
* \brief Multi-partite qudit ket
* \see ket template<char... Bits> qpp::operator "" _ket()
*
* Constructs the multi-partite qudit ket \f$|\mathrm{mask}\rangle\f$,
* all subsystem having equal dimension \a d.
* \a mask is a std::vector of non-negative integers, and
* each element in \a mask has to be strictly smaller than \a d.
*
* \param mask std::vector of non-negative integers
* \param d Subsystem dimensions
* \return Multi-partite qudit state vector, as a complex dynamic column vector
*/
inline ket mket(const std::vector<idx>& mask, idx d = 2) {
    idx N = mask.size();
    idx D = static_cast<idx>(std::llround(std::pow(d, N)));

    // EXCEPTION CHECKS

    // check zero size
    if (N == 0)
        throw exception::ZeroSize("qpp::mket()");

    // check valid dims
    if (d == 0)
        throw exception::DimsInvalid("qpp::mket()");

    // check mask is a valid vector
    for (idx i = 0; i < N; ++i)
        if (mask[i] >= d)
            throw exception::SubsysMismatchDims("qpp::mket()");
    // END EXCEPTION CHECKS

    ket result = ket::Zero(D);
    std::vector<idx> dims(N, d);
    idx pos = multiidx2n(mask, dims);
    result(pos) = 1;

    return result;
}

/**
* \brief Projector onto multi-partite qudit ket
* \see cmat template<char... Bits> qpp::operator "" _prj()
*
* Constructs the projector onto the multi-partite qudit ket
* \f$|\mathrm{mask}\rangle\f$,
* where \a mask is a std::vector of non-negative integers.
* Each element in \a mask has to be smaller than the corresponding element
* in \a dims.
*
* \param mask std::vector of non-negative integers
* \param dims Dimensions of the multi-partite system
* \return Projector onto multi-partite qudit state vector,
* as a complex dynamic matrix
*/
inline cmat mprj(const std::vector<idx>& mask, const std::vector<idx>& dims) {
    idx N = mask.size();

    idx D = std::accumulate(std::begin(dims), std::end(dims),
                            static_cast<idx>(1), std::multiplies<idx>());

    // EXCEPTION CHECKS

    // check zero size
    if (N == 0)
        throw exception::ZeroSize("qpp::mprj()");
    // check valid dims
    if (!internal::check_dims(dims))
        throw exception::DimsInvalid("qpp::mprj()");
    // check mask and dims have the same size
    if (mask.size() != dims.size())
        throw exception::SubsysMismatchDims("qpp::mprj()");
    // check mask is a valid vector
    for (idx i = 0; i < N; ++i)
        if (mask[i] >= dims[i])
            throw exception::SubsysMismatchDims("qpp::mprj()");
    // END EXCEPTION CHECKS

    cmat result = cmat::Zero(D, D);
    idx pos = multiidx2n(mask, dims);
    result(pos, pos) = 1;

    return result;
}

/**
* \brief Projector onto multi-partite qudit ket
* \see cmat template<char... Bits> qpp::operator "" _prj()
*
* Constructs the projector onto the multi-partite qudit ket
* \f$|\mathrm{mask}\rangle\f$,
* all subsystem having equal dimension \a d.
* \a mask is a std::vector of non-negative integers, and
* each element in \a mask has to be strictly smaller than \a d.
*
* \param mask std::vector of non-negative integers
* \param d Subsystem dimensions
* \return Projector onto multi-partite qudit state vector,
* as a complex dynamic matrix
*/
inline cmat mprj(const std::vector<idx>& mask, idx d = 2) {
    idx N = mask.size();
    idx D = static_cast<idx>(std::llround(std::pow(d, N)));

    // EXCEPTION CHECKS

    // check zero size
    if (N == 0)
        throw exception::ZeroSize("qpp::mprj()");

    // check valid dims
    if (d == 0)
        throw exception::DimsInvalid("qpp::mprj()");

    // check mask is a valid vector
    for (idx i = 0; i < N; ++i)
        if (mask[i] >= d)
            throw exception::SubsysMismatchDims("qpp::mprj()");
    // END EXCEPTION CHECKS

    cmat result = cmat::Zero(D, D);
    std::vector<idx> dims(N, d);
    idx pos = multiidx2n(mask, dims);
    result(pos, pos) = 1;

    return result;
}

/**
* \brief Computes the absolute values squared of an STL-like range
* of complex numbers

* \param first Iterator to the first element of the range
* \param last  Iterator to the last element of the range
* \return Real vector consisting of the range absolute values squared
*/
template <typename InputIterator>
std::vector<double> abssq(InputIterator first, InputIterator last) {
    std::vector<double> weights(std::distance(first, last));
    std::transform(first, last, std::begin(weights),
                   [](cplx z) -> double { return std::norm(z); });

    return weights;
}

/**
* \brief Computes the absolute values squared of an STL-like container
*
* \param c STL-like container
* \return Real vector consisting of the container's absolute values squared
*/
template <typename Container>
std::vector<double>
abssq(const Container& c,
      typename std::enable_if<is_iterable<Container>::value>::type* = nullptr)
// we need the std::enable_if to SFINAE out Eigen expressions
// that will otherwise match, instead of matching
// the overload below:

// template<typename Derived>
// abssq(const Eigen::MatrixBase<Derived>& A)
{
    return abssq(std::begin(c), std::end(c));
}

/**
* \brief Computes the absolute values squared of an Eigen expression

* \param A Eigen expression
* \return Real vector consisting of the absolute values squared
*/
template <typename Derived>
std::vector<double> abssq(const Eigen::MatrixBase<Derived>& A) {
    const dyn_mat<typename Derived::Scalar>& rA = A.derived();

    // EXCEPTION CHECKS

    // check zero-size
    if (!internal::check_nonzero_size(rA))
        throw exception::ZeroSize("qpp::abssq()");
    // END EXCEPTION CHECKS

    return abssq(rA.data(), rA.data() + rA.size());
}

/**
* \brief Element-wise sum of an STL-like range
*
* \param first Iterator to the first element of the range
* \param last  Iterator to the last element of the range
* \return Element-wise sum of the range,
* as a scalar over the same scalar field as the range
*/
template <typename InputIterator>
typename std::iterator_traits<InputIterator>::value_type
sum(InputIterator first, InputIterator last) {
    using value_type = typename std::iterator_traits<InputIterator>::value_type;

    return std::accumulate(first, last, static_cast<value_type>(0));
}

/**
* \brief Element-wise sum of the elements of an STL-like container
*
* \param c STL-like container
* \return Element-wise sum of the elements of the container,
* as a scalar over the same scalar field as the container
*/
template <typename Container>
typename Container::value_type
sum(const Container& c,
    typename std::enable_if<is_iterable<Container>::value>::type* = nullptr) {
    return sum(std::begin(c), std::end(c));
}

/**
* \brief Element-wise product of an STL-like range
*
* \param first Iterator to the first element of the range
* \param last  Iterator to the last element of the range
* \return Element-wise product of the range,
* as a scalar over the same scalar field as the range
*/
template <typename InputIterator>
typename std::iterator_traits<InputIterator>::value_type
prod(InputIterator first, InputIterator last) {
    using value_type = typename std::iterator_traits<InputIterator>::value_type;

    return std::accumulate(first, last, static_cast<value_type>(1),
                           std::multiplies<value_type>());
}

/**
* \brief Element-wise product of the elements of an STL-like container
*
* \param c STL-like container
* \return Element-wise product of the elements of the container,
* as a scalar over the same scalar field as the container
*/
template <typename Container>
typename Container::value_type
prod(const Container& c,
     typename std::enable_if<is_iterable<Container>::value>::type* = nullptr) {
    return prod(std::begin(c), std::end(c));
}

/**
* \brief Finds the pure state representation of a matrix
* proportional to a projector onto a pure state
*
* \note No purity check is done, the input state \a A must have rank one,
* otherwise the function returns the first non-zero eigenvector of \a A
*
* \param A Eigen expression, assumed to be proportional
* to a projector onto a pure state, i.e. \a A is assumed to have rank one
* \return The unique non-zero eigenvector of \a A (up to a phase),
* as a dynamic column vector over the same scalar field as \a A
*/
template <typename Derived>
dyn_col_vect<typename Derived::Scalar>
rho2pure(const Eigen::MatrixBase<Derived>& A) {
    const dyn_mat<typename Derived::Scalar>& rA = A.derived();

    // EXCEPTION CHECKS
    // check zero-size
    if (!internal::check_nonzero_size(rA))
        throw exception::ZeroSize("qpp::rho2pure()");
    // check square matrix
    if (!internal::check_square_mat(rA))
        throw exception::MatrixNotSquare("qpp::rho2pure()");
    // END EXCEPTION CHECKS

    dyn_col_vect<double> tmp_evals = hevals(rA);
    cmat tmp_evects = hevects(rA);
    dyn_col_vect<typename Derived::Scalar> result =
        dyn_col_vect<typename Derived::Scalar>::Zero(rA.rows());
    // find the non-zero eigenvector
    // there is only one, assuming the state is pure
    for (idx k = 0; k < static_cast<idx>(rA.rows()); ++k) {
        if (std::abs(tmp_evals(k)) > eps) {
            result = tmp_evects.col(k);
            break;
        }
    }

    return result;
}

/**
* \brief Constructs the complement of a subsystem vector
*
* \param subsys Subsystem vector
* \param N Total number of systems
* \return Complement of \a subsys with respect to the set
* \f$\{0, 1, \ldots, N - 1\}\f$
*/
template <typename T>
std::vector<T> complement(std::vector<T> subsys, idx N) {
    // EXCEPTION CHECKS

    if (N < subsys.size())
        throw exception::OutOfRange("qpp::complement()");
    // END EXCEPTION CHECKS

    std::vector<T> all(N);
    std::vector<T> subsys_bar(N - subsys.size());

    std::iota(std::begin(all), std::end(all), 0);
    std::sort(std::begin(subsys), std::end(subsys));
    std::set_difference(std::begin(all), std::end(all), std::begin(subsys),
                        std::end(subsys), std::begin(subsys_bar));

    return subsys_bar;
}

/**
* \brief Computes the 3-dimensional real Bloch vector
* corresponding to the qubit density matrix \a A
* \see qpp::bloch2rho()
*
* \note It is implicitly assumed that the density matrix is Hermitian
*
* \param A Eigen expression
* \return 3-dimensional Bloch vector
*/
template <typename Derived>
std::vector<double> rho2bloch(const Eigen::MatrixBase<Derived>& A) {
    const dyn_mat<typename Derived::Scalar>& rA = A.derived();

    // EXCEPTION CHECKS

    // check qubit matrix
    if (!internal::check_qubit_matrix(rA))
        throw exception::NotQubitMatrix("qpp::rho2bloch()");
    // END EXCEPTION CHECKS

    std::vector<double> result(3);
    cmat X(2, 2), Y(2, 2), Z(2, 2);
    X << 0, 1, 1, 0;
    Y << 0, -1_i, 1_i, 0;
    Z << 1, 0, 0, -1;
    result[0] = std::real(trace(rA * X));
    result[1] = std::real(trace(rA * Y));
    result[2] = std::real(trace(rA * Z));

    return result;
}

/**
* \brief Computes the density matrix corresponding to
* the 3-dimensional real Bloch vector \a r
* \see qpp::rho2bloch()
*
* \param r 3-dimensional real vector
* \return Qubit density matrix
*/
inline cmat bloch2rho(const std::vector<double>& r) {
    // EXCEPTION CHECKS

    // check 3-dimensional vector
    if (r.size() != 3)
        throw exception::CustomException("qpp::bloch2rho",
                                         "r is not a 3-dimensional vector!");
    // END EXCEPTION CHECKS

    cmat X(2, 2), Y(2, 2), Z(2, 2), Id2(2, 2);
    X << 0, 1, 1, 0;
    Y << 0, -1_i, 1_i, 0;
    Z << 1, 0, 0, -1;
    Id2 << 1, 0, 0, 1;

    return (Id2 + r[0] * X + r[1] * Y + r[2] * Z) / 2.;
}

inline namespace literals {
// Idea taken from http://techblog.altplus.co.jp/entry/2017/11/08/130921
/**
* \brief Multi-partite qubit ket user-defined literal
* \see qpp::mket()
*
* Constructs the multi-partite qubit ket \f$|\mathrm{Bits}\rangle\f$
*
* \tparam Bits String of binary numbers representing the qubit ket
* \return Multi-partite qubit ket, as a complex dynamic column vector
*/
template <char... Bits>
ket operator"" _ket() {
    constexpr idx n = sizeof...(Bits);
    constexpr char bits[n + 1] = {Bits..., '\0'};
    qpp::ket q = qpp::ket::Zero(std::pow(2, n));
    idx pos = 0;

    // EXCEPTION CHECKS

    // check valid multi-partite qubit state
    for (idx i = 0; i < n; ++i) {
        if (bits[i] != '0' && bits[i] != '1')
            throw exception::OutOfRange(R"xxx(qpp::operator "" _ket())xxx");
    }
    // END EXCEPTION CHECKS

    pos = std::stoi(bits, nullptr, 2);
    q(pos) = 1;

    return q;
}

/**
* \brief Multi-partite qubit bra user-defined literal
* \see qpp::mket() and qpp::adjoint()
*
* Constructs the multi-partite qubit bra \f$\langle\mathrm{Bits}|\f$
*
* \tparam Bits String of binary numbers representing the qubit bra
* \return Multi-partite qubit bra, as a complex dynamic row vector
*/
template <char... Bits>
bra operator"" _bra() {
    constexpr idx n = sizeof...(Bits);
    constexpr char bits[n + 1] = {Bits..., '\0'};
    qpp::bra q = qpp::ket::Zero(std::pow(2, n));
    idx pos = 0;

    // EXCEPTION CHECKS

    // check valid multi-partite qubit state
    for (idx i = 0; i < n; ++i) {
        if (bits[i] != '0' && bits[i] != '1')
            throw exception::OutOfRange(R"xxx(qpp::operator "" _bra())xxx");
    }
    // END EXCEPTION CHECKS

    pos = std::stoi(bits, nullptr, 2);
    q(pos) = 1;

    return q;
}

/**
* \brief Multi-partite qubit projector user-defined literal
* \see qpp::mprj()
*
* Constructs the multi-partite qubit projector
* \f$|\mathrm{Bits}\rangle\langle\mathrm{Bits}|\f$ (in the computational basis)
*
* \tparam Bits String of binary numbers representing the qubit state
* to project on
* \return Multi-partite qubit projector, as a complex dynamic matrix
*/
template <char... Bits>
cmat operator"" _prj() {
    constexpr idx n = sizeof...(Bits);
    constexpr char bits[n + 1] = {Bits..., '\0'};

    // EXCEPTION CHECKS

    // check valid multi-partite qubit state
    for (idx i = 0; i < n; ++i) {
        if (bits[i] != '0' && bits[i] != '1')
            throw exception::OutOfRange(R"xxx(qpp::operator "" _prj())xxx");
    }
    // END EXCEPTION CHECKS

    return kron(operator""_ket<Bits...>(), operator""_bra<Bits...>());
}
} /* inline namespace literals */

} /* namespace qpp */

#endif /* FUNCTIONS_H_ */
