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
* \file random.h
* \brief Randomness-related functions
*/

#ifndef RANDOM_H_
#define RANDOM_H_

namespace qpp {
/**
* \brief Generates a random real number uniformly distributed in
* the interval [a, b)
*
* \param a Beginning of the interval, belongs to it
* \param b End of the interval, does not belong to it
* \return Random real number (double) uniformly distributed in
* the interval [a, b)
*/
inline double rand(double a, double b) {
    // EXCEPTION CHECKS

    if (a >= b)
        throw exception::OutOfRange("qpp::rand()");
    // END EXCEPTION CHECKS

    std::uniform_real_distribution<> ud(a, b);

#ifdef NO_THREAD_LOCAL_
    return ud(RandomDevices::get_instance().get_prng());
#else
    return ud(RandomDevices::get_thread_local_instance().get_prng());
#endif // NO_THREAD_LOCAL_
}

/**
* \brief Generates a random big integer uniformly distributed in
* the interval [a, b]
*
* \note To avoid ambiguity with double qpp::rand(double, double) cast at
* least one of the arguments to qpp::bigint
*
* \param a Beginning of the interval, belongs to it
* \param b End of the interval, belongs to it
* \return Random big integer uniformly distributed in the interval [a, b]
*/
inline bigint rand(bigint a, bigint b) {
    // EXCEPTION CHECKS

    if (a > b)
        throw exception::OutOfRange("qpp::rand()");
    // END EXCEPTION CHECKS

    std::uniform_int_distribution<bigint> uid(a, b);

#ifdef NO_THREAD_LOCAL_
    return uid(RandomDevices::get_instance().get_prng());
#else
    return uid(RandomDevices::get_thread_local_instance().get_prng());
#endif // NO_THREAD_LOCAL_
}

/**
* \brief Generates a random index (idx) uniformly distributed in
* the interval [a, b]
*
* \param a Beginning of the interval, belongs to it
* \param b End of the interval, belongs to it
* \return Random index (idx) uniformly distributed in the interval [a, b]
*/
inline idx randidx(idx a = std::numeric_limits<idx>::min(),
                   idx b = std::numeric_limits<idx>::max()) {
    // EXCEPTION CHECKS

    if (a > b)
        throw exception::OutOfRange("qpp::randidx()");
    // END EXCEPTION CHECKS

    std::uniform_int_distribution<idx> uid(a, b);

#ifdef NO_THREAD_LOCAL_
    return uid(RandomDevices::get_instance().get_prng());
#else
    return uid(RandomDevices::get_thread_local_instance().get_prng());
#endif // NO_THREAD_LOCAL_
}

/**
* \brief Generates a random matrix with entries uniformly
* distributed in the interval [a, b)
*
* If complex, then both real and imaginary parts are uniformly distributed
* in [a, b)
*
* This is the generic version that always throws
* qpp::Exception::Type::UNDEFINED_TYPE. It is specialized only for
* qpp::dmat and qpp::cmat
*/
template <typename Derived>
Derived rand(idx rows, idx cols, double a = 0, double b = 1) {
    // silence -Wunused-parameter in clang++
    (void) rows;
    (void) cols;
    (void) a;
    (void) b;
    throw exception::UndefinedType("qpp::rand()");
}

/**
* \brief Generates a random real matrix with entries uniformly
* distributed in the interval [a, b),
* specialization for double matrices (qpp::dmat)
*
* The template parameter cannot be automatically deduced and
* must be explicitly provided
*
* Example:
* \code
* // generates a 3 x 3 random Eigen::MatrixXd,
* // with entries uniformly distributed in [-1,1)
* dmat mat = rand<dmat>(3, 3, -1, 1);
* \endcode
*
* \param rows Number of rows of the random generated matrix
* \param cols Number of columns of the random generated matrix
* \param a Beginning of the interval, belongs to it
* \param b End of the interval, does not belong to it
* \return Random real matrix
*/
template <>
inline dmat rand(idx rows, idx cols, double a, double b) {
    // EXCEPTION CHECKS

    if (rows == 0 || cols == 0)
        throw exception::ZeroSize("qpp::rand()");
    if (a >= b)
        throw exception::OutOfRange("qpp::rand()");
    // END EXCEPTION CHECKS

    return dmat::Zero(rows, cols).unaryExpr([a, b](double) {
        return rand(a, b);
    });
}

/**
* \brief Generates a random complex matrix with entries (both real and
* imaginary) uniformly distributed in the interval [a, b),
* specialization for complex matrices (qpp::cmat)
*
* The template parameter cannot be automatically deduced and
* must be explicitly provided
*
* Example:
* \code
* // generates a 3 x 3 random Eigen::MatrixXcd,
* // with entries (both real and imaginary) uniformly distributed in [-1,1)
* cmat mat = rand<cmat>(3, 3, -1, 1);
* \endcode
*
* \param rows Number of rows of the random generated matrix
* \param cols Number of columns of the random generated matrix
* \param a Beginning of the interval, belongs to it
* \param b End of the interval, does not belong to it
* \return Random complex matrix
*/
template <>
inline cmat rand(idx rows, idx cols, double a, double b) {
    // EXCEPTION CHECKS

    if (rows == 0 || cols == 0)
        throw exception::ZeroSize("qpp::rand()");
    if (a >= b)
        throw exception::OutOfRange("qpp::rand()");
    // END EXCEPTION CHECKS

    return rand<dmat>(rows, cols, a, b).cast<cplx>() +
           1_i * rand<dmat>(rows, cols, a, b).cast<cplx>();
}

/**
* \brief Generates a random matrix with entries normally
* distributed in N(mean, sigma)
*
* If complex, then both real and imaginary parts are normally distributed
* in N(mean, sigma)
*
* This is the generic version that always throws
* qpp::Exception::Type::UNDEFINED_TYPE. It is specialized only for
* qpp::dmat and qpp::cmat
*/
template <typename Derived>
Derived randn(idx rows, idx cols, double mean = 0, double sigma = 1) {
    // silence -Wunused-parameter in clang++
    (void) rows;
    (void) cols;
    (void) mean;
    (void) sigma;
    throw exception::UndefinedType("qpp::randn()");
}

/**
* \brief Generates a random real matrix with entries normally
* distributed in N(mean, sigma),
* specialization for double matrices (qpp::dmat)
*
* The template parameter cannot be automatically deduced and
* must be explicitly provided
*
* Example:
* \code
* // generates a 3 x 3 random Eigen::MatrixXd,
* // with entries normally distributed in N(0,2)
* dmat mat = randn<dmat>(3, 3, 0, 2);
* \endcode
*
* \param rows Number of rows of the random generated matrix
* \param cols Number of columns of the random generated matrix
* \param mean Mean
* \param sigma Standard deviation
* \return Random real matrix
*/
template <>
inline dmat randn(idx rows, idx cols, double mean, double sigma) {
    // EXCEPTION CHECKS

    if (rows == 0 || cols == 0)
        throw exception::ZeroSize("qpp::randn()");
    // END EXCEPTION CHECKS

    std::normal_distribution<> nd(mean, sigma);

    return dmat::Zero(rows, cols).unaryExpr([&nd](double) {
#ifdef NO_THREAD_LOCAL_
        return nd(RandomDevices::get_instance().get_prng());
#else
        return nd(RandomDevices::get_thread_local_instance().get_prng());
#endif // NO_THREAD_LOCAL_
    });
}

/**
* \brief Generates a random complex matrix with entries (both real and
* imaginary) normally distributed in N(mean, sigma),
* specialization for complex matrices (qpp::cmat)
*
* The template parameter cannot be automatically deduced and
* must be explicitly provided
*
* Example:
* \code
* // generates a 3 x 3 random Eigen::MatrixXcd,
* // with entries (both real and imaginary) normally distributed in N(0,2)
* cmat mat = randn<cmat>(3, 3, 0, 2);
* \endcode
*
* \param rows Number of rows of the random generated matrix
* \param cols Number of columns of the random generated matrix
* \param mean Mean
* \param sigma Standard deviation
* \return Random complex matrix
*/
template <>
inline cmat randn(idx rows, idx cols, double mean, double sigma) {
    // EXCEPTION CHECKS

    if (rows == 0 || cols == 0)
        throw exception::ZeroSize("qpp::randn()");
    // END EXCEPTION CHECKS

    return randn<dmat>(rows, cols, mean, sigma).cast<cplx>() +
           1_i * randn<dmat>(rows, cols, mean, sigma).cast<cplx>();
}

/**
* \brief Generates a random real number (double) normally distributed in
* N(mean, sigma)
*
* \param mean Mean
* \param sigma Standard deviation
* \return Random real number normally distributed in N(mean, sigma)
*/
inline double randn(double mean = 0, double sigma = 1) {
    std::normal_distribution<> nd(mean, sigma);

#ifdef NO_THREAD_LOCAL_
    return nd(RandomDevices::get_instance().get_prng());
#else
    return nd(RandomDevices::get_thread_local_instance().get_prng());
#endif // NO_THREAD_LOCAL_
}

/**
* \brief Generates a random unitary matrix
*
* \param D Dimension of the Hilbert space
* \return Random unitary
*/
inline cmat randU(idx D = 2)
// ~3 times slower than Toby Cubitt's MATLAB corresponding routine,
// because Eigen 3 QR algorithm is not parallelized
{
    // EXCEPTION CHECKS

    if (D == 0)
        throw exception::DimsInvalid("qpp::randU()");
    // END EXCEPTION CHECKS

    cmat X = 1 / std::sqrt(2.) * randn<cmat>(D, D);
    Eigen::HouseholderQR<cmat> qr(X);

    cmat Q = qr.householderQ();
    // phase correction so that the resultant matrix is
    // uniformly distributed according to the Haar measure

    Eigen::VectorXcd phases = (rand<dmat>(D, 1)).cast<cplx>();
    for (idx i = 0; i < static_cast<idx>(phases.rows()); ++i)
        phases(i) = std::exp(2 * pi * 1_i * phases(i));

    Q = Q * phases.asDiagonal();

    return Q;
}

/**
* \brief Generates a random isometry matrix
*
* \param Din Size of the input Hilbert space
* \param Dout Size of the output Hilbert space
* \return Random isometry matrix
*/
inline cmat randV(idx Din, idx Dout) {
    // EXCEPTION CHECKS

    if (Din == 0 || Dout == 0 || Din > Dout)
        throw exception::DimsInvalid("qpp::randV()");
    // END EXCEPTION CHECKS

    return randU(Dout).block(0, 0, Dout, Din);
}

/**
* \brief Generates a set of random Kraus operators
*
* \note The set of Kraus operators satisfy the closure condition
* \f$ \sum_i K_i^\dagger K_i = I\f$
*
* \param N Number of Kraus operators
* \param D Dimension of the Hilbert space
* \return Set of \a N Kraus operators satisfying the closure condition
*/
inline std::vector<cmat> randkraus(idx N, idx D = 2) {
    // EXCEPTION CHECKS

    if (N == 0)
        throw exception::OutOfRange("qpp::randkraus()");
    if (D == 0)
        throw exception::DimsInvalid("qpp::randkraus()");
    // END EXCEPTION CHECKS

    std::vector<cmat> result(N);
    for (idx i = 0; i < N; ++i)
        result[i] = cmat::Zero(D, D);

    cmat Fk(D, D);
    cmat U = randU(N * D);

#ifdef WITH_OPENMP_
#pragma omp parallel for collapse(3)
#endif // WITH_OPENMP_
    for (idx k = 0; k < N; ++k)
        for (idx a = 0; a < D; ++a)
            for (idx b = 0; b < D; ++b)
                result[k](a, b) = U(a * N + k, b * N);

    return result;
}

/**
* \brief Generates a random Hermitian matrix
*
* \param D Dimension of the Hilbert space
* \return Random Hermitian matrix
*/
inline cmat randH(idx D = 2) {
    // EXCEPTION CHECKS

    if (D == 0)
        throw exception::DimsInvalid("qpp::randH()");
    // END EXCEPTION CHECKS

    cmat H = 2 * rand<cmat>(D, D) - (1. + 1_i) * cmat::Ones(D, D);

    return H + adjoint(H);
}

/**
* \brief Generates a random normalized ket (pure state vector)
*
* \param D Dimension of the Hilbert space
* \return Random normalized ket
*/
inline ket randket(idx D = 2) {
    // EXCEPTION CHECKS

    if (D == 0)
        throw exception::DimsInvalid("qpp::randket()");
    // END EXCEPTION CHECKS

    /* slow
     ket kt = ket::Ones(D);
     ket result = static_cast<ket>(randU(D) * kt);

     return result;
     */

    ket kt = randn<cmat>(D, 1);

    return kt / norm(kt);
}

/**
* \brief Generates a random density matrix
*
* \param D Dimension of the Hilbert space
* \return Random density matrix
*/
inline cmat randrho(idx D = 2) {
    // EXCEPTION CHECKS

    if (D == 0)
        throw exception::DimsInvalid("qpp::randrho()");
    // END EXCEPTION CHECKS

    cmat result = 10 * randH(D);
    result = result * adjoint(result);

    return result / trace(result);
}

/**
* \brief Generates a random uniformly distributed permutation
*
* Uses Knuth shuffle method (as implemented by std::shuffle),
* so that all permutations are equally probable
*
* \param N Size of the permutation
* \return Random permutation of size \a N
*/
inline std::vector<idx> randperm(idx N) {
    // EXCEPTION CHECKS

    if (N == 0)
        throw exception::PermInvalid("qpp::randperm()");
    // END EXCEPTION CHECKS

    std::vector<idx> result(N);

    // fill in increasing order
    std::iota(std::begin(result), std::end(result), 0);
// shuffle
#ifdef NO_THREAD_LOCAL_
    std::shuffle(std::begin(result), std::end(result),
                 RandomDevices::get_instance().get_prng());
#else
    std::shuffle(std::begin(result), std::end(result),
                 RandomDevices::get_thread_local_instance().get_prng());

#endif // NO_THREAD_LOCAL_

    return result;
}

/**
* \brief Generates a random probability vector uniformly distributed over the
* probability simplex
*
* \param N Size of the probability vector
* \return Random probability vector
*/
inline std::vector<double> randprob(idx N) {
    // EXCEPTION CHECKS

    if (N == 0)
        throw exception::OutOfRange("qpp::randprob()");
    // END EXCEPTION CHECKS

    std::vector<double> result(N);

    // generate
    std::exponential_distribution<> ed(1);
    for (idx i = 0; i < N; ++i) {
#ifdef NO_THREAD_LOCAL_
        result[i] = ed(qpp::RandomDevices::get_instance().get_prng());
#else
        result[i] =
            ed(qpp::RandomDevices::get_thread_local_instance().get_prng());
#endif // NO_THREAD_LOCAL_
    }

    // normalize
    double sumprob = sum(result);
    for (idx i = 0; i < N; ++i)
        result[i] /= sumprob;

    return result;
}

} /* namespace qpp */

#endif /* RANDOM_H_ */
