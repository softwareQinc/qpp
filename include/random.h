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
* \file random.h
* \brief Randomness-related functions
*/

#ifndef RANDOM_H_
#define RANDOM_H_

// random matrices/states

namespace qpp
{

/**
* \brief Generates a random real number uniformly distributed in
* the interval [a, b)
* \param a Beginning of the interval, belongs to it
* \param b End of the interval, does not belong to it
* \return Random real number (double) uniformly distributed in
* the interval [a, b)
*/
inline double rand(double a = 0, double b = 1)
{
    std::uniform_real_distribution<> ud(a, b);

#ifndef _NO_THREAD_LOCAL
    return ud(RandomDevices::get_thread_local_instance()._rng);
#else
    return ud(RandomDevices::get_instance()._rng);
    #endif // _NO_THREAD_LOCAL
}

/**
* \brief Generates a random big integer uniformly distributed in
* the interval [a, b]
* \param a Beginning of the interval, belongs to it
* \param b End of the interval, belongs to it
* \return Random big integer uniformly distributed in
* the interval [a, b]
*/
inline bigint rand(bigint a = std::numeric_limits<bigint>::min(),
                   bigint b = std::numeric_limits<bigint>::max())
{
    std::uniform_int_distribution<bigint> uid(a, b);

#ifndef _NO_THREAD_LOCAL
    return uid(RandomDevices::get_thread_local_instance()._rng);
#else
    return uid(RandomDevices::get_instance()._rng);
    #endif // _NO_THREAD_LOCAL
}

/**
* \brief Generates a random non-negative big integer uniformly distributed in
* the interval [a, b]
* \param a Beginning of the interval, belongs to it
* \param b End of the interval, belongs to it
* \return Random non-negative big integer uniformly distributed in
* the interval [a, b]
*/
inline ubigint rand(ubigint a = std::numeric_limits<ubigint>::min(),
                    ubigint b = std::numeric_limits<ubigint>::max())
{
    std::uniform_int_distribution<ubigint> uid(a, b);

#ifndef _NO_THREAD_LOCAL
    return uid(RandomDevices::get_thread_local_instance()._rng);
#else
    return uid(RandomDevices::get_instance()._rng);
    #endif // _NO_THREAD_LOCAL
}

/**
* \brief Generates a random index (idx) uniformly distributed in
* the interval [a, b]
* \param a Beginning of the interval, belongs to it
* \param b End of the interval, belongs to it
* \return Random index (idx) uniformly distributed in
* the interval [a, b]
*/
inline idx randidx(idx a = std::numeric_limits<idx>::min(),
                   idx b = std::numeric_limits<idx>::max())
{
    std::uniform_int_distribution<idx> uid(a, b);

#ifndef _NO_THREAD_LOCAL
    return uid(RandomDevices::get_thread_local_instance()._rng);
#else
    return uid(RandomDevices::get_instance()._rng);
    #endif // _NO_THREAD_LOCAL
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
template<typename Derived>
Derived rand(idx rows, idx cols, double a = 0, double b = 1)
{
    // silence -Wunused-parameter in clang++
    (void) rows;
    (void) cols;
    (void) a;
    (void) b;
    throw Exception("qpp::rand()", Exception::Type::UNDEFINED_TYPE);
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
* auto mat = rand<dmat>(3, 3, -1, 1);
* \endcode
*
* \param rows Number of rows of the random generated matrix
* \param cols Number of columns of the random generated matrix
* \param a Beginning of the interval, belongs to it
* \param b End of the interval, does not belong to it
* \return Random real matrix
*/
template<>
inline dmat rand(idx rows, idx cols, double a, double b)
{
    if (rows == 0 || cols == 0)
        throw Exception("qpp::rand()", Exception::Type::ZERO_SIZE);

    return dmat::Zero(rows, cols).unaryExpr(
            [a, b](double)
            {
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
* auto mat = rand<cmat>(3, 3, -1, 1);
* \endcode
*
* \param rows Number of rows of the random generated matrix
* \param cols Number of columns of the random generated matrix
* \param a Beginning of the interval, belongs to it
* \param b End of the interval, does not belong to it
* \return Random complex matrix
*/
template<>
inline cmat rand(idx rows, idx cols, double a, double b)
{
    if (rows == 0 || cols == 0)
        throw Exception("qpp::rand()", Exception::Type::ZERO_SIZE);

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
template<typename Derived>
Derived randn(idx rows, idx cols, double mean = 0,
              double sigma = 1)
{
    // silence -Wunused-parameter in clang++
    (void) rows;
    (void) cols;
    (void) mean;
    (void) sigma;
    throw Exception("qpp::randn()", Exception::Type::UNDEFINED_TYPE);
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
* auto mat = randn<dmat>(3, 3, 0, 2);
* \endcode
*
* \param rows Number of rows of the random generated matrix
* \param cols Number of columns of the random generated matrix
* \param mean Mean
* \param sigma Standard deviation
* \return Random real matrix
*/
template<>
inline dmat randn(idx rows, idx cols,
                  double mean, double sigma)
{
    if (rows == 0 || cols == 0)
        throw Exception("qpp::randn()", Exception::Type::ZERO_SIZE);

    std::normal_distribution<> nd(mean, sigma);

    return dmat::Zero(rows, cols).unaryExpr(
            [&nd](double)
            {
#ifndef _NO_THREAD_LOCAL
                return nd(RandomDevices::get_thread_local_instance()._rng);
#else
                    return nd(RandomDevices::get_instance()._rng);
                    #endif // _NO_THREAD_LOCAL
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
* auto mat = randn<cmat>(3, 3, 0, 2);
* \endcode
*
* \param rows Number of rows of the random generated matrix
* \param cols Number of columns of the random generated matrix
* \param mean Mean
* \param sigma Standard deviation
* \return Random complex matrix
*/
template<>
inline cmat randn(idx rows, idx cols,
                  double mean, double sigma)
{
    if (rows == 0 || cols == 0)
        throw Exception("qpp::randn()", Exception::Type::ZERO_SIZE);

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
inline double randn(double mean = 0, double sigma = 1)
{
    std::normal_distribution<> nd(mean, sigma);

#ifndef _NO_THREAD_LOCAL
    return nd(RandomDevices::get_thread_local_instance()._rng);
#else
    return nd(RandomDevices::get_instance()._rng);
    #endif // _NO_THREAD_LOCAL
}

/**
* \brief Generates a random unitary matrix
*
* \param D Dimension of the Hilbert space
* \return Random unitary
*/
inline cmat randU(idx D)
// ~3 times slower than Toby Cubitt's MATLAB corresponding routine,
// because Eigen 3 QR algorithm is not parallelized
{
    if (D == 0)
        throw Exception("qpp::randU()", Exception::Type::DIMS_INVALID);

    cmat X(D, D);

    X = 1 / std::sqrt(2.) * randn < cmat > (D, D);
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
inline cmat randV(idx Din, idx Dout)
{
    if (Din == 0 || Dout == 0 || Din > Dout)
        throw Exception("qpp::randV()", Exception::Type::DIMS_INVALID);

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
inline std::vector <cmat> randkraus(idx N, idx D)
{
    if (N == 0)
        throw Exception("qpp::randkraus()", Exception::Type::OUT_OF_RANGE);
    if (D == 0)
        throw Exception("qpp::randkraus()", Exception::Type::DIMS_INVALID);

    std::vector <cmat> result(N);
    for (idx i = 0; i < N; ++i)
        result[i] = cmat::Zero(D, D);

    cmat Fk(D, D);
    cmat U = randU(N * D);

#pragma omp parallel for collapse(3)
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
inline cmat randH(idx D)
{
    if (D == 0)
        throw Exception("qpp::randH()", Exception::Type::DIMS_INVALID);

    cmat H = 2 * rand<cmat>(D, D) - (1. + 1_i) * cmat::Ones(D, D);

    return H + adjoint(H);
}

/**
* \brief Generates a random normalized ket (pure state vector)
*
* \param D Dimension of the Hilbert space
* \return Random normalized ket
*/
inline ket randket(idx D)
{
    if (D == 0)
        throw Exception("qpp::randket()", Exception::Type::DIMS_INVALID);

    /* slow
     ket kt = ket::Ones(D);
     ket result = static_cast<ket>(randU(D) * kt);

     return result;
     */

    ket kt = randn < cmat > (D, 1);

    return kt / norm(kt);
}

/**
* \brief Generates a random density matrix
*
* \param D Dimension of the Hilbert space
* \return Random density matrix
*/
inline cmat randrho(idx D)
{
    if (D == 0)
        throw Exception("qpp::randrho()", Exception::Type::DIMS_INVALID);

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
* \param n Size of the permutation
* \return Random permutation of size \a n
*/
inline std::vector <idx> randperm(idx n)
{
    if (n == 0)
        throw Exception("qpp::randperm()", Exception::Type::PERM_INVALID);

    std::vector <idx> result(n);

    // fill in increasing order
    std::iota(std::begin(result), std::end(result), 0);
    // shuffle
#ifndef _NO_THREAD_LOCAL
    std::shuffle(std::begin(result), std::end(result),
                 RandomDevices::get_thread_local_instance()._rng);
#else
    std::shuffle(std::begin(result), std::end(result),
                 RandomDevices::get_instance()._rng);
    #endif // _NO_THREAD_LOCAL

    return result;
}

} /* namespace qpp */

#endif /* RANDOM_H_ */
