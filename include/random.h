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

#ifndef INCLUDE_RANDOM_H_
#define INCLUDE_RANDOM_H_

// random matrices/states

namespace qpp
{

/**
* \brief Generates a random matrix with entries uniformly
* distributed in the interval [a, b)
*
* If complex, then both real and imaginary parts are uniformly distributed
* in [a, b)
*
* This is the generic version that always throws
* \a qpp::Exception::Type::UNDEFINED_TYPE. It is specialized only for
* \a qpp::dmat and \a qpp::cmat
*/
template<typename Derived>
Derived rand(std::size_t rows, std::size_t cols, double a = 0, double b = 1)
{
    throw Exception("qpp::rand()", Exception::Type::UNDEFINED_TYPE);
}

/**
* \brief Generates a random real matrix with entries uniformly
* distributed in the interval [a, b),
* specialization for double matrices (\a qpp::dmat)
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
inline dmat rand(std::size_t rows, std::size_t cols, double a, double b)
{
    if (rows == 0 || cols == 0)
        throw Exception("qpp::rand()", Exception::Type::ZERO_SIZE);

    return (0.5 * (b - a) * (dmat::Random(rows, cols) + dmat::Ones(rows, cols))
            + a * dmat::Ones(rows, cols));
}

/**
* \brief Generates a random complex matrix with entries (both real and
* imaginary) uniformly distributed in the interval [a, b),
* specialization for complex matrices (\a qpp::cmat)
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
inline cmat rand(std::size_t rows, std::size_t cols, double a, double b)
{
    if (rows == 0 || cols == 0)
        throw Exception("qpp::rand()", Exception::Type::ZERO_SIZE);

    return rand<dmat>(rows, cols, a, b).cast<cplx>() + 1_i * rand<dmat
    >(rows, cols, a, b).cast<cplx>();
}

/**
* \brief Generates a random real number uniformly distributed in
* the interval [a, b)
* \param a Beginning of the interval, belongs to it
* \param b End of the interval, does not belong to it
* \return Random real number (double) uniformly distributed in
* the interval [a, b)
*/
double rand(double a = 0, double b = 1)
{
    std::uniform_real_distribution<> ud(a, b);
    return ud(RandomDevices::get_instance()._rng);
}

/**
* \brief Generates a random integer (int) uniformly distributed in
* the interval [a, b]
* \param a Beginning of the interval, belongs to it
* \param b End of the interval, does not belong to it
* \return Random integer (int) uniformly distributed in
* the interval [a, b]
*/
int randint(int a = std::numeric_limits<int>::min(), int b =
std::numeric_limits<int>::max())
{
    std::uniform_int_distribution<int> ud(a, b);
    return ud(RandomDevices::get_instance()._rng);
}

/**
* \brief Generates a random matrix with entries normally
* distributed in N(mean, sigma)
*
* If complex, then both real and imaginary parts are normally distributed
* in N(mean, sigma)
*
* This is the generic version that always throws
* \a qpp::Exception::Type::UNDEFINED_TYPE. It is specialized only for
* \a qpp::dmat and \a qpp::cmat
*/
template<typename Derived>
Derived randn(std::size_t rows, std::size_t cols, double mean = 0,
        double sigma = 1)
{
    throw Exception("qpp::randn()", Exception::Type::UNDEFINED_TYPE);
}

/**
* \brief Generates a random real matrix with entries normally
* distributed in N(mean, sigma),
* specialization for double matrices (\a qpp::dmat)
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
inline dmat randn(std::size_t rows, std::size_t cols,
        double mean, double sigma)
{
    if (rows == 0 || cols == 0)
        throw Exception("qpp::randn()", Exception::Type::ZERO_SIZE);

    std::normal_distribution<> nd(mean, sigma);

    return dmat::Zero(rows, cols).unaryExpr([&nd](double)
    {
        return nd(RandomDevices::get_instance()._rng);
    });

}

/**
* \brief Generates a random complex matrix with entries (both real and
* imaginary) normally distributed in N(mean, sigma),
* specialization for complex matrices (\a qpp::cmat)
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
inline cmat randn(std::size_t rows, std::size_t cols,
        double mean, double sigma)
{
    if (rows == 0 || cols == 0)
        throw Exception("qpp::randn()", Exception::Type::ZERO_SIZE);

    return randn<dmat>(rows, cols, mean, sigma).cast<cplx>() + 1_i * randn
            <dmat>(rows, cols, mean, sigma).cast<cplx>();
}

/**
* \brief Generates a random real number (double) normally distributed in
* N(mean, sigma)
*
* \param mean Mean
* \param sigma Standard deviation
* \return Random real number normally distributed in N(mean, sigma)
*/
double randn(double mean = 0, double sigma = 1)
{
    std::normal_distribution<> nd(mean, sigma);
    return nd(RandomDevices::get_instance()._rng);
}

/**
* \brief Generates a random unitary matrix
*
* \param D Dimension of the Hilbert space
* \return Random unitary
*/
cmat randU(std::size_t D)
// ~3 times slower than Toby Cubitt's MATLAB corresponding routine,
// because 's QR algorithm is not parallelized
{
    if (D == 0)
        throw Exception("qpp::randU()", Exception::Type::DIMS_INVALID);

    cmat X(D, D);

    X = 1 / std::sqrt(2.) * randn < cmat > (D, D);
    Eigen::HouseholderQR<cmat> qr(X);

    cmat Q = qr.householderQ();
    // phase correction so that the resultant matrix is
    // uniformly distributed according to the Haar measure

    Eigen::VectorXcd phases = (rand < dmat > (D, 1)).cast<cplx>();
    for (std::size_t i = 0; i < static_cast<std::size_t>(phases.rows()); ++i)
        phases(i) = std::exp((cplx) (2 * pi * 1_i * phases(i)));

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
cmat randV(std::size_t Din, std::size_t Dout)
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
std::vector<cmat> randkraus(std::size_t N, std::size_t D)
{
    if (N == 0)
        throw Exception("qpp::randkraus()", Exception::Type::OUT_OF_RANGE);
    if (D == 0)
        throw Exception("qpp::randkraus()", Exception::Type::DIMS_INVALID);

    std::vector<cmat> result(N);
    for (std::size_t i = 0; i < N; ++i)
        result[i] = cmat::Zero(D, D);

    cmat Fk(D, D);
    cmat U = randU(N * D);

#pragma omp parallel for collapse(3)
    for (std::size_t k = 0; k < N; ++k)
        for (std::size_t a = 0; a < D; ++a)
            for (std::size_t b = 0; b < D; ++b)
                result[k](a, b) = U(a * N + k, b * N);

    return result;
}

/**
* \brief Generates a random Hermitian matrix
*
* \param D Dimension of the Hilbert space
* \return Random Hermitian matrix
*/
cmat randH(std::size_t D)
{
    if (D == 0)
        throw Exception("qpp::randH()", Exception::Type::DIMS_INVALID);

    cmat H = 2 * rand < cmat > (D, D) - (1. + 1_i) * cmat::Ones(D, D);

    return H + adjoint(H);
}

/**
* \brief Generates a random normalized ket (pure state vector)
*
* \param D Dimension of the Hilbert space
* \return Random normalized ket
*/
ket randket(std::size_t D)
{
    if (D == 0)
        throw Exception("qpp::randket()", Exception::Type::DIMS_INVALID);
    /* slow
     ket kt = ket::Ones(D);
     ket result = static_cast<ket>(randU(D) * kt);
     return result;
     */
    ket kt = static_cast<ket>(randn < cmat > (D, 1));
    return kt / norm(kt);
}

/**
* \brief Generates a random density matrix
*
* \param D Dimension of the Hilbert space
* \return Random density matrix
*/
cmat randrho(std::size_t D)
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
* Uses Knuth's shuffle method (as implemented by std::shuffle),
* so that all permutations are equally probable
*
* \param n Size of the permutation
* \return Random permutation of size \a n
*/
std::vector<std::size_t> randperm(std::size_t n)
{
    if (n == 0)
        throw Exception("qpp::randperm()", Exception::Type::PERM_INVALID);

    std::vector<std::size_t> result(n);

    // fill in increasing order
    std::iota(std::begin(result), std::end(result), 0);
    // shuffle
    std::shuffle(std::begin(result), std::end(result),
            RandomDevices::get_instance()._rng);

    return result;
}

} /* namespace qpp */

#endif /* INCLUDE_RANDOM_H_ */
