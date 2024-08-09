/*
 * This file is part of Quantum++.
 *
 * Copyright (c) 2017 - 2024 softwareQ Inc. All rights reserved.
 *
 * MIT License
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
 * \file qpp/random.hpp
 * \brief Randomness-related functions
 */

#ifndef QPP_RANDOM_HPP_
#define QPP_RANDOM_HPP_

#include <type_traits>

#include "qpp/constants.hpp"
#include "qpp/types.hpp"

#include "qpp/classes/exception.hpp"
#include "qpp/classes/random_devices.hpp"

namespace qpp {
/**
 * \brief Generates a random number in the interval [a, b] for integer types,
 * and in the interval [a, b) for floating-point types
 *
 * \tparam T Arithmetic type
 * \param a Beginning of the interval
 * \param b End of the interval
 * \return Random number uniformly distributed in the interval [a, b)/[a, b]
 */
template <typename T,
          typename std::enable_if_t<std::is_arithmetic_v<T>>* = nullptr>
T rand(T a, T b) {
    static_assert(std::is_integral_v<T> || std::is_floating_point_v<T>,
                  "Not defined for this type");

    auto& gen = RandomDevices::get_instance().get_prng();
    if constexpr (std::is_integral_v<T>) {
        if (a > b) {
            throw exception::OutOfRange("qpp::rand()", "a/b");
        }
        std::uniform_int_distribution<T> uid(a, b);
        return uid(gen);
    } else if constexpr (std::is_floating_point_v<T>) {
        if (a >= b) {
            throw exception::OutOfRange("qpp::rand()", "a/b");
        }
        std::uniform_real_distribution<realT> ud(a, b);
        return ud(gen);
    }
}

/**
 * \brief Generates a random real or complex matrix with entries uniformly
 * distributed in the interval [a, b)
 *
 * \note The template parameter cannot be automatically deduced and must be
 * explicitly provided. It is only specialized for qpp::rmat and qpp::cmat.
 *
 * \note If complex, then both real and imaginary parts are uniformly
 * distributed in [a, b)
 *
 * Example:
 * \code
 * // generates a 3 x 3 random Eigen::MatrixX over qpp::realT,
 * // with entries uniformly distributed in [-1,1)
 * rmat mat = rand<rmat>(3, 3, -1, 1);
 *
 * // generates a 3 x 3 random Eigen::MatrixX over qpp::cplx,
 * // with entries (both real and imaginary) uniformly distributed in [-1,1)
 * cmat mat = rand<cmat>(3, 3, -1, 1);
 * \endcode
 *
 * \tparam Derived Matrix type, must be either qpp::rmat or qpp::cmat
 * \param rows Number of rows of the randomly generated matrix
 * \param cols Number of columns of the randomly generated matrix
 * \param a Beginning of the interval, belongs to it
 * \param b End of the interval, does not belong to it
 * \return Random real (qpp::rmat specialization) or complex
 * (qpp::cmat specialization) matrix
 */
template <typename Derived,
          typename std::enable_if_t<!std::is_arithmetic_v<Derived>>* = nullptr>
Derived rand([[maybe_unused]] idx rows, [[maybe_unused]] idx cols,
             [[maybe_unused]] realT a = 0, [[maybe_unused]] realT b = 1) {
    throw exception::UndefinedType("qpp::rand()");
}

/// \cond DO_NOT_DOCUMENT
template <>
inline rmat rand(idx rows, idx cols, realT a, realT b) {
    // EXCEPTION CHECKS
    if (rows == 0 || cols == 0) {
        throw exception::ZeroSize("qpp::rand()", "rows/cols");
    }
    if (a >= b) {
        throw exception::OutOfRange("qpp::rand()", "a/b");
    }
    // END EXCEPTION CHECKS

    return rmat::Zero(rows, cols).unaryExpr([a, b](realT) {
        return rand(a, b);
    });
}

template <>
inline cmat rand(idx rows, idx cols, realT a, realT b) {
    // EXCEPTION CHECKS
    if (rows == 0 || cols == 0) {
        throw exception::ZeroSize("qpp::rand()", "rows/cols");
    }
    if (a >= b) {
        throw exception::OutOfRange("qpp::rand()", "a/b");
    }
    // END EXCEPTION CHECKS

    return rand<rmat>(rows, cols, a, b).cast<cplx>() +
           1_i * rand<rmat>(rows, cols, a, b).cast<cplx>();
}
/// \endcond

/**
 * \brief Generates a random real number normally distributed in
 * N(mean, sigma)
 *
 * \tparam T Arithmetic type
 * \param mean Mean
 * \param sigma Standard deviation
 * \return Random real number normally distributed in N(mean, sigma)
 */
template <typename T,
          typename std::enable_if_t<std::is_arithmetic_v<T>>* = nullptr>
T randn(T mean = 0, T sigma = 1) {
    std::normal_distribution<T> nd(mean, sigma);
    auto& gen = RandomDevices::get_instance().get_prng();

    return nd(gen);
}

/**
 * \brief Generates a random real or complex matrix with entries normally
 * distributed in N(mean, sigma).
 *
 * \note The template parameter cannot be automatically deduced and must be
 * explicitly provided. It is only specialized for qpp::rmat and qpp::cmat.
 *
 * \note If complex, then both real and imaginary parts are normally distributed
 * in N(mean, sigma)
 *
 * Example:
 * \code
 * // generates a 3 x 3 random Eigen::MatrixX over qpp::realT,
 * // with entries normally distributed in N(0,2)
 * rmat mat = randn<rmat>(3, 3, 0, 2);
 *
 * // generates a 3 x 3 random Eigen::MatrixX over qpp::cplx,
 * // with entries (both real and imaginary) normally distributed in N(0,2)
 * cmat mat = randn<cmat>(3, 3, 0, 2);
 * \endcode
 *
 * \tparam Derived Matrix type, must be either qpp::rmat or qpp::cmat
 * \param rows Number of rows of the randomly generated matrix
 * \param cols Number of columns of the randomly generated matrix
 * \param mean Mean
 * \param sigma Standard deviation
 * \return Random real (qpp::rmat specialization) or complex
 * (qpp::cmat specialization) matrix
 */
template <typename Derived,
          typename std::enable_if_t<!std::is_arithmetic_v<Derived>>* = nullptr>
Derived randn([[maybe_unused]] idx rows, [[maybe_unused]] idx cols,
              [[maybe_unused]] realT mean = 0,
              [[maybe_unused]] realT sigma = 1) {
    throw exception::UndefinedType("qpp::randn()");
}

/// \cond DO_NOT_DOCUMENT
template <>
inline rmat randn(idx rows, idx cols, realT mean, realT sigma) {
    // EXCEPTION CHECKS
    if (rows == 0 || cols == 0) {
        throw exception::ZeroSize("qpp::randn()", "rows/cols");
    }
    // END EXCEPTION CHECKS

    std::normal_distribution<realT> nd(mean, sigma);
    auto& gen = RandomDevices::get_instance().get_prng();

    return rmat::Zero(rows, cols).unaryExpr([&nd, &gen](realT) {
        return nd(gen);
    });
}

template <>
inline cmat randn(idx rows, idx cols, realT mean, realT sigma) {
    // EXCEPTION CHECKS
    if (rows == 0 || cols == 0) {
        throw exception::ZeroSize("qpp::randn()", "rows/cols");
    }
    // END EXCEPTION CHECKS

    return randn<rmat>(rows, cols, mean, sigma).cast<cplx>() +
           1_i * randn<rmat>(rows, cols, mean, sigma).cast<cplx>();
}
/// \endcond

/**
 * \brief Generates a random index (idx) uniformly distributed in the interval
 * [a, b]
 *
 * \param a Beginning of the interval, belongs to it
 * \param b End of the interval, belongs to it
 * \return Random index (idx) uniformly distributed in the interval [a, b]
 */
inline idx randidx(idx a = std::numeric_limits<idx>::min(),
                   idx b = std::numeric_limits<idx>::max()) {
    // EXCEPTION CHECKS
    if (a > b) {
        throw exception::OutOfRange("qpp::randidx()", "a/b");
    }
    // END EXCEPTION CHECKS

    std::uniform_int_distribution<idx> uid(a, b);
    auto& gen = RandomDevices::get_instance().get_prng();

    return uid(gen);
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
    if (D == 0) {
        throw exception::DimsInvalid("qpp::randU()", "D");
    }
    // END EXCEPTION CHECKS

    cmat X = 1 / std::sqrt(2.) * randn<cmat>(D, D);
    Eigen::HouseholderQR<cmat> qr(X);

    cmat Q = qr.householderQ();
    // phase correction so that the resultant matrix is
    // uniformly distributed according to the Haar measure

    dyn_col_vect<cplx> phases = (rand<rmat>(D, 1)).cast<cplx>();
    for (idx i = 0; i < static_cast<idx>(phases.rows()); ++i) {
        phases(i) =
            std::exp(static_cast<cplx::value_type>(2 * pi) * 1_i * phases(i));
    }

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
    if (Din == 0 || Dout == 0 || Din > Dout) {
        throw exception::DimsInvalid("qpp::randV()", "Din/Dout");
    }
    // END EXCEPTION CHECKS

    return randU(Dout).block(0, 0, Dout, Din);
}

/**
 * \brief Generates a set of random Kraus operators from an input space of
 * dimension \a Din to an output space of dimension \a Dout
 *
 * \note The set of Kraus operators satisfies the closure condition
 * \f$\sum_i K_i^\dagger K_i = I\f$. The Kraus operators can have their range
 * different from their domain (i.e., they can be rectangular matrices). The
 * dimensions \a Din, \a Dout and the number of Kraus operators \a N must be
 * chosen so that \f$D_{out}N/D_{in}\f$ is a positive integer.
 *
 * \param N Number of Kraus operators
 * \param Din Dimension of the input Hilbert space
 * \param Dout Dimension of the output Hilbert space
 * \return Set of \a N Kraus operators satisfying the closure condition
 */
[[qpp::parallel]] inline std::vector<cmat> randkraus(idx N, idx Din, idx Dout) {
    // EXCEPTION CHECKS
    if (N == 0) {
        throw exception::OutOfRange("qpp::randkraus()", "N");
    }
    if (Din == 0 || Dout == 0) {
        throw exception::DimsInvalid("qpp::randkraus()", "Din/Dout");
    }
    idx Din_env = (N * Dout) / Din; // dimension of the input environment
    if (Din_env * Din != Dout * N) {
        throw exception::DimsInvalid("qpp::randkraus()", "Din/Dout");
    }
    // END EXCEPTION CHECKS

    std::vector<cmat> result(N);
    for (idx i = 0; i < N; ++i) {
        result[i] = cmat::Zero(Dout, Din);
    }

    cmat Fk(Dout, Din);
    cmat U = randU(N * Dout);

#ifdef QPP_OPENMP
// NOLINTNEXTLINE
#pragma omp parallel for collapse(3)
#endif // QPP_OPENMP
    for (idx k = 0; k < N; ++k) {
        for (idx a = 0; a < Dout; ++a) {
            for (idx b = 0; b < Din; ++b) {
                result[k](a, b) = U(a * N + k, b * Din_env);
            }
        }
    }

    return result;
}

/**
 * \brief Generates a set of random Kraus operators from an input space of
 * dimension \a Din to an output space of dimension \a Dout
 *
 * \note The set of Kraus operators satisfies the closure condition
 * \f$\sum_i K_i^\dagger K_i = I\f$. The Kraus operators are assumed to have
 * their range equal to their domain (i.e., they are square matrices).
 *
 * \param N Number of Kraus operators
 * \param D Dimension of the Kraus input and output Hilbert space
 * \return Set of \a N Kraus operators satisfying the closure condition
 */
inline std::vector<cmat> randkraus(idx N, idx D = 2) {
    return randkraus(N, D, D);
}

/**
 * \brief Generates a random Hermitian matrix
 *
 * \param D Dimension of the Hilbert space
 * \return Random Hermitian matrix
 */
inline cmat randH(idx D = 2) {
    // EXCEPTION CHECKS
    if (D == 0) {
        throw exception::DimsInvalid("qpp::randH()", "D");
    }
    // END EXCEPTION CHECKS

    cmat H = 2 * rand<cmat>(D, D) -
             (static_cast<realT>(1.) + 1_i) * cmat::Ones(D, D);

    return H + H.adjoint();
}

/**
 * \brief Generates a random normalized ket (pure state vector)
 *
 * \param D Dimension of the Hilbert space
 * \return Random normalized ket
 */
inline ket randket(idx D = 2) {
    // EXCEPTION CHECKS
    if (D == 0) {
        throw exception::DimsInvalid("qpp::randket()", "D");
    }
    // END EXCEPTION CHECKS

    /* slow
     ket kt = ket::Ones(D);
     ket result = static_cast<ket>(randU(D) * kt);

     return result;
     */

    ket kt = randn<cmat>(D, 1);

    return kt / kt.norm();
}

/**
 * \brief Generates a random density matrix
 *
 * \param D Dimension of the Hilbert space
 * \return Random density matrix
 */
inline cmat randrho(idx D = 2) {
    // EXCEPTION CHECKS
    if (D == 0) {
        throw exception::DimsInvalid("qpp::randrho()", "D");
    }
    // END EXCEPTION CHECKS

    cmat result = 10 * randH(D);
    result = result * result.adjoint();

    return result / result.trace();
}

/**
 * \brief Generates a random uniformly distributed permutation
 *
 * Uses Knuth shuffle method (as implemented by std::shuffle), so that all
 * permutations are equally probable
 *
 * \param N Size of the permutation
 * \return Random permutation of size \a N
 */
inline std::vector<idx> randperm(idx N) {
    // EXCEPTION CHECKS
    if (N == 0) {
        throw exception::ZeroSize("qpp::randperm()", "N");
    }
    // END EXCEPTION CHECKS

    std::vector<idx> result(N);

    // fill in increasing order
    std::iota(result.begin(), result.end(), 0);

    // shuffle
    auto& gen = RandomDevices::get_instance().get_prng();
    std::shuffle(result.begin(), result.end(), gen);

    return result;
}

/**
 * \brief Generates a random probability vector uniformly distributed over the
 * probability simplex
 *
 * \param N Size of the probability vector
 * \return Random probability vector
 */
inline std::vector<realT> randprob(idx N) {
    // EXCEPTION CHECKS
    if (N == 0) {
        throw exception::ZeroSize("qpp::randprob()");
    }
    // END EXCEPTION CHECKS

    std::vector<realT> result(N);

    // generate
    std::exponential_distribution<realT> ed(1);
    auto& gen = RandomDevices::get_instance().get_prng();
    for (idx i = 0; i < N; ++i) {
        result[i] = ed(gen);
    }

    // normalize
    realT sumprob = std::accumulate(result.begin(), result.end(), 0.0);
    for (idx i = 0; i < N; ++i) {
        result[i] /= sumprob;
    }

    return result;
}

/**
 *\brief Generates a random boolean drawn from a Bernoulli-\f$p\f$ distribution
 *\note Outputs always false for \a p == 0, and always true for \a p == 1
 *
 * \param p Probability bias (0.5 by default)
 * \return Boolean drawn from a Bernoulli-\f$p\f$ distribution
 */
inline bool bernoulli(realT p = 0.5) {
    std::bernoulli_distribution bd(p);
    auto& gen = RandomDevices::get_instance().get_prng();

    return bd(gen);
}

} /* namespace qpp */

#endif /* QPP_RANDOM_HPP_ */
