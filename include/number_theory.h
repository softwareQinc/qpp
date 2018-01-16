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
* \file number_theory.h
* \brief Number theory functions
*/

#ifndef NUMBER_THEORY_H_
#define NUMBER_THEORY_H_

namespace qpp {
/**
* \brief Simple continued fraction expansion
* \see qpp::contfrac2x()
*
* \param x Real number
* \param N Maximum number of terms in the expansion
* \param cut Stop the expansion when the next term is greater than \a cut
* \return Integer vector containing the simple continued fraction expansion
* of \a x. If there are \a M less than \a N terms in the expansion,
* a shorter vector with \a M components is returned.
*/
inline std::vector<int> x2contfrac(double x, idx N, idx cut = 1e5) {
    // EXCEPTION CHECKS

    if (N == 0)
        throw exception::OutOfRange("qpp::x2contfrac()");
    // END EXCEPTION CHECKS

    std::vector<int> result;

    for (idx i = 0; i < N; ++i) {
        if (x > 0) {
            result.push_back(static_cast<int>(std::llround(std::floor(x))));
            x = 1 / (x - std::floor(x));
        } else // x < 0
        {
            result.push_back(static_cast<int>(std::llround(std::ceil(x))));
            x = 1 / (x - std::ceil(x));
        }
        if (!std::isfinite(x) || std::abs(x) > cut)
            return result;
    }

    return result;
}

/**
* \brief Real representation of a simple continued fraction
* \see qpp::x2contfrac()
*
* \note  If \a N is greater than the size of \a cf (by default it is), then all
* terms in \a cf are considered.
*
* \param cf Integer vector containing the simple continued fraction expansion
* \param N Number of terms considered in the continued fraction expansion.
*
* \return Real representation of the simple continued fraction
*/
inline double contfrac2x(const std::vector<int>& cf, idx N = idx(-1)) {
    // EXCEPTION CHECKS

    if (cf.size() == 0)
        throw exception::ZeroSize("qpp::contfrac2x()");

    if (N == 0)
        throw exception::OutOfRange("qpp::contfrac2x()");
    // END EXCEPTION CHECKS

    if (N > cf.size())
        N = cf.size();

    if (N == 1) // degenerate case, integer
        return cf[0];

    double tmp = 1. / cf[N - 1];
    for (idx i = N - 2; i != 0; --i) {
        tmp = 1. / (tmp + cf[i]);
    }

    return cf[0] + tmp;
}

/**
* \brief Greatest common divisor of two integers
* \see qpp::lcm()
*
* \param a Integer
* \param b Integer
* \return Greatest common divisor of \a a and \a b
*/
inline bigint gcd(bigint a, bigint b) {
    // EXCEPTION CHECKS

    if (a == 0 && b == 0)
        throw exception::OutOfRange("qpp::gcd()");
    // END EXCEPTION CHECKS

    if (a == 0 || b == 0)
        return (std::max(std::abs(a), std::abs(b)));

    bigint result = 1;
    while (b) {
        result = b;
        b = a % result;
        a = result;
    }

    return (result > 0) ? result : -result;
}

/**
* \brief Greatest common divisor of a list of integers
* \see qpp::lcm()
*
* \param as List of integers
* \return Greatest common divisor of all numbers in \a as
*/
inline bigint gcd(const std::vector<bigint>& as) {
    // EXCEPTION CHECKS

    if (as.size() == 0)
        throw exception::ZeroSize("qpp::gcd()");
    // END EXCEPTION CHECKS

    bigint result = as[0]; // convention: gcd({a}) = a
    for (idx i = 1; i < as.size(); ++i) {
        result = gcd(result, as[i]);
    }

    return (result > 0) ? result : -result;
}

/**
* \brief Least common multiple of two integers
* \see qpp::gcd()
*
* \param a Integer
* \param b Integer
* \return Least common multiple of \a a and \a b
*/
inline bigint lcm(bigint a, bigint b) {
    // EXCEPTION CHECKS

    if (a == 0 && b == 0)
        throw exception::OutOfRange("qpp::lcm()");
    // END EXCEPTION CHECKS

    bigint result = a * b / gcd(a, b);

    return (result > 0) ? result : -result;
}

/**
* \brief Least common multiple of a list of integers
* \see qpp::gcd()
*
* \param as List of integers
* \return Least common multiple of all numbers in \a as
*/
inline bigint lcm(const std::vector<bigint>& as) {
    // EXCEPTION CHECKS

    if (as.size() == 0)
        throw exception::ZeroSize("qpp::lcm()");

    if (as.size() == 1) // convention: lcm({a}) = a
        return as[0];

    if (std::find(std::begin(as), std::end(as), 0) != std::end(as))
        throw exception::OutOfRange("qpp::lcm()");
    // END EXCEPTION CHECKS

    bigint result = as[0]; // convention: lcm({n}) = a

    for (idx i = 1; i < as.size(); ++i) {
        result = lcm(result, as[i]);
    }

    return (result > 0) ? result : -result;
}

/**
* \brief Inverse permutation
*
* \param perm Permutation
* \return Inverse of the permutation \a perm
*/
inline std::vector<idx> invperm(const std::vector<idx>& perm) {
    // EXCEPTION CHECKS

    if (!internal::check_perm(perm))
        throw exception::PermInvalid("qpp::invperm()");
    // END EXCEPTION CHECKS

    // construct the inverse
    std::vector<idx> result(perm.size());
    for (idx i = 0; i < perm.size(); ++i)
        result[perm[i]] = i;

    return result;
}

/**
* \brief Compose permutations
*
* \param perm Permutation
* \param sigma Permutation
* \return Composition of the permutations \a perm \f$\circ\f$ \a sigma
*  = perm(sigma)
*/
inline std::vector<idx> compperm(const std::vector<idx>& perm,
                                 const std::vector<idx>& sigma) {
    // EXCEPTION CHECKS

    if (!internal::check_perm(perm))
        throw exception::PermInvalid("qpp::compperm()");
    if (!internal::check_perm(sigma))
        throw exception::PermInvalid("qpp::compperm()");
    if (perm.size() != sigma.size())
        throw exception::PermInvalid("qpp::compperm()");
    // END EXCEPTION CHECKS

    // construct the composition perm(sigma)
    std::vector<idx> result(perm.size());
    for (idx i = 0; i < perm.size(); ++i)
        result[i] = perm[sigma[i]];

    return result;
}

/**
* \brief Prime factor decomposition
*
* \note Runs in \f$\mathcal{O}(\sqrt{n})\f$ time complexity
*
* \param a Integer different from 0, 1 or -1
* \return Integer vector containing the factors
*/
inline std::vector<bigint> factors(bigint a) {
    // flip the sign if necessary
    a = std::abs(a);

    // EXCEPTION CHECKS

    if (a == 0 || a == 1)
        throw exception::OutOfRange("qpp::factors()");
    // END EXCEPTION CHECKS

    std::vector<bigint> result;
    bigint d = 2;

    while (a > 1) {
        while (a % d == 0) {
            result.push_back(d);
            a /= d;
        }
        ++d;
        if (d * d > a) // changes the running time from O(a) to O(sqrt(a))
        {
            if (a > 1) {
                result.push_back(a);
            }
            break;
        }
    }

    return result;
}

/**
 * \brief Modular multiplication without overflow
 *
 * Computes \f$ab\f$ \f$\mathrm{ mod }\f$ \f$p\f$ without overflow
 *
 * \param a Integer
 * \param b Integer
 * \param p Positive integer
 * \return \f$ab\f$ \f$\mathrm{ mod }\f$ \f$p\f$ avoiding overflow
 */
inline bigint modmul(bigint a, bigint b, bigint p) {
    using ubigint = unsigned long long int;

    // EXCEPTION CHECKS

    if (p < 1)
        throw exception::OutOfRange("qpp::modmul()");
    // END EXCEPTION CHECKS

    if (a == 0 || b == 0)
        return 0;

    ubigint ua, ub, up;

    bool is_positive = true;
    if (a < 0) {
        ua = -a;
        is_positive = false;
    } else
        ua = a;
    if (b < 0) {
        ub = -b;
        is_positive = false;
    } else
        ub = b;

    if (a < 0 && b < 0)
        is_positive = true;

    up = static_cast<ubigint>(p);
    ua %= up;
    ub %= up;

    // the code below is taken from
    // http://stackoverflow.com/a/18680280/3093378
    ubigint res = 0;
    ubigint temp_b;

    if (ub > ua)
        std::swap(ua, ub);

    /* only needed if un may be >= up */
    if (ub >= up) {
        if (up > std::numeric_limits<ubigint>::max() / 2u)
            ub -= up;
        else
            ub %= up;
    }

    while (ua != 0) {
        if (ua & 1) {
            /* add un to res, modulo p, without overflow */
            /* equiv to if (res + ub >= p), without overflow */
            if (ub >= up - res)
                res -= up;
            res += ub;
        }
        ua >>= 1;

        /* double b, modulo a */
        temp_b = ub;
        if (ub >= up - ub) /* equiv to if (2 * ub >= p), without overflow */
            temp_b -= up;
        ub += temp_b;
    }

    return is_positive ? static_cast<bigint>(res)
                       : static_cast<bigint>(p - res);
}

/**
* \brief Fast integer power modulo \a p based on
* the SQUARE-AND-MULTIPLY algorithm
*
* \note Uses qpp::modmul() that avoids overflows
*
* Computes \f$a^n\f$ \f$\mathrm{ mod }\f$ \f$p\f$
*
* \param a Non-negative integer
* \param n Non-negative integer
* \param p Strictly positive integer
* \return \f$a^n\f$ \f$\mathrm{ mod }\f$ \f$p\f$
*/
inline bigint modpow(bigint a, bigint n, bigint p) {
    // EXCEPTION CHECKS

    if (a < 0 || n < 0 || p < 1)
        throw exception::OutOfRange("qpp::modpow()");

    if (a == 0 && n == 0)
        throw exception::OutOfRange("qpp::modpow()");
    // END EXCEPTION CHECKS

    if (a == 0 && n > 0)
        return 0;

    if (n == 0 && p == 1)
        return 0;

    bigint result = 1;

    for (; n > 0; n /= 2) {
        if (n % 2)
            result = modmul(result, a, p) % p; // MULTIPLY
        a = modmul(a, a, p) % p;               // SQUARE
    }

    return result;
}

/**
 * \brief Extended greatest common divisor of two integers
 * \see qpp::gcd()
 *
 * \param a Integer
 * \param b Integer
 * \return Tuple of: 1. Integer \f$m\f$, 2. Integer \f$n\f$,
 * and 3. Non-negative integer \f$gcd(a, b)\f$ such that
 * \f$ma + nb = gcd(a, b)\f$
 */
inline std::tuple<bigint, bigint, bigint> egcd(bigint a, bigint b) {
    // EXCEPTION CHECKS

    if (a == 0 && b == 0)
        throw exception::OutOfRange("qpp::egcd()");
    // END EXCEPTION CHECKS

    bigint m, n, c, q, r;
    bigint m1 = 0, m2 = 1, n1 = 1, n2 = 0;

    while (b) {
        q = a / b, r = a - q * b;
        m = m2 - q * m1, n = n2 - q * n1;
        a = b, b = r;
        m2 = m1, m1 = m, n2 = n1, n1 = n;
    }
    c = a, m = m2, n = n2;

    // correct the signs
    if (c < 0) {
        m = -m;
        n = -n;
        c = -c;
    }

    return std::make_tuple(m, n, c);
}

/**
 * \brief Modular inverse of \a a mod \a p
 * \see qpp::egcd()
 *
 * \note \a a and \a p must be co-prime
 *
 * \param a Non-negative integer
 * \param p Non-negative integer
 * \return Modular inverse \f$a^{-1}\f$ \f$\textrm{ mod }\f$ \f$p\f$
 */
inline bigint modinv(bigint a, bigint p) {
    // EXCEPTION CHECKS

    if (a <= 0 || p <= 0)
        throw exception::OutOfRange("qpp::modinv()");

    bigint x, y;
    bigint gcd_ap;
    std::tie(x, y, gcd_ap) = egcd(p, a);

    if (gcd_ap != 1)
        throw exception::OutOfRange("qpp::modinv()");
    // END EXCEPTION CHECKS

    return (y > 0) ? y : y + p;
}

/**
 * \brief Primality test based on the Miller-Rabin's algorithm
 *
 * \param p Integer different from 0, 1 or -1
 * \param k Number of iterations. The probability of a
 * false positive is \f$2^{-k}\f$.
 * \return True if the number is (most-likely) prime, false otherwise
 */
inline bool isprime(bigint p, idx k = 80) {
    p = std::abs(p);

    // EXCEPTION CHECKS

    if (p < 2)
        throw exception::OutOfRange("qpp::isprime()");
    // END EXCEPTION CHECKS

    if (p == 2 || p == 3)
        return true;

    //    // perform a Fermat primality test
    bigint x = rand(2, p - 1);
    if (modpow(x, p - 1, p) != 1)
        return false;

    // compute u and r
    bigint u = 0, r = 1;

    // write n − 1 as 2^u * r
    for (bigint i = p - 1; i % 2 == 0; ++u, i /= 2)
        ;
    r = (p - 1) / static_cast<bigint>(std::llround(std::pow(2, u)));

    // repeat k times
    for (idx i = 0; i < k; ++i) {
        // pick a random integer a in the range [2, p − 2]
        bigint a = rand(2, p - 2);

        // set z = a^r mod p
        bigint z = modpow(a, r, p);

        if (z == 1 || z == p - 1)
            continue;

        // repeat u - 1 times
        bool jump = false;
        for (idx j = 0; j < static_cast<idx>(u); ++j) {
            z = (modmul(z, z, p)) % p;
            if (z == 1) {
                // composite
                return false;
            }
            if (z == p - 1) {
                jump = true;
                break;
            }
        }
        if (jump)
            continue;

        return false;
    }

    return true;
}

/**
* \brief Generates a random big prime uniformly distributed in
* the interval [a, b]
*
* \param a Beginning of the interval, belongs to it
* \param b End of the interval, belongs to it
* \param N Maximum number of candidates
* \return Random big integer uniformly distributed in
* the interval [a, b]
*/
// A std::optional<bigint> return type would have been awesome here!
inline bigint randprime(bigint a, bigint b, idx N = 1000) {
    // EXCEPTION CHECKS

    if (a > b)
        throw exception::OutOfRange("qpp::randprime()");
    // END EXCEPTION CHECKS

    idx i = 0;
    for (; i < N; ++i) {
        // select a candidate
        bigint candidate = rand(a, b);
        if (std::abs(candidate) < 2)
            continue;
        if (std::abs(candidate) == 2)
            return candidate;

        // perform a Fermat test
        bigint x = rand(2, candidate - 1);
        if (modpow(x, candidate - 1, candidate) != 1)
            continue; // candidate fails

        // passed the Fermat test, perform a Miller-Rabin test
        if (isprime(candidate))
            return candidate;
    }

    if (i == N)
        throw exception::CustomException("qpp::randprime()",
                                         "Prime not found!");

    return 0; // so we don't get a warning
}

} /* namespace qpp */

#endif /* NUMBER_THEORY_H_ */
