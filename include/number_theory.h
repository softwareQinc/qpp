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
* \file number_theory.h
* \brief Number theory functions
*/

#ifndef NUMBER_THEORY_H_
#define NUMBER_THEORY_H_

// number theory functions

namespace qpp
{
/**
* \brief Simple continued fraction expansion
* \see qpp::contfrac2x()
*
* \param x Real number
* \param n Number of terms in the expansion
* \param cut Stop the expansion when the next term is greater than \a cut
* \return Integer vector containing the simple continued fraction expansion
* of \a x. If there are \a m less than \a n terms in the expansion,
* a shorter vector with \a m components is returned.
*/
inline std::vector<int> x2contfrac(double x, idx n, idx cut = 1e5)
{
    // EXCEPTION CHECKS

    if (n == 0)
        throw Exception("qpp::x2contfrac()", Exception::Type::OUT_OF_RANGE);
    // END EXCEPTION CHECKS

    std::vector<int> result;

    for (idx i = 0; i < n; ++i)
    {
        result.push_back(std::llround(std::floor(x)));
        x = 1 / (x - std::floor(x));
        if (!std::isfinite(x) || x > cut)
            return result;
    }

    return result;
}

/**
* \brief Real representation of a simple continued fraction
* \see qpp::x2contfrac()
*
* \param cf Integer vector containing the simple continued fraction expansion
* \param n Number of terms considered in the continued fraction expansion.
* If \a n is greater than the size of \a cf,then all terms in \a cf
* are considered.
* \return Real representation of the simple continued fraction
*/
inline double contfrac2x(const std::vector<int>& cf, idx n)
{
    // EXCEPTION CHECKS

    if (cf.size() == 0)
        throw Exception("qpp::contfrac2x()", Exception::Type::ZERO_SIZE);

    if (n == 0)
        throw Exception("qpp::contfrac2x()", Exception::Type::OUT_OF_RANGE);
    // END EXCEPTION CHECKS

    if (n > cf.size())
        n = cf.size();

    if (n == 1) // degenerate case, integer
        return cf[0];

    double tmp = 1. / cf[n - 1];
    for (idx i = n - 2; i != 0; --i)
    {
        tmp = 1. / (tmp + cf[i]);
    }

    return cf[0] + tmp;
}

/**
* \brief Real representation of a simple continued fraction
* \see qpp::x2contfrac()
*
* \param cf Integer vector containing the simple continued fraction expansion
* \return Real representation of the simple continued fraction
*/
inline double contfrac2x(const std::vector<int>& cf)
{
    // EXCEPTION CHECKS

    if (cf.size() == 0)
        throw Exception("qpp::contfrac2x()", Exception::Type::ZERO_SIZE);
    // END EXCEPTION CHECKS

    if (cf.size() == 1) // degenerate case, integer
        return cf[0];

    double tmp = 1. / cf[cf.size() - 1];
    for (idx i = cf.size() - 2; i != 0; --i)
    {
        tmp = 1. / (tmp + cf[i]);
    }

    return cf[0] + tmp;
}

/**
* \brief Greatest common divisor of two integers
* \see qpp::lcm()
*
* \param m Integer
* \param n Integer
* \return Greatest common divisor of \a m and \a n
*/
inline bigint gcd(bigint m, bigint n)
{
    // EXCEPTION CHECKS

    if (m == 0 && n == 0)
        throw Exception("qpp::gcd()", Exception::Type::OUT_OF_RANGE);
    // END EXCEPTION CHECKS

    if (m == 0 || n == 0)
        return (std::max(m, n));

    bigint result = 1;
    while (n)
    {
        result = n;
        n = m % result;
        m = result;
    }

    return (result > 0) ?  result : -result;
}

/**
* \brief Greatest common divisor of a list of integers
* \see qpp::lcm()
*
* \param ns List of integers
* \return Greatest common divisor of all numbers in \a ns
*/
inline bigint gcd(const std::vector<bigint>& ns)
{
    // EXCEPTION CHECKS

    if (ns.size() == 0)
        throw Exception("qpp::gcd()", Exception::Type::ZERO_SIZE);
    // END EXCEPTION CHECKS

    bigint result = ns[0]; // convention: gcd({n}) = n
    for (idx i = 1; i < ns.size(); ++i)
    {
        result = gcd(result, ns[i]);
    }

    return  (result > 0) ?  result : -result;
}

/**
* \brief Least common multiple of two integers
* \see qpp::gcd()
*
* \param m Integer
* \param n Integer
* \return Least common multiple of \a m and \a n
*/
inline bigint lcm(bigint m, bigint n)
{
    if (m == 0 || n == 0)
        throw Exception("qpp::lcm()", Exception::Type::OUT_OF_RANGE);

    bigint result = m * n / gcd(m, n);

    return  (result > 0) ?  result : -result;
}

/**
* \brief Least common multiple of a list of integers
* \see qpp::gcd()
*
* \param ns List of integers
* \return Least common multiple of all numbers in \a ns
*/
inline bigint lcm(const std::vector<bigint>& ns)
{
    // EXCEPTION CHECKS

    if (ns.size() == 0)
        throw Exception("qpp::lcm()", Exception::Type::ZERO_SIZE);

    if (ns.size() == 1) // convention: lcm({n}) = n
        return ns[0];

    if (std::find(std::begin(ns), std::end(ns), 0) != std::end(ns))
        throw Exception("qpp::lcm()", Exception::Type::OUT_OF_RANGE);
    // END EXCEPTION CHECKS

    bigint prod =
            std::accumulate(std::begin(ns),
                            std::end(ns),
                            static_cast<bigint>(1),
                            std::multiplies<bigint>());

    bigint result = prod / gcd(ns);

    return  (result > 0) ?  result : -result;
}

/**
* \brief Inverse permutation
*
* \param perm Permutation
* \return Inverse of the permutation \a perm
*/
inline std::vector<idx> invperm(const std::vector<idx>& perm)
{
    // EXCEPTION CHECKS

    if (!internal::_check_perm(perm))
        throw Exception("qpp::invperm()", Exception::Type::PERM_INVALID);
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
                                 const std::vector<idx>& sigma)
{
    // EXCEPTION CHECKS

    if (!internal::_check_perm(perm))
        throw Exception("qpp::compperm()", Exception::Type::PERM_INVALID);
    if (!internal::_check_perm(sigma))
        throw Exception("qpp::compperm()", Exception::Type::PERM_INVALID);
    if (perm.size() != sigma.size())
        throw Exception("qpp::compperm()", Exception::Type::PERM_INVALID);
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
* \param n Integer different from 0, 1 or -1
* \return Integer vector containing the factors
*/
inline std::vector<bigint> factors(bigint n)
{
    if (n < 0)
        n = -n;

    // EXCEPTION CHECKS

    if (n == 0 || n == 1)
        throw Exception("qpp::factors()", Exception::Type::OUT_OF_RANGE);
    // END EXCEPTION CHECKS

    std::vector<bigint> result;
    bigint d = 2;

    while (n > 1)
    {
        while (n % d == 0)
        {
            result.push_back(d);
            n /= d;
        }
        ++d;
        if (d * d > n) // changes the running time from O(n) to O(sqrt(n))
        {
            if (n > 1)
            {
                result.push_back(n);
            }
            break;
        }
    }

    return result;
}

/**
* \brief Primality test
*
* \note Runs in \f$\mathcal{O}(\sqrt{n})\f$ time complexity. Use Miller-Rabin
* or something similar if you need high performance.
*
* \param n Integer different for 0, 1 or -1
* \return True if the number is prime, false otherwise
*/
inline bool isprime(bigint n)
{
    if (n < 0)
        n = -n;

    // EXCEPTION CHECKS

    if (n == 0 or n == 1)
        throw Exception("qpp::isprime()", Exception::Type::OUT_OF_RANGE);
    // END EXCEPTION CHECKS

    auto facts = factors(n);

    if (facts.size() == 1)
        return true;

    return false;
}

/**
* \brief Fast integer power modulo \a p based on
* the SQUARE-AND-MULTIPLY algorithm
*
* Computes \f$a^n \mathrm{ mod } p\f$
*
* \param a Non-negative integer
* \param n Non-negative integer
* \param p Strictly positive integer
* \return \f$a^n\mathrm{ mod }p\f$
*/
inline bigint modpow(bigint a, bigint n, bigint p)
{
    // EXCEPTION CHECKS

    if (p == 0 || (a == 0 && n == 0))
        throw Exception("qpp::modpow()",
                        Exception::Type::OUT_OF_RANGE);
    // END EXCEPTION CHECKS

    bigint result = 1;

    for (; n > 0; n /= 2)
    {
        if (n % 2)
            result = (result * a) % p; // MULTIPLY
        a = (a * a) % p; // SQUARE
    }

    return result;
}

/**
 * \brief Extended greatest common divisor of two integers
 * \see qpp::gcd()
 *
 * \param m Integer
 * \param n Integer
 * \return Tuple of: 1. Integer \f$a\f$, 2. Integer \f$b\f$,
 * and 3. Non-negative integer \f$gcd(a, b)\f$ such that
 * \f$am + bn = gcd(a, b)\f$
 */
inline std::tuple<bigint, bigint, bigint> egcd(bigint m, bigint n)
{
    // EXCEPTION CHECKS

    if (m == 0 && n == 0)
        throw Exception("qpp::egcd()", Exception::Type::OUT_OF_RANGE);

    // END EXCEPTION CHECKS

    bigint a, b, c;

    bigint q, r, a1, a2, b1, b2;

    a2 = 1, a1 = 0, b2 = 0, b1 = 1;
    while (n)
    {
        q = m / n, r = m - q * n;
        a = a2 - q * a1, b = b2 - q * b1;
        m = n, n = r;
        a2 = a1, a1 = a, b2 = b1, b1 = b;
    }
    c = m, a = a2, b = b2;

    // correct the signs
    if( c < 0)
    {
        a = -a;
        b = -b;
        c = -c;
    }

    return std::make_tuple(a, b, c);
}

/**
 * \brief Modular inverse of \a a mod \a p
 * \see qpp::egcd()
 *
 * \note \a a and \a p must be co-prime
 *
 * \param m Non-negative integer
 * \param p Non-negative integer
 * \return Modular inverse \f$a^{-1}\f$ of \a a mod \a p
 */
inline bigint modinv(bigint a, bigint p)
{
    // EXCEPTION CHECKS

    if (a == 0 || p == 0)
        throw Exception("qpp::modinv()", Exception::Type::OUT_OF_RANGE);

    bigint x, y;
    bigint gcd_ap;
    std::tie(x, y, gcd_ap) = egcd(p, a);

    if (gcd_ap != 1)
        throw Exception("qpp::modinv()", Exception::Type::OUT_OF_RANGE);

    // END EXCEPTION CHECKS

    return (y > 0) ? y : y + p;
}

} /* namespace qpp */

#endif /* NUMBER_THEORY_H_ */
