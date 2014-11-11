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

#ifndef INCLUDE_NUMBER_THEORY_H_
#define INCLUDE_NUMBER_THEORY_H_

// number theory functions

namespace qpp
{
/**
* \brief Simple continued fraction expansion
*
* \param x Real number
* \param n Number of terms in the expansion
* \param cut Stop the expansion when the next term is greater than \a cut
* \return Integer vector containing the simple continued fraction expansion of \a x.
* If there are \a m less than \a n terms in the expansion, a shorter vector with \a m components is returned.
*/
    std::vector<std::size_t> x2contfrac(double x, std::size_t n, std::size_t cut = 1e5)
    {
        if (n == 0)
            throw Exception("x2contfrac", Exception::Type::OUT_OF_RANGE);

        std::vector<std::size_t> result;

        for (std::size_t i = 0; i < n; ++i)
        {
            result.push_back(std::floor(x));
            x = 1. / (x - std::floor(x));
            if (!std::isfinite(x) || x > cut)
                return result;
        }
        return result;
    }

/**
* \brief Real representation of a simple continued fraction
*
* \param cf Integer vector containing the simple continued fraction expansion
* \param n Number of terms considered in the continued fraction expansion. If \a n is greater than the size of \a cf,
* then all terms in \a cf are considered.
* \return Real representation of the simple continued fraction
*/
    double contfrac2x(const std::vector<std::size_t> &cf, std::size_t n)
    {
        if (cf.size() == 0)
            throw Exception("contfrac2x", Exception::Type::ZERO_SIZE);

        if (n == 0)
            throw Exception("contfrac2x", Exception::Type::OUT_OF_RANGE);

        if (n > cf.size())
            n = cf.size();

        if (n == 1) // degenerate case, integer
            return cf[0];

        double tmp = 1. / cf[n - 1];
        for (int i = n - 2; i > 0; --i)
        {
            tmp = 1. / (tmp + cf[i]);
        }

        return cf[0] + tmp;
    }

/**
* \brief Real representation of a simple continued fraction
*
* \param cf Integer vector containing the simple continued fraction expansion
* \return Real representation of the simple continued fraction
*/
    double contfrac2x(const std::vector<std::size_t> &cf)
    {
        if (cf.size() == 0)
            throw Exception("contfrac2x", Exception::Type::ZERO_SIZE);

        if (cf.size() == 1) // degenerate case, integer
            return cf[0];

        double tmp = 1. / cf[cf.size() - 1];
        for (int i = cf.size() - 2; i > 0; --i)
        {
            tmp = 1. / (tmp + cf[i]);
        }

        return cf[0] + tmp;
    }

/**
* \brief Greatest common divisor of two non-negative integers
*
* \param m Non-negative integer
* \param n Non-negative integer
* \return Greatest common divisor of \a m and \a n
*/
    std::size_t gcd(std::size_t m, std::size_t n)
    {
        if (m == 0 || n == 0)
            return (std::max(m, n));

        std::size_t result = 1;
        while (n)
        {
            result = n;
            n = m % result;
            m = result;
        }

        return result;
    }

/**
* \brief Greatest common divisor of a list of non-negative integers
*
* \param ns List of non-negative integers
* \return Greatest common divisor of all numbers in \a ns
*/
    std::size_t gcd(const std::vector<std::size_t> &ns)
    {
        if (ns.size() == 0)
            throw Exception("gcd", Exception::Type::ZERO_SIZE);

        std::size_t result = ns[0]; // convention: gcd({n}) = n
        for (std::size_t i = 1; i < ns.size(); ++i)
        {
            result = gcd(result, ns[i]);
        }

        return result;
    }

/**
* \brief Least common multiple of two positive integers
*
* \param m Positive integer
* \param n Positive integer
* \return Least common multiple of \a m and \a n
*/
    std::size_t lcm(std::size_t m, std::size_t n)
    {
        if (m == 0 || n == 0)
            throw Exception("lcm", Exception::Type::OUT_OF_RANGE);

        return m * n / gcd(m, n);
    }

/**
* \brief Least common multiple of a list of positive integers
*
* \param ns List of positive integers
* \return Least common multiple of all numbers in \a ns
*/
    std::size_t lcm(const std::vector<std::size_t> &ns)
    {
        if (ns.size() == 0)
            throw Exception("lcm", Exception::Type::ZERO_SIZE);

        if (ns.size() == 1) // convention: lcm({n}) = n
            return ns[0];

        auto multiply = [](std::size_t x, std::size_t y) -> std::size_t
        {
            return x * y;
        };

        std::size_t prod = std::accumulate(std::begin(ns), std::end(ns), 1u,
                multiply);

        return prod / gcd(ns);
    }

} /* namespace qpp */

#endif /* INCLUDE_NUMBER_THEORY_H_ */