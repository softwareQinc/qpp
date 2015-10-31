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
* \file statistics.h
* \brief Statistics functions
*/

#ifndef STATISTICS_H_
#define STATISTICS_H_

// Collection of statistics-related functions

namespace qpp
{

/**
* \brief Uniform probability distribution vector
*
* \param N Size of the alphabet
* \return Real vector consisting of a uniform distribution of size \a N
*/
inline std::vector<double> uniform(idx N)
{
    // EXCEPTION CHECKS

    if (N == 0)
        throw Exception("qpp::uniform", Exception::Type::ZERO_SIZE);
    // END EXCEPTION CHECKS

    return std::vector<double>(N, 1. / N);
}

/**
* \brief Marginal distribution
*
* \param probXY Real matrix representing the joint probability distribution
* of \a X and \a Y in lexicographical order (\a X labels the rows, \a Y
* labels the columns)
* \return Real vector consisting of the marginal distribution of \a X
*/
inline std::vector<double> marginalX(const dmat& probXY)
{
    // EXCEPTION CHECKS

    if (!internal::_check_nonzero_size(probXY))
        throw Exception("qpp::marginalX", Exception::Type::ZERO_SIZE);
    // END EXCEPTION CHECKS

    std::vector<double> result(probXY.rows(), 0);
    for (idx i = 0; i < static_cast<idx>(probXY.rows()); ++i)
    {
        for (idx j = 0; j < static_cast<idx>(probXY.cols()); ++j)
        {
            result[i] += probXY(i, j);
        }
    }

    return result;
}

/**
* \brief Marginal distribution
*
* \param probXY Real matrix representing the joint probability distribution
* of \a X and \a Y in lexicographical order (\a X labels the rows, \a Y
* labels the columns)
* \return Real vector consisting of the marginal distribution of \a Y
*/
inline std::vector<double> marginalY(const dmat& probXY)
{
    // EXCEPTION CHECKS

    if (!internal::_check_nonzero_size(probXY))
        throw Exception("qpp::marginalY", Exception::Type::ZERO_SIZE);
    // END EXCEPTION CHECKS

    return marginalX(probXY.transpose());
}

/**
* \brief Average
*
* \param prob Real probability vector representing the probability
* distribution of \a X
*
* \param X Random variable values represented by an STL-like container 
* \return Average of \a X
*/
template<typename Container>
double avg(const std::vector<double>& prob, const Container& X)
{
    // EXCEPTION CHECKS

    if (!internal::_check_nonzero_size(prob))
        throw Exception("qpp:avg", Exception::Type::ZERO_SIZE);
    if (!internal::_check_matching_sizes(prob, X))
        throw Exception("qpp:avg", Exception::Type::SIZE_MISMATCH);
    // END EXCEPTION CHECKS

    double result = 0;
    for (idx i = 0; i < prob.size(); ++i)
        result += prob[i] * X[i];

    return result;
}

/**
* \brief Covariance
*
* \param probXY Real matrix representing the joint probability distribution
* of \a X and \a Y in lexicographical order (\a X labels the rows, \a Y
* labels the columns)
* \param X Random variable values represented by an STL-like container
* \param Y Random variable values represented by an STL-like container
* \return Covariance of \a X and \a Y
*/
template<typename Container>
double cov(const dmat& probXY,
           const Container& X,
           const Container& Y)
{
    // EXCEPTION CHECKS

    if (!internal::_check_nonzero_size(X) || !internal::_check_nonzero_size(Y))
        throw Exception("qpp:cov", Exception::Type::ZERO_SIZE);
    if (static_cast<idx>(probXY.rows()) != X.size() ||
        static_cast<idx>(probXY.cols()) != Y.size())
        throw Exception("qpp:cov", Exception::Type::SIZE_MISMATCH);
    // END EXCEPTION CHECKS

    std::vector<double> probX = marginalX(probXY); // marginals
    std::vector<double> probY = marginalY(probXY); // marginals

    double result = 0;
    for (idx i = 0; i < X.size(); ++i)
    {
        for (idx j = 0; j < Y.size(); ++j)
        {
            result += probXY(i, j) * X[i] * Y[j];
        }
    }

    return result - avg(probX, X) * avg(probY, Y);
}

/**
* \brief Variance
*
* \param prob Real probability vector representing the probability
* distribution of \a X
* \param X Random variable values represented by an STL-like container
* \return Variance of \a X
*/
template<typename Container>
double var(const std::vector<double>& prob, const Container& X)
{
    // EXCEPTION CHECKS

    if (!internal::_check_nonzero_size(prob))
        throw Exception("qpp:var", Exception::Type::ZERO_SIZE);
    if (!internal::_check_matching_sizes(prob, X))
        throw Exception("qpp:var", Exception::Type::SIZE_MISMATCH);
    // END EXCEPTION CHECKS

    Eigen::VectorXd diag(prob.size());
    for (idx i = 0; i < prob.size(); ++i)
        diag(i) = prob[i];
    dmat probXX = diag.asDiagonal();

    return cov(probXX, X, X);
}

/**
* \brief Standard deviation
*
* \param prob Real probability vector representing the probability
* distribution of \a X
* \param X Random variable values represented by an STL-like container
* \return Standard deviation of \a X
*/
template<typename Container>
double sigma(const std::vector<double>& prob, const Container& X)
{
    // EXCEPTION CHECKS

    if (!internal::_check_nonzero_size(prob))
        throw Exception("qpp:sigma", Exception::Type::ZERO_SIZE);
    if (!internal::_check_matching_sizes(prob, X))
        throw Exception("qpp:sigma", Exception::Type::SIZE_MISMATCH);
    // END EXCEPTION CHECKS

    return std::sqrt(var(prob, X));
}

/**
* \brief Correlation
*
* \param probXY Real matrix representing the joint probability distribution
* of \a X and \a Y in lexicographical order (\a X labels the rows, \a Y
* labels the columns)
* \param X Random variable values represented by an STL-like container
* \param Y Random variable values represented by an STL-like container
* \return Correlation of \a X and \a Y
*/
template<typename Container>
double cor(const dmat& probXY,
           const Container& X,
           const Container& Y)
{
    // EXCEPTION CHECKS

    if (!internal::_check_nonzero_size(X) || !internal::_check_nonzero_size(Y))
        throw Exception("qpp:cor", Exception::Type::ZERO_SIZE);
    if (static_cast<idx>(probXY.rows()) != X.size() ||
        static_cast<idx>(probXY.cols()) != Y.size())
        throw Exception("qpp:cor", Exception::Type::SIZE_MISMATCH);
    // END EXCEPTION CHECKS

    return cov(probXY, X, Y) / (sigma(marginalX(probXY), X) *
                                sigma(marginalX(probXY), Y));
}

} /* namespace qpp */

#endif /* STATISTICS_H_ */

