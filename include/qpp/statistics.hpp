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
 * \file statistics.hpp
 * \brief Statistics functions
 */

#ifndef QPP_STATISTICS_HPP_
#define QPP_STATISTICS_HPP_

#include <type_traits>
#include <vector>

#include "qpp/types.hpp"

#include "qpp/classes/exception.hpp"
#include "qpp/internal/util.hpp"

namespace qpp {
/**
 * \brief Uniform probability distribution vector
 *
 * \param N Size of the alphabet
 * \return Real vector consisting of a uniform distribution of size \a N
 */
inline std::vector<realT> uniform(idx N) {
    // EXCEPTION CHECKS

    if (N == 0) {
        throw exception::ZeroSize("qpp::uniform()", "N");
    }
    // END EXCEPTION CHECKS

    return {
        std::vector<realT>(N, static_cast<realT>(1.) / static_cast<realT>(N))};
}

/**
 * \brief Marginal distribution
 *
 * \param probXY Real matrix representing the joint probability distribution of
 * \a X and \a Y in lexicographical order (\a X labels the rows, \a Y labels the
 * columns)
 * \return Real vector consisting of the marginal distribution of \a X
 */
inline std::vector<realT> marginalX(const rmat& probXY) {
    // EXCEPTION CHECKS

    if (!internal::check_nonzero_size(probXY)) {
        throw exception::ZeroSize("qpp::marginalX()", "probXY");
    }
    // END EXCEPTION CHECKS

    std::vector<realT> result(probXY.rows(), 0);
    for (idx i = 0; i < static_cast<idx>(probXY.rows()); ++i) {
        for (idx j = 0; j < static_cast<idx>(probXY.cols()); ++j) {
            result[i] += probXY(i, j);
        }
    }

    return result;
}

/**
 * \brief Marginal distribution
 *
 * \param probXY Real matrix representing the joint probability distribution of
 * \a X and \a Y in lexicographical order (\a X labels the rows, \a Y labels the
 * columns)
 * \return Real vector consisting of the marginal distribution of \a Y
 */
inline std::vector<realT> marginalY(const rmat& probXY) {
    // EXCEPTION CHECKS

    if (!internal::check_nonzero_size(probXY)) {
        throw exception::ZeroSize("qpp::marginalY()", "probXY");
    }
    // END EXCEPTION CHECKS

    return marginalX(probXY.transpose());
}

/**
 * \brief Average
 *
 * \param prob Real probability vector representing the probability distribution
 * of \a X
 *
 * \param X Real random variable values represented by an STL-like container
 * \return Average of \a X
 */
template <typename Container>
realT avg(
    const std::vector<realT>& prob, const Container& X,
    typename std::enable_if<is_iterable<Container>::value>::type* = nullptr) {
    // EXCEPTION CHECKS

    if (!internal::check_nonzero_size(prob)) {
        throw exception::ZeroSize("qpp::avg()", "prob");
    }
    if (!internal::check_matching_sizes(prob, X)) {
        throw exception::SizeMismatch("qpp::avg()", "prob/X");
    }
    // END EXCEPTION CHECKS

    realT result = 0;
    for (idx i = 0; i < static_cast<idx>(prob.size()); ++i) {
        result += prob[i] * X[i];
    }

    return result;
}

/**
 * \brief Covariance
 *
 * \param probXY Real matrix representing the joint probability distribution of
 * \a X and \a Y in lexicographical order (\a X labels the rows, \a Y labels the
 * columns)
 * \param X Real random variable values represented by an STL-like container
 * \param Y Real random variable values represented by an STL-like container
 * \return Covariance of \a X and \a Y
 */
template <typename Container>
realT cov(
    const rmat& probXY, const Container& X, const Container& Y,
    typename std::enable_if<is_iterable<Container>::value>::type* = nullptr) {
    // EXCEPTION CHECKS

    if (!internal::check_nonzero_size(X)) {
        throw exception::ZeroSize("qpp::cov()", "X");
    }
    if (!internal::check_nonzero_size(Y)) {
        throw exception::ZeroSize("qpp::cov()", "Y");
    }
    if (static_cast<idx>(probXY.rows()) != static_cast<idx>(X.size())) {
        throw exception::SizeMismatch("qpp::cov()", "probXY/X");
    }
    if (static_cast<idx>(probXY.cols()) != static_cast<idx>(Y.size())) {
        throw exception::SizeMismatch("qpp::cov()", "probXY/Y");
    }
    // END EXCEPTION CHECKS

    std::vector<realT> probX = marginalX(probXY); // marginals
    std::vector<realT> probY = marginalY(probXY); // marginals

    realT result = 0;
    for (idx i = 0; i < static_cast<idx>(X.size()); ++i) {
        for (idx j = 0; j < static_cast<idx>(Y.size()); ++j) {
            result += probXY(i, j) * X[i] * Y[j];
        }
    }

    return result - avg(probX, X) * avg(probY, Y);
}

/**
 * \brief Variance
 *
 * \param prob Real probability vector representing the probability distribution
 * of \a X
 * \param X Real random variable values represented by an STL-like container
 * \return Variance of \a X
 */
template <typename Container>
realT var(
    const std::vector<realT>& prob, const Container& X,
    typename std::enable_if<is_iterable<Container>::value>::type* = nullptr) {
    // EXCEPTION CHECKS

    if (!internal::check_nonzero_size(prob)) {
        throw exception::ZeroSize("qpp::var()", "prob");
    }
    if (!internal::check_matching_sizes(prob, X)) {
        throw exception::SizeMismatch("qpp::var()", "prob/X");
    }
    // END EXCEPTION CHECKS

    dyn_col_vect<realT> diag(prob.size());
    for (idx i = 0; i < static_cast<idx>(prob.size()); ++i) {
        diag(i) = prob[i];
    }
    rmat probXX = diag.asDiagonal();

    return cov(probXX, X, X);
}

/**
 * \brief Standard deviation
 *
 * \param prob Real probability vector representing the probability distribution
 * of \a X
 * \param X Real random variable values represented by an STL-like container
 * \return Standard deviation of \a X
 */
template <typename Container>
realT sigma(
    const std::vector<realT>& prob, const Container& X,
    typename std::enable_if<is_iterable<Container>::value>::type* = nullptr) {
    // EXCEPTION CHECKS

    if (!internal::check_nonzero_size(prob)) {
        throw exception::ZeroSize("qpp::sigma()", "prob");
    }
    if (!internal::check_matching_sizes(prob, X)) {
        throw exception::SizeMismatch("qpp::sigma()", "prob/X");
    }
    // END EXCEPTION CHECKS

    return std::sqrt(var(prob, X));
}

/**
 * \brief Correlation
 *
 * \param probXY Real matrix representing the joint probability distribution of
 * \a X and \a Y in lexicographical order (\a X labels the rows, \a Y labels the
 * columns)
 * \param X Real random variable values represented by an STL-like container
 * \param Y Real random variable values represented by an STL-like container
 * \return Correlation of \a X and \a Y
 */
template <typename Container>
realT cor(
    const rmat& probXY, const Container& X, const Container& Y,
    typename std::enable_if<is_iterable<Container>::value>::type* = nullptr) {
    // EXCEPTION CHECKS

    if (!internal::check_nonzero_size(X)) {
        throw exception::ZeroSize("qpp::cor()", "X");
    }
    if (!internal::check_nonzero_size(Y)) {
        throw exception::ZeroSize("qpp::cor()", "Y");
    }
    if (static_cast<idx>(probXY.rows()) != X.size()) {
        throw exception::SizeMismatch("qpp::cor()", "probXY/X");
    }
    if (static_cast<idx>(probXY.cols()) != Y.size()) {
        throw exception::SizeMismatch("qpp::cor()", "probXY/Y");
    }
    // END EXCEPTION CHECKS

    return cov(probXY, X, Y) /
           (sigma(marginalX(probXY), X) * sigma(marginalY(probXY), Y));
}

} /* namespace qpp */

#endif /* QPP_STATISTICS_HPP_ */
