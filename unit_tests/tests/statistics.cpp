/*
 * Quantum++
 *
 * Copyright (c) 2013 - 2017 Vlad Gheorghiu (vgheorgh@gmail.com)
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

#include "gmock/gmock.h" // for matchers
#include "gtest/gtest.h"
#include "qpp.h"

using namespace qpp;

// Unit testing "statistics.h"

/******************************************************************************/
/// BEGIN template<typename Container> double qpp::avg(
///       const std::vector<double>& prob,
///       const Container& X,
///       typename std::enable_if<is_iterable<Container>::value>::type*
///       = nullptr)
TEST(qpp_avg, AllTests)
{
    // size 1
    std::vector<double> prob{1};
    std::vector<double> X{10};
    EXPECT_NEAR(10, qpp::avg(prob, X), 1e-7);

    // size 2
    prob = {0.5, 0.5};
    X = {10, 20};
    EXPECT_NEAR(15, qpp::avg(prob, X), 1e-7);

    // size 3
    prob = {0.7, 0.2, 0.1};
    X = {1, 2, 3};
    EXPECT_NEAR(1.4, qpp::avg(prob, X), 1e-7);

    // uniform probability distribution of size 10
    prob = uniform(10);
    X = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    EXPECT_NEAR(11 / 2., qpp::avg(prob, X), 1e-7);
}
/******************************************************************************/
/// BEGIN template<typename Container> double qpp::cor(
///       const dmat& probXY,
///       const Container& X,
///       const Container& Y,
///       typename std::enable_if<is_iterable<Container>::value>::type*
///       = nullptr)
TEST(qpp_cor, AllTests)
{
    // decoupled size 2
    dmat probX(2, 1), probY(1, 2);
    probX << 0.6, 0.4;
    probY << 0.8, 0.2;
    dmat probXY = kron(probX, probY);
    std::vector<double> X{1, 2};
    std::vector<double> Y{3, 4};
    EXPECT_NEAR(0, qpp::cor(probXY, X, Y), 1e-7);

    // fully correlated size 2
    probXY = dmat::Identity(2, 2) / 2.;
    EXPECT_NEAR(1, qpp::cor(probXY, X, Y), 1e-7);

    // random size 2 x 3
    idx NX = 2, NY = 3;
    probXY = dmat::Zero(NX, NY);
    probXY << 0.1, 0.2, 0.3, 0.05, 0.1, 0.25;
    X = std::vector<double>{1, 2};
    Y = std::vector<double>{3, 4, 5};
    double result = qpp::cor(probXY, X, Y);
    double expected = cov(probXY, X, Y) /
                      (sigma(marginalX(probXY), X) *
                       sigma(marginalY(probXY), Y));
    EXPECT_NEAR(expected, result, 1e-7);
    // symmetry
    EXPECT_NEAR(qpp::cor(probXY, X, Y),
                qpp::cor(transpose(probXY), Y, X), 1e-7);
}
/******************************************************************************/
/// BEGIN template<typename Container> double qpp::cov(
///       const dmat& probXY,
///       const Container& X,
///       const Container& Y,
///       typename std::enable_if<is_iterable<Container>::value>::type*
///       = nullptr)
TEST(qpp_cov, AllTests)
{
    // size 1
    dmat probXY(1, 1);
    std::vector<double> X{10};
    std::vector<double> Y{20};
    probXY << 1;
    EXPECT_NEAR(0, qpp::cov(probXY, X, Y), 1e-7);

    // decoupled size 2
    dmat probX(2, 1), probY(1, 2);
    probX << 0.6, 0.4;
    probY << 0.8, 0.2;
    probXY = kron(probX, probY);
    X = std::vector<double>{1, 2};
    Y = std::vector<double>{3, 4};
    EXPECT_NEAR(0, qpp::cov(probXY, X, Y), 1e-7);

    // fully correlated size 2
    probXY = dmat::Identity(2, 2) / 2.;
    X = std::vector<double>{1, 2};
    Y = std::vector<double>{3, 4};
    EXPECT_NEAR(0.25, qpp::cov(probXY, X, Y), 1e-7);

    // random size 2 x 3
    probXY = dmat::Zero(2, 3);
    probXY << 0.1, 0.2, 0.3, 0.05, 0.1, 0.25;
    X = std::vector<double>{1, 2};
    Y = std::vector<double>{3, 4, 5};
    double result = qpp::cov(probXY, X, Y);
    double expected = 0.04;

    EXPECT_NEAR(expected, result, 1e-7);
    // symmetry
    EXPECT_NEAR(qpp::cov(probXY, X, Y),
                qpp::cov(transpose(probXY), Y, X), 1e-7);
}
/******************************************************************************/
/// BEGIN inline std::vector<double> qpp::marginalX(const dmat& probXY)
TEST(qpp_marginalX, AllTests)
{
    // size 1
    dmat probXY(1, 1);
    probXY << 1;
    std::vector<double> result = marginalX(probXY);
    std::vector<double> expected{1};
    EXPECT_EQ(result.size(), expected.size());
    EXPECT_NEAR(1, sum(result.begin(), result.end()), 1e-7);
    for (idx i = 0; i < result.size(); ++i)
    {
        EXPECT_GE(result[i], 0);
        EXPECT_NEAR(result[i], expected[i], 1e-7);
    }

    // 1 x 3
    probXY = dmat::Zero(1, 3);
    probXY << 0.15, 0.3, 0.55;
    result = marginalX(probXY);
    expected = std::vector<double>{1};
    EXPECT_EQ(result.size(), expected.size());
    EXPECT_NEAR(1, sum(result.begin(), result.end()), 1e-7);
    for (idx i = 0; i < result.size(); ++i)
    {
        EXPECT_GE(result[i], 0);
        EXPECT_NEAR(result[i], expected[i], 1e-7);
    }

    // 3 x 1
    probXY = dmat::Zero(3, 1);
    probXY << 0.15, 0.3, 0.55;
    result = marginalX(probXY);
    expected = std::vector<double>{0.15, 0.3, 0.55};
    EXPECT_EQ(result.size(), expected.size());
    EXPECT_NEAR(1, sum(result.begin(), result.end()), 1e-7);
    for (idx i = 0; i < result.size(); ++i)
    {
        EXPECT_GE(result[i], 0);
        EXPECT_NEAR(result[i], expected[i], 1e-7);
    }

    // 3 x 2
    probXY = dmat::Zero(3, 2);
    probXY << 0.1, 0.2, 0.3, 0.05, 0.1, 0.25;
    result = marginalX(probXY);
    expected = std::vector<double>{0.3, 0.35, 0.35};
    EXPECT_EQ(result.size(), expected.size());
    EXPECT_NEAR(1, sum(result.begin(), result.end()), 1e-7);
    for (idx i = 0; i < result.size(); ++i)
    {
        EXPECT_GE(result[i], 0);
        EXPECT_NEAR(result[i], expected[i], 1e-7);
    }
}
/******************************************************************************/
/// BEGIN inline std::vector<double> qpp::marginalY(const dmat& probXY)
TEST(qpp_marginalY, AllTests)
{
    // size 1
    dmat probXY(1, 1);
    probXY << 1;
    std::vector<double> result = marginalY(probXY);
    std::vector<double> expected{1};
    EXPECT_EQ(result.size(), expected.size());
    EXPECT_NEAR(1, sum(result.begin(), result.end()), 1e-7);
    for (idx i = 0; i < result.size(); ++i)
    {
        EXPECT_GE(result[i], 0);
        EXPECT_NEAR(result[i], expected[i], 1e-7);
    }

    // 1 x 3
    probXY = dmat::Zero(1, 3);
    probXY << 0.15, 0.3, 0.55;
    result = marginalY(probXY);
    expected = std::vector<double>{0.15, 0.3, 0.55};
    EXPECT_EQ(result.size(), expected.size());
    EXPECT_NEAR(1, sum(result.begin(), result.end()), 1e-7);
    for (idx i = 0; i < result.size(); ++i)
    {
        EXPECT_GE(result[i], 0);
        EXPECT_NEAR(result[i], expected[i], 1e-7);
    }

    // 3 x 1
    probXY = dmat::Zero(3, 1);
    probXY << 0.15, 0.3, 0.55;
    result = marginalY(probXY);
    expected = std::vector<double>{1};
    EXPECT_EQ(result.size(), expected.size());
    EXPECT_NEAR(1, sum(result.begin(), result.end()), 1e-7);
    for (idx i = 0; i < result.size(); ++i)
    {
        EXPECT_GE(result[i], 0);
        EXPECT_NEAR(result[i], expected[i], 1e-7);
    }

    // 3 x 2
    probXY = dmat::Zero(3, 2);
    probXY << 0.1, 0.2, 0.3, 0.05, 0.15, 0.20;
    result = marginalY(probXY);
    expected = std::vector<double>{0.55, 0.45};
    EXPECT_EQ(result.size(), expected.size());
    EXPECT_NEAR(1, sum(result.begin(), result.end()), 1e-7);
    for (idx i = 0; i < result.size(); ++i)
    {
        EXPECT_GE(result[i], 0);
        EXPECT_NEAR(result[i], expected[i], 1e-7);
    }
}
/******************************************************************************/
/// BEGIN template<typename Container> double qpp::sigma(
///       const std::vector<double>& prob,
///       const Container& X,
///       typename std::enable_if<is_iterable<Container>::value>::type*
///       = nullptr)
TEST(qpp_sigma, AllTests)
{
    // size 1
    std::vector<double> prob{1};
    std::vector<double> X{10};
    EXPECT_NEAR(0, qpp::sigma(prob, X), 1e-7);

    // size 2
    prob = {0.5, 0.5};
    X = {10, 20};
    EXPECT_NEAR(5, qpp::sigma(prob, X), 1e-7);

    // size 3
    prob = {0.7, 0.2, 0.1};
    X = {1, 2, 3};
    EXPECT_NEAR(0.663324958071080, qpp::sigma(prob, X), 1e-7);

    // uniform probability distribution of size 10
    prob = uniform(10);
    X = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    EXPECT_NEAR(2.872281323269013, qpp::sigma(prob, X), 1e-7);
}
/******************************************************************************/
/// BEGIN inline std::vector<double> qpp::uniform(idx N)
TEST(qpp_uniform, AllTests)
{
    // size 1
    idx N = 1;
    EXPECT_THAT(qpp::uniform(N), testing::Each(1. / N));

    // size 2
    N = 2;
    EXPECT_THAT(qpp::uniform(N), testing::Each(1. / N));

    // size 10
    N = 10;
    EXPECT_THAT(qpp::uniform(N), testing::Each(1. / N));
}
/******************************************************************************/
/// BEGIN template<typename Container> double qpp::var(
///       const std::vector<double>& prob,
///       const Container& X,
///       typename std::enable_if<is_iterable<Container>::value>::type*
///       = nullptr)
TEST(qpp_var, AllTests)
{
    // size 1
    std::vector<double> prob{1};
    std::vector<double> X{10};
    EXPECT_NEAR(0, qpp::var(prob, X), 1e-7);

    // size 2
    prob = {0.5, 0.5};
    X = {10, 20};
    EXPECT_NEAR(25, qpp::var(prob, X), 1e-7);

    // size 3
    prob = {0.7, 0.2, 0.1};
    X = {1, 2, 3};
    EXPECT_NEAR(0.44, qpp::var(prob, X), 1e-7);

    // uniform probability distribution of size 10
    prob = uniform(10);
    X = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    EXPECT_NEAR(8.25, qpp::var(prob, X), 1e-7);
}
/******************************************************************************/
