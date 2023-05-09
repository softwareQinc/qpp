#include "gmock/gmock.h" // for matchers
#include "gtest/gtest.h"
#include "qpp.h"

using namespace qpp;

// Unit testing "statistics.hpp"

/******************************************************************************/
/// BEGIN template<typename Container> realT avg(
///       const std::vector<realT>& prob, const Container& X,
///       typename std::enable_if<is_iterable<Container>::value>::type*
///       = nullptr)
TEST(qpp_avg, AllTests) {
    // size 1
    std::vector<realT> prob{1};
    std::vector<realT> X{10};
    EXPECT_NEAR(10, avg(prob, X), 1e-5);

    // size 2
    prob = {0.5, 0.5};
    X = {10, 20};
    EXPECT_NEAR(15, avg(prob, X), 1e-5);

    // size 3
    prob = {0.7, 0.2, 0.1};
    X = {1, 2, 3};
    EXPECT_NEAR(1.4, avg(prob, X), 1e-5);

    // uniform probability distribution of size 10
    prob = uniform(10);
    X = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    EXPECT_NEAR(11 / 2., avg(prob, X), 1e-5);
}
/******************************************************************************/
/// BEGIN template<typename Container> realT cor(const rmat& probXY,
///       const Container& X, const Container& Y,
///       typename std::enable_if<is_iterable<Container>::value>::type*
///       = nullptr)
TEST(qpp_cor, AllTests) {
    // decoupled size 2
    rmat probX(2, 1), probY(1, 2);
    probX << 0.6, 0.4;
    probY << 0.8, 0.2;
    rmat probXY = kron(probX, probY);
    std::vector<realT> X{1, 2};
    std::vector<realT> Y{3, 4};
    EXPECT_NEAR(0, cor(probXY, X, Y), 1e-5);

    // fully correlated size 2
    probXY = rmat::Identity(2, 2) / 2.;
    EXPECT_NEAR(1, cor(probXY, X, Y), 1e-5);

    // random size 2 x 3
    idx NX = 2, NY = 3;
    probXY = rmat::Zero(NX, NY);
    probXY << 0.1, 0.2, 0.3, 0.05, 0.1, 0.25;
    X = std::vector<realT>{1, 2};
    Y = std::vector<realT>{3, 4, 5};
    realT result = cor(probXY, X, Y);
    realT expected = cov(probXY, X, Y) / (sigma(marginalX(probXY), X) *
                                          sigma(marginalY(probXY), Y));
    EXPECT_NEAR(expected, result, 1e-5);
    // symmetry
    EXPECT_NEAR(cor(probXY, X, Y), cor(transpose(probXY), Y, X), 1e-5);
}
/******************************************************************************/
/// BEGIN template<typename Container> realT cov(const rmat& probXY,
///       const Container& X, const Container& Y,
///       typename std::enable_if<is_iterable<Container>::value>::type*
///       = nullptr)
TEST(qpp_cov, AllTests) {
    // size 1
    rmat probXY(1, 1);
    std::vector<realT> X{10};
    std::vector<realT> Y{20};
    probXY << 1;
    EXPECT_NEAR(0, cov(probXY, X, Y), 1e-5);

    // decoupled size 2
    rmat probX(2, 1), probY(1, 2);
    probX << 0.6, 0.4;
    probY << 0.8, 0.2;
    probXY = kron(probX, probY);
    X = std::vector<realT>{1, 2};
    Y = std::vector<realT>{3, 4};
    EXPECT_NEAR(0, cov(probXY, X, Y), 1e-5);

    // fully correlated size 2
    probXY = rmat::Identity(2, 2) / 2.;
    X = std::vector<realT>{1, 2};
    Y = std::vector<realT>{3, 4};
    EXPECT_NEAR(0.25, cov(probXY, X, Y), 1e-5);

    // random size 2 x 3
    probXY = rmat::Zero(2, 3);
    probXY << 0.1, 0.2, 0.3, 0.05, 0.1, 0.25;
    X = std::vector<realT>{1, 2};
    Y = std::vector<realT>{3, 4, 5};
    realT result = cov(probXY, X, Y);
    realT expected = 0.04;

    EXPECT_NEAR(expected, result, 1e-5);
    // symmetry
    EXPECT_NEAR(cov(probXY, X, Y), cov(transpose(probXY), Y, X), 1e-5);
}
/******************************************************************************/
/// BEGIN inline std::vector<realT> marginalX(const rmat& probXY)
TEST(qpp_marginalX, AllTests) {
    // size 1
    rmat probXY(1, 1);
    probXY << 1;
    std::vector<realT> result = marginalX(probXY);
    std::vector<realT> expected{1};
    EXPECT_EQ(result.size(), expected.size());
    EXPECT_NEAR(1, sum(result.begin(), result.end()), 1e-5);
    for (idx i = 0; i < result.size(); ++i) {
        EXPECT_GE(result[i], 0);
        EXPECT_NEAR(result[i], expected[i], 1e-5);
    }

    // 1 x 3
    probXY = rmat::Zero(1, 3);
    probXY << 0.15, 0.3, 0.55;
    result = marginalX(probXY);
    expected = std::vector<realT>{1};
    EXPECT_EQ(result.size(), expected.size());
    EXPECT_NEAR(1, sum(result.begin(), result.end()), 1e-5);
    for (idx i = 0; i < result.size(); ++i) {
        EXPECT_GE(result[i], 0);
        EXPECT_NEAR(result[i], expected[i], 1e-5);
    }

    // 3 x 1
    probXY = rmat::Zero(3, 1);
    probXY << 0.15, 0.3, 0.55;
    result = marginalX(probXY);
    expected = std::vector<realT>{0.15, 0.3, 0.55};
    EXPECT_EQ(result.size(), expected.size());
    EXPECT_NEAR(1, sum(result.begin(), result.end()), 1e-5);
    for (idx i = 0; i < result.size(); ++i) {
        EXPECT_GE(result[i], 0);
        EXPECT_NEAR(result[i], expected[i], 1e-5);
    }

    // 3 x 2
    probXY = rmat::Zero(3, 2);
    probXY << 0.1, 0.2, 0.3, 0.05, 0.1, 0.25;
    result = marginalX(probXY);
    expected = std::vector<realT>{0.3, 0.35, 0.35};
    EXPECT_EQ(result.size(), expected.size());
    EXPECT_NEAR(1, sum(result.begin(), result.end()), 1e-5);
    for (idx i = 0; i < result.size(); ++i) {
        EXPECT_GE(result[i], 0);
        EXPECT_NEAR(result[i], expected[i], 1e-5);
    }
}
/******************************************************************************/
/// BEGIN inline std::vector<realT> marginalY(const rmat& probXY)
TEST(qpp_marginalY, AllTests) {
    // size 1
    rmat probXY(1, 1);
    probXY << 1;
    std::vector<realT> result = marginalY(probXY);
    std::vector<realT> expected{1};
    EXPECT_EQ(result.size(), expected.size());
    EXPECT_NEAR(1, sum(result.begin(), result.end()), 1e-5);
    for (idx i = 0; i < result.size(); ++i) {
        EXPECT_GE(result[i], 0);
        EXPECT_NEAR(result[i], expected[i], 1e-5);
    }

    // 1 x 3
    probXY = rmat::Zero(1, 3);
    probXY << 0.15, 0.3, 0.55;
    result = marginalY(probXY);
    expected = std::vector<realT>{0.15, 0.3, 0.55};
    EXPECT_EQ(result.size(), expected.size());
    EXPECT_NEAR(1, sum(result.begin(), result.end()), 1e-5);
    for (idx i = 0; i < result.size(); ++i) {
        EXPECT_GE(result[i], 0);
        EXPECT_NEAR(result[i], expected[i], 1e-5);
    }

    // 3 x 1
    probXY = rmat::Zero(3, 1);
    probXY << 0.15, 0.3, 0.55;
    result = marginalY(probXY);
    expected = std::vector<realT>{1};
    EXPECT_EQ(result.size(), expected.size());
    EXPECT_NEAR(1, sum(result.begin(), result.end()), 1e-5);
    for (idx i = 0; i < result.size(); ++i) {
        EXPECT_GE(result[i], 0);
        EXPECT_NEAR(result[i], expected[i], 1e-5);
    }

    // 3 x 2
    probXY = rmat::Zero(3, 2);
    probXY << 0.1, 0.2, 0.3, 0.05, 0.15, 0.20;
    result = marginalY(probXY);
    expected = std::vector<realT>{0.55, 0.45};
    EXPECT_EQ(result.size(), expected.size());
    EXPECT_NEAR(1, sum(result.begin(), result.end()), 1e-5);
    for (idx i = 0; i < result.size(); ++i) {
        EXPECT_GE(result[i], 0);
        EXPECT_NEAR(result[i], expected[i], 1e-5);
    }
}
/******************************************************************************/
/// BEGIN template<typename Container> realT sigma(
///       const std::vector<realT>& prob,
///       const Container& X,
///       typename std::enable_if<is_iterable<Container>::value>::type*
///       = nullptr)
TEST(qpp_sigma, AllTests) {
    // size 1
    std::vector<realT> prob{1};
    std::vector<realT> X{10};
    EXPECT_NEAR(0, sigma(prob, X), 1e-5);

    // size 2
    prob = {0.5, 0.5};
    X = {10, 20};
    EXPECT_NEAR(5, sigma(prob, X), 1e-5);

    // size 3
    prob = {0.7, 0.2, 0.1};
    X = {1, 2, 3};
    EXPECT_NEAR(0.663324958071080, sigma(prob, X), 1e-5);

    // uniform probability distribution of size 10
    prob = uniform(10);
    X = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    EXPECT_NEAR(2.872281323269013, sigma(prob, X), 1e-5);
}
/******************************************************************************/
/// BEGIN inline std::vector<realT> uniform(idx N)
TEST(qpp_uniform, AllTests) {
    // size 1
    idx N = 1;
    EXPECT_THAT(uniform(N), testing::Each(1. / N));

    // size 2
    N = 2;
    EXPECT_THAT(uniform(N), testing::Each(1. / N));

    // size 10
    N = 10;
    EXPECT_THAT(uniform(N), testing::Each(1. / N));
}
/******************************************************************************/
/// BEGIN template<typename Container> realT var(
///       const std::vector<realT>& prob, const Container& X,
///       typename std::enable_if<is_iterable<Container>::value>::type*
///       = nullptr)
TEST(qpp_var, AllTests) {
    // size 1
    std::vector<realT> prob{1};
    std::vector<realT> X{10};
    EXPECT_NEAR(0, var(prob, X), 1e-5);

    // size 2
    prob = {0.5, 0.5};
    X = {10, 20};
    EXPECT_NEAR(25, var(prob, X), 1e-5);

    // size 3
    prob = {0.7, 0.2, 0.1};
    X = {1, 2, 3};
    EXPECT_NEAR(0.44, var(prob, X), 1e-5);

    // uniform probability distribution of size 10
    prob = uniform(10);
    X = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    EXPECT_NEAR(8.25, var(prob, X), 1e-5);
}
/******************************************************************************/
