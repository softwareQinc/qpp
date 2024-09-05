#include <algorithm>

#include "gtest/gtest.h"

#include "qpp/qpp.hpp"

using namespace qpp;

// Unit testing "qpp/random.hpp"

/// BEGIN inline bool bernoulli(realT p = 0.5)
TEST(qpp_bernoulli, AllTests) {}

/// BEGIN template <typename T,
///       typename std::enable_if_t<std::is_arithmetic_v<T>>* = nullptr>
///       rand(T a, T b)
TEST(qpp_rand, Integer) {
    // 1 element equal boundaries
    bigint a = 42, b = a;
    EXPECT_EQ(a, rand(a, b));

    // 1000 elements between 0 and 1
    idx N = 1000;
    a = 0, b = 1;
    for (idx i = 0; i < N; ++i) {
        bigint n = rand(a, b);
        EXPECT_GE(n, a);
        EXPECT_LE(n, b);
    }

    // 10000 elements between -10 and 10
    N = 10000;
    a = -10, b = 10;
    realT average = 0;
    for (idx i = 0; i < N; ++i) {
        bigint n = rand(a, b);

        EXPECT_GE(n, a);
        EXPECT_LE(n, b);
        average += n;
    }
    average /= N;
    EXPECT_NEAR(0, average, 2e-1); // very likely
}
TEST(qpp_rand, Real) {
    // 1000 elements between 0 and 1
    idx N = 1000;
    realT a = 0, b = 1;
    for (idx i = 0; i < N; ++i) {
        realT n = rand(a, b);
        EXPECT_GE(n, a);
        EXPECT_LE(n, b);
    }

    // 10000 elements between -10 and 10
    N = 10000;
    a = -10, b = 10;
    realT average = 0;
    for (idx i = 0; i < N; ++i) {
        realT n = rand(a, b);
        EXPECT_GE(n, a);
        EXPECT_LE(n, b);
        average += n;
    }
    average /= N;
    EXPECT_NEAR(0, average, 2e-1); // very likely
}

/// BEGIN template<> inline cmat rand(idx rows, idx cols, realT a, realT b)
TEST(qpp_rand, ComplexMatrix) {
    // 1 x 1 matrix with elements between -1 and 1
    idx NA = 1, NB = 1;
    realT a = -1, b = 1;
    cmat A = rand<cmat>(NA, NB, a, b);
    for (idx i = 0; i < NA * NB; ++i) {
        cplx elem = A.data()[i];
        EXPECT_GE(elem.real(), a);
        EXPECT_GE(elem.imag(), a);
        EXPECT_LT(elem.real(), b);
        EXPECT_LT(elem.imag(), b);
    }

    // 1 x 10 matrix with elements between -10 and 10
    NA = 1, NB = 10;
    a = -10, b = 10;
    A = rand<cmat>(NA, NB, a, b);
    for (idx i = 0; i < NA * NB; ++i) {
        cplx elem = A.data()[i];
        EXPECT_GE(elem.real(), a);
        EXPECT_GE(elem.imag(), a);
        EXPECT_LT(elem.real(), b);
        EXPECT_LT(elem.imag(), b);
    }

    // 25 x 15 matrix with elements between -10 and 10
    NA = 25, NB = 15;
    a = -10, b = 10;
    A = rand<cmat>(NA, NB, a, b);
    cplx average = 0;
    for (idx i = 0; i < NA * NB; ++i) {
        cplx elem = A.data()[i];
        EXPECT_GE(elem.real(), a);
        EXPECT_GE(elem.imag(), a);
        EXPECT_LT(elem.real(), b);
        EXPECT_LT(elem.imag(), b);
        average += elem;
    }
    average /= static_cast<realT>(NA * NB);
    EXPECT_NEAR(0, average.real(), 1); // very likely
    EXPECT_NEAR(0, average.imag(), 1); // very likely
}

/// BEGIN template<> inline rmat rand(idx rows, idx cols, realT a, realT b)
TEST(qpp_rand, RealMatrix) {
    // 1 x 1 matrix with elements between -1 and 1
    idx NA = 1, NB = 1;
    realT a = -1, b = 1;
    rmat A = rand<rmat>(NA, NB, a, b);
    for (idx i = 0; i < NA * NB; ++i) {
        realT elem = A.data()[i];
        EXPECT_GE(elem, a);
        EXPECT_LT(elem, b);
    }

    // 1 x 10 matrix with elements between -10 and 10
    NA = 1, NB = 10;
    a = -10, b = 10;
    A = rand<rmat>(NA, NB, a, b);
    for (idx i = 0; i < NA * NB; ++i) {
        realT elem = A.data()[i];
        EXPECT_GE(elem, a);
        EXPECT_LT(elem, b);
    }

    // 25 x 15 matrix with elements between -10 and 10
    NA = 25, NB = 15;
    a = -10, b = 10;
    A = rand<rmat>(NA, NB, a, b);
    realT average = 0;
    for (idx i = 0; i < NA * NB; ++i) {
        realT elem = A.data()[i];
        EXPECT_GE(elem, a);
        EXPECT_LT(elem, b);
        average += elem;
    }
    average /= static_cast<realT>(NA * NB);
    EXPECT_NEAR(0, average, 1); // very likely
}

/// BEGIN inline cmat randH(idx D = 2)
TEST(qpp_randH, AllTests) {
    // D = 1 degenerate case
    idx D = 1;
    cmat A = randH(D);
    EXPECT_NEAR(0, norm(A - adjoint(A)), 1e-5);

    // D = 2
    D = 2;
    A = randH(D);
    EXPECT_NEAR(0, norm(A - adjoint(A)), 1e-5);

    // D = 10
    D = 10;
    A = randH(D);
    EXPECT_NEAR(0, norm(A - adjoint(A)), 1e-5);
}

/// BEGIN inline idx randidx(idx a = std::numeric_limits<idx>::min(),
///       idx b = std::numeric_limits<idx>::max())
TEST(qpp_randidx, AllTests) {
    // 1 element equal boundaries
    idx a = 42, b = a;
    EXPECT_EQ(a, randidx(a, b));

    // 1000 elements between 0 and 1
    idx N = 1000;
    a = 0, b = 1;
    for (idx i = 0; i < N; ++i) {
        idx n = randidx(a, b);
        EXPECT_GE(n, a);
        EXPECT_LE(n, b);
    }

    // 10000 elements between -10 and 10
    N = 10000;
    a = 0, b = 10;
    realT average = 0;
    for (idx i = 0; i < N; ++i) {
        idx n = randidx(a, b);
        EXPECT_GE(n, a);
        EXPECT_LE(n, b);
        average += n;
    }
    average /= N;
    EXPECT_NEAR(5, average, 1e-1); // very likely
}

/// BEGIN inline ket randket(idx D = 2)
TEST(qpp_randket, AllTests) {
    // D = 1 degenerate case
    idx D = 1;
    ket expected = ket::Ones(D);
    EXPECT_NEAR(0, norm(prj(expected) - prj(randket(D))), 1e-5);

    // D = 2, 10000 randomly generated kets, should average to zero
    D = 2;
    idx N = 10000; // number of runs
    ket avg_state = ket::Zero(D);
    for (idx i = 0; i < N; ++i) {
        avg_state += randket(D);
    }
    expected = ket::Zero(D);
    EXPECT_NEAR(0, norm(expected - avg_state / N), 2e-2);

    // D = 5, 10000 randomly generated kets, should average to zero
    D = 5;
    N = 10000; // number of runs
    avg_state = ket::Zero(D);
    for (idx i = 0; i < N; ++i) {
        avg_state += randket(D);
    }
    expected = ket::Zero(D);
    EXPECT_NEAR(0, norm(expected - avg_state / N), 2e-2);
}

/// BEGIN inline std::vector<cmat> randkraus(idx N, idx Din = 2, idx Dout)
TEST(qpp_randkraus, AllTests) {
    // D = 1, N = 1 degenerate case
    idx D = 1, N = 1;
    std::vector<cmat> Ks = randkraus(N, D);
    cmat closure = cmat::Zero(D, D);
    for (idx i = 0; i < N; ++i) {
        closure += adjoint(Ks[i]) * Ks[i];
    }
    EXPECT_NEAR(0, norm(closure - gt.Id(D)), 1e-5);

    // D = 1, N = 10 degenerate case
    D = 1, N = 10;
    Ks = randkraus(N, D);
    closure = cmat::Zero(D, D);
    for (idx i = 0; i < N; ++i) {
        closure += adjoint(Ks[i]) * Ks[i];
    }
    EXPECT_NEAR(0, norm(closure - gt.Id(D)), 1e-5);

    // D = 2, N = 1
    D = 2, N = 1;
    Ks = randkraus(N);
    closure = cmat::Zero(D, D);
    for (idx i = 0; i < N; ++i) {
        closure += adjoint(Ks[i]) * Ks[i];
    }
    EXPECT_NEAR(0, norm(closure - gt.Id(D)), 1e-5);

    // D = 2, N = 10
    D = 2, N = 10;
    Ks = randkraus(N);
    closure = cmat::Zero(D, D);
    for (idx i = 0; i < N; ++i) {
        closure += adjoint(Ks[i]) * Ks[i];
    }
    EXPECT_NEAR(0, norm(closure - gt.Id(D)), 1e-5);

    // D = 5, N = 5
    D = 5, N = 10;
    Ks = randkraus(N, D);
    closure = cmat::Zero(D, D);
    for (idx i = 0; i < N; ++i) {
        closure += adjoint(Ks[i]) * Ks[i];
    }
    EXPECT_NEAR(0, norm(closure - gt.Id(D)), 1e-5);

    // Din = 3, Dout = 5, N = 6
    idx Din = 3, Dout = 5;
    N = 6;
    Ks = randkraus(N, Din, Dout);
    closure = cmat::Zero(Din, Din);
    for (idx i = 0; i < N; ++i) {
        closure += adjoint(Ks[i]) * Ks[i];
    }
    EXPECT_NEAR(0, norm(closure - gt.Id(Din)), 1e-5);
}

/// BEGIN template <typename T,
//        typename std::enable_if_t<std::is_arithmetic_v<T>>* = nullptr>
//        T randn(T mean = 0, T sigma = 1)
TEST(qpp_randn, Real) {
    // 10000 elements with mean = 0 and sigma = 1
    realT mean = 0, sigma = 1;
    idx N = 10000;
    realT average = 0;
    for (idx i = 0; i < N; ++i) {
        realT n = randn(mean, sigma);
        average += n;
    }
    average /= N;
    EXPECT_NEAR(mean, average, 2.58 * sigma); // very likely (99%)

    // 10000 elements with mean = 1 and sigma = 2
    mean = 1, sigma = 2;
    N = 10000;
    average = 0;
    for (idx i = 0; i < N; ++i) {
        realT n = randn(mean, sigma);
        average += n;
    }
    average /= N;
    EXPECT_NEAR(mean, average, 2.58 * sigma); // very likely (99%)

    // 10000 elements with mean = -1 and sigma = 1e-5
    mean = 1, sigma = 1e-5;
    N = 10000;
    average = 0;
    for (idx i = 0; i < N; ++i) {
        realT n = randn(mean, sigma);
        average += n;
    }
    average /= N;
    EXPECT_NEAR(mean, average, 2.58 * sigma); // very likely (99%)
}

/// BEGIN template<> inline cmat randn(idx rows, idx cols, realT mean,
///       realT sigma)
TEST(qpp_randn, ComplexMatrix) {
    // 1 x 1 matrix with elements of mean = 0 and sigma = 0.1
    idx NA = 1, NB = 1;
    realT mean = 0, sigma = 0.1;
    cmat A = randn<cmat>(NA, NB, mean, sigma);
    cplx average = 0;
    for (idx i = 0; i < NA * NB; ++i) {
        cplx elem = A.data()[i];
        average += elem;
    }
    average /= static_cast<realT>(NA * NB);
    EXPECT_NEAR(mean, average.real(), 2.58 * sigma); // very likely (99%)
    EXPECT_NEAR(mean, average.imag(), 2.58 * sigma); // very likely (99%)

    // 1 x 10 matrix with elements of mean = -1 and sigma = 0.1
    NA = 1, NB = 10;
    mean = -1, sigma = 0.1;
    A = randn<cmat>(NA, NB, mean, sigma);
    average = 0;
    for (idx i = 0; i < NA * NB; ++i) {
        cplx elem = A.data()[i];
        average += elem;
    }
    average /= static_cast<realT>(NA * NB);
    EXPECT_NEAR(mean, average.real(), 2.58 * sigma); // very likely (99%)
    EXPECT_NEAR(mean, average.imag(), 2.58 * sigma); // very likely (99%)

    // 25 x 15 matrix with elements of mean = 10 and sigma = 1e-5
    NA = 25, NB = 15;
    mean = 10, sigma = 1e-5;
    A = randn<cmat>(NA, NB, mean, sigma);
    average = 0;
    for (idx i = 0; i < NA * NB; ++i) {
        cplx elem = A.data()[i];
        average += elem;
    }
    average /= static_cast<realT>(NA * NB);
    EXPECT_NEAR(mean, average.real(), 2.58 * sigma); // very likely (99%)
    EXPECT_NEAR(mean, average.imag(), 2.58 * sigma); // very likely (99%)
}

/// BEGIN template<> inline rmat randn(idx rows, idx cols, realT mean,
///       realT sigma)
TEST(qpp_randn, RealMatrix) {
    // 1 x 1 matrix with elements of mean = 0 and sigma = 0.1
    idx NA = 1, NB = 1;
    realT mean = 0, sigma = 0.1;
    rmat A = randn<rmat>(NA, NB, mean, sigma);
    realT average = 0;
    for (idx i = 0; i < NA * NB; ++i) {
        realT elem = A.data()[i];
        average += elem;
    }
    average /= static_cast<realT>(NA * NB);
    EXPECT_NEAR(mean, average, 2.58 * sigma); // very likely (99%)

    // 1 x 10 matrix with elements of mean = -1 and sigma = 0.1
    NA = 1, NB = 10;
    mean = -1, sigma = 0.1;
    A = randn<rmat>(NA, NB, mean, sigma);
    average = 0;
    for (idx i = 0; i < NA * NB; ++i) {
        realT elem = A.data()[i];
        average += elem;
    }
    average /= static_cast<realT>(NA * NB);
    EXPECT_NEAR(mean, average, 2.58 * sigma); // very likely (99%)

    // 25 x 15 matrix with elements of mean = 10 and sigma = 1e-5
    NA = 25, NB = 15;
    mean = 10, sigma = 1e-5;
    A = randn<rmat>(NA, NB, mean, sigma);
    average = 0;
    for (idx i = 0; i < NA * NB; ++i) {
        realT elem = A.data()[i];
        average += elem;
    }
    average /= static_cast<realT>(NA * NB);
    EXPECT_NEAR(mean, average, 2.58 * sigma); // very likely (99%)
}

/// BEGIN inline std::vector<idx> randperm(idx N)
TEST(qpp_randperm, AllTests) {
    idx N = 1;
    std::vector<idx> result = randperm(N);
    std::sort(std::begin(result), std::end(result));
    std::vector<idx> expected(N);
    std::iota(std::begin(expected), std::end(expected), 0u);
    EXPECT_TRUE(result == expected);

    N = 2;
    result = randperm(N);
    std::sort(std::begin(result), std::end(result));
    expected.resize(N);
    std::iota(std::begin(expected), std::end(expected), 0u);
    EXPECT_TRUE(result == expected);

    N = 20;
    result = randperm(N);
    expected.resize(N);
    std::iota(std::begin(expected), std::end(expected), 0u);
    EXPECT_FALSE(result == expected); // very very likely
    std::sort(std::begin(result), std::end(result));
    EXPECT_TRUE(result == expected);
}

/// BEGIN inline std::vector<realT> randprob(idx N)
TEST(qpp_randprob, AllTests) {
    idx N = 1;
    std::vector<realT> result = randprob(N);
    EXPECT_EQ(1u, result.size());
    EXPECT_NEAR(1, sum(result), 1e-5);

    N = 2;
    result = randprob(N);
    for (idx i = 0; i < N; ++i) {
        EXPECT_GE(result[i], 0);
    }
    EXPECT_EQ(2u, result.size());
    EXPECT_NEAR(1, sum(result), 1e-5);

    N = 10;
    result = randprob(N);
    for (idx i = 0; i < N; ++i) {
        EXPECT_GE(result[i], 0);
    }
    EXPECT_EQ(10u, result.size());
    EXPECT_NEAR(1, sum(result), 1e-5);
}

/// BEGIN inline cmat randrho(idx D = 2)
TEST(qpp_randrho, AllTests) {
    // D = 1 degenerate case
    idx D = 1;
    cmat rho = randrho(D);
    EXPECT_NEAR(0, norm(rho - adjoint(rho)), 1e-5);
    dyn_col_vect<realT> lambdas = hevals(rho);
    for (idx i = 0; i < static_cast<idx>(lambdas.size()); ++i) {
        EXPECT_GE(lambdas[i], 0);
        EXPECT_LE(lambdas[i], 1);
    }
    EXPECT_NEAR(1, trace(rho).real(), 1e-5);
    EXPECT_NEAR(0, trace(rho).imag(), 1e-5);

    // D = 2
    D = 2;
    rho = randrho(D);
    EXPECT_NEAR(0, norm(rho - adjoint(rho)), 1e-5);
    lambdas = hevals(rho);
    for (idx i = 0; i < static_cast<idx>(lambdas.size()); ++i) {
        EXPECT_GE(lambdas[i], 0);
        EXPECT_LE(lambdas[i], 1);
    }
    EXPECT_NEAR(1, trace(rho).real(), 1e-5);
    EXPECT_NEAR(0, trace(rho).imag(), 1e-5);

    // D = 5
    D = 5;
    rho = randrho(D);
    EXPECT_NEAR(0, norm(rho - adjoint(rho)), 1e-5);
    lambdas = hevals(rho);
    for (idx i = 0; i < static_cast<idx>(lambdas.size()); ++i) {
        EXPECT_GE(lambdas[i], 0);
        EXPECT_LE(lambdas[i], 1);
    }
    EXPECT_NEAR(1, trace(rho).real(), 1e-5);
    EXPECT_NEAR(0, trace(rho).imag(), 1e-5);
}

/// BEGIN inline cmat randU(idx D = 2)
TEST(qpp_randU, AllTests) {
    // D = 1 degenerate case
    idx D = 1;
    cmat U = randU(D);
    EXPECT_NEAR(0, norm(U * adjoint(U) - gt.Id(D)), 1e-5);

    // D = 2
    D = 2;
    U = randU();
    EXPECT_NEAR(0, norm(U * adjoint(U) - gt.Id(D)), 1e-5);

    idx N = 10000; // number of runs, average over random unitaries
    ket avg_state = ket::Zero(D);
    for (idx i = 0; i < N; ++i) {
        avg_state += randU() * st.z0;
    }
    ket expected = ket::Zero(D);
    EXPECT_NEAR(0, norm(expected - avg_state / N), 2e-2);

    // D = 10
    D = 10;
    U = randU(D);
    EXPECT_NEAR(0, norm(U * adjoint(U) - gt.Id(D)), 1e-5);
}

/// BEGIN inline cmat randV(idx Din, idx Dout)
TEST(qpp_randV, AllTests) {
    // Din = 1, Dout = 1 degenerate case
    idx Din = 1, Dout = 1;
    cmat V = randV(Din, Dout);
    EXPECT_NEAR(0, norm(adjoint(V) * V - gt.Id(Din)), 1e-5);

    // Din = 1, Dout = 10 degenerate case
    Din = 1, Dout = 10;
    V = randV(Din, Dout);
    EXPECT_NEAR(0, norm(adjoint(V) * V - gt.Id(Din)), 1e-5);

    // Din = 2, Dout = 2
    Din = 2, Dout = 2;
    V = randV(Din, Dout);
    EXPECT_NEAR(0, norm(adjoint(V) * V - gt.Id(Din)), 1e-5);

    // Din = 3, Dout = 5
    Din = 3, Dout = 5;
    V = randV(Din, Dout);
    EXPECT_NEAR(0, norm(adjoint(V) * V - gt.Id(Din)), 1e-5);
}
