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

#include <gtest/gtest.h>

#include <qpp.h>
// #include <MATLAB/matlab.h> // support for MATLAB
// #include <experimental/experimental.h> // support for experimental features

using namespace qpp;

// Write your unit tests here. Some examples are provided below.

// ********** qpp::sum() **********
TEST(qpp_sum_test, PositiveNumbers)
{
    std::vector<int> v{0, 1, 2, 3};
    EXPECT_EQ (6, qpp::sum(v.begin(), v.end()));
}

TEST(qpp_sum_test, NegativeNumbers)
{
    std::vector<int> v{0, -1, -2, -3};
    EXPECT_EQ (-6, qpp::sum(v.begin(), v.end()));
}

TEST(qpp_sum_test, MixedNumbers)
{
    std::vector<int> v{ -3, -2, -1, 0, 1, 2};
    EXPECT_EQ (-3, qpp::sum(v.begin(), v.end()));
}
// ********** END qpp::sum() **********

// ********** qpp::prod() **********
TEST(qpp_prod_test, PositiveNumbers)
{
    std::vector<int> v{1, 2, 3, 4};
    EXPECT_EQ (24, qpp::prod(v.begin(), v.end()));
}
// ********** END qpp::prod() **********

// ********** qpp::modpow() **********
TEST(qpp_modpow_test, PositiveNumbers)
{
    EXPECT_EQ (0, qpp::modpow(2, 3, 4));
    EXPECT_EQ (34, qpp::modpow(17, 176, 37));
    EXPECT_EQ (4042, qpp::modpow(178373, 9281623, 6217));
}
// ********** END qpp::modpow() **********

// ********** qpp::applyCTRL() **********
TEST(qpp_applyCTRL, EmptyControl)
{
    std::vector<idx> dims{2, 2, 2, 2};  // 3 qubits
    idx n = dims.size();                // total number of qudits
    idx D = prod(dims);                 // total dimension

    std::vector<idx> ctrl{2, 0};           // where we apply the control
    std::vector<idx> target{1, 3};   // target

    idx Dtarget = 1;                    // dimension of the target subsystems
    for(idx i = 0; i < target.size(); ++i)
    Dtarget *= dims[target[i]];         // compute it here

    // some random n qudit pure state
    ket psi = randket(D);

    cmat rho = psi * adjoint(psi); // the corresponding density matrix
    cmat U = randU(Dtarget);       // some random unitary on the target

    // applyCTRL on pure state
    ket A = applyCTRL(psi, U, ctrl, target, dims);

    // applyCTRL on density matrix
    cmat B = applyCTRL(rho, U, ctrl, target, dims);

    // result when using CTRL-U|psi><psi|CTRL-U^\dagger
    cmat result_psi = A * adjoint(A);
    // result when using CTRL-U(rho)CTRL-U^\dagger
    cmat result_rho = B;

    double res = norm(result_psi - result_rho);
    EXPECT_NEAR (0, res, 1e-8);
}

TEST(qpp_applyCTRL, NonEmptyControl)
{
    std::vector<idx> dims{2, 2, 2, 2};  // 3 qubits
    idx n = dims.size();                // total number of qudits
    idx D = prod(dims);                 // total dimension

    std::vector<idx> ctrl{};           // where we apply the control
    std::vector<idx> target{1, 0, 3};   // target

    idx Dtarget = 1;                    // dimension of the target subsystems
    for(idx i = 0; i < target.size(); ++i)
    Dtarget *= dims[target[i]];         // compute it here

    // some random n qudit pure state
    ket psi = randket(D);

    cmat rho = psi * adjoint(psi); // the corresponding density matrix
    cmat U = randU(Dtarget);       // some random unitary on the target

    // applyCTRL on pure state
    ket A = applyCTRL(psi, U, ctrl, target, dims);

    // applyCTRL on density matrix
    cmat B = applyCTRL(rho, U, ctrl, target, dims);

    // result when using CTRL-U|psi><psi|CTRL-U^\dagger
    cmat result_psi = A * adjoint(A);
    // result when using CTRL-U(rho)CTRL-U^\dagger
    cmat result_rho = B;

    double res = norm(result_psi - result_rho);
    EXPECT_NEAR (0, res, 1e-8);
}
// ********** END qpp::applyCTRL() **********

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
