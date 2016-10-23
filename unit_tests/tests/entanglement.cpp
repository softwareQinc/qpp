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

#include "gtest/gtest.h"
#include "qpp.h"

using namespace qpp;

// Unit testing "entanglement.h"

/******************************************************************************/
/// BEGIN template<typename Derived> double qpp::concurrence(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_concurrence, AllTests)
{
    // random qubit product state
    cmat rho1 = prj(kron(randU(), randU()) * mket({0, 0}));
    EXPECT_NEAR(0, qpp::concurrence(rho1), 1e-7);

    // random maximally entangled 2-qubit state
    cmat rho2 = prj(kron(randU(), randU()) * st.b00);
    EXPECT_NEAR(1, qpp::concurrence(rho2), 1e-7);

    // random 2-qubit state with Schmidt coefficients 0.8 and 0.6
    cmat rho3 = prj(kron(randU(), randU()) *
                    (0.8 * mket({0, 0}) + 0.6 * mket({1, 1})));
    EXPECT_NEAR(2 * 0.8 * 0.6, qpp::concurrence(rho3), 1e-7);
}
/******************************************************************************/
/// BEGIN template<typename Derived> double qpp::entanglement(
///       const Eigen::MatrixBase<Derived>& A, const std::vector<idx>& dims)
TEST(qpp_entanglement, AllTests)
{
    // random qutrit product state
    ket psi1 = kron(randU(3), randU(3)) * mket({0, 0}, 3);
    EXPECT_NEAR(0, qpp::entanglement(psi1, {3, 3}), 1e-7);

    // random 2-qutrit state with Schmidt coefficients 0.36, 0.09 and 0.01
    ket psi2 = kron(randU(3), randU(3)) * (0.6 * mket({0, 0}, 3) +
                                           0.3 * mket({1, 1}, 3) +
                                           0.1 * mket({2, 2}, 3));
    EXPECT_NEAR(-0.36 * std::log2(0.36) - 0.09 * std::log2(0.09) -
                0.01 * std::log2(0.01), qpp::entanglement(psi2, {3, 3}), 1e-7);

    // random maximally entangled 2-qutrit state
    ket psi3 = kron(randU(3), randU(3)) * st.mes(3);
    EXPECT_NEAR(std::log2(3), qpp::entanglement(psi3, {3, 3}), 1e-7);
}
/******************************************************************************/
/// BEGIN template<typename Derived> double qpp::entanglement(
///       const Eigen::MatrixBase<Derived>& A, idx d = 2)
TEST(qpp_entanglement_qubits, AllTests)
{
    // random qubit product state
    ket psi1 = kron(randU(), randU()) * mket({0, 0});
    EXPECT_NEAR(0, qpp::entanglement(psi1), 1e-7);

    // random maximally entangled 2-qubit state
    ket psi2 = kron(randU(), randU()) * st.b00;
    EXPECT_NEAR(1, qpp::entanglement(psi2), 1e-7);

    // random 2-qubit state with Schmidt coefficients 0.8 and 0.6
    ket psi3 = kron(randU(), randU()) *
               (0.8 * mket({0, 0}) + 0.6 * mket({1, 1}));
    EXPECT_NEAR(-0.64 * std::log2(0.64) - 0.36 * std::log2(0.36),
                qpp::entanglement(psi3), 1e-7);
}
/******************************************************************************/
/// BEGIN template<typename Derived> double qpp::gconcurrence(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_gconcurrence, AllTests)
{
    // random qubit product state
    ket psi1 = kron(randU(), randU()) * mket({0, 0});
    EXPECT_NEAR(0, qpp::gconcurrence(psi1), 1e-7);

    // random maximally entangled 2-qubit state
    ket psi2 = kron(randU(), randU()) * st.b00;
    EXPECT_NEAR(1, qpp::gconcurrence(psi2), 1e-7);

    // random 2-qubit state with Schmidt coefficients 0.8 and 0.6
    ket psi3 = kron(randU(), randU()) *
               (0.8 * mket({0, 0}) + 0.6 * mket({1, 1}));
    EXPECT_NEAR(2 * 0.8 * 0.6, qpp::gconcurrence(psi3), 1e-7);

    // random maximally entangled 2-qutrit state
    ket psi4 = kron(randU(3), randU(3)) * st.mes(3);
    EXPECT_NEAR(1, qpp::gconcurrence(psi4), 1e-7);

    // random 2-qutrit state with Schmidt coefficients 0.36, 0.09 and 0.01
    ket psi5 = kron(randU(3), randU(3)) * (0.6 * mket({0, 0}, 3) +
                                           0.3 * mket({1, 1}, 3) +
                                           0.1 * mket({2, 2}, 3));
    EXPECT_NEAR(3 * std::pow(0.36 * 0.09 * 0.01, 1. / 3),
                qpp::gconcurrence(psi5), 1e-7);

    // for qubits, expect the gconcurrence to be the same as concurrence
    ket psi6 = randket(4);
    EXPECT_NEAR(gconcurrence(psi6), concurrence(prj(psi6)), 1e-7);
}
/******************************************************************************/
/// BEGIN template<typename Derived> double qpp::lognegativity(
///       const Eigen::MatrixBase<Derived>& A, const std::vector<idx>& dims)
TEST(qpp_lognegativity, AllTests)
{
    // zero on product states (2 qutrits)
    cmat rho = randrho(3);
    cmat sigma = randrho(3);
    EXPECT_NEAR(0, qpp::lognegativity(kron(rho, sigma), {3, 3}), 1e-7);

    // additivity (2 ququads)
    rho = randrho(4);
    sigma = randrho(4);
    EXPECT_NEAR(qpp::lognegativity(syspermute(kron(rho, sigma), {0, 2, 1, 3}),
                                   {4, 4}),
                qpp::lognegativity(rho, {2, 2}) +
                qpp::lognegativity(sigma, {2, 2}), 1e-7);

    // random maximally entangled 2-qutrit state
    idx d = 3;
    rho = prj(kron(randU(d), randU(d)) * st.mes(d));
    EXPECT_NEAR(std::log2(d), qpp::lognegativity(rho, {d, d}), 1e-7);

    // random maximally entangled 2-ququad state
    d = 4;
    rho = prj(kron(randU(d), randU(d)) * st.mes(d));
    EXPECT_NEAR(std::log2(d), qpp::lognegativity(rho, {d, d}), 1e-7);

    // random maximally entangled state (d = 7)
    d = 7;
    rho = prj(kron(randU(d), randU(d)) * st.mes(d));
    EXPECT_NEAR(std::log2(d), qpp::lognegativity(rho, {d, d}), 1e-7);
}
/******************************************************************************/
/// BEGIN template<typename Derived> double qpp::lognegativity(
///       const Eigen::MatrixBase<Derived>& A, idx d = 2)
TEST(qpp_lognegativity_qubits, AllTests)
{
    // zero on product states (2 qubits)
    cmat rho = randrho();
    cmat sigma = randrho();
    EXPECT_NEAR(0, qpp::lognegativity(kron(rho, sigma)), 1e-7);

    // random 2-qubit state with Schmidt coefficients 0.8 and 0.6
    ket psi = kron(randU(), randU()) *
              (0.8 * mket({0, 0}) + 0.6 * mket({1, 1}));
    rho = prj(psi);
    EXPECT_NEAR(0.9708536, qpp::lognegativity(rho), 1e-7);

    // random maximally entangled 2-qubit state
    idx d = 2;
    rho = prj(kron(randU(d), randU(d)) * st.mes(d));
    EXPECT_NEAR(std::log2(d), qpp::lognegativity(rho), 1e-7);
}
/******************************************************************************/
/// BEGIN template<typename Derived> double qpp::negativity(
///       const Eigen::MatrixBase<Derived>& A, const std::vector<idx>& dims)
TEST(qpp_negativity, AllTests)
{
    // Must be (d - 1)/2 on MES

    // zero on product states (2 qutrits)
    cmat rho = randrho(3);
    cmat sigma = randrho(3);
    EXPECT_NEAR(0, qpp::negativity(kron(rho, sigma), {3, 3}), 1e-7);

    // convexity (10 ququads)
    idx N = 10;
    idx d = 4;
    std::vector<cmat> rhos(N);
    auto probs = randprob(N);
    rho = cmat::Zero(d * d, d * d);
    double sum_neg = 0;
    for (idx i = 0; i < N; ++i)
    {
        cmat rho_i = randrho(d * d);
        sum_neg += probs[i] * qpp::negativity(rho_i, {d, d});
        rho += probs[i] * rho_i;
    }
    EXPECT_LE(qpp::negativity(rho, {d, d}), sum_neg);

    // random maximally entangled 2-qutrit state
    d = 3;
    rho = prj(kron(randU(d), randU(d)) * st.mes(d));
    EXPECT_NEAR((d - 1) / 2., qpp::negativity(rho, {d, d}), 1e-7);

    // random maximally entangled 2-ququad state
    d = 4;
    rho = prj(kron(randU(d), randU(d)) * st.mes(d));
    EXPECT_NEAR((d - 1) / 2., qpp::negativity(rho, {d, d}), 1e-7);

    // random maximally entangled state (d = 7)
    d = 7;
    rho = prj(kron(randU(d), randU(d)) * st.mes(d));
    EXPECT_NEAR((d - 1) / 2., qpp::negativity(rho, {d, d}), 1e-7);
}
/******************************************************************************/
/// BEGIN template<typename Derived> double qpp::negativity(
///       const Eigen::MatrixBase<Derived>& A, idx d = 2)
TEST(qpp_negativity_qubits, AllTests)
{
    // Must be (d - 1)/2 on MES

    // zero on product states (2 qutrits)
    cmat rho = randrho();
    cmat sigma = randrho();
    EXPECT_NEAR(0, qpp::negativity(kron(rho, sigma)), 1e-7);

    // random 2-qubit state with Schmidt coefficients 0.8 and 0.6
    ket psi = kron(randU(), randU()) *
              (0.8 * mket({0, 0}) + 0.6 * mket({1, 1}));
    rho = prj(psi);
    EXPECT_NEAR(0.48, qpp::negativity(rho), 1e-7);

    // random maximally entangled 2-qubit state
    psi = kron(randU(), randU()) * st.mes(2);
    rho = prj(psi);
    EXPECT_NEAR(0.5, qpp::negativity(rho), 1e-7);
}
/******************************************************************************/
/// BEGIN template<typename Derived> cmat qpp::schmidtA(
///       const Eigen::MatrixBase<Derived>& A, const std::vector<idx>& dims)
TEST(qpp_schmidtA, AllTests)
{

}
/******************************************************************************/
/// BEGIN template<typename Derived> cmat qpp::schmidtA(
///       const Eigen::MatrixBase<Derived>& A, idx d = 2)
TEST(qpp_schmidtA_qubits, AllTests)
{

}
/******************************************************************************/
/// BEGIN template<typename Derived> cmat qpp::schmidtB(
///       const Eigen::MatrixBase<Derived>& A, const std::vector<idx>& dims)
TEST(qpp_schmidtB, AllTests)
{

}
/******************************************************************************/
/// BEGIN template<typename Derived> cmat qpp::schmidtB(
///       const Eigen::MatrixBase<Derived>& A, idx d = 2)
TEST(qpp_schmidtB_qubits, AllTests)
{

}
/******************************************************************************/
/// BEGIN template<typename Derived> dyn_col_vect<double> qpp::schmidtcoeffs(
///       const Eigen::MatrixBase<Derived>& A, const std::vector<idx>& dims)
TEST(qpp_schmidtcoeffs, AllTests)
{

}
/******************************************************************************/
/// BEGIN template<typename Derived> dyn_col_vect<double> qpp::schmidtcoeffs(
///       const Eigen::MatrixBase<Derived>& A, idx d = 2)
TEST(qpp_schmidtcoeffs_qubits, AllTests)
{

}
/******************************************************************************/
/// BEGIN template<typename Derived> std::vector<double> qpp::schmidtprobs(
///       const Eigen::MatrixBase<Derived>& A, const std::vector<idx>& dims)
TEST(qpp_schmidtprobs, AllTests)
{

}
/******************************************************************************/
/// BEGIN template<typename Derived> std::vector<double> qpp::schmidtprobs(
///       const Eigen::MatrixBase<Derived>& A, idx d = 2)
TEST(qpp_schmidtprobs_qubits, AllTests)
{

}
/******************************************************************************/
