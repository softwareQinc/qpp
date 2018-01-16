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

#include <cmath>
#include "gtest/gtest.h"
#include "qpp.h"

using namespace qpp;

// Unit testing "classes/states.h"

/******************************************************************************/
/// BEGIN ket qpp::States::jn(idx j, idx n, idx d = 2) const
TEST(qpp_States_jn, AllTests) {
    // qubits

    idx n = 1;
    EXPECT_NEAR(0, norm(qpp::st.z0 - qpp::st.jn(0, n)), 1e-7);
    EXPECT_NEAR(0, norm(qpp::st.z1 - qpp::st.jn(1, n)), 1e-7);

    n = 2;
    EXPECT_NEAR(0, norm(kron(st.z0, st.z0) - qpp::st.jn(0, n)), 1e-7);
    EXPECT_NEAR(0, norm(kron(st.z1, st.z1) - qpp::st.jn(1, n)), 1e-7);

    n = 4;
    EXPECT_NEAR(0, norm(kronpow(st.z0, n) - qpp::st.jn(0, n)), 1e-7);
    EXPECT_NEAR(0, norm(kronpow(st.z1, n) - qpp::st.jn(1, n)), 1e-7);

    // qudits
    idx d = 5;

    n = 1;
    EXPECT_NEAR(0, norm(mket({0}, d) - qpp::st.jn(0, n, d)), 1e-7);
    EXPECT_NEAR(0, norm(mket({1}, d) - qpp::st.jn(1, n, d)), 1e-7);
    EXPECT_NEAR(0, norm(mket({3}, d) - qpp::st.jn(3, n, d)), 1e-7);

    n = 4;
    EXPECT_NEAR(0, norm(kronpow(mket({0}, d), n) - qpp::st.jn(0, n, d)), 1e-7);
    EXPECT_NEAR(0, norm(kronpow(mket({1}, d), n) - qpp::st.jn(1, n, d)), 1e-7);
    EXPECT_NEAR(0, norm(kronpow(mket({4}, d), n) - qpp::st.jn(4, n, d)), 1e-7);
}
/******************************************************************************/
/// BEGIN ket qpp::States::mes(idx d = 2) const
TEST(qpp_States_mes, AllTests) {
    // d = 1 (number)
    ket ket_one(1);
    ket_one << 1;
    EXPECT_NEAR(0, norm(qpp::st.mes(1) - ket_one), 1e-7);

    // qubits
    EXPECT_NEAR(0, norm(qpp::st.mes() - st.b00), 1e-7);

    // qutrits
    idx d = 3;
    ket psi = mket({0, 0}, {d, d}) / std::sqrt(d);
    for (idx i = 1; i < d; ++i) {
        psi += mket({i, i}, {d, d});
    }
    ket mes_qutrit(d * d);
    mes_qutrit << 1, 0, 0, 0, 1, 0, 0, 0, 1;
    mes_qutrit /= std::sqrt(d);
    EXPECT_NEAR(0, norm(qpp::st.mes(d) - mes_qutrit), 1e-7);

    // ququads
    d = 4;
    ket mes_ququad(d * d);
    mes_ququad << 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1;
    mes_ququad /= std::sqrt(d);
    EXPECT_NEAR(0, norm(qpp::st.mes(d) - mes_ququad), 1e-7);
}
/******************************************************************************/
/// BEGIN ket qpp::States::minus(idx n) const
TEST(qpp_States_minus, AllTests) {
    idx n = 1;
    EXPECT_NEAR(0, norm(qpp::st.x1 - qpp::st.minus(n)), 1e-7);

    n = 2;
    EXPECT_NEAR(0, norm(kron(gt.H, gt.H) * mket({1, 1}) - qpp::st.minus(n)),
                1e-7);

    n = 4;
    EXPECT_NEAR(0, norm(kronpow(gt.H, n) * mket(std::vector<idx>(n, 1)) -
                        qpp::st.minus(n)),
                1e-7);
}
/******************************************************************************/
/// BEGIN ket qpp::States::one(idx n, idx d = 2) const
TEST(qpp_States_one, AllTests) {
    // qubits

    idx n = 1;
    EXPECT_NEAR(0, norm(qpp::st.z1 - qpp::st.one(n)), 1e-7);

    n = 2;
    EXPECT_NEAR(0, norm(kron(st.z1, st.z1) - qpp::st.one(n)), 1e-7);

    n = 4;
    EXPECT_NEAR(0, norm(kronpow(st.z1, n) - qpp::st.one(n)), 1e-7);

    // qudits
    idx d = 5;

    n = 1;
    EXPECT_NEAR(0, norm(mket({1}, d) - qpp::st.one(n, d)), 1e-7);

    n = 4;
    EXPECT_NEAR(0, norm(kronpow(mket({1}, d), n) - qpp::st.one(n, d)), 1e-7);
}
/******************************************************************************/
/// BEGIN ket qpp::States::plus(idx n) const
TEST(qpp_States_plus, AllTests) {
    idx n = 1;
    EXPECT_NEAR(0, norm(qpp::st.x0 - qpp::st.plus(n)), 1e-7);

    n = 2;
    EXPECT_NEAR(0, norm(kron(gt.H, gt.H) * mket({0, 0}) - qpp::st.plus(n)),
                1e-7);

    n = 4;
    EXPECT_NEAR(0, norm(kronpow(gt.H, n) * mket(std::vector<idx>(n, 0)) -
                        qpp::st.plus(n)),
                1e-7);
}
/******************************************************************************/
/// BEGIN ket qpp::States::zero(idx n, idx d = 2) const
TEST(qpp_States_zero, AllTests) {
    // qubits

    idx n = 1;
    EXPECT_NEAR(0, norm(qpp::st.z0 - qpp::st.zero(n)), 1e-7);

    n = 2;
    EXPECT_NEAR(0, norm(kron(st.z0, st.z0) - qpp::st.zero(n)), 1e-7);

    n = 4;
    EXPECT_NEAR(0, norm(kronpow(st.z0, n) - qpp::st.zero(n)), 1e-7);

    // qudits
    idx d = 5;

    n = 1;
    EXPECT_NEAR(0, norm(mket({0}, d) - qpp::st.zero(n, d)), 1e-7);

    n = 4;
    EXPECT_NEAR(0, norm(kronpow(mket({0}, d), n) - qpp::st.zero(n, d)), 1e-7);
}
/******************************************************************************/
