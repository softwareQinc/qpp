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

#include "gtest/gtest.h"
#include "qpp.h"

using namespace qpp;

// Unit testing "input_output.h"

/******************************************************************************/
/// BEGIN template<typename Derived> dyn_mat<typename Derived::Scalar>
///       qpp::load(const std::string& fname)
///
///       template<typename Derived> void qpp::save(
///       const Eigen::MatrixBase<Derived>& A, const std::string& fname)
TEST(qpp_load_save, Matrices) {
    // matrices,complex, real and integer

    // DA = 1, DB = 1 degenerate case
    idx DA = 1, DB = 1;
    cmat A = rand<cmat>(DA, DB);
    dmat B = rand<dmat>(DA, DB);
    dyn_mat<int> C = Eigen::MatrixXi::Random(DA, DB);

    qpp::save(A, "out.tmp");
    cmat loadA = qpp::load<cmat>("out.tmp");
    EXPECT_NEAR(0, norm(loadA - A), 1e-7);

    qpp::save(B, "out.tmp");
    dmat loadB = qpp::load<dmat>("out.tmp");
    EXPECT_NEAR(0, norm(loadB - B), 1e-7);

    qpp::save(C, "out.tmp");
    Eigen::MatrixXi loadC = qpp::load<Eigen::MatrixXi>("out.tmp");
    EXPECT_NEAR(0, norm(loadC - C), 1e-7);

    // DA = 1, DB = 10
    DA = 1, DB = 10;
    A = rand<cmat>(DA, DB);
    B = rand<dmat>(DA, DB);
    C = Eigen::MatrixXi::Random(DA, DB);

    qpp::save(A, "out.tmp");
    loadA = qpp::load<cmat>("out.tmp");
    EXPECT_NEAR(0, norm(loadA - A), 1e-7);

    qpp::save(B, "out.tmp");
    loadB = qpp::load<dmat>("out.tmp");
    EXPECT_NEAR(0, norm(loadB - B), 1e-7);

    qpp::save(C, "out.tmp");
    loadC = qpp::load<Eigen::MatrixXi>("out.tmp");
    EXPECT_NEAR(0, norm(loadC - C), 1e-7);

    // DA = 32, DB = 24
    DA = 32, DB = 24;
    A = rand<cmat>(DA, DB);
    B = rand<dmat>(DA, DB);
    C = Eigen::MatrixXi::Random(DA, DB);

    qpp::save(A, "out.tmp");
    loadA = qpp::load<cmat>("out.tmp");
    EXPECT_NEAR(0, norm(loadA - A), 1e-7);

    qpp::save(B, "out.tmp");
    loadB = qpp::load<dmat>("out.tmp");
    EXPECT_NEAR(0, norm(loadB - B), 1e-7);

    qpp::save(C, "out.tmp");
    loadC = qpp::load<Eigen::MatrixXi>("out.tmp");
    EXPECT_NEAR(0, norm(loadC - C), 1e-7);

    // expression
    A = rand<cmat>(5, 5);
    cmat expression = A * A + A;
    qpp::save(A * A + A, "out.tmp");
    cmat load_expression = qpp::load<cmat>("out.tmp");
    EXPECT_NEAR(0, norm(load_expression - expression), 1e-7);
}

TEST(qpp_load_save, Vectors) {
    // kets/row vectors, complex, real and integer

    // D = 1 degenerate case
    idx D = 1;
    ket A = randket(D);
    dyn_row_vect<double> B = dyn_row_vect<double>::Random(D);
    dyn_row_vect<int> C = dyn_row_vect<int>::Random(D);

    qpp::save(A, "out.tmp");
    ket loadA = qpp::load<ket>("out.tmp");
    EXPECT_NEAR(0, norm(loadA - A), 1e-7);

    qpp::save(B, "out.tmp");
    dyn_row_vect<double> loadB = qpp::load<dyn_row_vect<double>>("out.tmp");
    EXPECT_NEAR(0, norm(loadB - B), 1e-7);

    qpp::save(C, "out.tmp");
    dyn_row_vect<int> loadC = qpp::load<dyn_row_vect<int>>("out.tmp");
    EXPECT_NEAR(0, norm(loadC - C), 1e-7);

    // D = 32
    D = 32;
    A = randket(D);
    B = dyn_row_vect<double>::Random(D);
    C = dyn_row_vect<int>::Random(D);

    qpp::save(A, "out.tmp");
    loadA = qpp::load<ket>("out.tmp");
    EXPECT_NEAR(0, norm(loadA - A), 1e-7);

    qpp::save(B, "out.tmp");
    loadB = qpp::load<dyn_row_vect<double>>("out.tmp");
    EXPECT_NEAR(0, norm(loadB - B), 1e-7);

    qpp::save(C, "out.tmp");
    loadC = qpp::load<dyn_row_vect<int>>("out.tmp");
    EXPECT_NEAR(0, norm(loadC - C), 1e-7);

    // expression
    ket expression = 3. * A + A;
    qpp::save(3. * A + A, "out.tmp");
    ket load_expression = qpp::load<ket>("out.tmp");
    EXPECT_NEAR(0, norm(load_expression - expression), 1e-7);
}
/******************************************************************************/
