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

#ifdef WITH_MATLAB_

#include "gtest/gtest.h"
#include "qpp.h"

#include "MATLAB/matlab.h"

using namespace qpp;

// Unit testing "MATLAB/matlab.h"

/******************************************************************************/
/// BEGIN template<typename Derived> dyn_mat<typename Derived::Scalar>
///       qpp::loadMATLAB(const std::string& mat_file,
//                        const std::string& var_name)
///
///       template<typename Derived> dyn_mat<typename Derived::Scalar>
///       void qpp::saveMATLAB(const Eigen::MatrixBase <Derived>& A,
///                            const std::string& mat_file,
///                            const std::string& var_name,
///                       const std::string& mode)
TEST(qpp_MATLAB_load_save_MATLAB_Matrices, AllTests) {
    // matrices,complex, real and integer

    // DA = 1, DB = 1 degenerate case
    idx DA = 1, DB = 1;
    cmat A = rand<cmat>(DA, DB);
    dmat B = rand<dmat>(DA, DB);
    dyn_mat<int> C = Eigen::MatrixXi::Random(DA, DB);

    qpp::saveMATLAB(A, "out.mat", "A", "w");
    cmat loadA = qpp::loadMATLAB<cmat>("out.mat", "A");
    EXPECT_NEAR(0, norm(loadA - A), 1e-7);

    qpp::saveMATLAB(B, "out.mat", "B", "w");
    dmat loadB = qpp::loadMATLAB<dmat>("out.mat", "B");
    EXPECT_NEAR(0, norm(loadB - B), 1e-7);

    qpp::saveMATLAB(C, "out.mat", "C", "w");
    Eigen::MatrixXi loadC = qpp::loadMATLAB<Eigen::MatrixXi>("out.mat", "C");
    EXPECT_NEAR(0, norm(loadC - C), 1e-7);

    // DA = 1, DB = 10
    DA = 1, DB = 10;
    A = rand<cmat>(DA, DB);
    B = rand<dmat>(DA, DB);
    C = Eigen::MatrixXi::Random(DA, DB);

    qpp::saveMATLAB(A, "out.mat", "A", "w");
    loadA = qpp::loadMATLAB<cmat>("out.mat", "A");
    EXPECT_NEAR(0, norm(loadA - A), 1e-7);

    qpp::saveMATLAB(B, "out.mat", "B", "w");
    loadB = qpp::loadMATLAB<dmat>("out.mat", "B");
    EXPECT_NEAR(0, norm(loadB - B), 1e-7);

    qpp::saveMATLAB(C, "out.mat", "C", "w");
    loadC = qpp::loadMATLAB<Eigen::MatrixXi>("out.mat", "C");
    EXPECT_NEAR(0, norm(loadC - C), 1e-7);

    // DA = 32, DB = 24
    DA = 32, DB = 24;
    A = rand<cmat>(DA, DB);
    B = rand<dmat>(DA, DB);
    C = Eigen::MatrixXi::Random(DA, DB);

    qpp::saveMATLAB(A, "out.mat", "A", "w");
    loadA = qpp::loadMATLAB<cmat>("out.mat", "A");
    EXPECT_NEAR(0, norm(loadA - A), 1e-7);

    qpp::saveMATLAB(B, "out.mat", "B", "w");
    loadB = qpp::loadMATLAB<dmat>("out.mat", "B");
    EXPECT_NEAR(0, norm(loadB - B), 1e-7);

    qpp::saveMATLAB(C, "out.mat", "C", "w");
    loadC = qpp::loadMATLAB<Eigen::MatrixXi>("out.mat", "C");
    EXPECT_NEAR(0, norm(loadC - C), 1e-7);

    // expression
    A = rand<cmat>(5, 5);
    cmat expression = A * A + A;
    qpp::saveMATLAB(A * A + A, "out.mat", "expression", "w");
    cmat load_expression = qpp::loadMATLAB<cmat>("out.mat", "expression");
    EXPECT_NEAR(0, norm(load_expression - expression), 1e-7);
}

TEST(qpp_MATLAB_load_save_MATLAB_Vectors, AllTests) {
    // kets/row vectors, complex, real and integer

    // D = 1 degenerate case
    idx D = 1;
    ket A = randket(D);
    dyn_row_vect<double> B = dyn_row_vect<double>::Random(D);
    dyn_row_vect<int> C = dyn_row_vect<int>::Random(D);

    qpp::saveMATLAB(A, "out.mat", "A", "w");
    ket loadA = qpp::loadMATLAB<ket>("out.mat", "A");
    EXPECT_NEAR(0, norm(loadA - A), 1e-7);

    qpp::saveMATLAB(B, "out.mat", "B", "w");
    dyn_row_vect<double> loadB =
        qpp::loadMATLAB<dyn_row_vect<double>>("out.mat", "B");
    EXPECT_NEAR(0, norm(loadB - B), 1e-7);

    qpp::saveMATLAB(C, "out.mat", "C", "w");
    dyn_row_vect<int> loadC =
        qpp::loadMATLAB<dyn_row_vect<int>>("out.mat", "C");
    EXPECT_NEAR(0, norm(loadC - C), 1e-7);

    // D = 32
    D = 32;
    A = randket(D);
    B = dyn_row_vect<double>::Random(D);
    C = dyn_row_vect<int>::Random(D);

    qpp::saveMATLAB(A, "out.mat", "A", "w");
    loadA = qpp::loadMATLAB<ket>("out.mat", "A");
    EXPECT_NEAR(0, norm(loadA - A), 1e-7);

    qpp::saveMATLAB(B, "out.mat", "B", "w");
    loadB = qpp::loadMATLAB<dyn_row_vect<double>>("out.mat", "B");
    EXPECT_NEAR(0, norm(loadB - B), 1e-7);

    qpp::saveMATLAB(C, "out.mat", "C", "w");
    loadC = qpp::loadMATLAB<dyn_row_vect<int>>("out.mat", "C");
    EXPECT_NEAR(0, norm(loadC - C), 1e-7);

    // expression
    dyn_row_vect<int> expression = 3 * C + C;
    qpp::saveMATLAB(3 * C + C, "out.mat", "expression", "w");
    dyn_row_vect<int> load_expression =
        qpp::loadMATLAB<dyn_row_vect<int>>("out.mat", "expression");
    EXPECT_NEAR(0, norm(load_expression - expression), 1e-7);
}
/******************************************************************************/

#endif // WITH_MATLAB
