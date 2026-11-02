#include <filesystem>

#include "gtest/gtest.h"

#include "qpp/qpp.hpp"

#include "qpp/MATLAB/matlab.hpp"

using namespace qpp;

// Unit testing "qpp/MATLAB/matlab.hpp"

TEST(qpp_MATLAB_load_save_MATLAB, Matrix) {
    namespace fs = std::filesystem;
    fs::path mat_bin = fs::temp_directory_path() / "qpp_matlab_matrix.mat";

    // Matrices, complex, real and integer
    auto check_round_trip = [&](idx DA, idx DB) {
        SCOPED_TRACE(::testing::Message() << "DA = " << DA << ", DB = " << DB);

        cmat A = rand<cmat>(DA, DB);
        rmat B = rand<rmat>(DA, DB);
        dyn_mat<int> C(DA, DB);
        for (idx i = 0; i < DA; ++i) {
            for (idx j = 0; j < DB; ++j) {
                C(i, j) = static_cast<int>((i + j) % 21) - 10;
            }
        }

        save_MATLAB(A, mat_bin.string(), "A", "w");
        cmat loadA = load_MATLAB<cmat>(mat_bin.string(), "A");
        ASSERT_EQ(loadA.rows(), A.rows());
        ASSERT_EQ(loadA.cols(), A.cols());
        EXPECT_NEAR(0, norm(loadA - A), 1e-5);

        save_MATLAB(B, mat_bin.string(), "B", "w");
        rmat loadB = load_MATLAB<rmat>(mat_bin.string(), "B");
        ASSERT_EQ(loadB.rows(), B.rows());
        ASSERT_EQ(loadB.cols(), B.cols());
        EXPECT_NEAR(0, norm(loadB - B), 1e-5);

        save_MATLAB(C, mat_bin.string(), "C", "w");
        dyn_mat<int> loadC = load_MATLAB<dyn_mat<int>>(mat_bin.string(), "C");
        ASSERT_EQ(loadC.rows(), C.rows());
        ASSERT_EQ(loadC.cols(), C.cols());
        EXPECT_TRUE((loadC.array() == C.array()).all());
    };

    // DA = 1, DB = 1 degenerate case
    check_round_trip(1, 1);

    // DA = 1, DB = 10
    check_round_trip(1, 10);

    // DA = 32, DB = 24
    check_round_trip(32, 24);

    // Expression
    cmat A = rand<cmat>(5, 5);
    cmat expression = A * A + A;
    save_MATLAB(A * A + A, mat_bin.string(), "expression", "w");
    cmat load_expression = load_MATLAB<cmat>(mat_bin.string(), "expression");
    ASSERT_EQ(load_expression.rows(), expression.rows());
    ASSERT_EQ(load_expression.cols(), expression.cols());
    EXPECT_NEAR(0, norm(load_expression - expression), 1e-5);

    fs::remove(mat_bin);
}

TEST(qpp_MATLAB_load_save_MATLAB, Vector) {
    namespace fs = std::filesystem;
    fs::path mat_bin = fs::temp_directory_path() / "qpp_matlab_vector.mat";

    // Kets/row vectors, complex, real and integer
    auto check_round_trip = [&](idx D) {
        SCOPED_TRACE(::testing::Message() << "D = " << D);

        ket A = randket(D);
        dyn_row_vect<double> B = dyn_row_vect<double>::Random(D);
        dyn_row_vect<int> C(D);
        for (idx i = 0; i < D; ++i) {
            C(i) = static_cast<int>(i % 21) - 10;
        }

        save_MATLAB(A, mat_bin.string(), "A", "w");
        ket loadA = load_MATLAB<ket>(mat_bin.string(), "A");
        ASSERT_EQ(loadA.rows(), A.rows());
        ASSERT_EQ(loadA.cols(), A.cols());
        EXPECT_NEAR(0, norm(loadA - A), 1e-5);

        save_MATLAB(B, mat_bin.string(), "B", "w");
        dyn_row_vect<double> loadB =
            load_MATLAB<dyn_row_vect<double>>(mat_bin.string(), "B");
        ASSERT_EQ(loadB.rows(), B.rows());
        ASSERT_EQ(loadB.cols(), B.cols());
        EXPECT_NEAR(0, norm(loadB - B), 1e-5);

        save_MATLAB(C, mat_bin.string(), "C", "w");
        dyn_row_vect<int> loadC =
            load_MATLAB<dyn_row_vect<int>>(mat_bin.string(), "C");
        ASSERT_EQ(loadC.rows(), C.rows());
        ASSERT_EQ(loadC.cols(), C.cols());
        EXPECT_TRUE((loadC.array() == C.array()).all());
    };

    // D = 1 degenerate case
    check_round_trip(1);

    // D = 32
    check_round_trip(32);

    // Expression: bounded integer values prevent signed overflow
    dyn_row_vect<int> C(32);
    for (idx i = 0; i < 32; ++i) {
        C(i) = static_cast<int>(i) - 16;
    }
    dyn_row_vect<int> expression = 3 * C + C;
    save_MATLAB(3 * C + C, mat_bin.string(), "expression", "w");
    dyn_row_vect<int> load_expression =
        load_MATLAB<dyn_row_vect<int>>(mat_bin.string(), "expression");
    ASSERT_EQ(load_expression.rows(), expression.rows());
    ASSERT_EQ(load_expression.cols(), expression.cols());
    EXPECT_TRUE((load_expression.array() == expression.array()).all());

    fs::remove(mat_bin);
}
