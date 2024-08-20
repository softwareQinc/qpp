#include "gtest/gtest.h"

#include "qpp/qpp.hpp"

#include "qpp/MATLAB/matlab.hpp"

using namespace qpp;

// Unit testing "MATLAB/matlab.hpp"

/// BEGIN template <typename Derived> dyn_mat<typename Derived::Scalar>
///       load_MATLAB(const std::string& mat_file,
///                       const std::string& var_name)
///
///       template <typename Derived> dyn_mat<typename Derived::Scalar>
///       void save_MATLAB(const Eigen::MatrixBase <Derived>& A,
///                        const std::string& mat_file,
///                        const std::string& var_name,
///                        const std::string& mode)
TEST(qpp_MATLAB_load_save_MATLAB, Matrix) {
    // matrices, complex, real and integer

    // DA = 1, DB = 1 degenerate case
    idx DA = 1, DB = 1;
    cmat A = rand<cmat>(DA, DB);
    rmat B = rand<rmat>(DA, DB);
    dyn_mat<int> C = Eigen::MatrixXi::Random(DA, DB);

    save_MATLAB(A, "out.mat", "A", "w");
    cmat loadA = load_MATLAB<cmat>("out.mat", "A");
    EXPECT_NEAR(0, norm(loadA - A), 1e-5);

    save_MATLAB(B, "out.mat", "B", "w");
    rmat loadB = load_MATLAB<rmat>("out.mat", "B");
    EXPECT_NEAR(0, norm(loadB - B), 1e-5);

    save_MATLAB(C, "out.mat", "C", "w");
    Eigen::MatrixXi loadC = load_MATLAB<Eigen::MatrixXi>("out.mat", "C");
    EXPECT_NEAR(0, norm(loadC - C), 1e-5);

    // DA = 1, DB = 10
    DA = 1, DB = 10;
    A = rand<cmat>(DA, DB);
    B = rand<rmat>(DA, DB);
    C = Eigen::MatrixXi::Random(DA, DB);

    save_MATLAB(A, "out.mat", "A", "w");
    loadA = load_MATLAB<cmat>("out.mat", "A");
    EXPECT_NEAR(0, norm(loadA - A), 1e-5);

    save_MATLAB(B, "out.mat", "B", "w");
    loadB = load_MATLAB<rmat>("out.mat", "B");
    EXPECT_NEAR(0, norm(loadB - B), 1e-5);

    save_MATLAB(C, "out.mat", "C", "w");
    loadC = load_MATLAB<Eigen::MatrixXi>("out.mat", "C");
    EXPECT_NEAR(0, norm(loadC - C), 1e-5);

    // DA = 32, DB = 24
    DA = 32, DB = 24;
    A = rand<cmat>(DA, DB);
    B = rand<rmat>(DA, DB);
    C = Eigen::MatrixXi::Random(DA, DB);

    save_MATLAB(A, "out.mat", "A", "w");
    loadA = load_MATLAB<cmat>("out.mat", "A");
    EXPECT_NEAR(0, norm(loadA - A), 1e-5);

    save_MATLAB(B, "out.mat", "B", "w");
    loadB = load_MATLAB<rmat>("out.mat", "B");
    EXPECT_NEAR(0, norm(loadB - B), 1e-5);

    save_MATLAB(C, "out.mat", "C", "w");
    loadC = load_MATLAB<Eigen::MatrixXi>("out.mat", "C");
    EXPECT_NEAR(0, norm(loadC - C), 1e-5);

    // expression
    A = rand<cmat>(5, 5);
    cmat expression = A * A + A;
    save_MATLAB(A * A + A, "out.mat", "expression", "w");
    cmat load_expression = load_MATLAB<cmat>("out.mat", "expression");
    EXPECT_NEAR(0, norm(load_expression - expression), 1e-5);
}

TEST(qpp_MATLAB_load_save_MATLAB, Vector) {
    // kets/row vectors, complex, real and integer

    // D = 1 degenerate case
    idx D = 1;
    ket A = randket(D);
    dyn_row_vect<double> B = dyn_row_vect<double>::Random(D);
    dyn_row_vect<int> C = dyn_row_vect<int>::Random(D);

    save_MATLAB(A, "out.mat", "A", "w");
    ket loadA = load_MATLAB<ket>("out.mat", "A");
    EXPECT_NEAR(0, norm(loadA - A), 1e-5);

    save_MATLAB(B, "out.mat", "B", "w");
    dyn_row_vect<double> loadB =
        load_MATLAB<dyn_row_vect<double>>("out.mat", "B");
    EXPECT_NEAR(0, norm(loadB - B), 1e-5);

    save_MATLAB(C, "out.mat", "C", "w");
    dyn_row_vect<int> loadC = load_MATLAB<dyn_row_vect<int>>("out.mat", "C");
    EXPECT_NEAR(0, norm(loadC - C), 1e-5);

    // D = 32
    D = 32;
    A = randket(D);
    B = dyn_row_vect<double>::Random(D);
    C = dyn_row_vect<int>::Random(D);

    save_MATLAB(A, "out.mat", "A", "w");
    loadA = load_MATLAB<ket>("out.mat", "A");
    EXPECT_NEAR(0, norm(loadA - A), 1e-5);

    save_MATLAB(B, "out.mat", "B", "w");
    loadB = load_MATLAB<dyn_row_vect<double>>("out.mat", "B");
    EXPECT_NEAR(0, norm(loadB - B), 1e-5);

    save_MATLAB(C, "out.mat", "C", "w");
    loadC = load_MATLAB<dyn_row_vect<int>>("out.mat", "C");
    EXPECT_NEAR(0, norm(loadC - C), 1e-5);

    // expression
    dyn_row_vect<int> expression = 3 * C + C;
    save_MATLAB(3 * C + C, "out.mat", "expression", "w");
    dyn_row_vect<int> load_expression =
        load_MATLAB<dyn_row_vect<int>>("out.mat", "expression");
    EXPECT_NEAR(0, norm(load_expression - expression), 1e-5);
}
