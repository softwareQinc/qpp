#include "gtest/gtest.h"

#include "qpp/qpp.hpp"

#include "qpp/MATLAB/matlab.hpp"

using namespace qpp;

// Unit testing "qpp/MATLAB/matlab.hpp"

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
    namespace fs = std::filesystem;
    fs::path filename = "mat.bin";
    fs::path mat_bin = fs::temp_directory_path() / filename;

    // matrices, complex, real and integer

    // DA = 1, DB = 1 degenerate case
    idx DA = 1, DB = 1;
    cmat A = rand<cmat>(DA, DB);
    rmat B = rand<rmat>(DA, DB);
    dyn_mat<int> C = Eigen::MatrixXi::Random(DA, DB);

    save_MATLAB(A, mat_bin.string(), "A", "w");
    cmat loadA = load_MATLAB<cmat>(mat_bin.string(), "A");
    EXPECT_NEAR(0, norm(loadA - A), 1e-5);

    save_MATLAB(B, mat_bin.string(), "B", "w");
    rmat loadB = load_MATLAB<rmat>(mat_bin.string(), "B");
    EXPECT_NEAR(0, norm(loadB - B), 1e-5);

    save_MATLAB(C, mat_bin.string(), "C", "w");
    Eigen::MatrixXi loadC = load_MATLAB<Eigen::MatrixXi>(mat_bin.string(), "C");
    EXPECT_NEAR(0, norm(loadC - C), 1e-5);

    // DA = 1, DB = 10
    DA = 1, DB = 10;
    A = rand<cmat>(DA, DB);
    B = rand<rmat>(DA, DB);
    C = Eigen::MatrixXi::Random(DA, DB);

    save_MATLAB(A, mat_bin.string(), "A", "w");
    loadA = load_MATLAB<cmat>(mat_bin.string(), "A");
    EXPECT_NEAR(0, norm(loadA - A), 1e-5);

    save_MATLAB(B, mat_bin.string(), "B", "w");
    loadB = load_MATLAB<rmat>(mat_bin.string(), "B");
    EXPECT_NEAR(0, norm(loadB - B), 1e-5);

    save_MATLAB(C, mat_bin.string(), "C", "w");
    loadC = load_MATLAB<Eigen::MatrixXi>(mat_bin.string(), "C");
    EXPECT_NEAR(0, norm(loadC - C), 1e-5);

    // DA = 32, DB = 24
    DA = 32, DB = 24;
    A = rand<cmat>(DA, DB);
    B = rand<rmat>(DA, DB);
    C = Eigen::MatrixXi::Random(DA, DB);

    save_MATLAB(A, mat_bin.string(), "A", "w");
    loadA = load_MATLAB<cmat>(mat_bin.string(), "A");
    EXPECT_NEAR(0, norm(loadA - A), 1e-5);

    save_MATLAB(B, mat_bin.string(), "B", "w");
    loadB = load_MATLAB<rmat>(mat_bin.string(), "B");
    EXPECT_NEAR(0, norm(loadB - B), 1e-5);

    save_MATLAB(C, mat_bin.string(), "C", "w");
    loadC = load_MATLAB<Eigen::MatrixXi>(mat_bin.string(), "C");
    EXPECT_NEAR(0, norm(loadC - C), 1e-5);

    // expression
    A = rand<cmat>(5, 5);
    cmat expression = A * A + A;
    save_MATLAB(A * A + A, mat_bin.string(), "expression", "w");
    cmat load_expression = load_MATLAB<cmat>(mat_bin.string(), "expression");
    EXPECT_NEAR(0, norm(load_expression - expression), 1e-5);
}

TEST(qpp_MATLAB_load_save_MATLAB, Vector) {
    namespace fs = std::filesystem;
    fs::path filename = "mat.bin";
    fs::path mat_bin = fs::temp_directory_path() / filename;

    // kets/row vectors, complex, real and integer

    // D = 1 degenerate case
    idx D = 1;
    ket A = randket(D);
    dyn_row_vect<double> B = dyn_row_vect<double>::Random(D);
    dyn_row_vect<int> C = dyn_row_vect<int>::Random(D);

    save_MATLAB(A, mat_bin.string(), "A", "w");
    ket loadA = load_MATLAB<ket>(mat_bin.string(), "A");
    EXPECT_NEAR(0, norm(loadA - A), 1e-5);

    save_MATLAB(B, mat_bin.string(), "B", "w");
    dyn_row_vect<double> loadB =
        load_MATLAB<dyn_row_vect<double>>(mat_bin.string(), "B");
    EXPECT_NEAR(0, norm(loadB - B), 1e-5);

    save_MATLAB(C, mat_bin.string(), "C", "w");
    dyn_row_vect<int> loadC =
        load_MATLAB<dyn_row_vect<int>>(mat_bin.string(), "C");
    EXPECT_NEAR(0, norm(loadC - C), 1e-5);

    // D = 32
    D = 32;
    A = randket(D);
    B = dyn_row_vect<double>::Random(D);
    C = dyn_row_vect<int>::Random(D);

    save_MATLAB(A, mat_bin.string(), "A", "w");
    loadA = load_MATLAB<ket>(mat_bin.string(), "A");
    EXPECT_NEAR(0, norm(loadA - A), 1e-5);

    save_MATLAB(B, mat_bin.string(), "B", "w");
    loadB = load_MATLAB<dyn_row_vect<double>>(mat_bin.string(), "B");
    EXPECT_NEAR(0, norm(loadB - B), 1e-5);

    save_MATLAB(C, mat_bin.string(), "C", "w");
    loadC = load_MATLAB<dyn_row_vect<int>>(mat_bin.string(), "C");
    EXPECT_NEAR(0, norm(loadC - C), 1e-5);

    // expression
    dyn_row_vect<int> expression = 3 * C + C;
    save_MATLAB(3 * C + C, mat_bin.string(), "expression", "w");
    dyn_row_vect<int> load_expression =
        load_MATLAB<dyn_row_vect<int>>(mat_bin.string(), "expression");
    EXPECT_NEAR(0, norm(load_expression - expression), 1e-5);
}
