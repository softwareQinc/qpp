#include <fstream>

#include "gtest/gtest.h"

#include "qpp/qpp.hpp"

using namespace qpp;

// Unit testing "qpp/input_output.hpp"

/// BEGIN template <typename Derived> dyn_mat<typename Derived::Scalar>
///       load(std::istream& is)
///
///       template <typename Derived> void save(
///       const Eigen::MatrixBase<Derived>& A, std::ostream& os)
TEST(qpp_load_save, Matrix) {
    namespace fs = std::filesystem;
    fs::path filename = "mat.txt";
    fs::path mat_txt = fs::temp_directory_path() / filename;

    // matrices,complex, real and integer

    // DA = 1, DB = 1 degenerate case
    idx DA = 1, DB = 1;
    cmat A = rand<cmat>(DA, DB), loadA;
    rmat B = rand<rmat>(DA, DB), loadB;
    dyn_mat<int> C = dyn_mat<int>::Random(DA, DB), loadC;
    {
        std::ofstream fout(mat_txt.string());
        save(A, fout);
    }
    {
        std::ifstream fin(mat_txt.string());
        loadA = load<cmat>(fin);
    }
    EXPECT_NEAR(0, norm(loadA - A), 1e-5);

    {
        std::ofstream fout(mat_txt.string());
        save(B, fout);
    }
    {
        std::ifstream fin(mat_txt.string());
        loadB = load<rmat>(fin);
    }
    EXPECT_NEAR(0, norm(loadB - B), 1e-5);

    {
        std::ofstream fout(mat_txt.string());
        save(C, fout);
    }
    {
        std::ifstream fin(mat_txt.string());
        loadC = load<dyn_mat<int>>(fin);
    }
    EXPECT_NEAR(0, norm(loadC - C), 1e-5);

    // DA = 1, DB = 10
    DA = 1, DB = 10;
    A = rand<cmat>(DA, DB);
    B = rand<rmat>(DA, DB);
    C = dyn_mat<int>::Random(DA, DB);
    {
        std::ofstream fout(mat_txt.string());
        save(A, fout);
    }
    {
        std::ifstream fin(mat_txt.string());
        loadA = load<cmat>(fin);
    }
    EXPECT_NEAR(0, norm(loadA - A), 1e-5);

    {
        std::ofstream fout(mat_txt.string());
        save(B, fout);
    }
    {
        std::ifstream fin(mat_txt.string());
        loadB = load<rmat>(fin);
    }
    EXPECT_NEAR(0, norm(loadB - B), 1e-5);

    {
        std::ofstream fout(mat_txt.string());
        save(C, fout);
    }
    {
        std::ifstream fin(mat_txt.string());
        loadC = load<dyn_mat<int>>(fin);
    }
    EXPECT_NEAR(0, norm(loadC - C), 1e-5);

    // DA = 32, DB = 24
    DA = 32, DB = 24;
    A = rand<cmat>(DA, DB);
    B = rand<rmat>(DA, DB);
    C = dyn_mat<int>::Random(DA, DB);

    {
        std::ofstream fout(mat_txt.string());
        save(A, fout);
    }

    {
        std::ifstream fin(mat_txt.string());
        loadA = load<cmat>(fin);
    }
    EXPECT_NEAR(0, norm(loadA - A), 1e-5);

    {
        std::ofstream fout(mat_txt.string());
        save(B, fout);
    }
    {
        std::ifstream fin(mat_txt.string());
        loadB = load<rmat>(fin);
    }
    EXPECT_NEAR(0, norm(loadB - B), 1e-5);

    {
        std::ofstream fout(mat_txt.string());
        save(C, fout);
    }
    {
        std::ifstream fin(mat_txt.string());
        loadC = load<dyn_mat<int>>(fin);
    }
    EXPECT_NEAR(0, norm(loadC - C), 1e-5);

    // expression
    A = rand<cmat>(5, 5);
    cmat expression = A * A + A, load_expression;
    {
        std::ofstream fout(mat_txt.string());
        save(A * A + A, fout);
    }
    {
        std::ifstream fin(mat_txt.string());
        load_expression = load<cmat>(fin);
    }
    EXPECT_NEAR(0, norm(load_expression - expression), 1e-5);
}

TEST(qpp_load_save, Vector) {
    namespace fs = std::filesystem;
    fs::path filename = "mat.txt";
    fs::path mat_txt = fs::temp_directory_path() / filename;

    // kets/row vectors, complex, real and integer

    // D = 1 degenerate case
    idx D = 1;
    ket A = randket(D), loadA;
    dyn_row_vect<realT> B = dyn_row_vect<realT>::Random(D), loadB;
    dyn_row_vect<int> C = dyn_row_vect<int>::Random(D), loadC;

    {
        std::ofstream fout(mat_txt.string());
        save(A, fout);
    }
    {
        std::ifstream fin(mat_txt.string());
        loadA = load<ket>(fin);
    }
    EXPECT_NEAR(0, norm(loadA - A), 1e-5);

    {
        std::ofstream fout(mat_txt.string());
        save(B, fout);
    }
    {
        std::ifstream fin(mat_txt.string());
        loadB = load<dyn_row_vect<realT>>(fin);
    }
    EXPECT_NEAR(0, norm(loadB - B), 1e-5);

    {
        std::ofstream fout(mat_txt.string());
        save(C, fout);
    }
    {
        std::ifstream fin(mat_txt.string());
        loadC = load<dyn_row_vect<int>>(fin);
    }
    EXPECT_NEAR(0, norm(loadC - C), 1e-5);

    // D = 32
    D = 32;
    A = randket(D);
    B = dyn_row_vect<realT>::Random(D);
    C = dyn_row_vect<int>::Random(D);

    {
        std::ofstream fout(mat_txt.string());
        save(A, fout);
    }
    {
        std::ifstream fin(mat_txt.string());
        loadA = load<ket>(fin);
    }
    EXPECT_NEAR(0, norm(loadA - A), 1e-5);

    {
        std::ofstream fout(mat_txt.string());
        save(B, fout);
    }
    {
        std::ifstream fin(mat_txt.string());
        loadB = load<dyn_row_vect<realT>>(fin);
    }
    EXPECT_NEAR(0, norm(loadB - B), 1e-5);

    {
        std::ofstream fout(mat_txt.string());
        save(C, fout);
    }
    {
        std::ifstream fin(mat_txt.string());
        loadC = load<dyn_row_vect<int>>(fin);
    }
    EXPECT_NEAR(0, norm(loadC - C), 1e-5);

    // expression
    ket expression = 3. * A + A, load_expression;
    {
        std::ofstream fout(mat_txt.string());
        save(3. * A + A, fout);
    }
    {
        std::ifstream fin(mat_txt.string());
        load_expression = load<ket>(fin);
    }
    EXPECT_NEAR(0, norm(load_expression - expression), 1e-5);
}
