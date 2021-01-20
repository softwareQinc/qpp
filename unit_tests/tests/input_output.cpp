#include <fstream>
#include "gtest/gtest.h"
#include "qpp.h"

using namespace qpp;

// Unit testing "input_output.hpp"

#if (0)
/******************************************************************************/
/// BEGIN template <typename Derived> dyn_mat<typename Derived::Scalar>
///       qpp::load(std::istream& is)
///
///       template <typename Derived> void qpp::save(
///       const Eigen::MatrixBase<Derived>& A, std::ostream& os)
TEST(qpp_load_save, Matrix) {
    // matrices,complex, real and integer
    // DA = 1, DB = 1 degenerate case
    idx DA = 1, DB = 1;
    cmat A = rand<cmat>(DA, DB), loadA;
    dmat B = rand<dmat>(DA, DB), loadB;
    dyn_mat<int> C = dyn_mat<int>::Random(DA, DB), loadC;
    {
        std::ofstream fout("mat.dat", std::ios::out | std::ios::binary);
        qpp::save(A, fout);
    }
    {
        std::ifstream fin("mat.dat", std::ios::in | std::ios::binary);
        loadA = qpp::load<cmat>(fin);
    }
    EXPECT_NEAR(0, norm(loadA - A), 1e-7);

    {
        std::ofstream fout("mat.dat", std::ios::out | std::ios::binary);
        qpp::save(B, fout);
    }
    {
        std::ifstream fin("mat.dat", std::ios::in | std::ios::binary);
        loadB = qpp::load<dmat>(fin);
    }
    EXPECT_NEAR(0, norm(loadB - B), 1e-7);

    {
        std::ofstream fout("mat.dat", std::ios::out | std::ios::binary);
        qpp::save(C, fout);
    }
    {
        std::ifstream fin("mat.dat", std::ios::in | std::ios::binary);
        loadC = qpp::load<dyn_mat<int>>(fin);
    }
    EXPECT_NEAR(0, norm(loadC - C), 1e-7);

    // DA = 1, DB = 10
    DA = 1, DB = 10;
    A = rand<cmat>(DA, DB);
    B = rand<dmat>(DA, DB);
    C = dyn_mat<int>::Random(DA, DB);
    {
        std::ofstream fout("mat.dat", std::ios::out | std::ios::binary);
        qpp::save(A, fout);
    }
    {
        std::ifstream fin("mat.dat", std::ios::in | std::ios::binary);
        loadA = qpp::load<cmat>(fin);
    }
    EXPECT_NEAR(0, norm(loadA - A), 1e-7);

    {
        std::ofstream fout("mat.dat", std::ios::out | std::ios::binary);
        qpp::save(B, fout);
    }
    {
        std::ifstream fin("mat.dat", std::ios::in | std::ios::binary);
        loadB = qpp::load<dmat>(fin);
    }
    EXPECT_NEAR(0, norm(loadB - B), 1e-7);

    {
        std::ofstream fout("mat.dat", std::ios::out | std::ios::binary);
        qpp::save(C, fout);
    }
    {
        std::ifstream fin("mat.dat", std::ios::in | std::ios::binary);
        loadC = qpp::load<dyn_mat<int>>(fin);
    }
    EXPECT_NEAR(0, norm(loadC - C), 1e-7);

    // DA = 32, DB = 24
    DA = 32, DB = 24;
    A = rand<cmat>(DA, DB);
    B = rand<dmat>(DA, DB);
    C = dyn_mat<int>::Random(DA, DB);

    {
        std::ofstream fout("mat.dat", std::ios::out | std::ios::binary);
        qpp::save(A, fout);
    }

    {
        std::ifstream fin("mat.dat", std::ios::in | std::ios::binary);
        loadA = qpp::load<cmat>(fin);
    }
    EXPECT_NEAR(0, norm(loadA - A), 1e-7);

    {
        std::ofstream fout("mat.dat", std::ios::out | std::ios::binary);
        qpp::save(B, fout);
    }
    {
        std::ifstream fin("mat.dat", std::ios::in | std::ios::binary);
        loadB = qpp::load<dmat>(fin);
    }
    EXPECT_NEAR(0, norm(loadB - B), 1e-7);

    {
        std::ofstream fout("mat.dat", std::ios::out | std::ios::binary);
        qpp::save(C, fout);
    }
    {
        std::ifstream fin("mat.dat", std::ios::in | std::ios::binary);
        loadC = qpp::load<dyn_mat<int>>(fin);
    }
    EXPECT_NEAR(0, norm(loadC - C), 1e-7);

    // expression
    A = rand<cmat>(5, 5);
    {
        std::ofstream fout("mat.dat", std::ios::out | std::ios::binary);
        qpp::save(A * A + A, fout);
    }
    cmat load_expression;
    {
        std::ifstream fin("mat.dat", std::ios::in | std::ios::binary);
        load_expression = qpp::load<cmat>(fin);
    }
    EXPECT_NEAR(0, norm(load_expression - (A * A + A)), 1e-7);
}
/******************************************************************************/
TEST(qpp_load_save, Vector) {
    // kets/row vectors, complex, real and integer

    // D = 1 degenerate case
    idx D = 1;
    ket A = randket(D), loadA;
    dyn_row_vect<double> B = dyn_row_vect<double>::Random(D), loadB;
    dyn_row_vect<int> C = dyn_row_vect<int>::Random(D), loadC;

    {
        std::ofstream fout("mat.dat", std::ios::out | std::ios::binary);
        qpp::save(A, fout);
    }
    {
        std::ifstream fin("mat.dat", std::ios::in | std::ios::binary);
        loadA = qpp::load<ket>(fin);
    }
    EXPECT_NEAR(0, norm(loadA - A), 1e-7);

    {
        std::ofstream fout("mat.dat", std::ios::out | std::ios::binary);
        qpp::save(B, fout);
    }
    {
        std::ifstream fin("mat.dat", std::ios::in | std::ios::binary);
        loadB = qpp::load<dyn_row_vect<double>>(fin);
    }
    EXPECT_NEAR(0, norm(loadB - B), 1e-7);

    {
        std::ofstream fout("mat.dat", std::ios::out | std::ios::binary);
        qpp::save(C, fout);
    }
    {
        std::ifstream fin("mat.dat", std::ios::in | std::ios::binary);
        loadC = qpp::load<dyn_row_vect<int>>(fin);
    }
    EXPECT_NEAR(0, norm(loadC - C), 1e-7);

    // D = 32
    D = 32;
    A = randket(D);
    B = dyn_row_vect<double>::Random(D);
    C = dyn_row_vect<int>::Random(D);

    {
        std::ofstream fout("mat.dat", std::ios::out | std::ios::binary);
        qpp::save(A, fout);
    }
    {
        std::ifstream fin("mat.dat", std::ios::in | std::ios::binary);
        loadA = qpp::load<ket>(fin);
    }
    EXPECT_NEAR(0, norm(loadA - A), 1e-7);

    {
        std::ofstream fout("mat.dat", std::ios::out | std::ios::binary);
        qpp::save(B, fout);
    }
    {
        std::ifstream fin("mat.dat", std::ios::in | std::ios::binary);
        loadB = qpp::load<dyn_row_vect<double>>(fin);
    }
    EXPECT_NEAR(0, norm(loadB - B), 1e-7);

    {
        std::ofstream fout("mat.dat", std::ios::out | std::ios::binary);
        qpp::save(C, fout);
    }
    {
        std::ifstream fin("mat.dat", std::ios::in | std::ios::binary);
        loadC = qpp::load<dyn_row_vect<int>>(fin);
    }
    EXPECT_NEAR(0, norm(loadC - C), 1e-7);

    // expression
    {
        std::ofstream fout("mat.dat", std::ios::out | std::ios::binary);
        qpp::save(3. * A + A, fout);
    }
    ket load_expression;
    {
        std::ifstream fin("mat.dat", std::ios::in | std::ios::binary);
        load_expression = qpp::load<ket>(fin);
    }
    EXPECT_NEAR(0, norm(load_expression - (3. * A + A)), 1e-7);
}
#endif // if (0)

/******************************************************************************/
/// BEGIN template <typename Derived> dyn_mat<typename Derived::Scalar>
///       qpp::load(std::istream& is)
///
///       template <typename Derived> void qpp::save(
///       const Eigen::MatrixBase<Derived>& A, std::ostream& os)
TEST(qpp_load_save, Matrix) {
    // matrices,complex, real and integer

    // DA = 1, DB = 1 degenerate case
    idx DA = 1, DB = 1;
    cmat A = rand<cmat>(DA, DB), loadA;
    dmat B = rand<dmat>(DA, DB), loadB;
    dyn_mat<int> C = dyn_mat<int>::Random(DA, DB), loadC;
    {
        std::ofstream fout("mat.txt");
        qpp::save(A, fout);
    }
    {
        std::ifstream fin("mat.txt");
        loadA = qpp::load<cmat>(fin);
    }
    EXPECT_NEAR(0, norm(loadA - A), 1e-7);

    {
        std::ofstream fout("mat.txt");
        qpp::save(B, fout);
    }
    {
        std::ifstream fin("mat.txt");
        loadB = qpp::load<dmat>(fin);
    }
    EXPECT_NEAR(0, norm(loadB - B), 1e-7);

    {
        std::ofstream fout("mat.txt");
        qpp::save(C, fout);
    }
    {
        std::ifstream fin("mat.txt");
        loadC = qpp::load<dyn_mat<int>>(fin);
    }
    EXPECT_NEAR(0, norm(loadC - C), 1e-7);

    // DA = 1, DB = 10
    DA = 1, DB = 10;
    A = rand<cmat>(DA, DB);
    B = rand<dmat>(DA, DB);
    C = dyn_mat<int>::Random(DA, DB);
    {
        std::ofstream fout("mat.txt");
        qpp::save(A, fout);
    }
    {
        std::ifstream fin("mat.txt");
        loadA = qpp::load<cmat>(fin);
    }
    EXPECT_NEAR(0, norm(loadA - A), 1e-7);

    {
        std::ofstream fout("mat.txt");
        qpp::save(B, fout);
    }
    {
        std::ifstream fin("mat.txt");
        loadB = qpp::load<dmat>(fin);
    }
    EXPECT_NEAR(0, norm(loadB - B), 1e-7);

    {
        std::ofstream fout("mat.txt");
        qpp::save(C, fout);
    }
    {
        std::ifstream fin("mat.txt");
        loadC = qpp::load<dyn_mat<int>>(fin);
    }
    EXPECT_NEAR(0, norm(loadC - C), 1e-7);

    // DA = 32, DB = 24
    DA = 32, DB = 24;
    A = rand<cmat>(DA, DB);
    B = rand<dmat>(DA, DB);
    C = dyn_mat<int>::Random(DA, DB);

    {
        std::ofstream fout("mat.txt");
        qpp::save(A, fout);
    }

    {
        std::ifstream fin("mat.txt");
        loadA = qpp::load<cmat>(fin);
    }
    EXPECT_NEAR(0, norm(loadA - A), 1e-7);

    {
        std::ofstream fout("mat.txt");
        qpp::save(B, fout);
    }
    {
        std::ifstream fin("mat.txt");
        loadB = qpp::load<dmat>(fin);
    }
    EXPECT_NEAR(0, norm(loadB - B), 1e-7);

    {
        std::ofstream fout("mat.txt");
        qpp::save(C, fout);
    }
    {
        std::ifstream fin("mat.txt");
        loadC = qpp::load<dyn_mat<int>>(fin);
    }
    EXPECT_NEAR(0, norm(loadC - C), 1e-7);

    // expression
    A = rand<cmat>(5, 5);
    cmat expression = A * A + A, load_expression;
    {
        std::ofstream fout("mat.txt");
        qpp::save(A * A + A, fout);
    }
    {
        std::ifstream fin("mat.txt");
        load_expression = qpp::load<cmat>(fin);
    }
    EXPECT_NEAR(0, norm(load_expression - expression), 1e-7);
}
/******************************************************************************/
TEST(qpp_load_save, Vector) {
    // kets/row vectors, complex, real and integer

    // D = 1 degenerate case
    idx D = 1;
    ket A = randket(D), loadA;
    dyn_row_vect<double> B = dyn_row_vect<double>::Random(D), loadB;
    dyn_row_vect<int> C = dyn_row_vect<int>::Random(D), loadC;

    {
        std::ofstream fout("mat.txt");
        qpp::save(A, fout);
    }
    {
        std::ifstream fin("mat.txt");
        loadA = qpp::load<ket>(fin);
    }
    EXPECT_NEAR(0, norm(loadA - A), 1e-7);

    {
        std::ofstream fout("mat.txt");
        qpp::save(B, fout);
    }
    {
        std::ifstream fin("mat.txt");
        loadB = qpp::load<dyn_row_vect<double>>(fin);
    }
    EXPECT_NEAR(0, norm(loadB - B), 1e-7);

    {
        std::ofstream fout("mat.txt");
        qpp::save(C, fout);
    }
    {
        std::ifstream fin("mat.txt");
        loadC = qpp::load<dyn_row_vect<int>>(fin);
    }
    EXPECT_NEAR(0, norm(loadC - C), 1e-7);

    // D = 32
    D = 32;
    A = randket(D);
    B = dyn_row_vect<double>::Random(D);
    C = dyn_row_vect<int>::Random(D);

    {
        std::ofstream fout("mat.txt");
        qpp::save(A, fout);
    }
    {
        std::ifstream fin("mat.txt");
        loadA = qpp::load<ket>(fin);
    }
    EXPECT_NEAR(0, norm(loadA - A), 1e-7);

    {
        std::ofstream fout("mat.txt");
        qpp::save(B, fout);
    }
    {
        std::ifstream fin("mat.txt");
        loadB = qpp::load<dyn_row_vect<double>>(fin);
    }
    EXPECT_NEAR(0, norm(loadB - B), 1e-7);

    {
        std::ofstream fout("mat.txt");
        qpp::save(C, fout);
    }
    {
        std::ifstream fin("mat.txt");
        loadC = qpp::load<dyn_row_vect<int>>(fin);
    }
    EXPECT_NEAR(0, norm(loadC - C), 1e-7);

    // expression
    ket expression = 3. * A + A, load_expression;
    {
        std::ofstream fout("mat.txt");
        qpp::save(3. * A + A, fout);
    }
    {
        std::ifstream fin("mat.txt");
        load_expression = qpp::load<ket>(fin);
    }
    EXPECT_NEAR(0, norm(load_expression - expression), 1e-7);
}
/******************************************************************************/
