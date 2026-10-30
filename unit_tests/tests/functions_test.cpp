#include <cmath>
#include <vector>

#include "gtest/gtest.h"

#include "qpp/qpp.hpp"

using namespace qpp;

// Unit testing "qpp/functions.hpp"

/// BEGIN template <char... Bits> bra operator"" _bra()
TEST(qpp_literals_operator_bra, AllTests) {
    EXPECT_LT((0_bra - st.z0.adjoint()).norm(), 1e-12);
}

/// BEGIN template <char... Bits> ket operator"" _ket()
TEST(qpp_literals_operator_ket, AllTests) {
    EXPECT_LT((0_ket - st.z0).norm(), 1e-12);
}

/// BEGIN template <char... Bits> cmat operator"" _prj()
TEST(qpp_literals_operator_prj, AllTests) {
    EXPECT_LT((0_prj - prj(st.z0)).norm(), 1e-12);
}

/// BEGIN template <typename Derived> cmat absm(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_absm, AllTests) {
    // |A| = sqrt(A^dagger A), not the entrywise absolute value.
    cmat A(2, 2);
    A << 0., 2_i, 1., 0.;

    cmat expected = cmat::Zero(2, 2);
    expected(0, 0) = 1;
    expected(1, 1) = 2;

    EXPECT_LT((absm(A) - expected).norm(), 1e-12);
    EXPECT_LT((absm(A) * absm(A) - adjoint(A) * A).norm(), 1e-12);

    EXPECT_THROW(absm(cmat{}), exception::ZeroSize);
    EXPECT_THROW(absm(cmat::Zero(2, 3)), exception::MatrixNotSquare);
}

/// BEGIN template <typename Container> std::vector<realT> abssq(
///       const Container& c,
///       typename std::enable_if<is_iterable<Container>::value>::type*
///       = nullptr)
TEST(qpp_abssq, Container) {
    std::vector<cplx> a{3., 4_i, -2.};
    EXPECT_EQ((std::vector<realT>{9, 16, 4}), abssq(a));
}

/// BEGIN template <typename Derived> std::vector<realT> abssq(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_abssq, Matrix) {
    ket a(2);
    a << 3., 4_i;
    EXPECT_EQ((std::vector<realT>{9, 16}), abssq(a));
}

/// BEGIN template <typename InputIterator> std::vector<realT> abssq(
///       InputIterator first, InputIterator last)
TEST(qpp_abssq, Iterator) {
    std::vector<cplx> a{3., 4_i};
    EXPECT_EQ((std::vector<realT>{9, 16}), abssq(a.begin(), a.end()));
}

/// BEGIN template <typename Derived> dyn_mat<typename Derived::Scalar>
///       adjoint(const Eigen::MatrixBase<Derived>& A)
TEST(qpp_adjoint, AllTests) {
    cmat A(2, 2);
    A << 1., 2_i, 3., 4.;
    EXPECT_LT((adjoint(A) - A.adjoint()).norm(), 1e-12);
}

/// BEGIN template <typename Derived1, typename Derived2>
///       dyn_mat<typename Derived1::Scalar> anticomm(
///       const Eigen::MatrixBase<Derived1>& A,
///       const Eigen::MatrixBase<Derived2>& B)
TEST(qpp_anticomm, AllTests) {
    EXPECT_LT((anticomm(gt.X, gt.Z)).norm(), 1e-12);
    EXPECT_LT((anticomm(gt.X, gt.X) - 2 * cmat::Identity(2, 2)).norm(), 1e-12);
}

/// BEGIN inline cmat bloch2rho(const std::vector<realT>& r)
TEST(qpp_bloch2rho, AllTests) {
    EXPECT_LT((bloch2rho({0, 0, 1}) - prj(st.z0)).norm(), 1e-12);
    EXPECT_LT((bloch2rho({0, 0, 0}) - cmat::Identity(2, 2) / 2.).norm(), 1e-12);
}

/// BEGIN template <typename Derived1, typename Derived2>
///       dyn_mat<typename Derived1::Scalar> comm(
///       const Eigen::MatrixBase<Derived1>& A,
///       const Eigen::MatrixBase<Derived2>& B)
TEST(qpp_comm, AllTests) {
    EXPECT_LT(comm(gt.X, gt.X).norm(), 1e-12);
    EXPECT_LT((comm(gt.X, gt.Z) - (gt.X * gt.Z - gt.Z * gt.X)).norm(), 1e-12);
}

/// BEGIN inline std::vector<idx> complement(std::vector<idx> subsys, idx n)
TEST(qpp_complement, AllTests) {
    EXPECT_EQ((std::vector<idx>{0, 2, 4}), complement({3, 1}, 5));
    EXPECT_TRUE(complement({0, 1}, 2).empty());
}

/// BEGIN template <typename Derived> dyn_mat<typename Derived::Scalar>
///       conjugate(const Eigen::MatrixBase<Derived>& A)
TEST(qpp_conjugate, AllTests) {
    cmat A(2, 2);
    A << 1_i, 2., -3_i, 4.;
    EXPECT_LT((conjugate(A) - A.conjugate()).norm(), 1e-12);
}

/// BEGIN template <typename Derived> cmat cosm(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_cosm, AllTests) {
    EXPECT_LT((cosm(cmat::Zero(2, 2)) - cmat::Identity(2, 2)).norm(), 1e-12);
}

/// BEGIN template <typename OutputScalar, typename Derived>
///       dyn_mat <OutputScalar> cwise(const Eigen::MatrixBase<Derived>& A,
///       OutputScalar (*f)(typename Derived::Scalar))
TEST(qpp_cwise, AllTests) {
    rmat A = rmat::Constant(2, 2, 3);
    auto square = +[](realT x) -> realT { return x * x; };
    EXPECT_LT((cwise(A, square) - rmat::Constant(2, 2, 9)).norm(), 1e-12);
}

/// BEGIN template <typename Derived> typename Derived::Scalar det(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_det, AllTests) {
    cmat A(2, 2);
    A << 1., 2., 3., 4.;
    EXPECT_NEAR(-2, std::real(det(A)), 1e-12);
}

/// BEGIN template <typename Derived> dirac_t<typename Derived::Scalar>
///       dirac(const Eigen::MatrixBase<Derived>& A,
///       const std::vector<idx>& dims_rows, const std::vector<idx>& dims_cols)
TEST(qpp_dirac, Qudits) {
    auto d = dirac(mket({1, 2}, {2, 3}), std::vector<idx>{2, 3},
                   std::vector<idx>{1});
    EXPECT_EQ(1u, d.states.size());
}

/// BEGIN template <typename Derived> dirac_t<typename Derived::Scalar>
///       dirac(const Eigen::MatrixBase<Derived>& A, idx d = 2)
TEST(qpp_dirac, Qubits) {
    auto d = dirac(st.z0);
    EXPECT_EQ(1u, d.states.size());
}

/// BEGIN template <typename Derived> dyn_mat<typename Derived::Scalar>
///       dirsum(const std::initializer_list<Derived>& As)
TEST(qpp_dirsum, InitList) {
    cmat A = cmat::Identity(2, 2);
    cmat B = 2 * A;
    cmat C = dirsum({A, B});
    EXPECT_EQ(4, C.rows());
    EXPECT_LT((C.topLeftCorner(2, 2) - A).norm(), 1e-12);
    EXPECT_LT((C.bottomRightCorner(2, 2) - B).norm(), 1e-12);
}

/// BEGIN template <typename T, typename ... Args>
/// dyn_mat<typename T::Scalar> dirsum(const T& head, const Args& ... tail)
TEST(qpp_dirsum, Variadic) {
    cmat A = cmat::Identity(2, 2);
    cmat B = 2 * A;
    EXPECT_LT((dirsum(A, B) - dirsum({A, B})).norm(), 1e-12);
}

/// BEGIN template <typename Derived> dyn_mat<typename Derived::Scalar>
///       dirsum(const std::vector<Derived>& As)
TEST(qpp_dirsum, Vector) {
    cmat A = cmat::Identity(2, 2);
    EXPECT_LT((dirsum(std::vector<cmat>{A, A}) - cmat::Identity(4, 4)).norm(),
              1e-12);
}

/// BEGIN template <typename Derived> dyn_mat<typename Derived::Scalar>
///       dirsumpow(const Eigen::MatrixBase<Derived>& A, idx n)
TEST(qpp_dirsumpow, AllTests) {
    cmat A = 2 * cmat::Identity(2, 2);
    EXPECT_LT((dirsumpow(A, 3) - 2 * cmat::Identity(6, 6)).norm(), 1e-12);
}

/// BEGIN template <typename Derived> std::pair<dyn_col_vect<cplx>, cmat>
///       eig (const Eigen::MatrixBase<Derived>& A)
TEST(qpp_eig, AllTests) {
    auto [values, vectors] = eig(gt.Z);
    EXPECT_LT((gt.Z * vectors - vectors * values.asDiagonal()).norm(), 1e-10);
}

/// BEGIN template <typename Derived> dyn_col_vect <cplx> evals(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_evals, AllTests) {
    auto values = evals(gt.Z);
    EXPECT_EQ(2, values.size());
    EXPECT_NEAR(0, std::real(values.prod()) + 1, 1e-12);
}

/// BEGIN template <typename Derived> cmat evects(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_evects, AllTests) {
    cmat V = evects(gt.Z);
    EXPECT_LT((V.adjoint() * V - cmat::Identity(2, 2)).norm(), 1e-10);
}

/// BEGIN template <typename Derived> cmat expm(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_expm, AllTests) {
    EXPECT_LT((expm(cmat::Zero(2, 2)) - cmat::Identity(2, 2)).norm(), 1e-12);
}

/// BEGIN template <typename Derived1, typename Derived2> realT expect(
///       const Eigen::MatrixBase<Derived1>& H,
///       const Eigen::MatrixBase<Derived2>& state)
TEST(qpp_expect, PureStates) {
    ket psi = (st.z0 + st.z1).normalized();

    EXPECT_NEAR(expect(gt.Z, st.z0), 1, 1e-12);
    EXPECT_NEAR(expect(gt.Z, st.z1), -1, 1e-12);
    EXPECT_NEAR(expect(gt.X, psi), 1, 1e-12);
    EXPECT_NEAR(expect(gt.Y, psi), 0, 1e-12);
    EXPECT_NEAR(expect(gt.Z, psi), 0, 1e-12);
}

TEST(qpp_expect, MixedStates) {
    // pure state represented as a density matrix
    ket psi = (st.z0 + st.z1).normalized();
    cmat rho = psi * psi.adjoint();

    EXPECT_NEAR(expect(gt.X, rho), 1, 1e-12);
    EXPECT_NEAR(expect(gt.Y, rho), 0, 1e-12);
    EXPECT_NEAR(expect(gt.Z, rho), 0, 1e-12);

    // maximally mixed state
    cmat rho_mixed = cmat::Identity(2, 2) / 2.;

    EXPECT_NEAR(expect(gt.X, rho_mixed), 0, 1e-12);
    EXPECT_NEAR(expect(gt.Y, rho_mixed), 0, 1e-12);
    EXPECT_NEAR(expect(gt.Z, rho_mixed), 0, 1e-12);

    // equivalence between state vector and density matrix
    EXPECT_NEAR(expect(gt.X, psi), expect(gt.X, rho), 1e-12);
    EXPECT_NEAR(expect(gt.Y, psi), expect(gt.Y, rho), 1e-12);
    EXPECT_NEAR(expect(gt.Z, psi), expect(gt.Z, rho), 1e-12);
}

/// BEGIN template <typename Derived1, typename Derived2> realT fidelity(
///       const Eigen::MatrixBase<Derived1>& A,
///       const Eigen::MatrixBase<Derived2>& B)
TEST(qpp_fidelity, PureStates) {
    EXPECT_NEAR(1 / std::sqrt(2.), fidelity(st.z0, st.x0), 1e-5);
    EXPECT_NEAR(0, fidelity(st.z0, st.z1), 1e-5);
    EXPECT_NEAR(0, fidelity(st.b00, st.b11), 1e-5);
    // global phases do not matter
    EXPECT_NEAR(1, fidelity(st.GHZ, std::exp(1_i * 0.3) * st.GHZ), 1e-5);

    for (idx D : std::vector<idx>{2, 3, 8}) {
        ket psi = randket(D);
        ket phi = randket(D);
        realT F = std::abs(psi.dot(phi));
        EXPECT_NEAR(1, fidelity(psi, psi), 1e-5);
        EXPECT_NEAR(F, fidelity(psi, phi), 1e-5);
        EXPECT_NEAR(F, fidelity(phi, psi), 1e-5);
        // kets and their projectors are interchangeable
        EXPECT_NEAR(F, fidelity(psi, prj(phi)), 1e-5);
        EXPECT_NEAR(F, fidelity(prj(psi), phi), 1e-5);
        EXPECT_NEAR(F, fidelity(prj(psi), prj(phi)), 1e-5);
    }
}

TEST(qpp_fidelity, MixedStates) {
    for (idx D : std::vector<idx>{2, 3, 8}) {
        cmat rho = randrho(D);
        cmat sigma = randrho(D);
        realT F = fidelity(rho, sigma);
        EXPECT_NEAR(1, fidelity(rho, rho), 1e-5);
        EXPECT_NEAR(F, fidelity(sigma, rho), 1e-5); // symmetric
        EXPECT_GE(F, -1e-5);
        EXPECT_LE(F, 1 + 1e-5);
        // unitary invariance
        cmat U = randU(D);
        EXPECT_NEAR(F, fidelity(U * rho * adjoint(U), U * sigma * adjoint(U)),
                    1e-5);
        // definition, F = tr sqrt(sqrt(rho) sigma sqrt(rho))
        cmat sqrt_rho = sqrtm(rho);
        EXPECT_NEAR(F, std::real(trace(sqrtm(sqrt_rho * sigma * sqrt_rho))),
                    1e-5);
        // maximally mixed state vs. pure state, F = 1 / sqrt(D)
        cmat mixed = cmat::Identity(D, D) / static_cast<realT>(D);
        EXPECT_NEAR(1 / std::sqrt(static_cast<realT>(D)),
                    fidelity(mixed, randket(D)), 1e-5);
    }

    // commuting (diagonal) states, F = sum_i sqrt(p_i q_i)
    std::vector<realT> p = randprob(4);
    std::vector<realT> q = randprob(4);
    cmat rho = cmat::Zero(4, 4);
    cmat sigma = cmat::Zero(4, 4);
    realT F = 0;
    for (idx i = 0; i < 4; ++i) {
        rho(i, i) = p[i];
        sigma(i, i) = q[i];
        F += std::sqrt(p[i] * q[i]);
    }
    EXPECT_NEAR(F, fidelity(rho, sigma), 1e-5);

    // rank-deficient states
    EXPECT_NEAR(0, fidelity(prj(st.z0), prj(st.z1)), 1e-5);
    EXPECT_NEAR(1 / std::sqrt(2.), fidelity(prj(st.z0), prj(st.x0)), 1e-5);
}

TEST(qpp_fidelity, Exceptions) {
    EXPECT_THROW(fidelity(ket{}, st.z0), exception::ZeroSize);
    EXPECT_THROW(fidelity(st.z0, cmat{}), exception::ZeroSize);
    EXPECT_THROW(fidelity(cmat::Zero(2, 3), st.z0),
                 exception::MatrixNotSquareNorCvector);
    EXPECT_THROW(fidelity(st.z0, cmat::Zero(2, 3)),
                 exception::MatrixNotSquareNorCvector);
    EXPECT_THROW(fidelity(st.z0, st.b00), exception::DimsNotEqual);
    EXPECT_THROW(fidelity(randrho(2), randrho(3)), exception::DimsNotEqual);
}

/// BEGIN template <typename Derived> cmat funm(
///       const Eigen::MatrixBase<Derived>& A, cplx (*f)(const cplx&))
TEST(qpp_funm, AllTests) {
    auto square = +[](const cplx& z) -> cplx { return z * z; };
    EXPECT_LT((funm(gt.Z, square) - cmat::Identity(2, 2)).norm(), 1e-10);
}

/// BEGIN template <typename Derived> dyn_mat<typename Derived::Scalar>
///       grams(const std::initializer_list<Derived>& As)
TEST(qpp_grams, InitList) {
    ket a = st.z0;
    ket b = (st.z0 + st.z1).eval();
    cmat Q = grams({a, b});
    EXPECT_LT((Q.adjoint() * Q - cmat::Identity(2, 2)).norm(), 1e-10);
}

/// BEGIN template <typename Derived> dyn_mat<typename Derived::Scalar>
///       :grams(const Eigen::MatrixBase<Derived>& A)
TEST(qpp_grams, Matrix) {
    cmat A(2, 2);
    A << 1., 1., 0., 1.;
    cmat Q = grams(A);
    EXPECT_LT((Q.adjoint() * Q - cmat::Identity(2, 2)).norm(), 1e-10);
}

/// BEGIN template <typename Derived> dyn_mat<typename Derived::Scalar>
///       grams(const std::vector<Derived>& As)
TEST(qpp_grams, Vector) {
    ket a = st.z0;
    ket b = (st.z0 + st.z1).eval();
    cmat Q = grams(std::vector<ket>{a, b});
    EXPECT_LT((Q.adjoint() * Q - cmat::Identity(2, 2)).norm(), 1e-10);
}

/// BEGIN template <typename Derived> std::size_t hash_eigen(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_hash_eigen, AllTests) {
    cmat A = gt.X;
    cmat B = gt.X;
    EXPECT_EQ(hash_eigen(A), hash_eigen(B));
}
/// BEGIN template <typename Derived> std::pair<dyn_col_vect<realT>, cmat>
///       heig(const Eigen::MatrixBase<Derived>& A)
TEST(qpp_heig, AllTests) {
    auto [values, vectors] = heig(gt.Z);
    EXPECT_LT(
        (gt.Z * vectors - vectors * values.cast<cplx>().asDiagonal()).norm(),
        1e-10);
}

/// BEGIN template <typename Derived> dyn_col_vect<realT> hevals(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_hevals, AllTests) {
    auto values = hevals(gt.Z);
    EXPECT_NEAR(-1, values(0), 1e-12);
    EXPECT_NEAR(1, values(1), 1e-12);
}

/// BEGIN template <typename Derived> cmat hevects(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_hevects, AllTests) {
    cmat V = hevects(gt.X);
    EXPECT_LT((V.adjoint() * V - cmat::Identity(2, 2)).norm(), 1e-10);
}

/// BEGIN template <typename Derived> dyn_mat<typename Derived::Scalar>
///       inverse(const Eigen::MatrixBase<Derived>& A)
TEST(qpp_inverse, AllTests) {
    cmat A(2, 2);
    A << 1., 2., 3., 5.;
    EXPECT_LT((A * inverse(A) - cmat::Identity(2, 2)).norm(), 1e-10);
}

/// BEGIN template <typename Derived> dyn_mat<typename Derived::Scalar>
///       kron(const std::initializer_list<Derived>& As)
TEST(qpp_kron, InitList) {
    cmat A = gt.X;
    EXPECT_LT((kron({A, A}) - kron(A, A)).norm(), 1e-12);
}

/// BEGIN template <typename T, typename ... Args> dyn_mat<typename T::Scalar>
///       kron(const T& head, const Args& ... tail)
TEST(qpp_kron, Variadic) {
    EXPECT_LT((kron(gt.X, gt.Z) - kron(std::vector<cmat>{gt.X, gt.Z})).norm(),
              1e-12);
}

/// BEGIN template <typename Derived> dyn_mat<typename Derived::Scalar>
///       kron(const std::vector<Derived>& As)
TEST(qpp_kron, Vector) {
    EXPECT_LT((kron(std::vector<cmat>{gt.X, gt.Z}) - kron(gt.X, gt.Z)).norm(),
              1e-12);
}

/// BEGIN template <typename Derived> dyn_mat<typename Derived::Scalar>
///       kronpow(const Eigen::MatrixBase<Derived>& A, idx n)
TEST(qpp_kronpow, AllTests) {
    EXPECT_LT((kronpow(gt.X, 2) - kron(gt.X, gt.X)).norm(), 1e-12);
}

/// BEGIN template <typename Derived> typename Derived::Scalar logdet(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_logdet, AllTests) {
    cmat A = 2 * cmat::Identity(2, 2);
    EXPECT_NEAR(std::log(4.), std::real(logdet(A)), 1e-12);
}

/// BEGIN template <typename Derived> cmat logm(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_logm, AllTests) {
    EXPECT_LT(logm(cmat::Identity(2, 2)).norm(), 1e-10);
}

/// BEGIN inline ket mket(const std::vector<idx>& mask,
///       const std::vector<idx>& dims)
TEST(qpp_mket, Qudits) {
    EXPECT_LT((mket({1, 2}, {2, 3}) - mket({5}, 6)).norm(), 1e-12);
}

/// BEGIN inline ket mket(const std::vector<idx>& mask, idx d = 2)
TEST(qpp_mket, Qubits) { EXPECT_LT((mket({1, 0}) - 10_ket).norm(), 1e-12); }

/// BEGIN inline cmat mprj(const std::vector<idx>& mask,
///       const std::vector<idx>& dims)
TEST(qpp_mprj, Qudits) {
    EXPECT_LT((mprj({1, 2}, {2, 3}) - prj(mket({1, 2}, {2, 3}))).norm(), 1e-12);
}

/// BEGIN inline cmat mprj(const std::vector<idx>& mask, idx d = 2)
TEST(qpp_mprj, Qubits) {
    EXPECT_LT((mprj({1, 0}) - prj(10_ket)).norm(), 1e-12);
}

/// BEGIN template <typename V, typename T = V, typename U = T>
///       T multiidx2n(const std::vector<V>& midx, const std::vector<U>& dims)
TEST(qpp_multiidx2n, AllTests) {
    EXPECT_EQ(5u, multiidx2n(std::vector<idx>{1, 2}, std::vector<idx>{2, 3}));
}

/// BEGIN template <typename T, typename U = T, typename V = T> n2multiidx(T n,
///       const std::vector<U>& dims)
TEST(qpp_n2multiidx, AllTests) {
    EXPECT_EQ((std::vector<idx>{1, 2}),
              n2multiidx(idx{5}, std::vector<idx>{2, 3}));
}

/// BEGIN template <typename Derived> realT norm(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_norm, AllTests) {
    EXPECT_NEAR(5, norm((ket(2) << 3., 4.).finished()), 1e-12);
}

/// BEGIN template <typename Derived> dyn_mat<typename Derived::Scalar>
///       normalize(const Eigen::MatrixBase<Derived>& A)
TEST(qpp_normalize, AllTests) {
    ket A(2);
    A << 3., 4.;
    EXPECT_NEAR(1, norm(normalize(A)), 1e-12);
}

/// BEGIN template <typename Derived> dyn_mat<typename Derived::Scalar>
///       powm(const Eigen::MatrixBase<Derived>& A, idx n)
TEST(qpp_powm, AllTests) {
    EXPECT_LT((powm(gt.X, 2) - cmat::Identity(2, 2)).norm(), 1e-12);
}

/// BEGIN template <typename Derived> dyn_mat<typename Derived::Scalar>
///       prj(const Eigen::MatrixBase<Derived>& A)
TEST(qpp_prj, AllTests) { EXPECT_LT((prj(st.z0) - mprj({0})).norm(), 1e-12); }

/// BEGIN template <typename Derived> realT purity(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_purity, AllTests) {
    for (idx D : std::vector<idx>{2, 3, 8}) {
        // pure states
        ket psi = randket(D);
        EXPECT_NEAR(1, purity(psi), 1e-5);
        EXPECT_NEAR(1, purity(prj(psi)), 1e-5);
        // maximally mixed state
        cmat mixed = cmat::Identity(D, D) / static_cast<realT>(D);
        EXPECT_NEAR(1 / static_cast<realT>(D), purity(mixed), 1e-5);
        // random state, purity = tr(rho^2) in [1/D, 1]
        cmat rho = randrho(D);
        realT P = purity(rho);
        EXPECT_NEAR(std::real(trace(rho * rho)), P, 1e-5);
        EXPECT_GE(P, 1 / static_cast<realT>(D) - 1e-5);
        EXPECT_LE(P, 1 + 1e-5);
    }
    // real matrices
    rmat rho = rmat::Identity(2, 2) / 2;
    EXPECT_NEAR(0.5, purity(rho), 1e-5);

    EXPECT_THROW(purity(cmat{}), exception::ZeroSize);
    EXPECT_THROW(purity(cmat::Zero(2, 3)),
                 exception::MatrixNotSquareNorCvector);
}

/// BEGIN template <typename Container> typename Container::value_type
///       prod(const Container& c,
///       typename std::enable_if<is_iterable<Container>::value>::type*
///       = nullptr)
TEST(qpp_prod, Container) { EXPECT_EQ(24, prod(std::vector<int>{1, 2, 3, 4})); }

/// BEGIN template <typename T> T prod(const std::initializer_list<T>& Ts)
TEST(qpp_prod, InitList) { EXPECT_EQ(24, prod({1, 2, 3, 4})); }

/// BEGIN template <typename Derived> typename Derived::Scalar prod(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_prod, Matrix) {
    rmat A(2, 2);
    A << 1, 2, 3, 4;
    EXPECT_DOUBLE_EQ(24, prod(A));
}

/// BEGIN template <typename InputIterator, typename value_type =
///       std::decay_t<decltype(*std::declval<InputIterator>())>>
///       value_type prod(InputIterator first, InputIterator last)
TEST(qpp_prod, Iterator) {
    std::vector<int> v{1, 2, 3, 4};
    EXPECT_EQ(24, prod(v.begin(), v.end()));
}

/// BEGIN template <typename Derived> dyn_mat<typename Derived::Scalar>
///       reshape(const Eigen::MatrixBase<Derived>& A, idx rows, idx cols)
TEST(qpp_reshape, AllTests) {
    rmat A(2, 2);
    A << 1, 2, 3, 4;
    auto B = reshape(A, 1, 4);
    EXPECT_EQ(1, B.rows());
    EXPECT_EQ(4, B.cols());
    EXPECT_DOUBLE_EQ(1, B(0, 0));
}

/// BEGIN template <typename Derived> std::vector<realT> rho2bloch(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_rho2bloch, AllTests) {
    auto r = rho2bloch(prj(st.z0));
    ASSERT_EQ(3u, r.size());
    EXPECT_NEAR(1, r[2], 1e-12);
    EXPECT_NEAR(0, r[0], 1e-12);
}

/// BEGIN template <typename Derived> dyn_col_vect<typename Derived::Scalar>
///       rho2pure(const Eigen::MatrixBase<Derived>& A)
TEST(qpp_rho2pure, AllTests) {
    ket psi = rho2pure(prj(st.x0));
    EXPECT_NEAR(1, std::abs(st.x0.dot(psi)), 1e-10);
}

/// BEGIN template <typename Derived> realT schatten(
///       const Eigen::MatrixBase<Derived>& A, realT p)
TEST(qpp_schatten, AllTests) {
    EXPECT_NEAR(2, schatten(cmat::Identity(2, 2), 1), 1e-12);
    EXPECT_NEAR(std::sqrt(2.), schatten(cmat::Identity(2, 2), 2), 1e-12);
}

/// BEGIN template <typename Derived> cmat sinm(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_sinm, AllTests) { EXPECT_LT(sinm(cmat::Zero(2, 2)).norm(), 1e-12); }

/// BEGIN template <typename Derived> cmat spectralpowm(
///       const Eigen::MatrixBase<Derived>& A, const cplx z)
TEST(qpp_spectralpowm, AllTests) {
    EXPECT_LT((spectralpowm(2 * cmat::Identity(2, 2), cplx{0.5}) -
               std::sqrt(2.) * cmat::Identity(2, 2))
                  .norm(),
              1e-10);
}

/// BEGIN template <typename Derived> cmat sqrtm(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_sqrtm, AllTests) {
    EXPECT_LT(
        (sqrtm(4 * cmat::Identity(2, 2)) - 2 * cmat::Identity(2, 2)).norm(),
        1e-10);
}

/// BEGIN template <typename Container> typename Container::value_type
///       sum(const Container& c,
///       typename std::enable_if<is_iterable<Container>::value>::type*
///       = nullptr)
TEST(qpp_sum, Container) { EXPECT_EQ(10, sum(std::vector<int>{1, 2, 3, 4})); }

/// BEGIN template <typename T> T sum(const std::initializer_list<T>& Ts)
TEST(qpp_sum, InitList) { EXPECT_EQ(6, sum({0, 1, 2, 3})); }
/// BEGIN template <typename Derived> typename Derived::Scalar sum(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_sum, Matrix) {
    rmat A(2, 2);
    A << 1, 2, 3, 4;
    EXPECT_DOUBLE_EQ(10, sum(A));
}

/// BEGIN template <typename InputIterator, typename value_type =
///       std::decay_t<decltype(*std::declval<InputIterator>())>>
///      value_type sum(InputIterator first, InputIterator last)
TEST(qpp_sum, Iterator) {
    std::vector<int> v{0, 1, 2, 3};
    EXPECT_EQ(6, sum(v.begin(), v.end()));
}

/// BEGIN template <typename Derived> dyn_col_vect<realT> svals(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_svals, AllTests) {
    auto s = svals(2 * cmat::Identity(2, 2));
    ASSERT_EQ(2, s.size());
    EXPECT_NEAR(2, s(0), 1e-12);
}

/// BEGIN template <typename Derived>
///       std::tuple<cmat, dyn_col_vect<realT>, cmat> svd(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_svd, AllTests) {
    cmat A = gt.X;
    auto [U, s, V] = svd(A);
    EXPECT_LT((A - U * s.cast<cplx>().asDiagonal() * V.adjoint()).norm(),
              1e-10);
}

/// BEGIN template <typename Derived> cmat svdU(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_svdU, AllTests) {
    cmat U = svdU(gt.X);
    EXPECT_LT((U.adjoint() * U - cmat::Identity(2, 2)).norm(), 1e-10);
}

/// BEGIN template <typename Derived> cmat svdV(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_svdV, AllTests) {
    cmat V = svdV(gt.X);
    EXPECT_LT((V.adjoint() * V - cmat::Identity(2, 2)).norm(), 1e-10);
}

/// BEGIN template <typename Derived> typename Derived::Scalar trace(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_trace, AllTests) {
    EXPECT_NEAR(2, std::real(trace(cmat::Identity(2, 2))), 1e-12);
    EXPECT_NEAR(0, std::real(trace(gt.Z)), 1e-12);
}

/// BEGIN template <typename Derived1, typename Derived2> realT
///       trace_distance(const Eigen::MatrixBase<Derived1>& A,
///       const Eigen::MatrixBase<Derived2>& B)
TEST(qpp_trace_distance, AllTests) {
    // orthogonal and identical states
    EXPECT_NEAR(1, trace_distance(st.z0, st.z1), 1e-5);
    EXPECT_NEAR(1, trace_distance(prj(st.b00), prj(st.b01)), 1e-5);
    EXPECT_NEAR(0, trace_distance(st.x0, st.x0), 1e-5);
    EXPECT_NEAR(1 / std::sqrt(2.), trace_distance(st.z0, st.x0), 1e-5);

    for (idx D : std::vector<idx>{2, 3, 8}) {
        // pure states, T = sqrt(1 - F^2)
        ket psi = randket(D);
        ket phi = randket(D);
        realT F = fidelity(psi, phi);
        realT T = std::sqrt(1 - F * F);
        EXPECT_NEAR(T, trace_distance(psi, phi), 1e-5);
        EXPECT_NEAR(T, trace_distance(psi, prj(phi)), 1e-5);
        EXPECT_NEAR(T, trace_distance(prj(psi), phi), 1e-5);
        EXPECT_NEAR(T, trace_distance(prj(psi), prj(phi)), 1e-5);

        // mixed states, T = ||rho - sigma||_1 / 2, symmetric, and the
        // Fuchs-van de Graaf inequalities 1 - F <= T <= sqrt(1 - F^2)
        cmat rho = randrho(D);
        cmat sigma = randrho(D);
        T = trace_distance(rho, sigma);
        F = fidelity(rho, sigma);
        EXPECT_NEAR(schatten(rho - sigma, 1) / 2, T, 1e-5);
        EXPECT_NEAR(T, trace_distance(sigma, rho), 1e-5);
        EXPECT_LE(1 - F, T + 1e-5);
        EXPECT_LE(T, std::sqrt(1 - F * F) + 1e-5);
    }

    // unnormalized kets |a> are identified with |a><a| (not normalized)
    ket a = 2 * randket(3);
    ket b = randket(3);
    cmat rho_a = a * adjoint(a);
    cmat rho_b = b * adjoint(b);
    EXPECT_NEAR(trace_distance(rho_a, rho_b), trace_distance(a, b), 1e-5);
    EXPECT_NEAR(trace_distance(rho_a, rho_b), trace_distance(a, rho_b), 1e-5);

    // commuting (diagonal) states, T = sum_i |p_i - q_i| / 2
    std::vector<realT> p = randprob(4);
    std::vector<realT> q = randprob(4);
    cmat rho = cmat::Zero(4, 4);
    cmat sigma = cmat::Zero(4, 4);
    realT T = 0;
    for (idx i = 0; i < 4; ++i) {
        rho(i, i) = p[i];
        sigma(i, i) = q[i];
        T += std::abs(p[i] - q[i]) / 2;
    }
    EXPECT_NEAR(T, trace_distance(rho, sigma), 1e-5);

    EXPECT_THROW(trace_distance(ket{}, st.z0), exception::ZeroSize);
    EXPECT_THROW(trace_distance(cmat::Zero(2, 3), st.z0),
                 exception::MatrixNotSquareNorCvector);
    EXPECT_THROW(trace_distance(st.z0, st.GHZ), exception::DimsNotEqual);
}

/// BEGIN template <typename Derived> dyn_mat<typename Derived::Scalar>
///       transpose(const Eigen::MatrixBase<Derived>& A)
TEST(qpp_transpose, AllTests) {
    cmat A(2, 2);
    A << 1., 2_i, 3., 4.;
    EXPECT_LT((transpose(A) - A.transpose()).norm(), 1e-12);
}

/// BEGIN template <typename Derived> std::optional<std::vector<idx>>
///       zket2dits(const Eigen::MatrixBase<Derived>& psi,
///       const std::vector<idx>& dims, realT precision = 1e-12)
TEST(qpp_zket2dits, Qudits) {
    auto digits = zket2dits(mket({1, 2}, {2, 3}), std::vector<idx>{2, 3});
    ASSERT_TRUE(digits.has_value());
    EXPECT_EQ((std::vector<idx>{1, 2}), *digits);
}

/// BEGIN template <typename Derived> std::optional<std::vector<idx>>
///       zket2dits(const Eigen::MatrixBase<Derived>& psi, idx d = 2,
///       realT precision = 1e-12)
TEST(qpp_zket2dits, Qubits) {
    auto digits = zket2dits(10_ket);
    ASSERT_TRUE(digits.has_value());
    EXPECT_EQ((std::vector<idx>{1, 0}), *digits);
    EXPECT_FALSE(zket2dits(st.x0).has_value());
}
