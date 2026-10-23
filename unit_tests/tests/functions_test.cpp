#include <cmath>
#include <vector>

#include "gtest/gtest.h"

#include "qpp/qpp.hpp"

using namespace qpp;

// Unit testing "qpp/functions.hpp"

/// BEGIN template <char... Bits> bra operator"" _bra()
TEST(qpp_literals_operator_bra, AllTests) {}

/// BEGIN template <char... Bits> ket operator"" _ket()
TEST(qpp_literals_operator_ket, AllTests) {}

/// BEGIN template <char... Bits> cmat operator"" _prj()
TEST(qpp_literals_operator_prj, AllTests) {}

/// BEGIN template <typename Derived> cmat absm(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_absm, AllTests) {}

/// BEGIN template <typename Container> std::vector<realT> abssq(
///       const Container& c,
///       typename std::enable_if<is_iterable<Container>::value>::type*
///       = nullptr)
TEST(qpp_abssq, Container) {}

/// BEGIN template <typename Derived> std::vector<realT> abssq(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_abssq, Matrix) {}

/// BEGIN template <typename InputIterator> std::vector<realT> abssq(
///       InputIterator first, InputIterator last)
TEST(qpp_abssq, Iterator) {}

/// BEGIN template <typename Derived> dyn_mat<typename Derived::Scalar>
///       adjoint(const Eigen::MatrixBase<Derived>& A)
TEST(qpp_adjoint, AllTests) {}

/// BEGIN template <typename Derived1, typename Derived2>
///       dyn_mat<typename Derived1::Scalar> anticomm(
///       const Eigen::MatrixBase<Derived1>& A,
///       const Eigen::MatrixBase<Derived2>& B)
TEST(qpp_anticomm, AllTests) {}

/// BEGIN inline cmat bloch2rho(const std::vector<realT>& r)
TEST(qpp_bloch2rho, AllTests) {}

/// BEGIN template <typename Derived1, typename Derived2>
///       dyn_mat<typename Derived1::Scalar> comm(
///       const Eigen::MatrixBase<Derived1>& A,
///       const Eigen::MatrixBase<Derived2>& B)
TEST(qpp_comm, AllTests) {}

/// BEGIN inline std::vector<idx> complement(std::vector<idx> subsys, idx n)
TEST(qpp_complement, AllTests) {}

/// BEGIN template <typename Derived> dyn_mat<typename Derived::Scalar>
///       conjugate(const Eigen::MatrixBase<Derived>& A)
TEST(qpp_conjugate, AllTests) {}

/// BEGIN template <typename Derived> cmat cosm(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_cosm, AllTests) {}

/// BEGIN template <typename OutputScalar, typename Derived>
///       dyn_mat <OutputScalar> cwise(const Eigen::MatrixBase<Derived>& A,
///       OutputScalar (*f)(typename Derived::Scalar))
TEST(qpp_cwise, AllTests) {}

/// BEGIN template <typename Derived> typename Derived::Scalar det(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_det, AllTests) {}

/// BEGIN template <typename Derived> dirac_t<typename Derived::Scalar>
///       dirac(const Eigen::MatrixBase<Derived>& A,
///       const std::vector<idx>& dims_rows, const std::vector<idx>& dims_cols)
TEST(qpp_dirac, Qudits) {}

/// BEGIN template <typename Derived> dirac_t<typename Derived::Scalar>
///       dirac(const Eigen::MatrixBase<Derived>& A, idx d = 2)
TEST(qpp_dirac, Qubits) {}

/// BEGIN template <typename Derived> dyn_mat<typename Derived::Scalar>
///       dirsum(const std::initializer_list<Derived>& As)
TEST(qpp_dirsum, InitList) {}

/// BEGIN template <typename T, typename ... Args>
/// dyn_mat<typename T::Scalar> dirsum(const T& head, const Args& ... tail)
TEST(qpp_dirsum, Variadic) {}

/// BEGIN template <typename Derived> dyn_mat<typename Derived::Scalar>
///       dirsum(const std::vector<Derived>& As)
TEST(qpp_dirsum, Vector) {}

/// BEGIN template <typename Derived> dyn_mat<typename Derived::Scalar>
///       dirsumpow(const Eigen::MatrixBase<Derived>& A, idx n)
TEST(qpp_dirsumpow, AllTests) {}

/// BEGIN template <typename Derived> std::pair<dyn_col_vect<cplx>, cmat>
///       eig (const Eigen::MatrixBase<Derived>& A)
TEST(qpp_eig, AllTests) {}

/// BEGIN template <typename Derived> dyn_col_vect <cplx> evals(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_evals, AllTests) {}

/// BEGIN template <typename Derived> cmat evects(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_evects, AllTests) {}

/// BEGIN template <typename Derived> cmat expm(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_expm, AllTests) {}

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
TEST(qpp_funm, AllTests) {}

/// BEGIN template <typename Derived> dyn_mat<typename Derived::Scalar>
///       grams(const std::initializer_list<Derived>& As)
TEST(qpp_grams, InitList) {}

/// BEGIN template <typename Derived> dyn_mat<typename Derived::Scalar>
///       :grams(const Eigen::MatrixBase<Derived>& A)
TEST(qpp_grams, Matrix) {}

/// BEGIN template <typename Derived> dyn_mat<typename Derived::Scalar>
///       grams(const std::vector<Derived>& As)
TEST(qpp_grams, Vector) {}

/// BEGIN template <typename Derived> std::size_t hash_eigen(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_hash_eigen, AllTests) {}
/// BEGIN template <typename Derived> std::pair<dyn_col_vect<realT>, cmat>
///       heig(const Eigen::MatrixBase<Derived>& A)
TEST(qpp_heig, AllTests) {}

/// BEGIN template <typename Derived> dyn_col_vect<realT> hevals(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_hevals, AllTests) {}

/// BEGIN template <typename Derived> cmat hevects(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_hevects, AllTests) {}

/// BEGIN template <typename Derived> dyn_mat<typename Derived::Scalar>
///       inverse(const Eigen::MatrixBase<Derived>& A)
TEST(qpp_inverse, AllTests) {}

/// BEGIN template <typename Derived> dyn_mat<typename Derived::Scalar>
///       kron(const std::initializer_list<Derived>& As)
TEST(qpp_kron, InitList) {}

/// BEGIN template <typename T, typename ... Args> dyn_mat<typename T::Scalar>
///       kron(const T& head, const Args& ... tail)
TEST(qpp_kron, Variadic) {}

/// BEGIN template <typename Derived> dyn_mat<typename Derived::Scalar>
///       kron(const std::vector<Derived>& As)
TEST(qpp_kron, Vector) {}

/// BEGIN template <typename Derived> dyn_mat<typename Derived::Scalar>
///       kronpow(const Eigen::MatrixBase<Derived>& A, idx n)
TEST(qpp_kronpow, AllTests) {}

/// BEGIN template <typename Derived> typename Derived::Scalar logdet(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_logdet, AllTests) {}

/// BEGIN template <typename Derived> cmat logm(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_logm, AllTests) {}

/// BEGIN inline ket mket(const std::vector<idx>& mask,
///       const std::vector<idx>& dims)
TEST(qpp_mket, Qudits) {}

/// BEGIN inline ket mket(const std::vector<idx>& mask, idx d = 2)
TEST(qpp_mket, Qubits) {}

/// BEGIN inline cmat mprj(const std::vector<idx>& mask,
///       const std::vector<idx>& dims)
TEST(qpp_mprj, Qudits) {}

/// BEGIN inline cmat mprj(const std::vector<idx>& mask, idx d = 2)
TEST(qpp_mprj, Qubits) {}

/// BEGIN template <typename V, typename T = V, typename U = T>
///       T multiidx2n(const std::vector<V>& midx, const std::vector<U>& dims)
TEST(qpp_multiidx2n, AllTests) {}

/// BEGIN template <typename T, typename U = T, typename V = T> n2multiidx(T n,
///       const std::vector<U>& dims)
TEST(qpp_n2multiidx, AllTests) {}

/// BEGIN template <typename Derived> realT norm(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_norm, AllTests) {}

/// BEGIN template <typename Derived> dyn_mat<typename Derived::Scalar>
///       normalize(const Eigen::MatrixBase<Derived>& A)
TEST(qpp_normalize, AllTests) {}

/// BEGIN template <typename Derived> dyn_mat<typename Derived::Scalar>
///       powm(const Eigen::MatrixBase<Derived>& A, idx n)
TEST(qpp_powm, AllTests) {}

/// BEGIN template <typename Derived> dyn_mat<typename Derived::Scalar>
///       prj(const Eigen::MatrixBase<Derived>& A)
TEST(qpp_prj, AllTests) {}

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
TEST(qpp_prod, Container) {}

/// BEGIN template <typename T> T prod(const std::initializer_list<T>& Ts)
TEST(qpp_prod, InitList) { EXPECT_EQ(24, prod({1, 2, 3, 4})); }

/// BEGIN template <typename Derived> typename Derived::Scalar prod(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_prod, Matrix) {}

/// BEGIN template <typename InputIterator, typename value_type =
///       std::decay_t<decltype(*std::declval<InputIterator>())>>
///       value_type prod(InputIterator first, InputIterator last)
TEST(qpp_prod, Iterator) {
    std::vector<int> v{1, 2, 3, 4};
    EXPECT_EQ(24, prod(v.begin(), v.end()));
}

/// BEGIN template <typename Derived> dyn_mat<typename Derived::Scalar>
///       reshape(const Eigen::MatrixBase<Derived>& A, idx rows, idx cols)
TEST(qpp_reshape, AllTests) {}

/// BEGIN template <typename Derived> std::vector<realT> rho2bloch(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_rho2bloch, AllTests) {}

/// BEGIN template <typename Derived> dyn_col_vect<typename Derived::Scalar>
///       rho2pure(const Eigen::MatrixBase<Derived>& A)
TEST(qpp_rho2pure, AllTests) {}

/// BEGIN template <typename Derived> realT schatten(
///       const Eigen::MatrixBase<Derived>& A, realT p)
TEST(qpp_schatten, AllTests) {}

/// BEGIN template <typename Derived> cmat sinm(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_sinm, AllTests) {}

/// BEGIN template <typename Derived> cmat spectralpowm(
///       const Eigen::MatrixBase<Derived>& A, const cplx z)
TEST(qpp_spectralpowm, AllTests) {}

/// BEGIN template <typename Derived> cmat sqrtm(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_sqrtm, AllTests) {}

/// BEGIN template <typename Container> typename Container::value_type
///       sum(const Container& c,
///       typename std::enable_if<is_iterable<Container>::value>::type*
///       = nullptr)
TEST(qpp_sum, Container) {}

/// BEGIN template <typename T> T sum(const std::initializer_list<T>& Ts)
TEST(qpp_sum, InitList) { EXPECT_EQ(6, sum({0, 1, 2, 3})); }
/// BEGIN template <typename Derived> typename Derived::Scalar sum(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_sum, Matrix) {}

/// BEGIN template <typename InputIterator, typename value_type =
///       std::decay_t<decltype(*std::declval<InputIterator>())>>
///      value_type sum(InputIterator first, InputIterator last)
TEST(qpp_sum, Iterator) {
    std::vector<int> v{0, 1, 2, 3};
    EXPECT_EQ(6, sum(v.begin(), v.end()));
}

/// BEGIN template <typename Derived> dyn_col_vect<realT> svals(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_svals, AllTests) {}

/// BEGIN template <typename Derived>
///       std::tuple<cmat, dyn_col_vect<realT>, cmat> svd(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_svd, AllTests) {}

/// BEGIN template <typename Derived> cmat svdU(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_svdU, AllTests) {}

/// BEGIN template <typename Derived> cmat svdV(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_svdV, AllTests) {}

/// BEGIN template <typename Derived> typename Derived::Scalar trace(
///       const Eigen::MatrixBase<Derived>& A)
TEST(qpp_trace, AllTests) {}

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
TEST(qpp_transpose, AllTests) {}

/// BEGIN template <typename Derived> std::optional<std::vector<idx>>
///       zket2dits(const Eigen::MatrixBase<Derived>& psi,
///       const std::vector<idx>& dims, realT precision = 1e-12)
TEST(qpp_zket2dits, Qudits) {}

/// BEGIN template <typename Derived> std::optional<std::vector<idx>>
///       zket2dits(const Eigen::MatrixBase<Derived>& psi, idx d = 2,
///       realT precision = 1e-12)
TEST(qpp_zket2dits, Qubits) {}
