#include <limits>
#include <vector>

#include "gtest/gtest.h"

#include "qpp/qpp.hpp"

using namespace qpp;

// Unit testing "qpp/classes/qcircuit.hpp"

/// BEGIN QCircuit& QCircuit::add_dit(idx n = 1, idx i)
TEST(qpp_QCircuit_add_dit, SpecificPosition) {
    // qutrits
    QCircuit q3(2, 2, 3);
    q3.add_dit(2, 1);
    EXPECT_EQ(4u, q3.get_nc());
    // qubits
    QCircuit qc(2, 2);
    qc.add_dit(2, 1);
    EXPECT_EQ(4u, qc.get_nc());
    EXPECT_EQ((std::vector<idx>{0, 1, 2, 3}), qc.get_clean_dits());
}

/// BEGIN QCircuit& QCircuit::add_dit(idx n = 1)
TEST(qpp_QCircuit_add_dit, LastPosition) {
    // qutrits
    QCircuit q3(1, 1, 3);
    q3.add_dit(2);
    EXPECT_EQ(3u, q3.get_nc());
    // qubits
    QCircuit qc(2, 1);
    qc.add_dit(2);
    EXPECT_EQ(3u, qc.get_nc());
}

/// BEGIN QCircuit& QCircuit::add_qudit(idx n = 1, idx i)
TEST(qpp_QCircuit_add_qudit, SpecificPosition) {
    // qutrits
    QCircuit q3(2, 0, 3);
    q3.add_qudit(2, 1);
    EXPECT_EQ(4u, q3.get_nq());
    // qubits
    QCircuit qc(2);
    qc.add_qudit(2, 1);
    EXPECT_EQ(4u, qc.get_nq());
    EXPECT_EQ((std::vector<idx>{0, 1, 2, 3}), qc.get_clean_qudits());
}
///
/// BEGIN QCircuit& QCircuit::add_qudit(idx n = 1)
TEST(qpp_QCircuit_add_qudit, LastPosition) {
    // qutrits
    QCircuit q3(2, 0, 3);
    q3.add_qudit(2);
    EXPECT_EQ(4u, q3.get_nq());
    // qubits
    QCircuit qc(2);
    qc.add_qudit(3);
    EXPECT_EQ(5u, qc.get_nq());
}

/// BEGIN QCircuit& QCircuit::adjoint()
TEST(qpp_QCircuit_adjoint, AllTests) {
    // qutrits
    QCircuit q3(1, 0, 3);
    q3.gate(gt.Xd(3), 0);
    const auto orig3 = q3;
    q3.adjoint().adjoint();
    EXPECT_EQ(orig3, q3);
    // qubits
    QCircuit qc(1);
    qc.gate(gt.S, 0);
    const auto original = qc;
    qc.adjoint().adjoint();
    EXPECT_EQ(original, qc);

    QCircuit inverse(1);
    inverse.gate(gt.H, 0).gate(gt.S, 0);
    inverse.adjoint();
    QCircuit combined(1);
    combined.gate(gt.H, 0).gate(gt.S, 0);
    combined.compose_circuit(inverse, 0);
    QEngine engine{combined};
    engine.execute();
    EXPECT_NEAR(0, norm(engine.get_state() - st.z0), 1e-12);
}

/// BEGIN iterator QCircuit::begin()
TEST(qpp_QCircuit_begin, Iterator) {
    // qutrits
    QCircuit q3(1, 0, 3);
    q3.gate(gt.Xd(3), 0);
    EXPECT_NE(q3.begin(), q3.end());
    // qubits
    QCircuit qc(1);
    EXPECT_EQ(qc.begin(), qc.end());
    qc.gate(gt.X, 0);
    EXPECT_NE(qc.begin(), qc.end());
}

/// BEGIN const_iterator QCircuit::begin() const noexcept
TEST(qpp_QCircuit_begin, ConstIterator) {
    // qutrits
    QCircuit tmp3(1, 0, 3);
    tmp3.gate(gt.Xd(3), 0);
    const QCircuit& q3 = tmp3;
    EXPECT_NE(q3.begin(), q3.end());
    // qubits
    QCircuit qc(1);
    qc.gate(gt.X, 0);
    const QCircuit& ref = qc;
    EXPECT_NE(ref.begin(), ref.end());
}

/// BEGIN const_iterator QCircuit::cbegin() const noexcept
TEST(qpp_QCircuit_cbegin, AllTests) {
    // qutrits
    QCircuit q3(1, 0, 3);
    q3.gate(gt.Xd(3), 0);
    EXPECT_NE(q3.cbegin(), q3.cend());
    // qubits
    QCircuit qc(1);
    EXPECT_EQ(qc.cbegin(), qc.cend());
    qc.gate(gt.X, 0);
    EXPECT_NE(qc.cbegin(), qc.cend());
}

/// BEGIN QCircuit& QCircuit::cCTRL(const cmat& U,
///       const std::vector<idx>& ctrl_dits, const std::vector<idx>& target,
///       std::optional<std::vector<idx>> shift = std::nullopt,
///       std::optional<std::string> name = std::nullopt)
TEST(qpp_QCircuit_cCTRL, MultipleCtrlsMultipleTargets) {
    // qutrits
    QCircuit q3(3, 2, 3);
    q3.gate(gt.Xd(3), 0).gate(gt.Xd(3), 0).gate(gt.Xd(3), 1).gate(gt.Xd(3), 1);
    q3.measure({0, 1}, 0, false);
    q3.cCTRL(kron(gt.Xd(3), gt.Xd(3)), {0, 1}, std::vector<idx>{1, 2});
    QEngine e3{q3};
    e3.execute();
    EXPECT_NEAR(0, norm(e3.get_state() - mket({2, 1, 2}, 3)), 1e-12);
    // qubits
    QCircuit qc(3, 2);
    qc.gate(gt.X, 0).gate(gt.X, 1);
    qc.measure({0, 1}, 0, false);
    qc.cCTRL(gt.CNOT, {0, 1}, std::vector<idx>{1, 2});
    EXPECT_EQ(3u, qc.get_gate_count());
    QEngine engine{qc};
    engine.execute();
    EXPECT_EQ((std::vector<idx>{1, 1}), engine.get_dits());
    EXPECT_NEAR(0, norm(engine.get_state() - mket({1, 1, 1})), 1e-12);
}

/// BEGIN QCircuit& QCircuit::cCTRL(const cmat& U,
///       const std::vector<idx>& ctrl_dits, idx target,
///       std::optional<std::vector<idx>> shift = std::nullopt,
///       std::optional<std::string> name = std::nullopt)
TEST(qpp_QCircuit_cCTRL, MultipleCtrlsSingleTarget) {
    // qutrits
    QCircuit q3(2, 2, 3);
    q3.gate(gt.Xd(3), 0)
        .gate(gt.Xd(3), 0)
        .measure(0, 0, false)
        .measure(0, 1, false);
    q3.cCTRL(gt.Xd(3), std::vector<idx>{0, 1}, 1);
    QEngine e3{q3};
    e3.execute();
    EXPECT_NEAR(0, norm(e3.get_state() - mket({2, 2}, 3)), 1e-12);
    // qubits
    QCircuit qc(2, 2);
    qc.gate(gt.X, 0);
    qc.measure(0, 0, false).measure(0, 1, false);
    qc.cCTRL(gt.X, std::vector<idx>{0, 1}, 1);
    EXPECT_EQ(2u, qc.get_gate_count());
    QEngine engine{qc};
    engine.execute();
    EXPECT_EQ((std::vector<idx>{1, 1}), engine.get_dits());
    EXPECT_NEAR(0, norm(engine.get_state() - mket({1, 1})), 1e-12);
}

/// BEGIN QCircuit& QCircuit::cCTRL(const cmat& U, idx ctrl_dit,
///       const std::vector<idx>& target,
///       std::optional<idx> shift = std::nullopt,
///       std::optional<std::string> name = std::nullopt)
TEST(qpp_QCircuit_cCTRL, SingleCtrlMultipleTargets) {
    // qutrits
    QCircuit q3(3, 1, 3);
    q3.gate(gt.Xd(3), 0).gate(gt.Xd(3), 0).measure(0, 0, false);
    q3.cCTRL(kron(gt.Xd(3), gt.Xd(3)), 0, std::vector<idx>{1, 2});
    QEngine e3{q3};
    e3.execute();
    EXPECT_NEAR(0, norm(e3.get_state() - mket({2, 2, 2}, 3)), 1e-12);
    // qubits
    QCircuit qc(3, 1);
    qc.gate(gt.X, 0).gate(gt.X, 1).measure(0, 0, false);
    qc.cCTRL(gt.CNOT, 0, std::vector<idx>{1, 2});
    EXPECT_EQ(3u, qc.get_gate_count());
    QEngine engine{qc};
    engine.execute();
    EXPECT_EQ((std::vector<idx>{1}), engine.get_dits());
    EXPECT_NEAR(0, norm(engine.get_state() - mket({1, 1, 1})), 1e-12);
}

/// BEGIN QCircuit& QCircuit::cCTRL(const cmat& U, idx ctrl_dit, idx target,
///       std::optional<idx> shift = std::nullopt,
///       std::optional<std::string> name = std::nullopt)
TEST(qpp_QCircuit_cCTRL, SingleCtrlSingleTarget) {
    // qutrits
    QCircuit q3(2, 1, 3);
    q3.gate(gt.Xd(3), 0)
        .gate(gt.Xd(3), 0)
        .measure(0, 0, false)
        .cCTRL(gt.Xd(3), 0, 1);
    QEngine e3{q3};
    e3.execute();
    EXPECT_NEAR(0, norm(e3.get_state() - mket({2, 2}, 3)), 1e-12);
    // qubits
    QCircuit qc(2, 1);
    qc.gate(gt.X, 0).measure(0, 0, false).cCTRL(gt.X, 0, 1);
    EXPECT_EQ(2u, qc.get_gate_count());
    QEngine engine{qc};
    engine.execute();
    EXPECT_EQ((std::vector<idx>{1}), engine.get_dits());
    EXPECT_NEAR(0, norm(engine.get_state() - mket({1, 1})), 1e-12);
}

/// BEGIN QCircuit& QCircuit::cCTRL_fan(const cmat& U,
///       const std::vector<idx>& ctrl_dits, const std::vector<idx>& target,
///       std::optional<std::vector<idx>> shift = std::nullopt,
///       std::optional<std::string> name = std::nullopt)
TEST(qpp_QCircuit_cCTRL_fan, MultipleCtrlsMultipleTargets) {
    // qubits
    QCircuit qc(3, 2);
    qc.gate(gt.X, 0).measure(0, 0, false).measure(0, 1, false);
    qc.cCTRL_fan(gt.X, {0, 1}, {1, 2});
    EXPECT_EQ(3u, qc.get_gate_count());
    QEngine engine{qc};
    engine.execute();
    EXPECT_EQ((std::vector<idx>{1, 1}), engine.get_dits());
    EXPECT_NEAR(0, norm(engine.get_state() - mket({1, 1, 1})), 1e-12);

    // qutrits: both classical controls are 2, so Xd(3)^2 acts on each target
    QCircuit qc3(3, 2, 3);
    qc3.gate(gt.Xd(3), 0).gate(gt.Xd(3), 0);
    qc3.measure(0, 0, false).measure(0, 1, false);
    qc3.cCTRL_fan(gt.Xd(3), {0, 1}, {1, 2});
    QEngine engine3{qc3};
    engine3.execute();
    EXPECT_EQ((std::vector<idx>{2, 2}), engine3.get_dits());
    EXPECT_NEAR(0, norm(engine3.get_state() - mket({2, 2, 2}, 3)), 1e-12);
}

/// BEGIN QCircuit& QCircuit::cCTRL_fan(const cmat& U, idx ctrl_dit,
///       const std::vector<idx>& target,
///       std::optional<idx> shift = std::nullopt,
///       std::optional<std::string> name = std::nullopt)
TEST(qpp_QCircuit_cCTRL_fan, SingleCtrlMultipleTargets) {
    // qubits
    QCircuit qc(3, 1);
    qc.gate(gt.X, 0).measure(0, 0, false);
    qc.cCTRL_fan(gt.X, 0, {1, 2});
    EXPECT_EQ(3u, qc.get_gate_count());
    QEngine engine{qc};
    engine.execute();
    EXPECT_NEAR(0, norm(engine.get_state() - mket({1, 1, 1})), 1e-12);

    // qutrits
    QCircuit qc3(3, 1, 3);
    qc3.gate(gt.Xd(3), 0).gate(gt.Xd(3), 0);
    qc3.measure(0, 0, false);
    qc3.cCTRL_fan(gt.Xd(3), 0, {1, 2});
    QEngine engine3{qc3};
    engine3.execute();
    EXPECT_EQ((std::vector<idx>{2}), engine3.get_dits());
    EXPECT_NEAR(0, norm(engine3.get_state() - mket({2, 2, 2}, 3)), 1e-12);
}

/// BEGIN const_iterator QCircuit::cend() const noexcept
TEST(qpp_QCircuit_cend, AllTests) {
    // qutrits
    QCircuit q3(1, 0, 3);
    EXPECT_EQ(q3.cbegin(), q3.cend());
    q3.gate(gt.Xd(3), 0);
    EXPECT_NE(q3.cbegin(), q3.cend());
    // qubits
    QCircuit qc(1);
    EXPECT_EQ(qc.cbegin(), qc.cend());
    qc.gate(gt.X, 0);
    EXPECT_NE(qc.cbegin(), qc.cend());
}

/// BEGIN QCircuit& QCircuit::compose_circuit(QCircuit other, bigint pos_qudit,
///       std::optional<idx> pos_dit = std::nullopt)
TEST(qpp_QCircuit_compose_circuit, AllTests) {
    // qutrits
    QCircuit q3(2, 0, 3), o3(1, 0, 3);
    o3.gate(gt.Xd(3), 0);
    q3.compose_circuit(o3, 1);
    QEngine e3{q3};
    e3.execute();
    EXPECT_NEAR(0, norm(e3.get_state() - mket({0, 1}, 3)), 1e-12);
    // qubits
    QCircuit qc(2), other(1);
    other.gate(gt.X, 0);
    qc.compose_circuit(other, 1);
    EXPECT_EQ(2u, qc.get_nq());
    EXPECT_EQ(1u, qc.get_gate_count());
    QEngine engine{qc};
    engine.execute();
    EXPECT_NEAR(0, norm(engine.get_state() - mket({0, 1})), 1e-12);
}

/// BEGIN QCircuit& compose_CTRL_circuit(const std::vector<idx>& ctrl,
///       QCircuit qc_target, bigint pos_qudit,
///       std::optional<std::vector<idx>> shift = std::nullopt,
///       std::optional<idx> pos_dit = std::nullopt)
TEST(qpp_QCircuit_compose_CTRL_circuit, AllTests) {
    // qutrits
    QCircuit q3(2, 0, 3), o3(1, 0, 3);
    q3.gate(gt.Xd(3), 0).gate(gt.Xd(3), 0);
    o3.gate(gt.Xd(3), 0);
    q3.compose_CTRL_circuit({0}, o3, 1);
    QEngine e3{q3};
    e3.execute();
    EXPECT_NEAR(0, norm(e3.get_state() - mket({2, 2}, 3)), 1e-12);
    // qubits
    QCircuit qc(2), other(1);
    qc.gate(gt.X, 0);
    other.gate(gt.X, 0);
    qc.compose_CTRL_circuit({0}, other, 1);
    EXPECT_EQ(2u, qc.get_nq());
    EXPECT_EQ(2u, qc.get_gate_count());
    QEngine engine{qc};
    engine.execute();
    EXPECT_NEAR(0, norm(engine.get_state() - mket({1, 1})), 1e-12);
}

/// BEGIN QCircuit& QCircuit::compress()
TEST(qpp_QCircuit_compress, AllTests) {
    // qutrits
    QCircuit q3(3, 0, 3);
    q3.gate(gt.Xd(3), 1).compress();
    EXPECT_EQ(1u, q3.get_nq());
    QEngine e3{q3};
    e3.execute();
    EXPECT_NEAR(0, norm(e3.get_state() - mket({1}, 3)), 1e-12);
    // qubits
    QCircuit qc(3, 2);
    qc.gate(gt.X, 1).compress();
    EXPECT_EQ(1u, qc.get_nq());
    EXPECT_EQ(1u, qc.get_gate_count());
    QEngine engine{qc};
    engine.execute();
    EXPECT_NEAR(0, norm(engine.get_state() - st.z1), 1e-12);
}

/// BEGIN QCircuit& QCircuit::cond_if(cond_pred_t pred)
TEST(qpp_QCircuit_cond_if, AllTests) {
    // qutrits
    QCircuit q3(2, 1, 3);
    q3.gate(gt.Xd(3), 0).gate(gt.Xd(3), 0).measure(0, 0, false);
    q3.cond_if([](auto dits) { return dits[0] == 2; })
        .gate(gt.Xd(3), 1)
        .cond_end();
    QEngine e3{q3};
    e3.execute();
    EXPECT_NEAR(0, norm(e3.get_state() - mket({2, 1}, 3)), 1e-12);
    // qubits
    QCircuit qc(2, 1);
    qc.measure(0, 0, false);
    auto pred = [](auto dits) { return dits[0] == 1; };
    qc.cond_if(pred);
    qc.gate(gt.X, 1);
    qc.cond_end();
    EXPECT_TRUE(qc.has_runtime_steps());
    EXPECT_TRUE(qc.validate_conditionals());
    QEngine engine{qc};
    engine.execute();
    EXPECT_NEAR(0, norm(engine.get_state() - mket({0, 0})), 1e-12);

    QCircuit true_branch(2, 1);
    true_branch.gate(gt.X, 0).measure(0, 0, false);
    true_branch.cond_if(pred).gate(gt.X, 1).cond_end();
    QEngine true_engine{true_branch};
    true_engine.execute();
    EXPECT_EQ((std::vector<idx>{1}), true_engine.get_dits());
    EXPECT_NEAR(0, norm(true_engine.get_state() - mket({1, 1})), 1e-12);
}

/// BEGIN QCircuit& QCircuit::cond_else()
TEST(qpp_QCircuit_cond_else, AllTests) {
    // qutrits
    QCircuit q3(2, 1, 3);
    q3.measure(0, 0, false);
    q3.cond_if([](auto dits) { return dits[0] == 2; })
        .gate(gt.Xd(3), 0)
        .cond_else()
        .gate(gt.Xd(3), 1)
        .cond_end();
    QEngine e3{q3};
    e3.execute();
    EXPECT_NEAR(0, norm(e3.get_state() - mket({0, 1}, 3)), 1e-12);
    // qubits
    QCircuit qc(3, 1);
    qc.measure(0, 0, false);
    auto pred = [](auto dits) { return dits[0] == 1; };
    qc.cond_if(pred);
    qc.gate(gt.X, 1);
    qc.cond_else();
    qc.gate(gt.X, 2);
    qc.cond_end();
    EXPECT_TRUE(qc.has_runtime_steps());
    EXPECT_TRUE(qc.validate_conditionals());
    QEngine false_engine{qc};
    false_engine.execute();
    EXPECT_NEAR(0, norm(false_engine.get_state() - mket({0, 0, 1})), 1e-12);

    QCircuit true_branch(3, 1);
    true_branch.gate(gt.X, 0).measure(0, 0, false);
    true_branch.cond_if(pred).gate(gt.X, 1);
    true_branch.cond_else().gate(gt.X, 2).cond_end();
    QEngine true_engine{true_branch};
    true_engine.execute();
    EXPECT_NEAR(0, norm(true_engine.get_state() - mket({1, 1, 0})), 1e-12);
}

/// BEGIN QCircuit& QCircuit::cond_while(cond_pred_t pred)
TEST(qpp_QCircuit_cond_while, AllTests) {
    // qutrits
    QCircuit q3(1, 1, 3);
    q3.cond_while([](auto dits) { return dits[0] < 2; })
        .gate(gt.Xd(3), 0)
        .set_dits_runtime([](auto dits) { ++dits[0]; })
        .cond_end();
    QEngine e3{q3};
    e3.execute();
    EXPECT_NEAR(0, norm(e3.get_state() - mket({2}, 3)), 1e-12);
    // qubits
    QCircuit qc(2, 2);
    auto pred = [](auto dits) { return !(dits[0] == 1 && dits[1] == 1); };
    qc.cond_while(pred);
    qc.reset({0, 1});
    qc.gate_fan(gt.H, {0, 1});
    qc.measure({0, 1}, 0, false);
    qc.cond_end();
    EXPECT_TRUE(qc.has_runtime_steps());
    EXPECT_TRUE(qc.validate_conditionals());
    QCircuit terminating(1, 1);
    terminating.cond_while([](auto dits) { return dits[0] == 0; });
    terminating.gate(gt.X, 0);
    terminating.set_dits_runtime([](auto dits) { dits[0] = 1; });
    terminating.cond_end();
    QEngine engine{terminating};
    engine.execute();
    EXPECT_EQ((std::vector<idx>{1}), engine.get_dits());
    EXPECT_NEAR(0, norm(engine.get_state() - st.z1), 1e-12);
}

/// BEGIN QCircuit& QCircuit::cond_end()
TEST(qpp_QCircuit_cond_end, AllTests) {
    // qutrits
    QCircuit q3(1, 1, 3);
    q3.cond_if([](auto dits) { return dits[0] == 0; }).cond_end();
    EXPECT_TRUE(q3.validate_conditionals());
    // qubits
    QCircuit qc(1);
    EXPECT_THROW(qc.cond_end(), exception::InvalidConditional);
    qc.cond_if([](auto) { return true; }).cond_end();
    EXPECT_TRUE(qc.validate_conditionals());
}

/// BEGIN QCircuit& QCircuit::couple_circuit_left(QCircuit other,
///       const std::vector<idx>& target,
///       std::optional<idx> pos_dit = std::nullopt)
TEST(qpp_QCircuit_couple_circuit_left, AllTests) {
    // qutrits
    QCircuit q3(2, 0, 3), o3(1, 0, 3);
    o3.gate(gt.Xd(3), 0);
    q3.couple_circuit_left(o3, {1});
    QEngine e3{q3};
    e3.execute();
    EXPECT_NEAR(0, norm(e3.get_state() - mket({0, 1}, 3)), 1e-12);
    // qubits
    QCircuit qc(2), other(1);
    other.gate(gt.X, 0);
    qc.couple_circuit_left(other, {1});
    EXPECT_EQ(2u, qc.get_nq());
    EXPECT_EQ(1u, qc.get_gate_count());
    QEngine engine{qc};
    engine.execute();
    EXPECT_NEAR(0, norm(engine.get_state() - mket({0, 1})), 1e-12);
}

/// BEGIN QCircuit& QCircuit::couple_circuit_right(QCircuit other,
///       const std::vector<idx>& target,
///       std::optional<idx> pos_dit = std::nullopt)
TEST(qpp_QCircuit_couple_circuit_right, AllTests) {
    // qutrits
    QCircuit q3(2, 0, 3), o3(1, 0, 3);
    o3.gate(gt.Xd(3), 0);
    q3.couple_circuit_right(o3, {1});
    QEngine e3{q3};
    e3.execute();
    EXPECT_NEAR(0, norm(e3.get_state() - mket({0, 1}, 3)), 1e-12);
    // qubits
    QCircuit qc(2), other(1);
    other.gate(gt.X, 0);
    qc.couple_circuit_right(other, {1});
    EXPECT_EQ(2u, qc.get_nq());
    EXPECT_EQ(1u, qc.get_gate_count());
    QEngine engine{qc};
    engine.execute();
    EXPECT_NEAR(0, norm(engine.get_state() - mket({0, 1})), 1e-12);
}

/// BEGIN QCircuit& QCircuit::CTRL(const cmat& U,
///       const std::vector<idx>& ctrl, const std::vector<idx>& target,
///       std::optional<std::vector<idx>> shift = std::nullopt,
///       std::optional<std::string> name = std::nullopt)
TEST(qpp_QCircuit_CTRL, MultipleCtrlsMultipleTargets) {
    // qutrits
    QCircuit q3(4, 0, 3);
    q3.gate(gt.Xd(3), 0).gate(gt.Xd(3), 0).gate(gt.Xd(3), 1).gate(gt.Xd(3), 1);
    q3.CTRL(kron(gt.Xd(3), gt.Xd(3)), {0, 1}, std::vector<idx>{2, 3});
    QEngine e3{q3};
    e3.execute();
    EXPECT_NEAR(0, norm(e3.get_state() - mket({2, 2, 2, 2}, 3)), 1e-12);
    // qubits
    QCircuit qc(4);
    qc.gate(gt.X, 0).gate(gt.X, 1).gate(gt.X, 2);
    qc.CTRL(gt.CNOT, {0, 1}, std::vector<idx>{2, 3});
    EXPECT_EQ(4u, qc.get_gate_count());
    QEngine engine{qc};
    engine.execute();
    EXPECT_NEAR(0, norm(engine.get_state() - mket({1, 1, 1, 1})), 1e-12);
}

/// BEGIN QCircuit& QCircuit::CTRL(const cmat& U,
///       const std::vector<idx>& ctrl, idx target,
///       std::optional<std::vector<idx>> shift = std::nullopt,
///       std::optional<std::string> name = std::nullopt)
TEST(qpp_QCircuit_CTRL, MultipleCtrlsSingleTarget) {
    // qutrits
    QCircuit q3(3, 0, 3);
    q3.gate(gt.Xd(3), 0).gate(gt.Xd(3), 0).gate(gt.Xd(3), 1).gate(gt.Xd(3), 1);
    q3.CTRL(gt.Xd(3), std::vector<idx>{0, 1}, 2);
    QEngine e3{q3};
    e3.execute();
    EXPECT_NEAR(0, norm(e3.get_state() - mket({2, 2, 2}, 3)), 1e-12);
    // qubits
    QCircuit qc(3);
    qc.gate(gt.X, 0).gate(gt.X, 1);
    qc.CTRL(gt.X, std::vector<idx>{0, 1}, 2);
    EXPECT_EQ(3u, qc.get_gate_count());
    QEngine engine{qc};
    engine.execute();
    EXPECT_NEAR(0, norm(engine.get_state() - mket({1, 1, 1})), 1e-12);
}

/// BEGIN QCircuit& QCircuit::CTRL(const cmat& U, idx ctrl,
///       const std::vector<idx>& target,
///       std::optional<idx> shift = std::nullopt,
///       std::optional<std::string> name = std::nullopt)
TEST(qpp_QCircuit_CTRL, SingleCtrlMultipleTargets) {
    // qutrits
    QCircuit q3(3, 0, 3);
    q3.gate(gt.Xd(3), 0).gate(gt.Xd(3), 0);
    q3.CTRL(kron(gt.Xd(3), gt.Xd(3)), 0, std::vector<idx>{1, 2});
    QEngine e3{q3};
    e3.execute();
    EXPECT_NEAR(0, norm(e3.get_state() - mket({2, 2, 2}, 3)), 1e-12);
    // qubits
    QCircuit qc(3);
    qc.gate(gt.X, 0).gate(gt.X, 1);
    qc.CTRL(gt.CNOT, 0, std::vector<idx>{1, 2});
    EXPECT_EQ(3u, qc.get_gate_count());
    QEngine engine{qc};
    engine.execute();
    EXPECT_NEAR(0, norm(engine.get_state() - mket({1, 1, 1})), 1e-12);
}

/// BEGIN QCircuit& QCircuit::CTRL(const cmat& U, idx ctrl, idx target,
///       std::optional<idx> shift = std::nullopt,
///       std::optional<std::string> name = std::nullopt)
TEST(qpp_QCircuit_CTRL, SingleCtrlSingleTarget) {
    // qutrits
    QCircuit q3(2, 0, 3);
    q3.gate(gt.Xd(3), 0).gate(gt.Xd(3), 0).CTRL(gt.Xd(3), 0, 1);
    QEngine e3{q3};
    e3.execute();
    EXPECT_NEAR(0, norm(e3.get_state() - mket({2, 2}, 3)), 1e-12);
    // qubits
    QCircuit qc(2);
    qc.gate(gt.X, 0).CTRL(gt.X, 0, 1);
    EXPECT_EQ(2u, qc.get_gate_count());
    QEngine engine{qc};
    engine.execute();
    EXPECT_NEAR(0, norm(engine.get_state() - mket({1, 1})), 1e-12);
}

/// BEGIN QCircuit& QCircuit::CTRL_fan(const cmat& U,
///       const std::vector<idx>& ctrl, const std::vector<idx>& target,
///       std::optional<std::vector<idx>> shift = std::nullopt,
///       std::optional<std::string> name = std::nullopt)
TEST(qpp_QCircuit_CTRL_fan, MultipleCtrlsMultipleTargets) {
    // qubits
    QCircuit qc(4);
    qc.gate(gt.X, 0).gate(gt.X, 1);
    qc.CTRL_fan(gt.X, {0, 1}, {2, 3});
    EXPECT_EQ(4u, qc.get_gate_count());
    QEngine engine{qc};
    engine.execute();
    std::cout << disp(dirac(engine.get_state())) << '\n';
    EXPECT_NEAR(0, norm(engine.get_state() - mket({1, 1, 1, 1})), 1e-12);

    // qutrits: controls in |2> apply Xd(3)^2 to each target
    QCircuit qc3(4, 0, 3);
    qc3.gate(gt.Xd(3), 0).gate(gt.Xd(3), 0);
    qc3.gate(gt.Xd(3), 1).gate(gt.Xd(3), 1);
    qc3.CTRL_fan(gt.Xd(3), {0, 1}, {2, 3});
    QEngine engine3{qc3};
    engine3.execute();
    EXPECT_NEAR(0, norm(engine3.get_state() - mket({2, 2, 2, 2}, 3)), 1e-12);
}

/// BEGIN QCircuit& QCircuit::CTRL_fan(const cmat& U, idx ctrl,
///       const std::vector<idx>& target,
///       std::optional<idx> shift = std::nullopt,
///       std::optional<std::string> name = std::nullopt)
TEST(qpp_QCircuit_CTRL_fan, SingleCtrlMultipleTargets) {
    // qubits
    QCircuit qc(3);
    qc.gate(gt.X, 0).CTRL_fan(gt.X, 0, {1, 2});
    EXPECT_EQ(3u, qc.get_gate_count());
    QEngine engine{qc};
    engine.execute();
    std::cout << disp(dirac(engine.get_state())) << '\n';
    EXPECT_NEAR(0, norm(engine.get_state() - mket({1, 1, 1})), 1e-12);

    // qutrits
    QCircuit qc3(3, 0, 3);
    qc3.gate(gt.Xd(3), 0).gate(gt.Xd(3), 0);
    qc3.CTRL_fan(gt.Xd(3), 0, {1, 2});
    QEngine engine3{qc3};
    engine3.execute();
    EXPECT_NEAR(0, norm(engine3.get_state() - mket({2, 2, 2}, 3)), 1e-12);
}

/// BEGIN QCircuit& QCircuit::discard(const std::vector<idx>& target,
///       std::optional<std::string> name = "discard")
TEST(qpp_QCircuit_discard, MultipleTargets) {
    // qutrits
    QCircuit q3(3, 0, 3);
    q3.gate(gt.Xd(3), 1).discard({0, 2});
    QEngine e3{q3};
    e3.execute();
    EXPECT_NEAR(0, norm(e3.get_state() - mket({1}, 3)), 1e-12);
    // qubits
    QCircuit qc(3);
    qc.discard({0, 2});
    EXPECT_TRUE(qc.removes_qudits());
    EXPECT_EQ(1u, qc.get_step_count());
}

/// BEGIN QCircuit& QCircuit::discard(idx target,
///       std::optional<std::string> name = "discard"")
TEST(qpp_QCircuit_discard, SingleTarget) {
    // qutrits
    QCircuit q3(2, 0, 3);
    q3.gate(gt.Xd(3), 0).discard(1);
    QEngine e3{q3};
    e3.execute();
    EXPECT_NEAR(0, norm(e3.get_state() - mket({1}, 3)), 1e-12);
    // qubits
    QCircuit qc(2);
    qc.gate(gt.X, 0).discard(1);
    EXPECT_TRUE(qc.removes_qudits());
    EXPECT_EQ(2u, qc.get_step_count());
    QEngine engine{qc};
    engine.execute();
    EXPECT_NEAR(0, norm(engine.get_state() - st.z1), 1e-12);
}

/// BEGIN iterator QCircuit::end()
TEST(qpp_QCircuit_end, Iterator) {
    // qutrits
    QCircuit q3(1, 0, 3);
    EXPECT_EQ(q3.begin(), q3.end());
    q3.gate(gt.Xd(3), 0);
    EXPECT_NE(q3.begin(), q3.end());
    // qubits
    QCircuit qc(1);
    EXPECT_EQ(qc.begin(), qc.end());
    qc.gate(gt.X, 0);
    EXPECT_NE(qc.begin(), qc.end());
}

/// BEGIN const_iterator QCircuit::end() const noexcept
TEST(qpp_QCircuit_end, ConstIterator) {
    // qutrits
    QCircuit m3(1, 0, 3);
    m3.gate(gt.Xd(3), 0);
    const QCircuit& q3 = m3;
    EXPECT_NE(q3.begin(), q3.end());
    // qubits
    QCircuit qc(1);
    qc.gate(gt.X, 0);
    const QCircuit& ref = qc;
    EXPECT_NE(ref.begin(), ref.end());
}

/// BEGIN QCircuit& QCircuit::gate(const cmat& U, idx i, idx j, idx k,
///       std::optional<std::string> name = std::nullopt)
TEST(qpp_QCircuit_gate, ThreeQudits) {
    // qutrits
    QCircuit q3(3, 0, 3);
    q3.gate(cmat::Identity(27, 27), 0, 1, 2);
    QEngine e3{q3};
    e3.execute();
    EXPECT_NEAR(0, norm(e3.get_state() - mket({0, 0, 0}, 3)), 1e-12);
    // qubits
    QCircuit qc(3);
    qc.gate(cmat::Identity(8, 8), 0, 1, 2);
    EXPECT_EQ(1u, qc.get_gate_count());
    QEngine engine{qc};
    engine.execute();
    EXPECT_NEAR(0, norm(engine.get_state() - mket({0, 0, 0})), 1e-12);
}

/// BEGIN QCircuit& QCircuit::gate(const cmat& U, idx i, idx j,
///       std::optional<std::string> name = std::nullopt)
TEST(qpp_QCircuit_gate, TwoQudits) {
    // qutrits
    QCircuit q3(2, 0, 3);
    q3.gate(kron(gt.Xd(3), gt.Xd(3)), 0, 1);
    QEngine e3{q3};
    e3.execute();
    EXPECT_NEAR(0, norm(e3.get_state() - mket({1, 1}, 3)), 1e-12);
    // qubits
    QCircuit qc(2);
    qc.gate(gt.X, 0).gate(gt.CNOT, 0, 1);
    EXPECT_EQ(2u, qc.get_gate_count());
    QEngine engine{qc};
    engine.execute();
    EXPECT_NEAR(0, norm(engine.get_state() - mket({1, 1})), 1e-12);
}

/// BEGIN QCircuit& QCircuit::gate(const cmat& U, idx i,
///       std::optional<std::string> name = std::nullopt)
TEST(qpp_QCircuit_gate, SingleQudit) {
    // qutrits
    QCircuit q3(2, 0, 3);
    q3.gate(gt.Xd(3), 1);
    QEngine e3{q3};
    e3.execute();
    EXPECT_NEAR(0, norm(e3.get_state() - mket({0, 1}, 3)), 1e-12);
    // qubits
    QCircuit qc(2);
    qc.gate(gt.X, 1);
    EXPECT_EQ(1u, qc.get_gate_count());
    QEngine engine{qc};
    engine.execute();
    EXPECT_NEAR(0, norm(engine.get_state() - mket({0, 1})), 1e-12);
}

/// BEGIN QCircuit& QCircuit::gate(const cmat& U,
///       const std::vector<idx>& target,
///       std::optional<std::string> name = std::nullopt)
TEST(qpp_QCircuit_gate, JointQudits) {
    // qutrits
    QCircuit q3(3, 0, 3);
    q3.gate(kron(gt.Xd(3), gt.Xd(3)), {0, 2});
    QEngine e3{q3};
    e3.execute();
    EXPECT_NEAR(0, norm(e3.get_state() - mket({1, 0, 1}, 3)), 1e-12);
    // qubits
    QCircuit qc(3);
    qc.gate(gt.X, 0).gate(gt.CNOT, {0, 2});
    EXPECT_EQ(2u, qc.get_gate_count());
    QEngine engine{qc};
    engine.execute();
    EXPECT_NEAR(0, norm(engine.get_state() - mket({1, 0, 1})), 1e-12);
}

/// BEGIN QCircuit& QCircuit::gate_fan(const cmat& U,
///       const std::vector<idx>& target,
///       std::optional<std::string> name = std::nullopt)
TEST(qpp_QCircuit_gate_fan, SpecificQudits) {
    // qutrits
    QCircuit q3(3, 0, 3);
    q3.gate_fan(gt.Xd(3), {0, 2});
    QEngine e3{q3};
    e3.execute();
    EXPECT_NEAR(0, norm(e3.get_state() - mket({1, 0, 1}, 3)), 1e-12);
    // qubits
    QCircuit qc(3);
    qc.gate_fan(gt.X, {0, 2});
    EXPECT_EQ(2u, qc.get_gate_count());
    QEngine engine{qc};
    engine.execute();
    EXPECT_NEAR(0, norm(engine.get_state() - mket({1, 0, 1})), 1e-12);
}

/// BEGIN QCircuit& QCircuit::gate_fan(const cmat& U,
///       std::optional<std::string> name = std::nullopt)
TEST(qpp_QCircuit_gate_fan, AllQudits) {
    // qutrits
    QCircuit q3(3, 0, 3);
    q3.gate_fan(gt.Xd(3));
    QEngine e3{q3};
    e3.execute();
    EXPECT_NEAR(0, norm(e3.get_state() - mket({1, 1, 1}, 3)), 1e-12);
    // qubits
    QCircuit qc(3);
    qc.gate_fan(gt.X);
    EXPECT_EQ(3u, qc.get_gate_count());
    QEngine engine{qc};
    engine.execute();
    EXPECT_NEAR(0, norm(engine.get_state() - mket({1, 1, 1})), 1e-12);
}

/// BEGIN std::vector<idx> QCircuit::get_clean_dits() const
TEST(qpp_QCircuit_get_clean_dits, AllTests) {
    // qutrits
    QCircuit q3(1, 2, 3);
    q3.measure(0, 1);
    EXPECT_EQ((std::vector<idx>{0}), q3.get_clean_dits());
    // qubits
    QCircuit qc(2, 2);
    EXPECT_EQ((std::vector<idx>{0, 1}), qc.get_clean_dits());
    qc.measure(0, 1);
    EXPECT_EQ((std::vector<idx>{0}), qc.get_clean_dits());
}

/// BEGIN std::vector<idx> QCircuit::get_clean_qudits() const
TEST(qpp_QCircuit_get_clean_qudits, AllTests) {
    // qutrits
    QCircuit q3(2, 0, 3);
    q3.gate(gt.Xd(3), 1);
    EXPECT_EQ((std::vector<idx>{0}), q3.get_clean_qudits());
    // qubits
    QCircuit qc(2);
    EXPECT_EQ((std::vector<idx>{0, 1}), qc.get_clean_qudits());
    qc.gate(gt.X, 1);
    EXPECT_EQ((std::vector<idx>{0}), qc.get_clean_qudits());
}

/// BEGIN const std::unordered_map<std::size_t, cmat>&
///       QCircuit::get_cmat_hash_tbl() const noexcept
TEST(qpp_QCircuit_get_cmat_hash_tbl, AllTests) {
    // qutrits
    QCircuit q3(1, 0, 3);
    q3.gate(gt.Xd(3), 0);
    EXPECT_FALSE(q3.get_cmat_hash_tbl().empty());
    // qubits
    QCircuit qc(1);
    EXPECT_TRUE(qc.get_cmat_hash_tbl().empty());
    qc.gate(gt.X, 0);
    EXPECT_FALSE(qc.get_cmat_hash_tbl().empty());
}

/// BEGIN idx QCircuit::get_d() const noexcept
TEST(qpp_QCircuit_get_d, AllTests) {
    // qutrits
    EXPECT_EQ(3u, QCircuit(2, 0, 3).get_d());
    EXPECT_EQ(5u, QCircuit(1, 0, 5).get_d());
    // qubits
    EXPECT_EQ(2u, QCircuit(1).get_d());
    EXPECT_EQ(3u, QCircuit(2, 0, 3).get_d());
}

/// BEGIN idx QCircuit::get_depth() const
TEST(qpp_QCircuit_get_depth, AllTests) {
    // qutrits
    QCircuit q3(1, 1, 3);
    q3.gate(gt.Xd(3), 0).measure(0, 0);
    EXPECT_EQ(2u, q3.get_depth());
    // qubits
    QCircuit qc(1, 1);
    EXPECT_EQ(0u, qc.get_depth());
    qc.gate(gt.X, 0).measure(0, 0);
    EXPECT_EQ(2u, qc.get_depth());
}

/// BEGIN std::vector<idx> QCircuit::get_dirty_dits() const
TEST(qpp_QCircuit_get_dirty_dits, AllTests) {
    // qutrits
    QCircuit q3(1, 2, 3);
    q3.measure(0, 1);
    EXPECT_EQ((std::vector<idx>{1}), q3.get_dirty_dits());
    // qubits
    QCircuit qc(2, 2);
    EXPECT_TRUE(qc.get_dirty_dits().empty());
    qc.measure(0, 1);
    EXPECT_EQ((std::vector<idx>{1}), qc.get_dirty_dits());
}

/// BEGIN std::vector<idx> QCircuit::get_dirty_qudits() const
TEST(qpp_QCircuit_get_dirty_qudits, AllTests) {
    // qutrits
    QCircuit q3(2, 0, 3);
    q3.gate(gt.Xd(3), 1);
    EXPECT_EQ((std::vector<idx>{1}), q3.get_dirty_qudits());
    // qubits
    QCircuit qc(2);
    EXPECT_TRUE(qc.get_dirty_qudits().empty());
    qc.gate(gt.X, 1);
    EXPECT_EQ((std::vector<idx>{1}), qc.get_dirty_qudits());
}

/// BEGIN idx QCircuit::get_gate_count(std::optional<cmat> U = std::nullopt)
///       const
TEST(qpp_QCircuit_get_gate_count, TotalGateCount) {
    // qutrits
    QCircuit q3(2, 0, 3);
    q3.gate(gt.Xd(3), 0).gate(gt.Zd(3), 1);
    EXPECT_EQ(2u, q3.get_gate_count());
    // qubits
    QCircuit qc(2);
    EXPECT_EQ(0u, qc.get_gate_count());
    qc.gate(gt.X, 0).gate(gt.Z, 1);
    EXPECT_EQ(2u, qc.get_gate_count());
    QEngine engine{qc};
    engine.execute();
    EXPECT_NEAR(0, norm(engine.get_state() - mket({1, 0})), 1e-12);
}
TEST(qpp_QCircuit_get_gate_count, SpecificGateCount) {
    // qutrits
    QCircuit q3(2, 0, 3);
    q3.gate(gt.Xd(3), 0).gate(gt.Xd(3), 1).gate(gt.Zd(3), 0);
    EXPECT_EQ(2u, q3.get_gate_count(gt.Xd(3)));
    EXPECT_EQ(1u, q3.get_gate_count(gt.Zd(3)));
    // qubits
    QCircuit qc(2);
    qc.gate(gt.X, 0).gate(gt.X, 1).gate(gt.Z, 0);
    EXPECT_EQ(2u, qc.get_gate_count(gt.X));
    EXPECT_EQ(1u, qc.get_gate_count(gt.Z));
    QEngine engine{qc};
    engine.execute();
    EXPECT_NEAR(0, norm(engine.get_state() + mket({1, 1})), 1e-12);
}

/// BEGIN idx QCircuit::get_gate_depth(std::optional<cmat> U = std::nullopt)
///       const
TEST(qpp_QCircuit_get_gate_depth, TotalGateDepth) {
    // qutrits
    QCircuit q3(1, 0, 3);
    q3.gate(gt.Xd(3), 0).gate(gt.Zd(3), 0);
    EXPECT_EQ(2u, q3.get_gate_depth());
    // qubits
    QCircuit qc(1);
    qc.gate(gt.X, 0).gate(gt.Z, 0);
    EXPECT_EQ(2u, qc.get_gate_depth());
}
TEST(qpp_QCircuit_get_gate_depth, SpecificGateDepth) {
    // qutrits
    QCircuit q3(1, 0, 3);
    q3.gate(gt.Xd(3), 0).gate(gt.Zd(3), 0).gate(gt.Xd(3), 0);
    EXPECT_EQ(2u, q3.get_gate_depth(gt.Xd(3)));
    // qubits
    QCircuit qc(1);
    qc.gate(gt.X, 0).gate(gt.Z, 0).gate(gt.X, 0);
    EXPECT_EQ(2u, qc.get_gate_depth(gt.X));
    EXPECT_EQ(1u, qc.get_gate_depth(gt.Z));
}

/// BEGIN std::vector<idx> QCircuit::get_measured_d() const
TEST(qpp_QCircuit_get_measured_d, AllTests) {
    // qutrits
    QCircuit q3(2, 1, 3);
    q3.measure(0, 0);
    EXPECT_EQ((std::vector<idx>{0}), q3.get_measured_d());
    // qubits
    QCircuit qc(2, 1);
    qc.measure(0, 0);
    EXPECT_EQ((std::vector<idx>{0}), qc.get_measured_d());
}

/// BEGIN std::vector<idx> QCircuit::get_measured_nd() const
TEST(qpp_QCircuit_get_measured_nd, AllTests) {
    // qutrits
    QCircuit q3(2, 1, 3);
    q3.measure(1, 0, false);
    EXPECT_EQ((std::vector<idx>{1}), q3.get_measured_nd());
    // qubits
    QCircuit qc(2, 1);
    qc.measure(1, 0, false);
    EXPECT_EQ((std::vector<idx>{1}), qc.get_measured_nd());
}

/// BEGIN idx QCircuit::get_measurement_count(std::optional<cmat> V =
///       std::nullopt) const
TEST(qpp_QCircuit_get_measurement_count, TotalMeasurementCount) {
    // qutrits
    QCircuit q3(2, 2, 3);
    q3.measure(0, 0).measure(1, 1);
    EXPECT_EQ(2u, q3.get_measurement_count());
    // qubits
    QCircuit qc(2, 2);
    qc.measure(0, 0).measure(1, 1);
    EXPECT_EQ(2u, qc.get_measurement_count());
}
TEST(qpp_QCircuit_get_measurement_count, SpecificMeasurementCount) {
    // qutrits
    QCircuit q3(2, 2, 3);
    q3.measureV(gt.Fd(3), 0, 0).measure(1, 1);
    EXPECT_EQ(1u, q3.get_measurement_count(gt.Fd(3)));
    // qubits
    QCircuit qc(2, 2);
    qc.measureV(gt.H, 0, 0).measure(1, 1);
    EXPECT_EQ(1u, qc.get_measurement_count(gt.H));
}

/// BEGIN idx QCircuit::get_measurement_depth(std::optional<cmat> V =
///       std::nullopt) const
TEST(qpp_QCircuit_get_measurement_depth, TotalMeasurementDepth) {
    // qutrits
    QCircuit q3(2, 2, 3);
    q3.measure(0, 0).measure(1, 1);
    EXPECT_EQ(1u, q3.get_measurement_depth());
    // qubits
    QCircuit qc(2, 2);
    qc.measure(0, 0).measure(1, 1);
    EXPECT_EQ(1u, qc.get_measurement_depth());
}
TEST(qpp_QCircuit_get_measurement_depth, SpecificMeasurementDepth) {
    // qutrits
    QCircuit q3(1, 1, 3);
    q3.measureV(gt.Fd(3), 0, 0, false).measureV(gt.Fd(3), 0, 0, false);
    EXPECT_EQ(2u, q3.get_measurement_depth(gt.Fd(3)));
    // qubits
    QCircuit qc(1, 1);
    qc.measureV(gt.H, 0, 0, false).measureV(gt.H, 0, 0, false);
    EXPECT_EQ(2u, qc.get_measurement_depth(gt.H));
}

/// BEGIN std::vector<idx> QCircuit::get_measurement_dits() const
TEST(qpp_QCircuit_get_measurement_dits, AllTests) {
    // qutrits
    QCircuit q3(1, 2, 3);
    q3.measure(0, 1);
    EXPECT_EQ((std::vector<idx>{1}), q3.get_measurement_dits());
    // qubits
    QCircuit qc(2, 2);
    qc.measure(0, 1);
    EXPECT_EQ((std::vector<idx>{1}), qc.get_measurement_dits());
}

/// BEGIN std::optional<std::string> QCircuit::get_name() const
TEST(qpp_QCircuit_get_name, AllTests) {
    // qutrits
    QCircuit q3(1, 0, 3);
    q3.set_name("qutrit");
    ASSERT_TRUE(q3.get_name());
    EXPECT_EQ("qutrit", *q3.get_name());
    // qubits
    QCircuit qc(1);
    EXPECT_FALSE(qc.get_name().has_value());
    qc.set_name("demo");
    ASSERT_TRUE(qc.get_name());
    EXPECT_EQ("demo", *qc.get_name());
}

/// BEGIN idx QCircuit::get_nc() const noexcept
TEST(qpp_QCircuit_get_nc, AllTests) {
    // qutrits
    EXPECT_EQ(3u, QCircuit(1, 3, 3).get_nc());
    // qubits
    EXPECT_EQ(0u, QCircuit(2).get_nc());
    EXPECT_EQ(3u, QCircuit(2, 3).get_nc());
}

/// BEGIN std::vector<idx> QCircuit::get_non_measured_d() const
TEST(qpp_QCircuit_get_non_measured_d, AllTests) {
    // qutrits
    QCircuit q3(2, 1, 3);
    q3.measure(0, 0);
    EXPECT_EQ((std::vector<idx>{1}), q3.get_non_measured_d());
    // qubits
    QCircuit qc(2, 1);
    qc.measure(0, 0);
    EXPECT_EQ((std::vector<idx>{1}), qc.get_non_measured_d());
}

/// BEGIN idx QCircuit::get_nop_count() const
TEST(qpp_QCircuit_get_nop_count, AllTests) {
    // qutrits
    QCircuit q3(1, 0, 3);
    q3.nop().nop();
    EXPECT_EQ(2u, q3.get_nop_count());
    // qubits
    QCircuit qc(1);
    qc.nop().nop();
    EXPECT_EQ(2u, qc.get_nop_count());
}

/// BEGIN idx QCircuit::get_nq() const noexcept
TEST(qpp_QCircuit_get_nq, AllTests) {
    // qutrits
    EXPECT_EQ(4u, QCircuit(4, 0, 3).get_nq());
    // qubits
    EXPECT_EQ(1u, QCircuit().get_nq());
    EXPECT_EQ(4u, QCircuit(4).get_nq());
}

/// BEGIN QCircuit::Resources QCircuit::get_resources() const
TEST(qpp_QCircuit_get_resources, AllTests) {
    // qutrits
    QCircuit q3(1, 1, 3);
    q3.gate(gt.Xd(3), 0).measure(0, 0);
    const auto r3 = q3.get_resources();
    EXPECT_EQ(1u, r3.gate_count);
    EXPECT_EQ(1u, r3.measurement_count);
    // qubits
    QCircuit qc(1, 1);
    qc.gate(gt.X, 0).measure(0, 0);
    const auto resources = qc.get_resources();
    EXPECT_EQ(1u, resources.gate_count);
    EXPECT_EQ(1u, resources.measurement_count);
}

/// BEGIN idx QCircuit::get_step_count() const noexcept
TEST(qpp_QCircuit_get_step_count, AllTests) {
    // qutrits
    QCircuit q3(1, 0, 3);
    q3.gate(gt.Xd(3), 0).nop();
    EXPECT_EQ(2u, q3.get_step_count());
    // qubits
    QCircuit qc(1);
    EXPECT_EQ(0u, qc.get_step_count());
    qc.gate(gt.X, 0).nop();
    EXPECT_EQ(2u, qc.get_step_count());
}

/// BEGIN bool QCircuit::has_measurements() const noexcept
TEST(qpp_QCircuit_has_measurements, AllTests) {
    // qutrits
    QCircuit q3(1, 1, 3);
    EXPECT_FALSE(q3.has_measurements());
    q3.measure(0, 0);
    EXPECT_TRUE(q3.has_measurements());
    // qubits
    QCircuit qc(1, 1);
    EXPECT_FALSE(qc.has_measurements());
    qc.measure(0, 0);
    EXPECT_TRUE(qc.has_measurements());
}

/// BEGIN bool QCircuit::has_runtime_steps() const noexcept
TEST(qpp_QCircuit_has_runtime_steps, AllTests) {
    // qutrits
    QCircuit q3(1, 1, 3);
    q3.cond_if([](auto dits) { return dits[0] == 2; }).cond_end();
    EXPECT_TRUE(q3.has_runtime_steps());
    // qubits
    QCircuit qc(1);
    EXPECT_FALSE(qc.has_runtime_steps());
    qc.cond_if([](auto) { return true; }).cond_end();
    EXPECT_TRUE(qc.has_runtime_steps());
}

/// BEGIN inline static bool QCircuit::is_cCTRL(const GateStep& gate_step)
TEST(qpp_QCircuit_is_cCTRL, AllTests) {
    // qutrits
    QCircuit q3(2, 1, 3);
    q3.cCTRL(gt.Xd(3), 0, 1);
    EXPECT_TRUE(internal::is_cCTRL(q3.begin()));
    // qubits
    QCircuit qc(2, 1);
    qc.cCTRL(gt.X, 0, 1);
    EXPECT_TRUE(internal::is_cCTRL(qc.begin()));
    const bool is_quantum_control = std::visit(
        [](const auto& step) {
            using Step = std::decay_t<decltype(step)>;
            if constexpr (std::is_same_v<Step, internal::QCircuitGateStep>) {
                return QCircuit::is_CTRL(step);
            }
            return false;
        },
        qc.begin()->get_step());
    EXPECT_FALSE(is_quantum_control);
}

/// BEGIN bool QCircuit::is_clean_dit(idx i) const
TEST(qpp_QCircuit_is_clean_dit, AllTests) {
    // qutrits
    QCircuit q3(1, 2, 3);
    q3.measure(0, 0);
    EXPECT_FALSE(q3.is_clean_dit(0));
    EXPECT_TRUE(q3.is_clean_dit(1));
    // qubits
    QCircuit qc(1, 2);
    EXPECT_TRUE(qc.is_clean_dit(0));
    qc.measure(0, 0);
    EXPECT_FALSE(qc.is_clean_dit(0));
    EXPECT_TRUE(qc.is_clean_dit(1));
}

/// BEGIN bool QCircuit::is_clean_qudit(idx i) const
TEST(qpp_QCircuit_is_clean_qudit, AllTests) {
    // qutrits
    QCircuit q3(2, 0, 3);
    q3.gate(gt.Xd(3), 0);
    EXPECT_FALSE(q3.is_clean_qudit(0));
    EXPECT_TRUE(q3.is_clean_qudit(1));
    // qubits
    QCircuit qc(2);
    EXPECT_TRUE(qc.is_clean_qudit(0));
    qc.gate(gt.X, 0);
    EXPECT_FALSE(qc.is_clean_qudit(0));
    EXPECT_TRUE(qc.is_clean_qudit(1));
}

/// BEGIN inline static bool QCircuit::is_CTRL(const GateStep& gate_step)
TEST(qpp_QCircuit_is_CTRL, AllTests) {
    // qutrits
    QCircuit q3(2, 0, 3);
    q3.CTRL(gt.Xd(3), 0, 1);
    EXPECT_FALSE(internal::is_cCTRL(q3.begin()));
    // qubits
    QCircuit qc(2);
    qc.CTRL(gt.X, 0, 1);
    const bool is_quantum_control = std::visit(
        [](const auto& step) {
            using Step = std::decay_t<decltype(step)>;
            if constexpr (std::is_same_v<Step, internal::QCircuitGateStep>) {
                return QCircuit::is_CTRL(step);
            }
            return false;
        },
        qc.begin()->get_step());
    EXPECT_TRUE(is_quantum_control);
    EXPECT_FALSE(internal::is_cCTRL(qc.begin()));
}

/// BEGIN bool QCircuit::is_measurement_dit(idx i) const
TEST(qpp_QCircuit_is_measurement_dit, AllTests) {
    // qutrits
    QCircuit q3(1, 2, 3);
    q3.measure(0, 1);
    EXPECT_TRUE(q3.is_measurement_dit(1));
    // qubits
    QCircuit qc(1, 2);
    EXPECT_FALSE(qc.is_measurement_dit(1));
    qc.measure(0, 1);
    EXPECT_TRUE(qc.is_measurement_dit(1));
}

/// BEGIN QCircuit& QCircuit::kron(QCircuit qc)
TEST(qpp_QCircuit_kron, AllTests) {
    // qutrits
    QCircuit q3(1, 0, 3), o3(1, 0, 3);
    q3.gate(gt.Xd(3), 0);
    o3.gate(gt.Xd(3), 0);
    q3.kron(o3);
    QEngine e3{q3};
    e3.execute();
    EXPECT_NEAR(0, norm(e3.get_state() - mket({1, 1}, 3)), 1e-12);
    // qubits
    QCircuit qc(1), other(2);
    qc.gate(gt.X, 0);
    other.gate(gt.X, 1);
    qc.kron(other);
    EXPECT_EQ(3u, qc.get_nq());
    EXPECT_EQ(2u, qc.get_gate_count());
    QEngine engine{qc};
    engine.execute();
    EXPECT_NEAR(0, norm(engine.get_state() - mket({1, 0, 1})), 1e-12);
}

/// BEGIN QCircuit& QCircuit::measureV(const cmat& V,
///       const std::vector<idx>& target, idx c_reg, bool destructive = true,
///       std::optional<std::string> name = std::nullopt)
TEST(qpp_QCircuit_measureV, JointMultipleTargets) {
    // qutrits
    QCircuit q3(2, 2, 3);
    q3.gate_fan(gt.Xd(3)).measureV(cmat::Identity(9, 9), {0, 1}, 0, false);
    QEngine e3{q3};
    e3.execute();
    EXPECT_EQ(4u, e3.get_dits()[0]);
    EXPECT_NEAR(0, norm(e3.get_state() - mket({1, 1}, 3)), 1e-12);
    // qubits
    QCircuit qc(2, 2);
    qc.gate_fan(gt.X).measureV(cmat::Identity(4, 4), {0, 1}, 0, false);
    EXPECT_EQ(1u, qc.get_measurement_count());
    QEngine engine{qc};
    engine.execute();
    EXPECT_EQ(3u, engine.get_dits()[0]);
    EXPECT_NEAR(0, norm(engine.get_state() - mket({1, 1})), 1e-12);
}

/// BEGIN QCircuit& QCircuit::measureV(const cmat& V, idx target, idx c_reg,
///       bool destructive = true,
///       std::optional<std::string> name = std::nullopt)
TEST(qpp_QCircuit_measureV, SingleTarget) {
    // qutrits
    QCircuit q3(1, 1, 3);
    q3.gate(gt.Fd(3), 0).measureV(gt.Fd(3), 0, 0, false);
    QEngine e3{q3};
    e3.execute();
    EXPECT_EQ(0u, e3.get_dits()[0]);
    EXPECT_NEAR(0, norm(e3.get_state() - gt.Fd(3) * mket({0}, 3)), 1e-12);
    // qubits
    QCircuit qc(1, 1);
    qc.gate(gt.H, 0).measureV(gt.H, 0, 0, false);
    EXPECT_EQ(1u, qc.get_measurement_count());
    QEngine engine{qc};
    engine.execute();
    EXPECT_EQ((std::vector<idx>{0}), engine.get_dits());
    EXPECT_NEAR(0, norm(engine.get_state() - st.x0), 1e-12);
}

/// BEGIN QCircuit& QCircuit::measure(const std::vector<idx>& target,
///       idx c_reg = 0, bool destructive = true,
///       std::optional<std::string> name = "mZ")
TEST(qpp_QCircuit_measure, MultipleTargets) {
    // qutrits
    QCircuit q3(2, 2, 3);
    q3.gate(gt.Xd(3), 0).gate(gt.Xd(3), 0).measure({0, 1}, 0, false);
    QEngine e3{q3};
    e3.execute();
    EXPECT_EQ((std::vector<idx>{2, 0}), e3.get_dits());
    EXPECT_NEAR(0, norm(e3.get_state() - mket({2, 0}, 3)), 1e-12);
    // qubits
    QCircuit qc(2, 2);
    qc.gate(gt.X, 0).measure({0, 1}, 0, false);
    EXPECT_EQ(1u, qc.get_measurement_count());
    EXPECT_EQ(2u, qc.get_measured_nd().size());
    QEngine engine{qc};
    engine.execute();
    EXPECT_EQ((std::vector<idx>{1, 0}), engine.get_dits());
    EXPECT_NEAR(0, norm(engine.get_state() - mket({1, 0})), 1e-12);
}

/// BEGIN QCircuit& QCircuit::measure(idx target, idx c_reg,
///       bool destructive = true, std::optional<std::string> name = "mZ")
TEST(qpp_QCircuit_measure, SingleTarget) {
    // qutrits
    QCircuit q3(2, 1, 3);
    q3.gate(gt.Xd(3), 1).gate(gt.Xd(3), 1).measure(1, 0, false);
    QEngine e3{q3};
    e3.execute();
    EXPECT_EQ((std::vector<idx>{2}), e3.get_dits());
    // qubits
    QCircuit qc(2, 1);
    qc.gate(gt.X, 1).measure(1, 0, false);
    EXPECT_EQ((std::vector<idx>{1}), qc.get_measured_nd());
    QEngine engine{qc};
    engine.execute();
    EXPECT_EQ((std::vector<idx>{1}), engine.get_dits());
    EXPECT_NEAR(0, norm(engine.get_state() - mket({0, 1})), 1e-12);
}

/// BEGIN QCircuit& QCircuit::measure_all(idx c_reg = 0,
///       bool destructive = true, std::optional<std::string> name = "mZ")
TEST(qpp_QCircuit_measure_all, AllTests) {
    // qutrits
    QCircuit q3(2, 2, 3);
    q3.gate(gt.Xd(3), 0)
        .gate(gt.Xd(3), 1)
        .gate(gt.Xd(3), 1)
        .measure_all(0, false);
    QEngine e3{q3};
    e3.execute();
    EXPECT_EQ((std::vector<idx>{1, 2}), e3.get_dits());
    // qubits
    QCircuit qc(3, 3);
    qc.gate(gt.X, 0).gate(gt.X, 2).measure_all(0, false);
    EXPECT_TRUE(qc.has_measurements());
    EXPECT_EQ(3u, qc.get_measured_nd().size());
    QEngine engine{qc};
    engine.execute();
    EXPECT_EQ((std::vector<idx>{1, 0, 1}), engine.get_dits());
    EXPECT_NEAR(0, norm(engine.get_state() - mket({1, 0, 1})), 1e-12);
}

/// BEGIN QCircuit& QCircuit::nop()
TEST(qpp_QCircuit_nop, AllTests) {
    // qutrits
    QCircuit q3(1, 0, 3);
    q3.nop();
    QEngine e3{q3};
    e3.execute();
    EXPECT_NEAR(0, norm(e3.get_state() - mket({0}, 3)), 1e-12);
    // qubits
    QCircuit qc(1);
    qc.nop();
    EXPECT_EQ(1u, qc.get_nop_count());
    EXPECT_EQ(1u, qc.get_step_count());
}

/// BEGIN bool QCircuit::operator!=(const QCircuit& rhs) const noexcept
TEST(qpp_QCircuit_operator_noneq, AllTests) {
    // qutrits
    QCircuit a3(1, 0, 3), b3(1, 0, 3);
    b3.gate(gt.Xd(3), 0);
    EXPECT_TRUE(a3 != b3);
    // qubits
    QCircuit a(1), b(1);
    EXPECT_FALSE(a != b);
    b.gate(gt.X, 0);
    EXPECT_TRUE(a != b);
}

/// BEGIN bool QCircuit::operator==(const QCircuit& rhs) const noexcept
TEST(qpp_QCircuit_operator_eq, AllTests) {
    // qutrits
    QCircuit a3(1, 0, 3), b3(1, 0, 3);
    a3.gate(gt.Xd(3), 0);
    b3.gate(gt.Xd(3), 0);
    EXPECT_TRUE(a3 == b3);
    // qubits
    QCircuit a(1), b(1);
    EXPECT_TRUE(a == b);
    a.gate(gt.X, 0);
    EXPECT_FALSE(a == b);
    b.gate(gt.X, 0);
    EXPECT_TRUE(a == b);
}

/// BEGIN QCircuit& QCircuit::post_select(const std::vector<idx>& target,
///       const std::vector<idx>& ps_vals, idx c_reg, bool destructive = true,
///       std::optional<std::string> name = "pZ")
TEST(qpp_QCircuit_post_select, MultipleTargets) {
    // qutrits
    QCircuit q3(2, 2, 3);
    q3.gate(gt.Xd(3), 0).gate(gt.Xd(3), 0).gate(gt.Xd(3), 1);
    q3.post_select({0, 1}, {2, 1}, 0, false);
    QEngine e3{q3};
    e3.execute();
    EXPECT_EQ((std::vector<idx>{2, 1}), e3.get_dits());
    EXPECT_NEAR(0, norm(e3.get_state() - mket({2, 1}, 3)), 1e-12);
    // qubits
    QCircuit qc{2, 2};
    qc.gate_fan(gt.H, {0, 1});
    qc.post_select({0, 1}, {0, 1}, 0, false);

    QEngine qe{qc};
    qe.set_ensure_post_selection(true); // enforce post-selection
    qe.execute(2);

    auto dits = qe.get_dits();
    ket state = qe.get_state();

    std::vector<idx> expected_dits{0, 1};
    ket expected_state = 01_ket;

    ASSERT_EQ(dits, expected_dits);
    ASSERT_EQ(state, expected_state);
}

/// BEGIN QCircuit& QCircuit::post_select(idx target, idx ps_val, idx c_reg,
///       bool destructive = true, std::optional<std::string> name = "pZ")
TEST(qpp_QCircuit_post_select, SingleTarget) {
    // qutrits
    QCircuit q3(1, 1, 3);
    q3.gate(gt.Xd(3), 0).gate(gt.Xd(3), 0).post_select(0, 2, 0, false);
    QEngine e3{q3};
    e3.execute();
    EXPECT_EQ(2u, e3.get_dits()[0]);
    EXPECT_NEAR(0, norm(e3.get_state() - mket({2}, 3)), 1e-12);
    // qubits
    QCircuit qc{2, 2};
    qc.gate_fan(gt.H, {0, 1});
    qc.post_select(0, 1, 0, false);
    qc.post_select(1, 1, 1, false);

    QEngine qe{qc};
    qe.set_ensure_post_selection(true); // enforce post-selection
    qe.execute(2);

    auto dits = qe.get_dits();
    ket state = qe.get_state();

    std::vector<idx> expected_dits{1, 1};
    ket expected_state = 11_ket;

    ASSERT_EQ(dits, expected_dits);
    ASSERT_EQ(state, expected_state);
}

/// BEGIN QCircuit& QCircuit::post_selectV(const cmat& V,
///       const std::vector<idx>& target, idx ps_val, idx c_reg,
///       bool destructive = true,
///       std::optional<std::string> name = std::nullopt)
TEST(qpp_QCircuit_post_selectV, MultipleTargets) {
    // qutrits
    QCircuit q3(2, 2, 3);
    q3.gate_fan(gt.Fd(3));
    q3.post_selectV(kron(gt.Fd(3), gt.Fd(3)), {0, 1}, 0, 0, false);
    QEngine e3{q3};
    e3.execute();
    EXPECT_EQ(0u, e3.get_dits()[0]);
    EXPECT_NEAR(
        0, norm(e3.get_state() - kron(gt.Fd(3), gt.Fd(3)) * mket({0, 0}, 3)),
        1e-11);
    // qubits
    QCircuit qc{2, 2};
    cmat HH_basis = kron(gt.H, gt.H);
    qc.post_selectV(HH_basis, {0, 1}, 2, 0, false);

    QEngine qe{qc};
    qe.set_ensure_post_selection(true); // enforce post-selection
    qe.execute(2);

    auto dits = qe.get_dits();
    ket state = qe.get_state();

    std::vector<idx> expected_dits{2, 0};
    ket expected_state = kron(st.minus(), st.plus());

    ASSERT_EQ(dits, expected_dits);
    ASSERT_NEAR(0, norm(state - expected_state), 1e-5);
}

/// BEGIN QCircuit& QCircuit::post_selectV(const cmat& V, idx target,
///       idx ps_val, idx c_reg, bool destructive = true,
///       std::optional<std::string> name = std::nullopt)
TEST(qpp_QCircuit_post_selectV, SingleTarget) {
    // qutrits
    QCircuit q3(1, 1, 3);
    q3.gate(gt.Fd(3), 0).post_selectV(gt.Fd(3), 0, 0, 0, false);
    QEngine e3{q3};
    e3.execute();
    EXPECT_EQ(0u, e3.get_dits()[0]);
    // qubits
    QCircuit qc{2, 2};
    qc.post_selectV(gt.H, 0, 1, 0, false);
    qc.post_selectV(gt.H, 1, 1, 1, false);

    QEngine qe{qc};
    qe.set_ensure_post_selection(true); // enforce post-selection
    qe.execute(2);

    auto dits = qe.get_dits();
    ket state = qe.get_state();

    std::vector<idx> expected_dits{1, 1};
    ket expected_state = kron(st.minus(), st.minus());

    ASSERT_EQ(dits, expected_dits);
    ASSERT_NEAR(0, norm(state - expected_state), 1e-5);
}

/// BEGIN QCircuit& QCircuit::QFT(bool swap = true)
TEST(qpp_QCircuit_QFT, AllQudits) {
    // qutrits
    QCircuit q3(2, 0, 3);
    q3.QFT();
    QEngine e3{q3};
    e3.execute();
    EXPECT_NEAR(0, norm(e3.get_state() - ket::Ones(9) / 3.), 1e-11);
    // qubits
    QCircuit qc(2);
    qc.QFT();
    EXPECT_GT(qc.get_gate_count(), 0u);
    EXPECT_EQ(2u, qc.get_nq());
    QEngine engine{qc};
    engine.execute();
    const ket expected = ket::Ones(4) / 2.;
    EXPECT_NEAR(0, norm(engine.get_state() - expected), 1e-12);
}

/// BEGIN QCircuit& QCircuit::QFT(const std::vector<idx>& target,
///       bool swap = true)
TEST(qpp_QCircuit_QFT, SpecificQudits) {
    // qutrits
    QCircuit q3(3, 0, 3);
    q3.QFT({0, 2});
    QEngine e3{q3};
    e3.execute();
    EXPECT_NEAR(
        0,
        norm(e3.get_state() - kron(kron(gt.Fd(3) * mket({0}, 3), mket({0}, 3)),
                                   gt.Fd(3) * mket({0}, 3))),
        1e-11);
    // qubits
    QCircuit qc(3);
    qc.QFT({0, 2});
    EXPECT_GT(qc.get_gate_count(), 0u);
    EXPECT_TRUE(qc.is_clean_qudit(1));
    QEngine engine{qc};
    engine.execute();
    const ket expected = kron(kron(st.x0, st.z0), st.x0);
    EXPECT_NEAR(0, norm(engine.get_state() - expected), 1e-12);
}

/// BEGIN QCircuit& QCircuit::remove_clean_dit(idx target)
TEST(qpp_QCircuit_remove_clean_dit, AllTests) {
    // qutrits
    QCircuit q3(1, 3, 3);
    q3.remove_clean_dit(1);
    EXPECT_EQ(2u, q3.get_nc());
    // qubits
    QCircuit qc(2, 3);
    qc.remove_clean_dit(1);
    EXPECT_EQ(2u, qc.get_nc());
}

/// BEGIN QCircuit& QCircuit::remove_clean_dits(std::vector<idx> target)
TEST(qpp_QCircuit_remove_clean_dits, AllTests) {
    // qutrits
    QCircuit q3(1, 4, 3);
    q3.remove_clean_dits({0, 2});
    EXPECT_EQ(2u, q3.get_nc());
    // qubits
    QCircuit qc(2, 4);
    qc.remove_clean_dits({0, 2});
    EXPECT_EQ(2u, qc.get_nc());
}

/// BEGIN QCircuit& QCircuit::remove_clean_qudit(idx target)
TEST(qpp_QCircuit_remove_clean_qudit, AllTests) {
    // qutrits
    QCircuit q3(3, 0, 3);
    q3.remove_clean_qudit(1);
    EXPECT_EQ(2u, q3.get_nq());
    // qubits
    QCircuit qc(3);
    qc.remove_clean_qudit(1);
    EXPECT_EQ(2u, qc.get_nq());
}

/// BEGIN QCircuit& QCircuit::remove_clean_qudits(std::vector<idx> target)
TEST(qpp_QCircuit_remove_clean_qudits, AllTests) {
    // qutrits
    QCircuit q3(4, 0, 3);
    q3.remove_clean_qudits({0, 2});
    EXPECT_EQ(2u, q3.get_nq());
    // qubits
    QCircuit qc(4);
    qc.remove_clean_qudits({0, 2});
    EXPECT_EQ(2u, qc.get_nq());
}

/// BEGIN bool QCircuit::removes_qudits() const noexcept
TEST(qpp_QCircuit_removes_qudits, AllTests) {
    // qutrits
    QCircuit q3(1, 0, 3);
    EXPECT_FALSE(q3.removes_qudits());
    q3.discard(0);
    EXPECT_TRUE(q3.removes_qudits());
    // qubits
    QCircuit qc(2);
    EXPECT_FALSE(qc.removes_qudits());
    qc.discard(0);
    EXPECT_TRUE(qc.removes_qudits());
}
/// BEGIN QCircuit& QCircuit::replicate(idx n)
TEST(qpp_QCircuit_replicate, AllTests) {
    // qutrits
    QCircuit q3(1, 0, 3);
    q3.gate(gt.Xd(3), 0).replicate(2);
    QEngine e3{q3};
    e3.execute();
    EXPECT_NEAR(0, norm(e3.get_state() - mket({2}, 3)), 1e-12);
    // qubits
    QCircuit qc(1);
    qc.gate(gt.X, 0).replicate(3);
    EXPECT_EQ(3u, qc.get_gate_count());
    QEngine engine{qc};
    engine.execute();
    EXPECT_NEAR(0, norm(engine.get_state() - st.z1), 1e-12);
}

/// BEGIN QCircuit& QCircuit::reset(const std::vector<idx>& target,
///       std::optional<std::string> name = "reset")
TEST(qpp_QCircuit_reset, MultipleTargets) {
    // qutrits
    QCircuit q3(3, 0, 3);
    q3.gate_fan(gt.Xd(3)).reset({0, 2});
    QEngine e3{q3};
    e3.execute();
    EXPECT_NEAR(0, norm(e3.get_state() - mket({0, 1, 0}, 3)), 1e-12);
    // qubits
    QCircuit qc(3);
    qc.gate_fan(gt.X).reset({0, 2});
    EXPECT_EQ(2u, qc.get_step_count());
    QEngine engine{qc};
    engine.execute();
    EXPECT_NEAR(0, norm(engine.get_state() - mket({0, 1, 0})), 1e-12);
}

/// BEGIN QCircuit& QCircuit::reset(idx target,
///       std::optional<std::string> name = "reset")
TEST(qpp_QCircuit_reset, SingleTarget) {
    // qutrits
    QCircuit q3(2, 0, 3);
    q3.gate_fan(gt.Xd(3)).reset(1);
    QEngine e3{q3};
    e3.execute();
    EXPECT_NEAR(0, norm(e3.get_state() - mket({1, 0}, 3)), 1e-12);
    // qubits
    QCircuit qc(2);
    qc.gate_fan(gt.X).reset(1);
    EXPECT_EQ(2u, qc.get_step_count());
    QEngine engine{qc};
    engine.execute();
    EXPECT_NEAR(0, norm(engine.get_state() - mket({1, 0})), 1e-12);
}

/// BEGIN QCircuit& QCircuit::set_dits_runtime(mutable_dits_functor_t functor)
TEST(qpp_QCircuit_set_dits_runtime, AllTests) {
    // qutrits
    QCircuit q3(1, 1, 3);
    q3.set_dits_runtime([](auto dits) { dits[0] = 2; });
    QEngine e3{q3};
    e3.execute();
    EXPECT_EQ(2u, e3.get_dits()[0]);
    // qubits
    QCircuit qc(1, 1);
    qc.set_dits_runtime([](auto dits) { dits[0] = 1; });
    EXPECT_TRUE(qc.has_runtime_steps());
}

/// BEGIN QCircuit& QCircuit::set_name(const std::string& name)
TEST(qpp_QCircuit_set_name, AllTests) {
    // qutrits
    QCircuit q3(1, 0, 3);
    q3.set_name("qutrit");
    ASSERT_TRUE(q3.get_name());
    EXPECT_EQ("qutrit", *q3.get_name());
    // qubits
    QCircuit qc(1);
    qc.set_name("example");
    ASSERT_TRUE(qc.get_name());
    EXPECT_EQ("example", *qc.get_name());
}

/// BEGIN QCircuit& QCircuit::TFQ(bool swap = true)
TEST(qpp_QCircuit_TFQ, AllQudits) {
    // qutrits
    QCircuit q3(2, 0, 3);
    q3.gate(gt.Xd(3), 0).QFT().TFQ();
    QEngine e3{q3};
    e3.execute();
    EXPECT_NEAR(0, norm(e3.get_state() - mket({1, 0}, 3)), 1e-11);
    // qubits
    QCircuit qc(2);
    qc.gate(gt.X, 0).gate(gt.X, 1).QFT().TFQ();
    EXPECT_GT(qc.get_gate_count(), 2u);
    QEngine engine{qc};
    engine.execute();
    EXPECT_NEAR(0, norm(engine.get_state() - mket({1, 1})), 1e-12);
}

/// BEGIN QCircuit& QCircuit::TFQ(const std::vector<idx>& target,
///       bool swap = true)
TEST(qpp_QCircuit_TFQ, SpecificQudits) {
    // qutrits
    QCircuit q3(3, 0, 3);
    q3.gate(gt.Xd(3), 2).QFT({0, 2}).TFQ({0, 2});
    QEngine e3{q3};
    e3.execute();
    EXPECT_NEAR(0, norm(e3.get_state() - mket({0, 0, 1}, 3)), 1e-11);
    // qubits
    QCircuit qc(3);
    qc.gate(gt.X, 0).gate(gt.X, 2).QFT({0, 2}).TFQ({0, 2});
    EXPECT_GT(qc.get_gate_count(), 2u);
    EXPECT_TRUE(qc.is_clean_qudit(1));
    QEngine engine{qc};
    engine.execute();
    EXPECT_NEAR(0, norm(engine.get_state() - mket({1, 0, 1})), 1e-12);
}

/// BEGIN std::string QCircuit::to_JSON(
///       bool enclosed_in_curly_brackets = true)
///       const override
TEST(qpp_QCircuit_to_JSON, AllTests) {
    // qutrits
    QCircuit q3(1, 0, 3);
    q3.gate(gt.Xd(3), 0);
    EXPECT_NE(std::string::npos, q3.to_JSON().find("3"));
    // qubits
    QCircuit qc(1);
    qc.gate(gt.X, 0);
    const auto json = qc.to_JSON();
    EXPECT_FALSE(json.empty());
    EXPECT_EQ('{', json.front());
    EXPECT_EQ('}', json.back());
}

/// BEGIN bool QCircuit::validate_conditionals() const
TEST(qpp_QCircuit_validate_conditionals, AllTests) {
    // qutrits
    QCircuit q3(1, 1, 3);
    q3.cond_if([](auto dits) { return dits[0] == 2; }).cond_end();
    EXPECT_TRUE(q3.validate_conditionals());
    // qubits
    QCircuit qc(1);
    EXPECT_TRUE(qc.validate_conditionals());
    qc.cond_if([](auto) { return true; });
    EXPECT_FALSE(qc.validate_conditionals());
    qc.cond_end();
    EXPECT_TRUE(qc.validate_conditionals());
}

/// BEGIN bool QCircuit::was_measured_d(idx i) const
TEST(qpp_QCircuit_was_measured_d, AllTests) {
    // qutrits
    QCircuit q3(2, 1, 3);
    q3.measure(0, 0);
    EXPECT_TRUE(q3.was_measured_d(0));
    // qubits
    QCircuit qc(2, 1);
    EXPECT_FALSE(qc.was_measured_d(0));
    qc.measure(0, 0);
    EXPECT_TRUE(qc.was_measured_d(0));
}

/// BEGIN bool QCircuit::was_measured_nd(idx i) const
TEST(qpp_QCircuit_was_measured_nd, AllTests) {
    // qutrits
    QCircuit q3(2, 1, 3);
    q3.measure(0, 0, false);
    EXPECT_TRUE(q3.was_measured_nd(0));
    // qubits
    QCircuit qc(2, 1);
    EXPECT_FALSE(qc.was_measured_nd(1));
    qc.measure(1, 0, false);
    EXPECT_TRUE(qc.was_measured_nd(1));
}

// free functions

/// BEGIN inline QCircuit adjoint(QCircuit qc,
///       std::optional<std::string> name = std::nullopt)
TEST(qpp_adjoint, QCircuitAllTests) {
    // qutrits
    QCircuit q3(1, 0, 3);
    q3.gate(gt.Xd(3), 0);
    EXPECT_EQ(q3, qpp::adjoint(qpp::adjoint(q3)));
    // qubits
    QCircuit qc(1);
    qc.gate(gt.S, 0);
    EXPECT_EQ(qc, qpp::adjoint(qpp::adjoint(qc)));
}

/// BEGIN inline QCircuit compose_circuit(QCircuit qc1, const QCircuit& qc2,
///       bigint pos_qudit, std::optional<std::string> name = std::nullopt,
///       std::optional<idx> pos_dit = std::nullopt)
TEST(qpp_compose_circuit, AllTests) {
    // qutrits
    QCircuit a3(2, 0, 3), b3(1, 0, 3);
    b3.gate(gt.Xd(3), 0);
    auto c3 = qpp::compose_circuit(a3, b3, 1);
    QEngine e3{c3};
    e3.execute();
    EXPECT_NEAR(0, norm(e3.get_state() - mket({0, 1}, 3)), 1e-12);
    // qubits
    QCircuit a(2), b(1);
    b.gate(gt.X, 0);
    const auto result = qpp::compose_circuit(a, b, 1);
    EXPECT_EQ(2u, result.get_nq());
    EXPECT_EQ(1u, result.get_gate_count());
    QEngine engine{result};
    engine.execute();
    EXPECT_NEAR(0, norm(engine.get_state() - mket({0, 1})), 1e-12);
}

/// BEGIN QCircuit& compose_CTRL_circuit(QCircuit qc_ctrl,
///       const std::vector<idx>& ctrl, const QCircuit& qc_target,
///       bigint pos_qudit,
///       std::optional<std::vector<idx>> shift = std::nullopt,
///       bigint pos_qudit, std::optional<idx> pos_dit = std::nullopt,
///       std::optional<std::string> name = std::nullopt)
TEST(qpp_compose_CTRL_circuit, AllTests) {
    QCircuit qc_ctrl{5, 3, 2, "Control circuit"};
    qc_ctrl.gate_fan(gt.H);

    QCircuit qc_target{6, 4, 2, "Target circuit"};
    qc_target.gate_fan(gt.Id2, {1, 2});
    qc_target.gate(gt.X, 0);
    qc_target.gate(gt.Y, 1);
    qc_target.gate(gt.Z, 2);
    qc_target.CTRL(gt.X, 1, 2);
    qc_target.CTRL(gt.CNOT, {0, 2}, {1, 3}, std::vector<idx>{1, 0});

    std::vector<idx> ctrl = {1, 2}, shift = {0, 1};

    // compose at the end
    QCircuit result1 =
        QCircuit{qc_ctrl}.compose_CTRL_circuit(ctrl, qc_target, 5, shift);

    // compose at the end, free function
    QCircuit result1a =
        compose_CTRL_circuit(qc_ctrl, ctrl, qc_target, 5, shift);

    QCircuit expected_result1{11, 7};
    expected_result1.gate_fan(gt.H, {0, 1, 2, 3, 4});
    expected_result1.CTRL_fan(gt.Id2, {1, 2}, {6, 7}, std::vector<idx>{0, 1});
    expected_result1.CTRL(gt.X, {1, 2}, 5, std::vector<idx>{0, 1});
    expected_result1.CTRL(gt.Y, {1, 2}, 6, std::vector<idx>{0, 1});
    expected_result1.CTRL(gt.Z, {1, 2}, 7, std::vector<idx>{0, 1});
    expected_result1.CTRL(gt.X, {1, 2, 6}, 7, std::vector<idx>{0, 1, 0});
    expected_result1.CTRL(gt.CNOT, {1, 2, 5, 7}, {6, 8},
                          std::vector<idx>{0, 1, 1, 0});

    EXPECT_EQ(result1, result1a);
    EXPECT_EQ(expected_result1, result1);
    EXPECT_EQ(expected_result1.get_nq(), 11);
    EXPECT_EQ(expected_result1.get_nc(), 7);

    // compose at the middle
    ctrl = {1, 2}, shift = {1, 1};
    QCircuit result2 =
        QCircuit{qc_ctrl}.compose_CTRL_circuit(ctrl, qc_target, 3, shift);

    // compose at the middle, free function
    QCircuit result2a =
        compose_CTRL_circuit(qc_ctrl, ctrl, qc_target, 3, shift);

    QCircuit expected_result2{9, 7};
    expected_result2.gate_fan(gt.H, {0, 1, 2, 3, 4});
    expected_result2.CTRL_fan(gt.Id2, {1, 2}, {4, 5}, std::vector<idx>{1, 1});
    expected_result2.CTRL(gt.X, {1, 2}, 3, std::vector<idx>{1, 1});
    expected_result2.CTRL(gt.Y, {1, 2}, 4, std::vector<idx>{1, 1});
    expected_result2.CTRL(gt.Z, {1, 2}, 5, std::vector<idx>{1, 1});
    expected_result2.CTRL(gt.X, {1, 2, 4}, 5, std::vector<idx>{1, 1, 0});
    expected_result2.CTRL(gt.CNOT, {1, 2, 3, 5}, {4, 6},
                          std::vector<idx>{1, 1, 1, 0});

    EXPECT_EQ(result2, result2a);
    EXPECT_EQ(expected_result2, result2);
    EXPECT_EQ(expected_result2.get_nq(), 9);
    EXPECT_EQ(expected_result2.get_nc(), 7);
}

/// BEGIN inline QCircuit couple_circuit_left(QCircuit qc1, const QCircuit& qc2,
///       const std::vector<idx>& target,
///       std::optional<std::string> name = std::nullopt,
///       std::optional<idx> pos_dit = std::nullopt)
TEST(qpp_couple_circuit_left, AllTests) {
    // qutrits
    QCircuit a3(1, 0, 3), b3(1, 0, 3);
    b3.gate(gt.Xd(3), 0);
    auto c3 = qpp::couple_circuit_left(a3, b3, {0});
    QEngine e3{c3};
    e3.execute();
    EXPECT_NEAR(0, norm(e3.get_state() - mket({1}, 3)), 1e-12);
    // qubits
    QCircuit a(1), b(1);
    b.gate(gt.X, 0);
    const auto result = qpp::couple_circuit_left(a, b, {0});
    EXPECT_EQ(1u, result.get_nq());
    EXPECT_EQ(1u, result.get_gate_count());
    QEngine engine{result};
    engine.execute();
    EXPECT_NEAR(0, norm(engine.get_state() - st.z1), 1e-12);
}

/// BEGIN inline QCircuit couple_circuit_right(QCircuit qc1,
///       const QCircuit& qc2, const std::vector<idx>& target,
///       std::optional<std::string> name = std::nullopt,
///       std::optional<idx> pos_dit = std::nullopt)
TEST(qpp_couple_circuit_right, AllTests) {
    // qutrits
    QCircuit a3(1, 0, 3), b3(1, 0, 3);
    b3.gate(gt.Xd(3), 0);
    auto c3 = qpp::couple_circuit_right(a3, b3, {0});
    QEngine e3{c3};
    e3.execute();
    EXPECT_NEAR(0, norm(e3.get_state() - mket({1}, 3)), 1e-12);
    // qubits
    QCircuit a(1), b(1);
    b.gate(gt.X, 0);
    const auto result = qpp::couple_circuit_right(a, b, {0});
    EXPECT_EQ(1u, result.get_nq());
    EXPECT_EQ(1u, result.get_gate_count());
    QEngine engine{result};
    engine.execute();
    EXPECT_NEAR(0, norm(engine.get_state() - st.z1), 1e-12);
}

/// BEGIN inline QCircuit kron(QCircuit qc1, const QCircuit& qc2,
///       std::optional<std::string> name = std::nullopt)
TEST(qpp_kron, QCircuitAllTests) {
    // qutrits
    QCircuit a3(1, 0, 3), b3(1, 0, 3);
    a3.gate(gt.Xd(3), 0);
    b3.gate(gt.Xd(3), 0);
    auto c3 = qpp::kron(a3, b3);
    QEngine e3{c3};
    e3.execute();
    EXPECT_NEAR(0, norm(e3.get_state() - mket({1, 1}, 3)), 1e-12);
    // qubits
    QCircuit a(1), b(2);
    a.gate(gt.X, 0);
    b.gate(gt.X, 1);
    const auto result = qpp::kron(a, b);
    EXPECT_EQ(3u, result.get_nq());
    EXPECT_EQ(2u, result.get_gate_count());
    QEngine engine{result};
    engine.execute();
    EXPECT_NEAR(0, norm(engine.get_state() - mket({1, 0, 1})), 1e-12);
}

/// BEGIN inline QCircuit qpe_circuit(cmat U, qpp::idx n,
///       bool omit_measurements = true, idx d = 2,
///       std::optional<std::string> name = "qpe")
TEST(qpp_qpe_circuit, AllTests) {
    // qutrits
    auto q3 = qpe_circuit(gt.Zd(3), 2, true, 3);
    QEngine e3{q3};
    e3.execute();
    EXPECT_NEAR(0, norm(e3.get_state() - mket({0, 0, 0}, 3)), 1e-10);
    // qubits
    const auto qc = qpe_circuit(gt.Z, 2);
    EXPECT_EQ(3u, qc.get_nq());
    EXPECT_GT(qc.get_gate_count(), 0u);
    QEngine engine{qc};
    engine.execute();
    EXPECT_NEAR(0, norm(engine.get_state() - mket({0, 0, 0})), 1e-12);
}

/// BEGIN inline QCircuit state_prep_circuit(const ket& psi, idx d = 2,
///       bool decompose = true,
///       std::optional<std::string> name = "state_prep")
TEST(qpp_state_prep_circuit, Qubits) {
    for (idx n = 1; n <= 8; ++n) {
        idx D = idx{1} << n;
        ket psi = randket(D);
        QCircuit qc = state_prep_circuit(psi);
        EXPECT_EQ(n, qc.get_nq());
        EXPECT_EQ(0, qc.get_nc());
        ket result = QEngine{qc}.execute().get_state();
        EXPECT_NEAR(0, norm(result - psi), 1e-5);
    }
}

TEST(qpp_state_prep_circuit, KnownStates) {
    std::vector<ket> states{st.z0,      st.z1,      st.x0,           st.x1,
                            st.y0,      st.y1,      st.b00,          st.b01,
                            st.b10,     st.b11,     st.GHZ,          st.W,
                            st.plus(4), st.zero(3), mket({1, 0, 1}), st.one(5)};
    for (auto&& psi : states) {
        QCircuit qc = state_prep_circuit(psi);
        ket result = QEngine{qc}.execute().get_state();
        EXPECT_NEAR(0, norm(result - psi), 1e-5);
    }

    // global phases are reproduced exactly
    ket psi = std::exp(1_i * 0.7) * st.GHZ;
    QCircuit qc = state_prep_circuit(psi);
    ket result = QEngine{qc}.execute().get_state();
    EXPECT_NEAR(0, norm(result - psi), 1e-5);
}

TEST(qpp_state_prep_circuit, RealAmplitudes) {
    for (idx n = 1; n <= 7; ++n) {
        idx D = idx{1} << n;
        ket psi = rand<rmat>(D, 1, -1, 1).cast<cplx>();
        psi.normalize();
        QCircuit qc = state_prep_circuit(psi);
        ket result = QEngine{qc}.execute().get_state();
        EXPECT_NEAR(0, norm(result - psi), 1e-5);
        EXPECT_LE(qc.get_gate_count(gt.CNOT), D - 2);
        // only RY rotations and CNOTs
        for (auto&& elem : qc) {
            auto step = elem.get_step();
            ASSERT_TRUE(
                std::holds_alternative<internal::QCircuitGateStep>(step));
            auto name =
                std::get<internal::QCircuitGateStep>(step).name_.value();
            EXPECT_TRUE(name == "RY" || name == "CNOT");
        }
    }
}

TEST(qpp_state_prep_circuit, GateCountBound) {
    for (idx n = 1; n <= 7; ++n) {
        idx D = idx{1} << n;
        ket psi = randket(D);
        QCircuit qc = state_prep_circuit(psi);
        EXPECT_LE(qc.get_gate_count(gt.CNOT), 2 * D - 4);
    }
}

TEST(qpp_state_prep_circuit, BasisStates) {
    for (idx n = 1; n <= 5; ++n) {
        idx D = idx{1} << n;
        for (idx x = 0; x < D; ++x) {
            ket psi = ket::Zero(D);
            psi[x] = 1;
            QCircuit qc = state_prep_circuit(psi);
            EXPECT_EQ(0, qc.get_gate_count(gt.CNOT));
            ket result = QEngine{qc}.execute().get_state();
            EXPECT_NEAR(0, norm(result - psi), 1e-5);

            // with a phase
            psi[x] = 1_i;
            qc = state_prep_circuit(psi);
            EXPECT_EQ(0, qc.get_gate_count(gt.CNOT));
            result = QEngine{qc}.execute().get_state();
            EXPECT_NEAR(0, norm(result - psi), 1e-5);
        }
    }

    // product states require no CNOTs
    ket psi = kron(randket(), randket(), randket(), randket());
    QCircuit qc = state_prep_circuit(psi);
    EXPECT_EQ(0, qc.get_gate_count(gt.CNOT));
    ket result = QEngine{qc}.execute().get_state();
    EXPECT_NEAR(0, norm(result - psi), 1e-5);
}

TEST(qpp_state_prep_circuit, DecomposedMatchesJoint) {
    for (idx n = 1; n <= 6; ++n) {
        idx D = idx{1} << n;
        ket psi = randket(D);
        QCircuit qc_decomposed = state_prep_circuit(psi);
        QCircuit qc_joint = state_prep_circuit(psi, 2, false);
        ket result_decomposed = QEngine{qc_decomposed}.execute().get_state();
        ket result_joint = QEngine{qc_joint}.execute().get_state();
        EXPECT_NEAR(0, norm(result_decomposed - result_joint), 1e-5);
        EXPECT_NEAR(0, norm(result_joint - psi), 1e-5);
    }
}

TEST(qpp_state_prep_circuit, Adjoint) {
    for (idx n = 1; n <= 6; ++n) {
        idx D = idx{1} << n;
        ket psi = randket(D);
        QCircuit qc = adjoint(state_prep_circuit(psi));
        ket result = QEngine{qc}.reset(psi).execute().get_state();
        EXPECT_NEAR(0, norm(result - st.zero(n)), 1e-5);
    }
}

TEST(qpp_state_prep_circuit, DensityEngine) {
    for (idx n = 1; n <= 5; ++n) {
        idx D = idx{1} << n;
        ket psi = randket(D);
        QCircuit qc = state_prep_circuit(psi);
        cmat result = QDensityEngine{qc}.execute().get_state();
        EXPECT_NEAR(0, norm(result - prj(psi)), 1e-5);
    }
}

TEST(qpp_state_prep_circuit, Qudits) {
    for (idx d : std::vector<idx>{3, 5}) {
        for (idx n = 1; n <= 3; ++n) {
            idx D = internal::ipow_rounded(d, n);
            ket psi = randket(D);
            QCircuit qc = state_prep_circuit(psi, d);
            EXPECT_EQ(n, qc.get_nq());
            EXPECT_EQ(d, qc.get_d());
            ket result = QEngine{qc}.execute().get_state();
            EXPECT_NEAR(0, norm(result - psi), 1e-5);

            // computational basis state
            ket basis = st.jn(d - 1, n, d);
            qc = state_prep_circuit(basis, d);
            result = QEngine{qc}.execute().get_state();
            EXPECT_NEAR(0, norm(result - basis), 1e-5);
        }
    }
}

TEST(qpp_state_prep_circuit, NearBasisStates) {
    // small perturbations of computational basis states are not lost
    realT tol = 100 * std::numeric_limits<realT>::epsilon();
    for (idx d : std::vector<idx>{2, 3}) {
        idx D = internal::ipow_rounded(d, idx{3});
        ket psi = st.zero(3, d) + 1e-9 * randket(D);
        psi.normalize();
        QCircuit qc_decomposed = state_prep_circuit(psi, d);
        QCircuit qc_joint = state_prep_circuit(psi, d, false);
        ket result_decomposed = QEngine{qc_decomposed}.execute().get_state();
        ket result_joint = QEngine{qc_joint}.execute().get_state();
        EXPECT_NEAR(0, norm(result_decomposed - psi), tol);
        EXPECT_NEAR(0, norm(result_joint - psi), tol);
    }
}

TEST(qpp_state_prep_circuit, Exceptions) {
    EXPECT_THROW(state_prep_circuit(ket{}), exception::ZeroSize);
    EXPECT_THROW(state_prep_circuit(st.z0, 1), exception::DimsInvalid);
    EXPECT_THROW(state_prep_circuit(randket(6)),
                 exception::DimsMismatchCvector);
    EXPECT_THROW(state_prep_circuit(randket(1)),
                 exception::DimsMismatchCvector);
    EXPECT_THROW(state_prep_circuit(randket(4), 3),
                 exception::DimsMismatchCvector);
    ket psi = 2 * randket(4);
    EXPECT_THROW(state_prep_circuit(psi), exception::CustomException);
}

/// BEGIN inline QCircuit replicate(QCircuit qc, idx n,
///       std::optional<std::string> name = std::nullopt)
TEST(qpp_replicate, AllTests) {
    // qutrits
    QCircuit q3(1, 0, 3);
    q3.gate(gt.Xd(3), 0);
    auto c3 = qpp::replicate(q3, 2);
    QEngine e3{c3};
    e3.execute();
    EXPECT_NEAR(0, norm(e3.get_state() - mket({2}, 3)), 1e-12);
    // qubits
    QCircuit qc(1);
    qc.gate(gt.X, 0);
    const auto repeated = qpp::replicate(qc, 4);
    EXPECT_EQ(4u, repeated.get_gate_count());
    QEngine engine{repeated};
    engine.execute();
    EXPECT_NEAR(0, norm(engine.get_state() - st.z0), 1e-12);
}

/// BEGIN inline QCircuit random_circuit_count(idx nq, idx d, idx gate_count,
///       std::optional<realT> p_two,
///       std::optional<cmat> with_respect_to_gate = std::nullopt,
///       std::optional<std::vector<cmat>> one_qudit_gate_set = std::nullopt,
///       std::optional<std::vector<cmat>> two_qudit_gate_set = std::nullopt,
///       std::optional<std::vector<std::string>> one_qudit_gate_names =
///           std::nullopt,
///       std::optional<std::vector<std::string>> two_qudit_gate_names =
///           std::nullopt)
TEST(qpp_random_circuit_count, AllTests) {
    // qutrits
    auto q3 = random_circuit_count(3, 3, 8, 0.5, std::nullopt,
                                   std::vector<cmat>{gt.Xd(3), gt.Zd(3)},
                                   std::vector<cmat>{kron(gt.Xd(3), gt.Xd(3))});
    EXPECT_EQ(3u, q3.get_d());
    EXPECT_EQ(8u, q3.get_gate_count());
    // qubits
    const auto qc = random_circuit_count(3, 2, 12, 0.5);
    EXPECT_EQ(3u, qc.get_nq());
    EXPECT_EQ(12u, qc.get_gate_count());
}

/// BEGIN inline QCircuit random_circuit_depth(idx nq, idx d, idx gate_depth,
///       std::optional<realT> p_two,
///       std::optional<cmat> with_respect_to_gate = std::nullopt,
///       std::optional<std::vector<cmat>> one_qudit_gate_set = std::nullopt,
///       std::optional<std::vector<cmat>> two_qudit_gate_set = std::nullopt,
///       std::optional<std::vector<std::string>> one_qudit_gate_names =
///           std::nullopt,
///       std::optional<std::vector<std::string>> two_qudit_gate_names =
///           std::nullopt)
TEST(qpp_random_circuit_depth, AllTests) {
    // qutrits
    auto q3 = random_circuit_depth(3, 3, 3, 0.5, std::nullopt,
                                   std::vector<cmat>{gt.Xd(3), gt.Zd(3)},
                                   std::vector<cmat>{kron(gt.Xd(3), gt.Xd(3))});
    EXPECT_EQ(3u, q3.get_d());
    EXPECT_EQ(3u, q3.get_gate_depth());
    // qubits
    const auto qc = random_circuit_depth(3, 2, 4, 0.5);
    EXPECT_EQ(3u, qc.get_nq());
    EXPECT_EQ(4u, qc.get_gate_depth());
}
