#include "gtest/gtest.h"

#include "qpp/qpp.hpp"

using namespace qpp;

// Unit testing "qpp/classes/qengine.hpp"

/// BEGIN QBaseEngine& QEngineT::execute(
///       typename QCircuitTraits<QCT>::iterator_type& it)
TEST(qpp_QEngineT_execute, Iterator) {
    // qubits
    QCircuit qc(1);
    qc.gate(gt.X, 0).gate(gt.X, 0);
    QEngine engine{qc};
    auto it = qc.begin();
    engine.execute(it);
    EXPECT_NEAR(0, norm(engine.get_state() - st.z1), 1e-12);
    engine.execute(++it);
    EXPECT_NEAR(0, norm(engine.get_state() - st.z0), 1e-12);

    // qutrits
    QCircuit qc3(1, 0, 3);
    qc3.gate(gt.Xd(3), 0).gate(gt.Xd(3), 0);
    QEngine engine3{qc3};
    auto it3 = qc3.begin();
    engine3.execute(it3);
    EXPECT_NEAR(0, norm(engine3.get_state() - mket({1}, 3)), 1e-12);
    engine3.execute(++it3);
    EXPECT_NEAR(0, norm(engine3.get_state() - mket({2}, 3)), 1e-12);
}

/// BEGIN QEngineT& QEngineT::execute(
///       const typename QCircuitTraits<QCircuit>::value_type& elem) override
TEST(qpp_QEngineT_execute, ValueType) {
    // qubits
    QCircuit qc(1);
    qc.gate(gt.X, 0);
    QEngine engine{qc};
    const auto elem = *qc.begin();
    engine.execute(elem);
    EXPECT_NEAR(0, norm(engine.get_state() - st.z1), 1e-12);

    // qutrits
    QCircuit qc3(1, 0, 3);
    qc3.gate(gt.Xd(3), 0);
    QEngine engine3{qc3};
    const auto elem3 = *qc3.begin();
    engine3.execute(elem3);
    EXPECT_NEAR(0, norm(engine3.get_state() - mket({1}, 3)), 1e-12);
}

/// BEGIN QEngineT& QEngineT::execute(idx reps = 1) override
TEST(qpp_QEngineT_execute, AllCircuitWithRepetitions) {
    // qubits
    QCircuit qc(2, 2);
    qc.gate(gt.X, 0).measure_all(0, false);
    QEngine engine{qc};
    engine.execute(5);
    EXPECT_EQ((std::vector<idx>{1, 0}), engine.get_dits());
    EXPECT_NEAR(0, norm(engine.get_state() - mket({1, 0})), 1e-12);
    EXPECT_EQ(5u, engine.get_stats().get_num_reps());
    EXPECT_EQ(1u, engine.get_stats().get_num_outcomes());
    EXPECT_THROW(engine.execute(0), exception::OutOfRange);

    // qutrits
    QCircuit qc3(1, 1, 3);
    qc3.gate(gt.Xd(3), 0).gate(gt.Xd(3), 0).measure(0, 0, false);
    QEngine engine3{qc3};
    engine3.execute(7);
    EXPECT_EQ((std::vector<idx>{2}), engine3.get_dits());
    EXPECT_NEAR(0, norm(engine3.get_state() - mket({2}, 3)), 1e-12);
    EXPECT_EQ(7u, engine3.get_stats().get_num_reps());
}

/// BEGIN const QCT& QBaseEngine::get_circuit() const& noexcept
TEST(qpp_QEngineT_get_circuit, Lvalue) {
    // qubits
    QCircuit qc(1);
    qc.gate(gt.X, 0);
    const QEngine engine{qc};
    EXPECT_EQ(qc, engine.get_circuit());

    // qutrits
    QCircuit qc3(1, 0, 3);
    qc3.gate(gt.Xd(3), 0);
    const QEngine engine3{qc3};
    EXPECT_EQ(qc3, engine3.get_circuit());
}

/// BEGIN QCT QBaseEngine::get_circuit() const&& noexcept
TEST(qpp_QEngineT_get_circuit, Rvalue) {
    // qubits
    QCircuit qc(1);
    qc.gate(gt.X, 0);
    const QEngine engine{qc};
    EXPECT_EQ(qc, std::move(engine).get_circuit());

    // qutrits
    QCircuit qc3(1, 0, 3);
    qc3.gate(gt.Xd(3), 0);
    const QEngine engine3{qc3};
    EXPECT_EQ(qc3, std::move(engine3).get_circuit());
}

/// BEGIN idx QEngineT::get_dit(idx i) const
TEST(qpp_QEngineT_get_dit, AllTests) {
    // qubits
    QCircuit qc(1, 2);
    QEngine engine{qc};
    engine.set_dits({1, 0});
    EXPECT_EQ(1u, engine.get_dit(0));
    EXPECT_EQ(0u, engine.get_dit(1));
    EXPECT_THROW(engine.get_dit(2), exception::OutOfRange);

    // qutrits
    QCircuit qc3(1, 2, 3);
    QEngine engine3{qc3};
    engine3.set_dits({2, 1});
    EXPECT_EQ(2u, engine3.get_dit(0));
    EXPECT_EQ(1u, engine3.get_dit(1));
}

/// BEGIN std::vector<idx> QEngineT::get_dits() const
TEST(qpp_QEngineT_get_dits, AllTests) {
    // qubits
    QCircuit qc(1, 2);
    QEngine engine{qc};
    EXPECT_EQ((std::vector<idx>{0, 0}), engine.get_dits());
    engine.set_dits({1, 0});
    EXPECT_EQ((std::vector<idx>{1, 0}), engine.get_dits());

    // qutrits
    QCircuit qc3(1, 2, 3);
    QEngine engine3{qc3};
    engine3.set_dits({2, 1});
    EXPECT_EQ((std::vector<idx>{2, 1}), engine3.get_dits());
}

/// BEGIN bool QEngineT::get_ensure_post_selection() const
TEST(qpp_QEngineT_get_ensure_post_selection, AllTests) {
    QCircuit qc(1);
    QEngine engine{qc};
    EXPECT_FALSE(engine.get_ensure_post_selection());
    engine.set_ensure_post_selection(true);
    EXPECT_TRUE(engine.get_ensure_post_selection());
    engine.set_ensure_post_selection(false);
    EXPECT_FALSE(engine.get_ensure_post_selection());
}

/// BEGIN idx QEngineT::get_max_post_selection_reps() const
TEST(qpp_QEngineT_get_max_post_selection_reps, AllTests) {
    QCircuit qc(1);
    QEngine engine{qc};
    engine.set_max_post_selection_reps(17);
    EXPECT_EQ(17u, engine.get_max_post_selection_reps());
    engine.set_max_post_selection_reps(4);
    EXPECT_EQ(4u, engine.get_max_post_selection_reps());
}

/// BEGIN std::vector<idx> QEngineT::get_measured_d() const
TEST(qpp_QEngineT_get_measured_d, AllTests) {
    // qubits
    QCircuit qc(3, 1);
    qc.measure(1, 0);
    QEngine engine{qc};
    engine.execute();
    EXPECT_EQ((std::vector<idx>{1}), engine.get_measured_d());

    // qutrits
    QCircuit qc3(3, 1, 3);
    qc3.measure(2, 0);
    QEngine engine3{qc3};
    engine3.execute();
    EXPECT_EQ((std::vector<idx>{2}), engine3.get_measured_d());
}

/// BEGIN std::vector<idx> QEngineT::get_non_measured_d() const
TEST(qpp_QEngineT_get_non_measured_d, AllTests) {
    // qubits
    QCircuit qc(3, 1);
    qc.measure(1, 0);
    QEngine engine{qc};
    engine.execute();
    EXPECT_EQ((std::vector<idx>{0, 2}), engine.get_non_measured_d());

    // qutrits
    QCircuit qc3(3, 1, 3);
    qc3.measure(0, 0);
    QEngine engine3{qc3};
    engine3.execute();
    EXPECT_EQ((std::vector<idx>{1, 2}), engine3.get_non_measured_d());
}

/// BEGIN std::vector<realT> QEngineT::get_probs() const
TEST(qpp_QEngineT_get_probs, AllTests) {
    // qubits
    QCircuit qc(1, 1);
    qc.gate(gt.X, 0).measure(0, 0, false);
    QEngine engine{qc};
    engine.execute();
    ASSERT_EQ(1u, engine.get_probs().size());
    EXPECT_NEAR(1, engine.get_probs()[0], 1e-12);

    // qutrits
    QCircuit qc3(1, 1, 3);
    qc3.gate(gt.Xd(3), 0).measure(0, 0, false);
    QEngine engine3{qc3};
    engine3.execute();
    ASSERT_EQ(1u, engine3.get_probs().size());
    EXPECT_NEAR(1, engine3.get_probs()[0], 1e-12);
}

/// BEGIN T QEngineT::get_state() const override
TEST(qpp_QEngineT_get_state, AllTests) {
    // qubits
    QCircuit qc(2);
    qc.gate(gt.X, 1);
    QEngine engine{qc};
    EXPECT_NEAR(0, norm(engine.get_state() - mket({0, 0})), 1e-12);
    engine.execute();
    EXPECT_NEAR(0, norm(engine.get_state() - mket({0, 1})), 1e-12);
    QDensityEngine mixed{qc};
    mixed.execute();
    EXPECT_NEAR(0, norm(mixed.get_state() - mprj({0, 1})), 1e-12);

    // qutrits
    QCircuit qc3(2, 0, 3);
    qc3.gate(gt.Xd(3), 1).gate(gt.Xd(3), 1);
    QEngine engine3{qc3};
    engine3.execute();
    EXPECT_NEAR(0, norm(engine3.get_state() - mket({0, 2}, 3)), 1e-12);
}

/// BEGIN internal::QEngineStatistics QEngineT::get_stats() const
TEST(qpp_QEngineT_get_stats, AllTests) {
    QCircuit qc(1, 1);
    qc.gate(gt.X, 0).measure(0, 0);
    QEngine engine{qc};
    engine.execute(10);
    const auto stats = engine.get_stats();
    EXPECT_EQ(10u, stats.get_num_reps());
    EXPECT_EQ(1u, stats.get_num_outcomes());
    EXPECT_EQ(10u, stats.data().at(std::vector<idx>{1}));
}

/// BEGIN bool QEngineT::post_select_ok() const
TEST(qpp_QEngineT_post_select_ok, AllTests) {
    // qubits
    QCircuit qc(1, 1);
    qc.gate(gt.X, 0).post_select(0, 1, 0, false);
    QEngine engine{qc};
    engine.execute();
    EXPECT_TRUE(engine.post_select_ok());

    // qutrits
    QCircuit qc3(1, 1, 3);
    qc3.gate(gt.Xd(3), 0).gate(gt.Xd(3), 0);
    qc3.post_select(0, 2, 0, false);
    QEngine engine3{qc3};
    engine3.execute();
    EXPECT_TRUE(engine3.post_select_ok());
}

/// BEGIN QEngineT& QEngineT::reset(std::optional<T> qstate = std::nullopt,
///       bool reset_stats = true)
TEST(qpp_QEngineT_reset, AllTests) {
    auto circuit = qpp::QCircuit{2, 1};
    circuit.measure(0, 0);

    auto engine = qpp::QEngine{circuit};
    auto const psi = qpp::randket(4).eval();
    engine.reset().set_state(psi).execute();

    EXPECT_NO_THROW(engine.reset().set_state(psi).execute(2));
}

/// BEGIN QEngineT& QEngineT::reset_stats()
TEST(qpp_QEngineT_reset_stats, AllTests) {
    QCircuit qc(1, 1);
    qc.measure(0, 0);
    QEngine engine{qc};
    engine.execute(5);
    EXPECT_EQ(5u, engine.get_stats().get_num_reps());
    engine.reset_stats();
    EXPECT_EQ(0u, engine.get_stats().get_num_reps());
    EXPECT_EQ(0u, engine.get_stats().get_num_outcomes());
}

/// BEGIN QEngineT& QEngineT::set_dit(idx i, idx value)
TEST(qpp_QEngineT_set_dit, AllTests) {
    // qubits
    QCircuit qc(1, 2);
    QEngine engine{qc};
    engine.set_dit(1, 1);
    EXPECT_EQ((std::vector<idx>{0, 1}), engine.get_dits());
    EXPECT_THROW(engine.set_dit(2, 1), exception::OutOfRange);

    // qutrits
    QCircuit qc3(1, 2, 3);
    QEngine engine3{qc3};
    engine3.set_dit(0, 2);
    EXPECT_EQ((std::vector<idx>{2, 0}), engine3.get_dits());
}

/// BEGIN QEngineT& QEngineT::set_dits(std::vector<idx> dits)
TEST(qpp_QEngineT_set_dits, AllTests) {
    // qubits
    QCircuit qc(1, 2);
    QEngine engine{qc};
    engine.set_dits({1, 0});
    EXPECT_EQ((std::vector<idx>{1, 0}), engine.get_dits());
    EXPECT_THROW(engine.set_dits({1}), exception::SizeMismatch);

    // qutrits
    QCircuit qc3(1, 2, 3);
    QEngine engine3{qc3};
    engine3.set_dits({2, 1});
    EXPECT_EQ((std::vector<idx>{2, 1}), engine3.get_dits());
}

/// BEGIN QEngineT& QEngineT::set_ensure_post_selection(bool val)
TEST(qpp_QEngineT_set_ensure_post_selection, AllTests) {
    QCircuit qc(1);
    QEngine engine{qc};
    engine.set_ensure_post_selection(true);
    EXPECT_TRUE(engine.get_ensure_post_selection());
    engine.set_ensure_post_selection(false);
    EXPECT_FALSE(engine.get_ensure_post_selection());
}

/// BEGIN QEngineT& QEngineT::set_max_post_selection_reps(
///       idx max_post_selection_reps)
TEST(qpp_QEngineT_set_max_post_selection_reps, AllTests) {
    QCircuit qc(1);
    QEngine engine{qc};
    engine.set_max_post_selection_reps(23);
    EXPECT_EQ(23u, engine.get_max_post_selection_reps());
}

/// BEGIN QEngineT& QEngineT::set_state(const T& state) override
TEST(qpp_QEngineT_set_state, AllTests) {
    // qubits
    QCircuit qc(2);
    QEngine engine{qc};
    engine.set_state(mket({1, 0}));
    EXPECT_NEAR(0, norm(engine.get_state() - mket({1, 0})), 1e-12);
    EXPECT_THROW(engine.set_state(st.z0), exception::DimsNotEqual);
    QCircuit mixed_qc(1);
    QDensityEngine mixed{mixed_qc};
    mixed.set_state(prj(st.x0));
    EXPECT_NEAR(0, norm(mixed.get_state() - prj(st.x0)), 1e-12);

    // qutrits
    QCircuit qc3(2, 0, 3);
    QEngine engine3{qc3};
    engine3.set_state(mket({2, 1}, 3));
    EXPECT_NEAR(0, norm(engine3.get_state() - mket({2, 1}, 3)), 1e-12);
}

/// BEGIN std::string QEngineT::to_JSON(
///       bool enclosed_in_curly_brackets = true) const override
TEST(qpp_QEngineT_to_JSON, AllTests) {
    QCircuit qc(1);
    QEngine engine{qc};
    const auto json = engine.to_JSON();
    ASSERT_FALSE(json.empty());
    EXPECT_EQ('{', json.front());
    EXPECT_EQ('}', json.back());
    EXPECT_EQ(json.substr(1, json.size() - 2), engine.to_JSON(false));
}

/// BEGIN  std::string QEngineT::traits_get_name() const override
TEST(qpp_QEngineT_traits_get_name, AllTests) {
    QCircuit qc(1);
    EXPECT_EQ("QEngine", QEngine{qc}.traits_get_name());
    EXPECT_EQ("QDensityEngine", QDensityEngine{qc}.traits_get_name());
}

/// BEGIN bool QEngineT::traits_is_noisy() const override
TEST(qpp_QEngineT_traits_is_noisy, AllTests) {
    QCircuit qc(1);
    EXPECT_FALSE(QEngine{qc}.traits_is_noisy());
    EXPECT_FALSE(QDensityEngine{qc}.traits_is_noisy());
}

/// BEGIN bool QEngineT::traits_is_pure() const override
TEST(qpp_QEngineT_traits_is_pure, AllTests) {
    QCircuit qc(1);
    EXPECT_TRUE(QEngine{qc}.traits_is_pure());
    EXPECT_FALSE(QDensityEngine{qc}.traits_is_pure());
}

/// BEGIN bool QEngineT::was_measured_d(idx i) const
TEST(qpp_QEngineT_was_measured_d, AllTests) {
    // qubits
    QCircuit qc(2, 1);
    qc.measure(0, 0);
    QEngine engine{qc};
    EXPECT_FALSE(engine.was_measured_d(0));
    engine.execute();
    EXPECT_TRUE(engine.was_measured_d(0));
    EXPECT_FALSE(engine.was_measured_d(1));

    // qutrits
    QCircuit qc3(2, 1, 3);
    qc3.measure(1, 0);
    QEngine engine3{qc3};
    engine3.execute();
    EXPECT_TRUE(engine3.was_measured_d(1));
    EXPECT_FALSE(engine3.was_measured_d(0));
}
