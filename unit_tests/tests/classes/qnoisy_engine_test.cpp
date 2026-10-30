#include "gtest/gtest.h"

#include "qpp/qpp.hpp"

using namespace qpp;

// Unit testing "qpp/classes/qnoisy_engine.hpp"

/// BEGIN QBaseEngine& QNoisyEngineT::execute(
///       typename QCircuitTraits<QCT>::iterator_type& it)
TEST(qpp_QNoisyEngineT_execute, Iterator) {
    // qubits
    QCircuit qc(1);
    qc.gate(gt.X, 0).gate(gt.X, 0);
    QubitBitFlipNoise noise(0);
    QNoisyEngine engine{qc, noise};
    auto it = qc.begin();
    engine.execute(it);
    EXPECT_NEAR(0, norm(engine.get_state() - st.z1), 1e-12);
    engine.execute(++it);
    EXPECT_NEAR(0, norm(engine.get_state() - st.z0), 1e-12);

    // deterministic noise: X before an identity gate
    QCircuit noisy_qc(1);
    noisy_qc.gate(gt.Id2, 0);
    QubitBitFlipNoise always_flip(1);
    QNoisyEngine noisy_engine{noisy_qc, always_flip};
    auto noisy_it = noisy_qc.begin();
    noisy_engine.execute(noisy_it);
    EXPECT_NEAR(0, norm(noisy_engine.get_state() - st.z1), 1e-12);

    // qutrits
    QCircuit qc3(1, 0, 3);
    qc3.gate(gt.Xd(3), 0).gate(gt.Xd(3), 0);
    QuditDepolarizingNoise noise3(0, 3);
    QNoisyEngine engine3{qc3, noise3};
    auto it3 = qc3.begin();
    engine3.execute(it3);
    EXPECT_NEAR(0, norm(engine3.get_state() - mket({1}, 3)), 1e-12);
    engine3.execute(++it3);
    EXPECT_NEAR(0, norm(engine3.get_state() - mket({2}, 3)), 1e-12);
}

/// BEGIN QNoisyEngineT& QNoisyEngineT::execute(
///       const typename QCircuitTraits<QCircuit>::value_type& elem) override
TEST(qpp_QNoisyEngineT_execute, ValueType) {
    // qubits
    QCircuit qc(1);
    qc.gate(gt.X, 0);
    QubitBitFlipNoise noise(0);
    QNoisyEngine engine{qc, noise};
    const auto elem = *qc.begin();
    engine.execute(elem);
    EXPECT_NEAR(0, norm(engine.get_state() - st.z1), 1e-12);
    QDensityNoisyEngine density_engine{qc, noise};
    density_engine.execute(elem);
    EXPECT_NEAR(0, norm(density_engine.get_state() - prj(st.z1)), 1e-12);

    // qutrits
    QCircuit qc3(1, 0, 3);
    qc3.gate(gt.Xd(3), 0);
    QuditDepolarizingNoise noise3(0, 3);
    QNoisyEngine engine3{qc3, noise3};
    const auto elem3 = *qc3.begin();
    engine3.execute(elem3);
    EXPECT_NEAR(0, norm(engine3.get_state() - mket({1}, 3)), 1e-12);
    QDensityNoisyEngine density_engine3{qc3, noise3};
    density_engine3.execute(elem3);
    EXPECT_NEAR(0, norm(density_engine3.get_state() - mprj({1}, 3)), 1e-12);
}

/// BEGIN QNoisyEngineT& QNoisyEngineT::execute(idx reps = 1) override
TEST(qpp_QNoisyEngineT_execute, AllCircuitWithRepetitions) {
    // qubits
    QCircuit qc(2, 2);
    qc.gate(gt.X, 0).measure_all(0, false);
    QubitBitFlipNoise noise(0);
    QNoisyEngine engine{qc, noise};
    engine.execute(5);
    EXPECT_EQ((std::vector<idx>{1, 0}), engine.get_dits());
    EXPECT_NEAR(0, norm(engine.get_state() - mket({1, 0})), 1e-12);
    EXPECT_EQ(5u, engine.get_stats().get_num_reps());
    EXPECT_THROW(engine.execute(0), exception::OutOfRange);

    // qutrits
    QCircuit qc3(1, 1, 3);
    qc3.gate(gt.Xd(3), 0).gate(gt.Xd(3), 0).measure(0, 0, false);
    QuditDepolarizingNoise noise3(0, 3);
    QNoisyEngine engine3{qc3, noise3};
    engine3.execute(4);
    EXPECT_EQ((std::vector<idx>{2}), engine3.get_dits());
    EXPECT_NEAR(0, norm(engine3.get_state() - mket({2}, 3)), 1e-12);
    EXPECT_EQ(4u, engine3.get_stats().get_num_reps());
}

/// BEGIN std::vector<std::vector<idx>>
///       QNoisyEngineT::get_noise_results() const
TEST(qpp_QNoisyEngineT_get_noise_results, AllTests) {
    // qubits: noise is applied before every non-measurement step
    QCircuit qc(2);
    qc.gate(gt.X, 0).gate(gt.X, 1);
    QubitBitFlipNoise noise(0);
    QNoisyEngine engine{qc, noise};
    engine.execute();
    const auto results = engine.get_noise_results();
    ASSERT_EQ(2u, results.size());
    EXPECT_EQ((std::vector<idx>{0, 0}), results[0]);
    EXPECT_EQ((std::vector<idx>{0, 0}), results[1]);

    // qutrits
    QCircuit qc3(2, 0, 3);
    qc3.gate(gt.Xd(3), 0);
    QuditDepolarizingNoise noise3(0, 3);
    QNoisyEngine engine3{qc3, noise3};
    engine3.execute();
    const auto results3 = engine3.get_noise_results();
    ASSERT_EQ(1u, results3.size());
    EXPECT_EQ((std::vector<idx>{0, 0}), results3[0]);
}

/// BEGIN QEngineT& QNoisyEngineT::reset(std::optional<T> qstate = std::nullopt,
///       bool reset_stats = true)
TEST(qpp_QNoisyEngineT_reset, AllTests) {
    // qubits
    QCircuit qc(1, 1);
    qc.gate(gt.X, 0).measure(0, 0, false);
    QubitBitFlipNoise noise(0);
    QNoisyEngine engine{qc, noise};
    engine.execute(3);
    EXPECT_EQ(3u, engine.get_stats().get_num_reps());
    engine.reset();
    EXPECT_NEAR(0, norm(engine.get_state() - st.z0), 1e-12);
    EXPECT_EQ((std::vector<idx>{0}), engine.get_dits());
    EXPECT_EQ(0u, engine.get_stats().get_num_reps());
    EXPECT_TRUE(engine.get_noise_results().empty());
    engine.reset(st.z1, false);
    EXPECT_NEAR(0, norm(engine.get_state() - st.z1), 1e-12);

    // qutrits
    QCircuit qc3(1, 0, 3);
    qc3.gate(gt.Xd(3), 0);
    QuditDepolarizingNoise noise3(0, 3);
    QNoisyEngine engine3{qc3, noise3};
    engine3.execute();
    engine3.reset(mket({2}, 3));
    EXPECT_NEAR(0, norm(engine3.get_state() - mket({2}, 3)), 1e-12);
    EXPECT_TRUE(engine3.get_noise_results().empty());
}

/// BEGIN std::string QNoisyEngineT::traits_get_name() const override
TEST(qpp_QNoisyEngineT_traits_get_name, AllTests) {
    QCircuit qc(1);
    QubitBitFlipNoise noise(0);
    QNoisyEngine engine{qc, noise};
    QKetNoisyEngine ket_engine{qc, noise};
    QDensityNoisyEngine density_engine{qc, noise};
    EXPECT_EQ("QNoisyEngine", engine.traits_get_name());
    EXPECT_EQ("QKetNoisyEngine", ket_engine.traits_get_name());
    EXPECT_EQ("QDensityNoisyEngine", density_engine.traits_get_name());

    QCircuit qc3(1, 0, 3);
    QuditDepolarizingNoise noise3(0, 3);
    QNoisyEngine engine3{qc3, noise3};
    EXPECT_EQ("QNoisyEngine", engine3.traits_get_name());
}

/// BEGIN bool QNoisyEngineT::traits_is_noisy() const override
TEST(qpp_QNoisyEngineT_traits_is_noisy, AllTests) {
    QCircuit qc(1);
    QubitBitFlipNoise noise(0);
    QNoisyEngine engine{qc, noise};
    QDensityNoisyEngine density_engine{qc, noise};
    EXPECT_TRUE(engine.traits_is_noisy());
    EXPECT_TRUE(density_engine.traits_is_noisy());

    QCircuit qc3(1, 0, 3);
    QuditDepolarizingNoise noise3(0, 3);
    QNoisyEngine engine3{qc3, noise3};
    EXPECT_TRUE(engine3.traits_is_noisy());
}

/// BEGIN bool QNoisyEngineT::traits_is_pure() const override
TEST(qpp_QNoisyEngineT_traits_is_pure, AllTests) {
    QCircuit qc(1);
    QubitBitFlipNoise noise(0);
    QNoisyEngine engine{qc, noise};
    QKetNoisyEngine ket_engine{qc, noise};
    QDensityNoisyEngine density_engine{qc, noise};
    EXPECT_TRUE(engine.traits_is_pure());
    EXPECT_TRUE(ket_engine.traits_is_pure());
    EXPECT_FALSE(density_engine.traits_is_pure());

    QCircuit qc3(1, 0, 3);
    QuditDepolarizingNoise noise3(0, 3);
    QNoisyEngine engine3{qc3, noise3};
    EXPECT_TRUE(engine3.traits_is_pure());
}
