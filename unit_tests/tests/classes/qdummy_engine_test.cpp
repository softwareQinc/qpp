#include "gtest/gtest.h"

#include "qpp/qpp.hpp"

using namespace qpp;

// Unit testing "qpp/classes/qdummy_engine.hpp"

/// BEGIN QBaseEngine& QDummyEngineT::execute(
///       typename QCircuitTraits<QCT>::iterator_type& it)
TEST(qpp_QDummyEngine_execute, Iterator) {
    // qubits
    QCircuit qc(2);
    qc.gate(gt.X, 0);
    QKetDummyEngine engine{qc};
    auto it = qc.begin();
    EXPECT_EQ(&engine, &engine.execute(it));
    EXPECT_EQ(qc.begin(), it);
    EXPECT_EQ(0, engine.get_state().size());

    // qutrits
    QCircuit qc3(2, 0, 3);
    qc3.gate(gt.Xd(3), 1);
    QKetDummyEngine engine3{qc3};
    auto it3 = qc3.begin();
    EXPECT_EQ(&engine3, &engine3.execute(it3));
    EXPECT_EQ(qc3.begin(), it3);
    EXPECT_EQ(0, engine3.get_state().size());
}

/// BEGIN QBaseEngine& QDummyEngine::execute(
///       const typename QCircuitTraits<QCircuit>::value_type& elem) override
TEST(qpp_QDummyEngine_execute, ValueType) {
    // qubits
    QCircuit qc(1);
    qc.gate(gt.X, 0);
    QKetDummyEngine engine{qc};
    const auto& elem = *qc.begin();
    EXPECT_EQ(&engine, &engine.execute(elem));
    EXPECT_EQ(0, engine.get_state().size());

    // qutrits
    QCircuit qc3(1, 0, 3);
    qc3.gate(gt.Xd(3), 0);
    QKetDummyEngine engine3{qc3};
    const auto& elem3 = *qc3.begin();
    EXPECT_EQ(&engine3, &engine3.execute(elem3));
    EXPECT_EQ(0, engine3.get_state().size());
}

/// BEGIN QBaseEngine& QDummyEngine::execute(idx reps = 1) override
TEST(qpp_QDummyEngine_execute, AllCircuitWithRepetitions) {
    // qubits
    QCircuit qc(1);
    qc.gate(gt.X, 0);
    QKetDummyEngine engine{qc};
    EXPECT_EQ(&engine, &engine.execute());
    EXPECT_EQ(&engine, &engine.execute(0));
    EXPECT_EQ(&engine, &engine.execute(5));
    EXPECT_EQ(0, engine.get_state().size());

    // qutrits
    QCircuit qc3(1, 0, 3);
    qc3.gate(gt.Xd(3), 0);
    QKetDummyEngine engine3{qc3};
    EXPECT_EQ(&engine3, &engine3.execute(3));
    EXPECT_EQ(0, engine3.get_state().size());
}

/// BEGIN const QCT& QDummyEngine::get_circuit() const& noexcept
TEST(qpp_QDummyEngine_get_circuit, Lvalue) {
    // qubits
    QCircuit qc(2);
    qc.gate(gt.X, 0);
    const QKetDummyEngine engine{qc};
    EXPECT_EQ(&qc, &engine.get_circuit());
    EXPECT_EQ(qc, engine.get_circuit());

    // qutrits
    QCircuit qc3(2, 0, 3);
    qc3.gate(gt.Xd(3), 1);
    const QKetDummyEngine engine3{qc3};
    EXPECT_EQ(&qc3, &engine3.get_circuit());
    EXPECT_EQ(qc3, engine3.get_circuit());
}

/// BEGIN QCT QDummyEngine::get_circuit() const&& noexcept
TEST(qpp_QDummyEngine_get_circuit, Rvalue) {
    // qubits
    QCircuit qc(2);
    qc.gate(gt.X, 0);
    const QKetDummyEngine engine{qc};
    const QCircuit copy = std::move(engine).get_circuit();
    EXPECT_EQ(qc, copy);
    EXPECT_NE(&qc, &copy);

    // qutrits
    QCircuit qc3(2, 0, 3);
    qc3.gate(gt.Xd(3), 1);
    const QKetDummyEngine engine3{qc3};
    const QCircuit copy3 = std::move(engine3).get_circuit();
    EXPECT_EQ(qc3, copy3);
    EXPECT_NE(&qc3, &copy3);
}

/// BEGIN T QDummyEngine::get_state() const override
TEST(qpp_QDummyEngine_get_state, AllTests) {
    // qubits
    QCircuit qc(2);
    qc.gate(gt.X, 0);
    QKetDummyEngine pure{qc};
    QDensityDummyEngine mixed{qc};
    EXPECT_EQ(0, pure.get_state().size());
    EXPECT_EQ(0, mixed.get_state().size());
    pure.execute();
    mixed.execute();
    EXPECT_EQ(0, pure.get_state().size());
    EXPECT_EQ(0, mixed.get_state().size());

    // qutrits
    QCircuit qc3(1, 0, 3);
    qc3.gate(gt.Xd(3), 0);
    QKetDummyEngine pure3{qc3};
    QDensityDummyEngine mixed3{qc3};
    EXPECT_EQ(0, pure3.get_state().size());
    EXPECT_EQ(0, mixed3.get_state().size());
}

/// BEGIN QBaseEngine& QDummyEngine::set_state(const T& state) override
TEST(qpp_QDummyEngine_set_state, AllTests) {
    // qubits
    QCircuit qc(1);
    QKetDummyEngine pure{qc};
    QDensityDummyEngine mixed{qc};
    EXPECT_EQ(&pure, &pure.set_state(st.z1));
    EXPECT_EQ(&mixed, &mixed.set_state(prj(st.z1)));
    EXPECT_EQ(0, pure.get_state().size());
    EXPECT_EQ(0, mixed.get_state().size());

    // qutrits
    QCircuit qc3(1, 0, 3);
    QKetDummyEngine pure3{qc3};
    QDensityDummyEngine mixed3{qc3};
    const ket psi3 = mket({2}, 3);
    EXPECT_EQ(&pure3, &pure3.set_state(psi3));
    EXPECT_EQ(&mixed3, &mixed3.set_state(prj(psi3)));
    EXPECT_EQ(0, pure3.get_state().size());
    EXPECT_EQ(0, mixed3.get_state().size());
}

/// BEGIN std::string QDummyEngine::to_JSON(
///       bool enclosed_in_curly_brackets = true) const override
TEST(qpp_QDummyEngine_to_JSON, AllTests) {
    // qubits
    QCircuit qc(2, 1);
    QKetDummyEngine pure{qc};
    const std::string json = pure.to_JSON();
    ASSERT_GE(json.size(), 2u);
    EXPECT_EQ('{', json.front());
    EXPECT_EQ('}', json.back());
    EXPECT_NE(std::string::npos, json.find("\"name\": \"QKetDummyEngine\""));
    EXPECT_NE(std::string::npos, json.find("\"is_noisy\": false"));
    EXPECT_NE(std::string::npos, json.find("\"is_pure\": true"));
    EXPECT_NE(std::string::npos, json.find("\"nq\": 2"));
    EXPECT_NE(std::string::npos, json.find("\"nc\": 1"));
    EXPECT_NE(std::string::npos, json.find("\"d\": 2"));
    EXPECT_EQ(json.substr(1, json.size() - 2), pure.to_JSON(false));

    // qutrits
    QCircuit qc3(1, 0, 3);
    QDensityDummyEngine mixed3{qc3};
    const std::string json3 = mixed3.to_JSON();
    EXPECT_NE(std::string::npos,
              json3.find("\"name\": \"QDensityDummyEngine\""));
    EXPECT_NE(std::string::npos, json3.find("\"is_pure\": false"));
    EXPECT_NE(std::string::npos, json3.find("\"d\": 3"));
    EXPECT_EQ(json3.substr(1, json3.size() - 2), mixed3.to_JSON(false));
}

/// BEGIN std::string QDummyEngine::traits_get_name() const override
TEST(qpp_QDummyEngine_traits_get_name, AllTests) {
    // qubits
    QCircuit qc(1);
    QDummyEngine<ket, QCircuit> base{qc};
    QKetDummyEngine pure{qc};
    QDensityDummyEngine mixed{qc};
    EXPECT_EQ("QDummyEngine", base.traits_get_name());
    EXPECT_EQ("QKetDummyEngine", pure.traits_get_name());
    EXPECT_EQ("QDensityDummyEngine", mixed.traits_get_name());

    // qutrits
    QCircuit qc3(1, 0, 3);
    QKetDummyEngine pure3{qc3};
    QDensityDummyEngine mixed3{qc3};
    EXPECT_EQ("QKetDummyEngine", pure3.traits_get_name());
    EXPECT_EQ("QDensityDummyEngine", mixed3.traits_get_name());
}

/// BEGIN bool QDummyEngine::traits_is_noisy() const override
TEST(qpp_QDummyEngine_traits_is_noisy, AllTests) {
    // qubits
    QCircuit qc(1);
    QDummyEngine<ket, QCircuit> base{qc};
    QKetDummyEngine pure{qc};
    QDensityDummyEngine mixed{qc};
    EXPECT_FALSE(base.traits_is_noisy());
    EXPECT_FALSE(pure.traits_is_noisy());
    EXPECT_FALSE(mixed.traits_is_noisy());

    // qutrits
    QCircuit qc3(1, 0, 3);
    QKetDummyEngine pure3{qc3};
    QDensityDummyEngine mixed3{qc3};
    EXPECT_FALSE(pure3.traits_is_noisy());
    EXPECT_FALSE(mixed3.traits_is_noisy());
}

/// BEGIN bool QDummyEngine::traits_is_pure() const override
TEST(qpp_QDummyEngine_traits_is_pure, AllTests) {
    // qubits
    QCircuit qc(1);
    QDummyEngine<ket, QCircuit> base_pure{qc};
    QDummyEngine<cmat, QCircuit> base_mixed{qc};
    QKetDummyEngine pure{qc};
    QDensityDummyEngine mixed{qc};
    EXPECT_TRUE(base_pure.traits_is_pure());
    EXPECT_FALSE(base_mixed.traits_is_pure());
    EXPECT_TRUE(pure.traits_is_pure());
    EXPECT_FALSE(mixed.traits_is_pure());

    // qutrits
    QCircuit qc3(1, 0, 3);
    QKetDummyEngine pure3{qc3};
    QDensityDummyEngine mixed3{qc3};
    EXPECT_TRUE(pure3.traits_is_pure());
    EXPECT_FALSE(mixed3.traits_is_pure());
}
