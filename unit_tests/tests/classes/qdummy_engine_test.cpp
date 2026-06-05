#include "gtest/gtest.h"

#include "qpp/qpp.hpp"

using namespace qpp;

// Unit testing "qpp/classes/qdummy_engine.hpp"

/// BEGIN QBaseEngine& QDummyEngineT::execute(
///       typename QCircuitTraits<QCT>::iterator_type& it)
TEST(qpp_QDummyEngine_execute, Iterator) {}

/// BEGIN QBaseEngine& QDummyEngine::execute(
///       const typename QCircuitTraits<QCircuit>::value_type& elem) override
TEST(qpp_QDummyEngine_execute, ValueType) {}

/// BEGIN QBaseEngine& QDummyEngine::execute(idx reps = 1) override
TEST(qpp_QDummyEngine_execute, AllCircuitWithRepetitions) {}

/// BEGIN const QCT& QDummyEngine::get_circuit() const& noexcept
TEST(qpp_QDummyEngine_get_circuit, Lvalue) {}

/// BEGIN QCT QDummyEngine::get_circuit() const&& noexcept
TEST(qpp_QDummyEngine_get_circuit, Rvalue) {}

/// BEGIN T QDummyEngine::get_state() const override
TEST(qpp_QDummyEngine_get_state, AllTests) {}

/// BEGIN QBaseEngine& QDummyEngine::set_state(const T& state) override
TEST(qpp_QDummyEngine_set_state, AllTests) {}

/// BEGIN std::string QDummyEngine::to_JSON(
///       bool enclosed_in_curly_brackets = true) const override
TEST(qpp_QDummyEngine_to_JSON, AllTests) {}

/// BEGIN std::string QDummyEngine::traits_get_name() const override
TEST(qpp_QDummyEngine_traits_get_name, AllTests) {}

/// BEGIN bool QDummyEngine::traits_is_noisy() const override
TEST(qpp_QDummyEngine_traits_is_noisy, AllTests) {}

/// BEGIN bool QDummyEngine::traits_is_pure() const override
TEST(qpp_QDummyEngine_traits_is_pure, AllTests) {}
