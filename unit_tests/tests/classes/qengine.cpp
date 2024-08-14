#include "gtest/gtest.h"

#include "qpp/qpp.h"

using namespace qpp;

// Unit testing "classes/qengine.hpp"

/// BEGIN QBaseEngine& QEngineT::execute(
///       const typename QCircuitTraits<QCT>::iterator_type& it)
TEST(qpp_QEngineT_execute, Iterator) {}

/// BEGIN QEngineT& QEngineT::execute(
///       const typename QCircuitTraits<QCircuit>::value_type& elem) override
TEST(qpp_QEngineT_execute, ValueType) {}

/// BEGIN QEngineT& QEngineT::execute(idx reps = 1) override
TEST(qpp_QEngineT_execute, AllCircuitWithRepetitions) {}

/// BEGIN const QCT& QBaseEngine::get_circuit() const& noexcept
TEST(qpp_QEngineT_get_circuit, Lvalue) {}

/// BEGIN QCT QBaseEngine::get_circuit() const&& noexcept
TEST(qpp_QEngineT_get_circuit, Rvalue) {}

/// BEGIN idx QEngineT::get_dit(idx i) const
TEST(qpp_QEngineT_get_dit, AllTests) {}

/// BEGIN std::vector<idx> QEngineT::get_dits() const
TEST(qpp_QEngineT_get_dits, AllTests) {}

/// BEGIN bool QEngineT::get_measured_destructively() const
TEST(qpp_QEngineT_get_measured_destructively, AllTests) {}

/// BEGIN std::vector<idx> QEngineT::get_non_measured_destructively() const
TEST(qpp_QEngineT_get_non_measured_destructively, AllTests) {}

/// BEGIN std::vector<realT> QEngineT::get_probs() const
TEST(qpp_QEngineT_get_probs, AllTests) {}

/// BEGIN T QEngineT::get_state() const override
TEST(qpp_QEngineT_get_state, AllTests) {}

/// BEGIN internal::QEngineStatistics QEngineT::get_stats() const
TEST(qpp_QEngineT_get_stats, AllTests) {}

/// BEGIN bool QEngineT::post_select_ok() const
TEST(qpp_QEngineT_post_select_ok, AllTests) {}

/// BEGIN QEngineT& QEngineT::reset(bool reset_stats = true,
///       bool ensure_post_selection = false)
TEST(qpp_QEngineT_reset, AllTests) {}

/// BEGIN QEngineT& QEngineT::reset_stats()
TEST(qpp_QEngineT_reset_stats, AllTests) {}

/// BEGIN QEngineT& QEngineT::set_dit(idx i, idx value)
TEST(qpp_QEngineT_set_dit, AllTests) {}

/// BEGIN QEngineT& QEngineT::set_dits(std::vector<idx> dits)
TEST(qpp_QEngineT_set_dits, AllTests) {}

/// BEGIN QEngineT& QEngineT::set_state(const T& state) override
TEST(qpp_QEngineT_set_state, AllTests) {}

/// BEGIN std::string QEngineT::to_JSON(
///       bool enclosed_in_curly_brackets = true) const override
TEST(qpp_QEngineT_to_JSON, AllTests) {}

/// BEGIN  std::string QEngineT::traits_get_name() const override
TEST(qpp_QEngineT_traits_get_name, AllTests) {}

/// BEGIN bool QEngineT::traits_is_noisy() const override
TEST(qpp_QEngineT_traits_is_noisy, AllTests) {}

/// BEGIN bool QEngineT::traits_is_pure() const override
TEST(qpp_QEngineT_traits_is_pure, AllTests) {}

/// BEGIN std::vector<idx> QEngineT::was_measured_destructively(idx i) const
TEST(qpp_QEngineT_was_measured_destructively, AllTests) {}
