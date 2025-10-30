#include "gtest/gtest.h"

#include "qpp/qpp.hpp"

using namespace qpp;

// Unit testing "qpp/classes/qengine.hpp"

/// BEGIN QBaseEngine& QEngineT::execute(
///       typename QCircuitTraits<QCT>::iterator_type& it)
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

/// BEGIN bool QEngineT::get_ensure_post_selection() const
TEST(qpp_QEngineT_get_ensure_post_selection, AllTests) {}

/// BEGIN idx QEngineT::get_max_post_selection_reps() const
TEST(qpp_QEngineT_get_max_post_selection_reps, AllTests) {}

/// BEGIN std::vector<idx> QEngineT::get_measured_d() const
TEST(qpp_QEngineT_get_measured_d, AllTests) {}

/// BEGIN std::vector<idx> QEngineT::get_non_measured_d() const
TEST(qpp_QEngineT_get_non_measured_d, AllTests) {}

/// BEGIN std::vector<realT> QEngineT::get_probs() const
TEST(qpp_QEngineT_get_probs, AllTests) {}

/// BEGIN T QEngineT::get_state() const override
TEST(qpp_QEngineT_get_state, AllTests) {}

/// BEGIN internal::QEngineStatistics QEngineT::get_stats() const
TEST(qpp_QEngineT_get_stats, AllTests) {}

/// BEGIN bool QEngineT::post_select_ok() const
TEST(qpp_QEngineT_post_select_ok, AllTests) {}

/// BEGIN QEngineT& QEngineT::reset(bool reset_stats = true)
TEST(qpp_QEngineT_reset, AllTests) {
    auto circuit = qpp::QCircuit{2, 1};
    circuit.measure(0, 0);

    auto engine = qpp::QEngine{circuit};
    auto const psi = qpp::randket(4).eval();
    engine.reset().set_state(psi).execute();

    EXPECT_NO_THROW(engine.reset().set_state(psi).execute());
}

/// BEGIN QEngineT& QEngineT::reset_stats()
TEST(qpp_QEngineT_reset_stats, AllTests) {}

/// BEGIN QEngineT& QEngineT::set_dit(idx i, idx value)
TEST(qpp_QEngineT_set_dit, AllTests) {}

/// BEGIN QEngineT& QEngineT::set_dits(std::vector<idx> dits)
TEST(qpp_QEngineT_set_dits, AllTests) {}

/// BEGIN QEngineT& QEngineT::set_ensure_post_selection(bool val)
TEST(qpp_QEngineT_set_ensure_post_selection, AllTests) {}

/// BEGIN QEngineT& QEngineT::set_max_post_selection_reps(
///       idx max_post_selection_reps)
TEST(qpp_QEngineT_set_max_post_selection_reps, AllTests) {}

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

/// BEGIN bool QEngineT::was_measured_d(idx i) const
TEST(qpp_QEngineT_was_measured_d, AllTests) {}
