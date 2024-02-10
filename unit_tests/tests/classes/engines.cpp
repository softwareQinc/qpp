#include "gtest/gtest.h"

#include "qpp/qpp.h"

using namespace qpp;

// Unit testing "classes/engines.hpp"

/// BEGIN QEngineT& QEngineT::execute(const QCircuit::iterator& it)
TEST(qpp_QEngineT_execute, Iterator) {}

/// BEGIN virtual QEngineT& QEngineT::execute(
///       const QCircuit::iterator::value_type& elem) override
TEST(qpp_QEngineT_execute, ValueType) {}

/// BEGIN QEngineT& QEngineT::execute(idx reps = 1) override
TEST(qpp_QEngineT_execute, AllCircuitWithRepetitions) {}

/// BEGIN const QEngineT& QEngineT::get_circuit() const& noexcept
TEST(qpp_QEngineT_get_circuit, Lvalue) {}

/// BEGIN QCircuit get_circuit() const&& noexcept
TEST(qpp_QEngineT_get_circuit, Rvalue) {}

/// BEGIN idx QEngineT::get_dit(idx i) const
TEST(qpp_QEngineT_get_dit, AllTests) {}

/// BEGIN std::vector<idx> QEngineT::get_dits() const
TEST(qpp_QEngineT_get_dits, AllTests) {}

/// BEGIN bool QEngineT::get_measured(idx i) const
TEST(qpp_QEngineT_get_measured, SpecificQudit) {}

/// BEGIN std::vector<idx> QEngineT::get_non_measured() const
TEST(qpp_QEngineT_get_non_measured, AllTests) {}

/// BEGIN std::vector<realT> QEngineT::get_probs() const
TEST(qpp_QEngineT_get_probs, AllTests) {}

/// BEGIN ket QEngineT::get_state() const
TEST(qpp_QEngineT_get_state, AllTests) {}

/// BEGIN std::map<std::string, idx, internal::EqualSameSizeStringDits>
///       QEngineT::get_stats() const
TEST(qpp_QEngineT_get_stats, AllTests) {}

/// BEGIN QEngineT& QEngineT::reset(bool reset_stats = true)
TEST(qpp_QEngineT_reset, AllTests) {}

/// BEGIN QEngineT& QEngineT::reset_stats()
TEST(qpp_QEngineT_reset_stats, AllTests) {}

/// BEGIN QEngineT& QEngineT::set_dit(idx i, idx value)
TEST(qpp_QEngineT_set_dit, AllTests) {}

/// BEGIN QEngineT& QEngineT::set_dits(std::vector<idx> dits)
TEST(qpp_QEngineT_set_dits, AllTests) {}

/// BEGIN QEngineT& QEngineT::set_psi(const ket& psi)
TEST(qpp_QEngineT_set_state, AllTests) {}

/// BEGIN std::string QEngineT::to_JSON(
///       bool enclosed_in_curly_brackets = true) const override
TEST(qpp_QEngineT_to_JSON, AllTests) {}

/// BEGIN virtual std::string QEngineT::traits_get_name() const
TEST(qpp_QEngineT_traits_get_name, AllTests) {}

/// BEGIN virtual bool QEngineT::traits_is_noisy() const
TEST(qpp_QEngineT_traits_is_noisy, AllTests) {}

/// BEGIN virtual bool QEngineT::traits_is_pure() const
TEST(qpp_QEngineT_traits_is_pure, AllTests) {}

/// BEGIN std::vector<idx> QEngineT::was_measured() const
TEST(qpp_QEngineT_was_measured, AllQudits) {}

/// BEGIN QNoisyEngineT& QNoisyEngineT::execute(const QCircuit::iterator& it)
TEST(qpp_QNoisyEngineT_execute, Iterator) {}

/// BEGIN QNoisyEngineT& QNoisyEngineT::execute(
///       const QCircuit::iterator::value_type& elem) override
TEST(qpp_QNoisyEngineT_execute, ValueType) {}

/// BEGIN QEngineT& QNoisyEngineT::execute(idx reps = 1) override
TEST(qpp_QNoisyEngineT_execute, AllCircuitWithRepetitions) {}

/// BEGIN std::vector<std::vector<idx>>
///       QNoisyEngineT::get_noise_results() const
TEST(qpp_QNoisyEngineT_get_noise_results, AllTests) {}

/// BEGIN QEngineT& QNoisyEngineT::reset(bool reset_stats = true)
TEST(qpp_QNoisyEngineT_reset, AllTests) {}

/// BEGIN virtual std::string QNoisyEngineT::traits_get_name() const
TEST(qpp_QNoisyEngineT_traits_get_name, AllTests) {}

/// BEGIN virtual bool QNoisyEngineT::traits_is_noisy() const
TEST(qpp_QNoisyEngineT_traits_is_noisy, AllTests) {}

/// BEGIN virtual bool QNoisyEngineT::traits_is_pure() const
TEST(qpp_QNoisyEngineT_traits_is_pure, AllTests) {}
