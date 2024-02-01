#include <vector>

#include "gtest/gtest.h"

#include "qpp/qpp.h"

using namespace qpp;

// Unit testing "classes/engines.hpp"

/// BEGIN QEngine& QEngine::execute(const QCircuit::iterator& it)
TEST(qpp_QEngine_execute, Iterator) {}

/// BEGIN virtual QEngine& QEngine::execute(
///       const QCircuit::iterator::value_type& elem)
TEST(qpp_QEngine_execute, ValueType) {}

/// BEGIN QEngine& QEngine::execute(idx reps = 1, bool try_sampling = true)
TEST(qpp_QEngine_execute, AllCircuitWithRepetitions) {}

/// BEGIN const QEngine& QEngine::get_circuit() const& noexcept
TEST(qpp_QEngine_get_circuit, Lvalue) {}

/// BEGIN QCircuit get_circuit() const&& noexcept
TEST(qpp_QEngine_get_circuit, Rvalue) {}

/// BEGIN idx QEngine::get_dit(idx i) const
TEST(qpp_QEngine_get_dit, AllTests) {}

/// BEGIN std::vector<idx> QEngine::get_dits() const
TEST(qpp_QEngine_get_dits, AllTests) {}

/// BEGIN bool QEngine::get_measured(idx i) const
TEST(qpp_QEngine_get_measured, SpecificQudit) {}

/// BEGIN std::vector<idx> QEngine::get_non_measured() const
TEST(qpp_QEngine_get_non_measured, AllTests) {}

/// BEGIN std::vector<realT> QEngine::get_probs() const
TEST(qpp_QEngine_get_probs, AllTests) {}

/// BEGIN ket QEngine::get_state() const
TEST(qpp_QEngine_get_state, AllTests) {}

/// BEGIN std::map<std::string, idx, internal::EqualSameSizeStringDits>
///       QEngine::get_stats() const
TEST(qpp_QEngine_get_stats, AllTests) {}

/// BEGIN QEngine& QEngine::reset(bool reset_stats = true)
TEST(qpp_QEngine_reset, AllTests) {}

/// BEGIN QEngine& QEngine::reset_stats()
TEST(qpp_QEngine_reset_stats, AllTests) {}

/// BEGIN QEngine& QEngine::set_dit(idx i, idx value)
TEST(qpp_QEngine_set_dit, AllTests) {}

/// BEGIN QEngine& QEngine::set_dits(std::vector<idx> dits)
TEST(qpp_QEngine_set_dits, AllTests) {}

/// BEGIN QEngine& QEngine::set_psi(const ket& psi)
TEST(qpp_QEngine_set_psi, AllTests) {}

/// BEGIN std::string QEngine::to_JSON(
///       bool enclosed_in_curly_brackets = true) const override
TEST(qpp_QEngine_to_JSON, AllTests) {}

/// BEGIN virtual std::string QEngine::traits_get_name() const
TEST(qpp_QEngine_traits_get_name, AllTests) {}

/// BEGIN virtual bool QEngine::traits_is_mixed() const
TEST(qpp_QEngine_traits_is_mixed, AllTests) {}

/// BEGIN virtual bool QEngine::traits_is_noisy() const
TEST(qpp_QEngine_traits_is_noisy, AllTests) {}

/// BEGIN virtual bool QEngine::traits_is_pure() const
TEST(qpp_QEngine_traits_is_pure, AllTests) {}

/// BEGIN std::vector<idx> QEngine::was_measured() const
TEST(qpp_QEngine_was_measured, AllQudits) {}

/// BEGIN QNoisyEngine& QNoisyEngine::execute(
///       const QCircuit::iterator::value_type& elem) override
TEST(qpp_QNoisyEngine_execute, ValueType) {}

/// BEGIN QEngine& QNoisyEngine::execute(idx reps = 1, bool = true)
TEST(qpp_QNoisyEngine_execute, AllCircuitWithRepetitions) {}

/// BEGIN std::vector<std::vector<idx>>
///       QNoisyEngine::get_noise_results() const
TEST(qpp_QNoisyEngine_get_noise_results, AllTests) {}

/// BEGIN QEngine& QNoisyEngine::reset(bool reset_stats = true)
TEST(qpp_QNoisyEngine_reset, AllTests) {}

/// BEGIN virtual std::string QNoisyEngine::traits_get_name() const
TEST(qpp_QNoisyEngine_traits_get_name, AllTests) {}

/// BEGIN virtual bool QNoisyEngine::traits_is_mixed() const
TEST(qpp_QNoisyEngine_traits_is_mixed, AllTests) {}

/// BEGIN virtual bool QNoisyEngine::traits_is_noisy() const
TEST(qpp_QNoisyEngine_traits_is_noisy, AllTests) {}

/// BEGIN virtual bool QNoisyEngine::traits_is_pure() const
TEST(qpp_QNoisyEngine_traits_is_pure, AllTests) {}
