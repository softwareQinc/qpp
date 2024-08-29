#include "gtest/gtest.h"

#include "qpp/qpp.hpp"

using namespace qpp;

// Unit testing "qpp/classes/qnoisy_engine.hpp"

/// BEGIN QBaseEngine& QNoisyEngineT::execute(
///       const typename QCircuitTraits<QCT>::iterator_type& it)
TEST(qpp_QNoisyEngineT_execute, Iterator) {}

/// BEGIN QNoisyEngineT& QNoisyEngineT::execute(
///       const typename QCircuitTraits<QCircuit>::value_type& elem) override
TEST(qpp_QNoisyEngineT_execute, ValueType) {}

/// BEGIN QNoisyEngineT& QNoisyEngineT::execute(idx reps = 1) override
TEST(qpp_QNoisyEngineT_execute, AllCircuitWithRepetitions) {}

/// BEGIN std::vector<std::vector<idx>>
///       QNoisyEngineT::get_noise_results() const
TEST(qpp_QNoisyEngineT_get_noise_results, AllTests) {}

/// BEGIN QEngineT& QNoisyEngineT::reset(bool reset_stats = true,
///       bool ensure_post_selection = false, idx max_post_selection_reps =
///       std::numeric_limits<idx>::max())
TEST(qpp_QNoisyEngineT_reset, AllTests) {}

/// BEGIN std::string QNoisyEngineT::traits_get_name() const override
TEST(qpp_QNoisyEngineT_traits_get_name, AllTests) {}

/// BEGIN bool QNoisyEngineT::traits_is_noisy() const override
TEST(qpp_QNoisyEngineT_traits_is_noisy, AllTests) {}

/// BEGIN bool QNoisyEngineT::traits_is_pure() const override
TEST(qpp_QNoisyEngineT_traits_is_pure, AllTests) {}
