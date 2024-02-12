#include "gtest/gtest.h"

#include "qpp/qpp.h"

using namespace qpp;

// Unit testing "classes/qnoisy_engine.hpp"

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
