#include <vector>

#include "gtest/gtest.h"

#include "qpp/qpp.hpp"

using namespace qpp;

// Unit testing "qpp/classes/noise.hpp"

/// BEGIN idx NoiseBase::get_d() const noexcept
TEST(qpp_NoiseBase_get_d, AllTests) {}

/// BEGIN std::vector<cmat> NoiseBase::get_Ks() const
TEST(qpp_NoiseBase_get_Ks, AllTests) {}

/// BEGIN idx NoiseBase::get_last_idx() const
TEST(qpp_NoiseBase_get_last_idx, AllTests) {}

/// BEGIN cmat NoiseBase::get_last_K() const
TEST(qpp_NoiseBase_get_last_K, AllTests) {}

/// BEGIN realT NoiseBase::get_last_p() const
TEST(qpp_NoiseBase_get_last_p, AllTests) {}

/// BEGIN std::vector<realT> NoiseBase::get_probs() const
TEST(qpp_NoiseBase_get_probs, AllTests) {}

/// BEGIN virtual cmat NoiseBase::operator()(const cmat& state) const
TEST(qpp_NoiseBase_functor, AllQudits) {}

/// BEGIN virtual cmat NoiseBase::operator()(const cmat& state, idx target)
///       const
TEST(qpp_NoiseBase_functor, SpecificQudit) {}

/// BEGIN virtual cmat NoiseBase::operator()(const cmat& state,
///       const std::vector<idx>& target) const
TEST(qpp_NoiseBase_functor, CorrelatedNoise) {}

/// BEGIN explicit QubitAmplitudeDampingNoise::QubitAmplitudeDampingNoise(
///       realT gamma)
TEST(qpp_QubitAmplitudeDampingNoise_QubitAmplitudeDampingNoise, AllTests) {}

/// BEGIN explicit QubitBitFlipNoise::QubitBitFlipNoise(realT p)
TEST(qpp_QubitBitFlipNoise_QubitBitFlipNoise, AllTests) {}

/// BEGIN explicit QubitBitPhaseFlipNoise::QubitBitPhaseFlipNoise(realT p)
TEST(qpp_QubitBitPhaseFlipNoise_QubitBitPhaseFlipNoise, AllTests) {}

/// BEGIN explicit QubitDepolarizingNoise::QubitDepolarizingNoise(realT p)
TEST(qpp_QubitDepolarizingNoise_QubitDepolarizingNoise, AllTests) {}

/// BEGIN explicit QubitPhaseFlipNoise::QubitPhaseFlipNoise(realT p)
TEST(qpp_QubitPhaseFlipNoise_QubitPhaseFlipNoise, AllTests) {}

/// BEGIN explicit QuditDepolarizingNoise::QuditDepolarizingNoise(realT p,
///       idx d)
TEST(qpp_QuditDepolarizingNoise_QuditDepolarizingNoise, AllTests) {}
