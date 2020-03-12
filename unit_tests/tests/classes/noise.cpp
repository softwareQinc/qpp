#include <vector>
#include "gtest/gtest.h"
#include "qpp.h"

using namespace qpp;

// Unit testing "classes/noise.h"

/******************************************************************************/
/// BEGIN idx qpp::NoiseBase::get_d() const noexcept
TEST(qpp_NoiseBase_get_d, AllTests) {}
/******************************************************************************/
/// BEGIN std::vector<cmat> qpp::NoiseBase::get_Ks() const
TEST(qpp_NoiseBase_get_Ks, AllTests) {}
/******************************************************************************/
/// BEGIN idx qpp::NoiseBase::get_last_idx() const
TEST(qpp_NoiseBase_get_last_idx, AllTests) {}
/******************************************************************************/
/// BEGIN cmat qpp::NoiseBase::get_last_K() const
TEST(qpp_NoiseBase_get_last_K, AllTests) {}
/******************************************************************************/
/// BEGIN double qpp::NoiseBase::get_last_p() const
TEST(qpp_NoiseBase_get_last_p, AllTests) {}
/******************************************************************************/
/// BEGIN std::vector<double> qpp::NoiseBase::get_probs() const
TEST(qpp_NoiseBase_get_probs, AllTests) {}
/******************************************************************************/
/// BEGIN virtual cmat qpp::NoiseBase::operator()(const cmat& state) const
TEST(qpp_NoiseBase_functor, AllQudits) {}
/******************************************************************************/
/// BEGIN virtual cmat qpp::NoiseBase::operator()(const cmat& state, idx target)
///       const
TEST(qpp_NoiseBase_functor, SpecificQudit) {}
/******************************************************************************/
/// BEGIN virtual cmat qpp::NoiseBase::operator()(const cmat& state,
///       const std::vector<idx>& target) const
TEST(qpp_NoiseBase_functor, CorrelatedNoise) {}
/******************************************************************************/

/******************************************************************************/
/// BEGIN explicit qpp::QubitAmplitudeDampingNoise::QubitAmplitudeDampingNoise(
///       double gamma)
TEST(qpp_QubitAmplitudeDampingNoise_QubitAmplitudeDampingNoise, AllTests) {}
/******************************************************************************/

/******************************************************************************/
/// BEGIN explicit qpp::QubitBitFlipNoise::QubitBitFlipNoise(double p)
TEST(qpp_QubitBitFlipNoise_QubitBitFlipNoise, AllTests) {}
/******************************************************************************/

/******************************************************************************/
/// BEGIN explicit qpp::QubitBitPhaseFlipNoise::QubitBitPhaseFlipNoise(double p)
TEST(qpp_QubitBitPhaseFlipNoise_QubitBitPhaseFlipNoise, AllTests) {}
/******************************************************************************/

/******************************************************************************/
/// BEGIN explicit qpp::QubitDepolarizingNoise::QubitDepolarizingNoise(double p)
TEST(qpp_QubitDepolarizingNoise_QubitDepolarizingNoise, AllTests) {}
/******************************************************************************/

/******************************************************************************/
/// BEGIN explicit qpp::QubitPhaseFlipNoise::QubitPhaseFlipNoise(double p)
TEST(qpp_QubitPhaseFlipNoise_QubitPhaseFlipNoise, AllTests) {}
/******************************************************************************/

/******************************************************************************/
/// BEGIN explicit qpp::QuditDepolarizingNoise::QuditDepolarizingNoise(double p,
///       idx d)
TEST(qpp_QuditDepolarizingNoise_QuditDepolarizingNoise, AllTests) {}
/******************************************************************************/
