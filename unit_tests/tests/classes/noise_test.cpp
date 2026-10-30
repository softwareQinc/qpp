#include <vector>

#include "gtest/gtest.h"

#include "qpp/qpp.hpp"

using namespace qpp;

// Unit testing "qpp/classes/noise.hpp"

/// BEGIN idx NoiseBase::get_d() const noexcept
TEST(qpp_NoiseBase_get_d, AllTests) {
    const QubitBitFlipNoise qubit(0.25);
    const QuditDepolarizingNoise qutrit(0, 3);
    EXPECT_EQ(2u, qubit.get_d());
    EXPECT_EQ(3u, qutrit.get_d());
}

/// BEGIN std::vector<cmat> NoiseBase::get_Ks() const
TEST(qpp_NoiseBase_get_Ks, AllTests) {
    const QubitBitFlipNoise noise(0.25);
    const auto Ks = noise.get_Ks();
    ASSERT_EQ(2u, Ks.size());
    EXPECT_LT((Ks[0] - gt.Id2).norm(), 1e-12);
    EXPECT_LT((Ks[1] - gt.X).norm(), 1e-12);
}

/// BEGIN idx NoiseBase::get_last_idx() const
TEST(qpp_NoiseBase_get_last_idx, AllTests) {
    QubitBitFlipNoise noise(1);
    EXPECT_THROW(noise.get_last_idx(), exception::CustomException);
    noise(prj(st.z0));
    EXPECT_EQ(1u, noise.get_last_idx());
}

/// BEGIN cmat NoiseBase::get_last_K() const
TEST(qpp_NoiseBase_get_last_K, AllTests) {
    QubitPhaseFlipNoise noise(1);
    EXPECT_THROW(noise.get_last_K(), exception::CustomException);
    noise(prj(st.x0));
    EXPECT_LT((noise.get_last_K() - gt.Z).norm(), 1e-12);
}

/// BEGIN realT NoiseBase::get_last_p() const
TEST(qpp_NoiseBase_get_last_p, AllTests) {
    QubitBitFlipNoise noise(1);
    EXPECT_THROW(noise.get_last_p(), exception::CustomException);
    noise(prj(st.z0));
    EXPECT_NEAR(1, noise.get_last_p(), 1e-12);
}

/// BEGIN std::vector<realT> NoiseBase::get_probs() const
TEST(qpp_NoiseBase_get_probs, AllTests) {
    QubitBitFlipNoise independent(0.3);
    EXPECT_EQ((std::vector<realT>{0.7, 0.3}), independent.get_probs());

    QubitAmplitudeDampingNoise dependent(0.36);
    EXPECT_THROW(dependent.get_probs(), exception::CustomException);
    dependent(prj(st.z1));
    const auto probs = dependent.get_probs();
    ASSERT_EQ(2u, probs.size());
    EXPECT_NEAR(0.36, probs[0], 1e-12);
    EXPECT_NEAR(0.64, probs[1], 1e-12);
}

/// BEGIN virtual cmat NoiseBase::operator()(const cmat& state) const
TEST(qpp_NoiseBase_functor, AllQudits) {
    QubitBitFlipNoise noise(1);
    EXPECT_LT((noise(prj(st.z0)) - prj(st.z1)).norm(), 1e-12);
    EXPECT_LT((noise(prj(st.z1)) - prj(st.z0)).norm(), 1e-12);

    QubitPhaseFlipNoise identity(0);
    EXPECT_LT((identity(prj(st.x0)) - prj(st.x0)).norm(), 1e-12);
}

/// BEGIN virtual cmat NoiseBase::operator()(const cmat& state, idx target)
///       const
TEST(qpp_NoiseBase_functor, SpecificQudit) {
    QubitBitFlipNoise noise(1);
    EXPECT_LT((noise(mprj({0, 0}), 0) - mprj({1, 0})).norm(), 1e-12);
    EXPECT_LT((noise(mprj({0, 0}), 1) - mprj({0, 1})).norm(), 1e-12);
}

/// BEGIN virtual cmat NoiseBase::operator()(const cmat& state,
///       const std::vector<idx>& target) const
TEST(qpp_NoiseBase_functor, CorrelatedNoise) {
    NoiseBase<NoiseType::StateIndependent> noise({gt.X}, {1});
    EXPECT_LT((noise(mprj({0, 0}), std::vector<idx>{1}) - mprj({0, 1})).norm(),
              1e-12);
    EXPECT_EQ(0u, noise.get_last_idx());
}

/// BEGIN explicit QubitAmplitudeDampingNoise::QubitAmplitudeDampingNoise(
///       realT gamma)
TEST(qpp_QubitAmplitudeDampingNoise_QubitAmplitudeDampingNoise, AllTests) {
    QubitAmplitudeDampingNoise noise(0.36);
    const auto Ks = noise.get_Ks();
    ASSERT_EQ(2u, Ks.size());
    EXPECT_NEAR(0.6, std::real(Ks[0](1, 1)), 1e-12);
    EXPECT_NEAR(0.8, std::real(Ks[1](0, 1)), 1e-12);
    EXPECT_LT((adjoint(Ks[0]) * Ks[0] + adjoint(Ks[1]) * Ks[1] - gt.Id2).norm(),
              1e-12);
    EXPECT_THROW(QubitAmplitudeDampingNoise(-0.1), exception::OutOfRange);
    EXPECT_THROW(QubitAmplitudeDampingNoise(1.1), exception::OutOfRange);
}

/// BEGIN explicit QubitBitFlipNoise::QubitBitFlipNoise(realT p)
TEST(qpp_QubitBitFlipNoise_QubitBitFlipNoise, AllTests) {
    QubitBitFlipNoise noise(1);
    EXPECT_EQ((std::vector<realT>{0, 1}), noise.get_probs());
    EXPECT_LT((noise(prj(st.z0)) - prj(st.z1)).norm(), 1e-12);
    EXPECT_THROW(QubitBitFlipNoise(-0.1), exception::OutOfRange);
    EXPECT_THROW(QubitBitFlipNoise(1.1), exception::OutOfRange);
}

/// BEGIN explicit QubitBitPhaseFlipNoise::QubitBitPhaseFlipNoise(realT p)
TEST(qpp_QubitBitPhaseFlipNoise_QubitBitPhaseFlipNoise, AllTests) {
    QubitBitPhaseFlipNoise noise(1);
    EXPECT_EQ((std::vector<realT>{0, 1}), noise.get_probs());
    EXPECT_LT((noise(prj(st.z0)) - prj(st.z1)).norm(), 1e-12);
    EXPECT_THROW(QubitBitPhaseFlipNoise(-0.1), exception::OutOfRange);
    EXPECT_THROW(QubitBitPhaseFlipNoise(1.1), exception::OutOfRange);
}

/// BEGIN explicit QubitDepolarizingNoise::QubitDepolarizingNoise(realT p)
TEST(qpp_QubitDepolarizingNoise_QubitDepolarizingNoise, AllTests) {
    QubitDepolarizingNoise noise(0.6);
    const auto probs = noise.get_probs();
    ASSERT_EQ(4u, probs.size());
    EXPECT_NEAR(0.4, probs[0], 1e-12);
    for (idx i = 1; i < 4; ++i) {
        EXPECT_NEAR(0.2, probs[i], 1e-12);
    }
    EXPECT_EQ(4u, noise.get_Ks().size());
    EXPECT_LT((QubitDepolarizingNoise(0)(prj(st.z0)) - prj(st.z0)).norm(),
              1e-12);
    EXPECT_THROW(QubitDepolarizingNoise(-0.1), exception::OutOfRange);
    EXPECT_THROW(QubitDepolarizingNoise(1.1), exception::OutOfRange);
}

/// BEGIN explicit QubitPhaseFlipNoise::QubitPhaseFlipNoise(realT p)
TEST(qpp_QubitPhaseFlipNoise_QubitPhaseFlipNoise, AllTests) {
    QubitPhaseFlipNoise noise(1);
    EXPECT_EQ((std::vector<realT>{0, 1}), noise.get_probs());
    EXPECT_LT((noise(prj(st.x0)) - prj(st.x1)).norm(), 1e-12);
    EXPECT_THROW(QubitPhaseFlipNoise(-0.1), exception::OutOfRange);
    EXPECT_THROW(QubitPhaseFlipNoise(1.1), exception::OutOfRange);
}

/// BEGIN explicit QuditDepolarizingNoise::QuditDepolarizingNoise(realT p,
///       idx d)
TEST(qpp_QuditDepolarizingNoise_QuditDepolarizingNoise, AllTests) {
    QuditDepolarizingNoise noise(0, 3);
    EXPECT_EQ(3u, noise.get_d());
    EXPECT_EQ(9u, noise.get_Ks().size());
    EXPECT_EQ(9u, noise.get_probs().size());
    EXPECT_NEAR(1, noise.get_probs()[0], 1e-12);
    EXPECT_LT((noise(mprj({2}, 3)) - mprj({2}, 3)).norm(), 1e-12);
    EXPECT_THROW(QuditDepolarizingNoise(-0.1, 3), exception::OutOfRange);
    EXPECT_THROW(QuditDepolarizingNoise(1.1, 3), exception::OutOfRange);
}
