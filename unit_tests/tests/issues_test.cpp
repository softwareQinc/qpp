#include "gtest/gtest.h"

#include "qpp/qpp.hpp"

using namespace qpp;

// Unit testing reported issues and pull requests

TEST(qpp_PR, 110) {
    auto const state = (0_ket + 1_ket).normalized().eval();
    auto const res = measure(state, gt.Id2);
    auto const& resulting_states = std::get<measure_idx::st>(res);

    EXPECT_EQ(resulting_states[0][0], static_cast<realT>(1.));
    EXPECT_EQ(resulting_states[1][1], static_cast<realT>(1.));
}

TEST(qpp_PR, 113) {
    auto circuit = QCircuit{1, 1};
    circuit.measureV(gt.Id2, 0, 0);

    auto engine = QEngine{circuit};
    EXPECT_NO_THROW(engine.execute());
}

TEST(qpp_issue, 76) {
    // https://github.com/softwareQinc/qpp/issues/76
    ket psi(16);
    psi << 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1;
    psi /= std::sqrt(static_cast<realT>(6));

    auto qc = state_prep_circuit(psi);
    auto result = QEngine{qc}.execute().get_state();
    EXPECT_NEAR(0, norm(result - psi), 1e-5);
}
