#include "gtest/gtest.h"

#include "qpp/qpp.h"

using namespace qpp;

// Unit testing reported issues and pull requests

/******************************************************************************/
TEST(qpp_PR, 110) {
    auto const state = (0_ket + 1_ket).normalized().eval();
    auto const res = measure(state, gt.Id2);
    auto const& resulting_states = std::get<ST>(res);

    EXPECT_EQ(resulting_states[0][0], static_cast<realT>(1.));
    EXPECT_EQ(resulting_states[1][1], static_cast<realT>(1.));
}
/******************************************************************************/
TEST(qpp_PR, 113) {
    auto circuit = QCircuit{1, 1};
    circuit.measureV(gt.Id2, 0, 0);

    auto engine = QEngine{circuit};
    EXPECT_NO_THROW(engine.execute());
}
/******************************************************************************/
