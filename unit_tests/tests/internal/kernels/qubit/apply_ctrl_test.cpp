#include "gmock/gmock.h" // for matchers
#include "gtest/gtest.h"

#include "qpp/qpp.hpp"

#include "qpp/internal/kernels/qubit/apply.hpp"

using namespace qpp;

// Unit testing "qpp/internal/kernels/qubit/apply_ctrl_test.hpp"

/// BEGIN template <typename Derived1, typename Derived2> void
/// apply_ctrl_psi_1q_inplace(...)
TEST(qpp_apply_ctrl_psi_1q_inplace, AllTests) {
    using namespace qpp;

    // ------------------------------------------------------------------------
    // Test 1: Standard CNOT (1 Control, 1 Target)
    // ------------------------------------------------------------------------
    {
        idx n = 2;                 // 2 qubits
        idx i = 1;                 // Target qubit 1 (LSB)
        std::vector<idx> ctrl{0};  // Control qubit 0 (MSB)
        std::vector<idx> shift{0}; // Trigger when control is |1>

        // State |10> -> Big-Endian index 2
        ket state = ket::Zero(4);
        state(2) = 1.0;

        // Gate A = X
        cmat A = gt.X;

        // Apply controlled gate. Expected outcome: |11> -> index 3
        internal::kernels::qubit::apply_ctrl_psi_1q_inplace(state, A, ctrl, i,
                                                            shift, n);

        EXPECT_NEAR(state(0).real(), 0.0, 1e-5);
        EXPECT_NEAR(state(1).real(), 0.0, 1e-5);
        EXPECT_NEAR(state(2).real(), 0.0, 1e-5);
        EXPECT_NEAR(state(3).real(), 1.0, 1e-5);
    }

    // ------------------------------------------------------------------------
    // Test 2: Control with Shift (Triggers when Control is |0>)
    // ------------------------------------------------------------------------
    {
        idx n = 2;
        idx i = 1;                 // Target qubit 1
        std::vector<idx> ctrl{0};  // Control qubit 0
        std::vector<idx> shift{1}; // shift=1 means trigger when control is |0>

        // State |00> -> Big-Endian index 0
        ket state = ket::Zero(4);
        state(0) = 1.0;

        cmat A = gt.X;

        // Apply shift-controlled X gate. Expected outcome: |01> -> index 1
        internal::kernels::qubit::apply_ctrl_psi_1q_inplace(state, A, ctrl, i,
                                                            shift, n);

        EXPECT_NEAR(state(0).real(), 0.0, 1e-5);
        EXPECT_NEAR(state(1).real(), 1.0, 1e-5);
        EXPECT_NEAR(state(2).real(), 0.0, 1e-5);
        EXPECT_NEAR(state(3).real(), 0.0, 1e-5);
    }

    // ------------------------------------------------------------------------
    // Test 3: Multiple Control Qubits on a Random Ket State
    // ------------------------------------------------------------------------
    {
        idx n = 3;                    // 3 qubits (size 8)
        idx i = 2;                    // Target qubit 2 (LSB)
        std::vector<idx> ctrl{0, 1};  // Control qubits 0 and 1
        std::vector<idx> shift{0, 0}; // Trigger when both controls are |11>

        ket state = randket(8);
        ket original_state = state;
        cmat A = randU(2); // Random 1-qubit unitary gate

        // Apply Toffoli-like custom unitary transformation
        internal::kernels::qubit::apply_ctrl_psi_1q_inplace(state, A, ctrl, i,
                                                            shift, n);

        ket expected_state = ket::Zero(8);

        // Big-Endian: Control condition |11x> matches indices 6 (|110>) and 7
        // (|111>) Indices 0-5 remain completely unchanged
        expected_state.segment(0, 6) = original_state.segment(0, 6);

        // Target subspace (indices 6 and 7) gets multiplied by the 1-qubit gate
        // A
        Eigen::Vector2cd target_subspace = original_state.segment(6, 2);
        expected_state.segment(6, 2) = A * target_subspace;

        // Verify the entire state vector matches within error bounds
        EXPECT_NEAR(norm(state - expected_state), 0.0, 1e-5);
    }
}
