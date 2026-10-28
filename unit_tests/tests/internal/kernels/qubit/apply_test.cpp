#include "gmock/gmock.h" // for matchers
#include "gtest/gtest.h"

#include "qpp/qpp.hpp"

#include "qpp/internal/kernels/qubit/apply.hpp"

using namespace qpp;

// Unit testing "qpp/internal/kernels/qubit/apply_test.hpp"

/// BEGIN template <typename Derived1, typename Derived2> void
/// apply_psi_1q_inplace(Eigen::MatrixBase<Derived1>& state, const
/// Eigen::MatrixBase<Derived2>& A, idx i, idx n)
TEST(qpp_apply_psi_1q_inplace, AllTests) {
    using namespace qpp;

    // 2 qubits (size 4), applying Hadamard to qubit 0
    {
        idx n = 2;
        idx i = 0; // Target qubit 0

        // State |00> = [1, 0, 0, 0]^T
        ket state = ket::Zero(4);
        state(0) = 1.0;

        // Hadamard gate matrix
        cmat H(2, 2);
        H << 1.0 / std::sqrt(2), 1.0 / std::sqrt(2), 1.0 / std::sqrt(2),
            -1.0 / std::sqrt(2);

        internal::kernels::qubit::apply_psi_1q_inplace(state, H, i, n);

        // Expected output: 1/sqrt(2) * (|00> + |10>) -> [1/sqrt(2), 0,
        // 1/sqrt(2), 0]^T
        EXPECT_NEAR(state(0).real(), 1.0 / std::sqrt(2), 1e-5);
        EXPECT_NEAR(state(1).real(), 0.0, 1e-5);
        EXPECT_NEAR(state(2).real(), 1.0 / std::sqrt(2), 1e-5);
        EXPECT_NEAR(state(3).real(), 0.0, 1e-5);
    }

    // 2 qubits (size 4), applying Hadamard to qubit 1
    {
        idx n = 2;
        idx i = 1; // Target qubit 1

        // State |00> = [1, 0, 0, 0]^T
        ket state = ket::Zero(4);
        state(0) = 1.0;

        cmat H(2, 2);
        H << 1.0 / std::sqrt(2), 1.0 / std::sqrt(2), 1.0 / std::sqrt(2),
            -1.0 / std::sqrt(2);

        internal::kernels::qubit::apply_psi_1q_inplace(state, H, i, n);

        // Expected output: 1/sqrt(2) * (|00> + |01>) -> [1/sqrt(2), 1/sqrt(2),
        // 0, 0]^T
        EXPECT_NEAR(state(0).real(), 1.0 / std::sqrt(2), 1e-5);
        EXPECT_NEAR(state(1).real(), 1.0 / std::sqrt(2), 1e-5);
        EXPECT_NEAR(state(2).real(), 0.0, 1e-5);
        EXPECT_NEAR(state(3).real(), 0.0, 1e-5);
    }
}

/// BEGIN template <typename Derived1, typename Derived2> void
/// apply_psi_1q_inplace(Eigen::MatrixBase<Derived1>& state, const
/// Eigen::MatrixBase<Derived2>& A, idx i, idx n)
TEST(qpp_apply_psi_1q_inplace, TeleportationRandomKet) {
    using namespace qpp;

    // Total qubits = 3 (0: Alice's input, 1: Alice's half of EPR, 2: Bob's half
    // of EPR)
    idx n = 3;

    // 1. Generate a random single-qubit state to teleport: |\psi> = \alpha|0> +
    // \beta|1>
    ket psi = randket(2);

    // 2. Initialize the full 3-qubit system state to |000>
    ket state = ket::Zero(8);

    // Embed the random state into qubit 0: |\psi> \otimes |00>
    state(0) = psi(0); // |000> coefficient
    state(4) = psi(1); // |100> coefficient (Big-Endian: qubit 0 is the MSB)

    // Define standard gates
    cmat H(2, 2);
    H << 1.0 / std::sqrt(2), 1.0 / std::sqrt(2), 1.0 / std::sqrt(2),
        -1.0 / std::sqrt(2);

    cmat X(2, 2);
    X << 0.0, 1.0, 1.0, 0.0;

    cmat Z(2, 2);
    Z << 1.0, 0.0, 0.0, -1.0;

    // 3. Create EPR pair between qubit 1 and qubit 2
    // Apply H on qubit 1
    internal::kernels::qubit::apply_psi_1q_inplace(state, H, 1, n);
    // Apply CNOT from qubit 1 to qubit 2 manually via inplace 1-qubit gates on
    // subspaces (Equivalent to applying X on qubit 2 if qubit 1 is |1>)
    {
        // In Big-Endian (qubit 0, qubit 1, qubit 2), qubit 1 is set when
        // indices have bit 1 equal to 1. We can isolate the 2x2 subspace of
        // qubit 2 where qubit 1 is |1>
        for (idx q0 : {0, 4}) { // Qubit 0 can be 0 or 1
            idx base = q0 + 2;  // Qubit 1 is |1> -> binary xx1x -> +2
            // Extract the 2-element vector for qubit 2
            Eigen::Vector2cd sub_state;
            sub_state(0) = state(base + 0); // |q0 1 0>
            sub_state(1) = state(base + 1); // |q0 1 1>

            sub_state = X * sub_state;

            state(base + 0) = sub_state(0);
            state(base + 1) = sub_state(1);
        }
    }

    // 4. Alice applies CNOT from qubit 0 to qubit 1
    {
        // Apply X on qubit 1 if qubit 0 is |1>
        idx base = 4; // Qubit 0 is |1> -> binary 1xxx -> +4
        for (idx q2 : {0, 1}) {
            idx idx0 = base + q2;     // |1 0 q2>
            idx idx1 = base + 2 + q2; // |1 1 q2>
            Eigen::Vector2cd sub_state;
            sub_state(0) = state(idx0);
            sub_state(1) = state(idx1);

            sub_state = X * sub_state;

            state(idx0) = sub_state(0);
            state(idx1) = sub_state(1);
        }
    }

    // 5. Alice applies H on qubit 0
    internal::kernels::qubit::apply_psi_1q_inplace(state, H, 0, n);

    // 6. Projective measurement or classic controlled corrections.
    // Instead of sampling, we conditionally apply corrections to Bob's qubit
    // (2) based on the state of Alice's qubits (0 and 1) to reconstruct the
    // state.

    // If qubit 1 is |1>, apply X on qubit 2
    for (idx q0 : {0, 4}) {
        idx base = q0 + 2;
        Eigen::Vector2cd sub_state;
        sub_state(0) = state(base + 0);
        sub_state(1) = state(base + 1);
        sub_state = X * sub_state;
        state(base + 0) = sub_state(0);
        state(base + 1) = sub_state(1);
    }

    // If qubit 0 is |1>, apply Z on qubit 2
    {
        idx base = 4;
        for (idx q1 : {0, 2}) {
            Eigen::Vector2cd sub_state;
            sub_state(0) = state(base + q1 + 0);
            sub_state(1) = state(base + q1 + 1);
            sub_state = Z * sub_state;
            state(base + q1 + 0) = sub_state(0);
            state(base + q1 + 1) = sub_state(1);
        }
    }

    // 7. Verify that Bob's qubit (2) now holds the original random state |\psi>
    // The full state is now normalized to 1/2 across Alice's 4 possible
    // measurement outcomes. Each outcome branch |00>, |01>, |10>, |11> for
    // qubits 0 and 1 should match \psi on qubit 2.
    for (idx q0 : {0, 4}) {
        for (idx q1 : {0, 2}) {
            idx base = q0 + q1;
            // Bob's qubit amplitudes extracted from the branch
            cplx b0 = state(base + 0);
            cplx b1 = state(base + 1);

            // Re-normalize this specific branch (factor of 2 accounts for
            // 1/sqrt(4))
            b0 *= 2.0;
            b1 *= 2.0;

            EXPECT_NEAR(std::abs(b0), std::abs(psi(0)), 1e-5);
            EXPECT_NEAR(std::abs(b1), std::abs(psi(1)), 1e-5);
        }
    }
}

/// BEGIN template <typename Derived1, typename Derived2> void
/// apply_psi_1q_inplace(Eigen::MatrixBase<Derived1>& state, const
/// Eigen::MatrixBase<Derived2>& A, idx i, idx n)
TEST(qpp_apply_psi_1q_inplace, TeleportationStandardProtocol) {
    using namespace qpp;

    // Total qubits = 3 (0: Alice's input, 1: Alice's EPR half, 2: Bob's EPR
    // half)
    idx n = 3;

    // 1. Generate a random single-qubit state to teleport
    ket psi_a = randket();

    // 2. Setup the global input state: kron(psi_a, phi_AB) where phi_AB = |00>
    // + |11>
    ket input_aAB = ket::Zero(8);
    // |psi_a> kron (|00> + |11>) / sqrt(2)
    // Big-endian indexing for elements:
    // |000> = psi_a(0) / sqrt(2), |011> = psi_a(0) / sqrt(2)
    // |100> = psi_a(1) / sqrt(2), |111> = psi_a(1) / sqrt(2)
    double inv_sqrt2 = 1.0 / std::sqrt(2.0);
    input_aAB(0) = psi_a(0) * inv_sqrt2; // |000>
    input_aAB(3) = psi_a(0) * inv_sqrt2; // |011>
    input_aAB(4) = psi_a(1) * inv_sqrt2; // |100>
    input_aAB(7) = psi_a(1) * inv_sqrt2; // |111>

    // 3. Apply a CNOT on qubits 'AB' (control: 0, target: 1)
    // Instead of using applyCTRL_inplace, we isolate target qubit 1 where
    // control qubit 0 is |1> Qubit 0 is |1> at indices 4, 5, 6, 7. We apply
    // gt.X to qubit 1 inside this subspace.
    {
        // For indices 4-7, qubit 1 corresponds to bit 1 (from right, 0-indexed:
        // 2^1 = position 1) Subspace pairs for qubit 1: (4, 6) [|100> and
        // |110>] and (5, 7) [|101> and |111>]
        for (idx q2 : {0, 1}) {
            idx idx0 = 4 + q2;     // |1 0 q2>
            idx idx1 = 4 + 2 + q2; // |1 1 q2>

            cplx tmp = input_aAB(idx0);
            input_aAB(idx0) = input_aAB(idx1);
            input_aAB(idx1) = tmp;
        }
    }

    // 4. Apply an H on qubit 'a' (qubit 0) using the tested internal kernel
    internal::kernels::qubit::apply_psi_1q_inplace(input_aAB, gt.H, 0, n);

    // 5. Instead of random non-deterministic measurement via measure_seq,
    // extract and verify all 4 possible outcomes (z, x) for Alice's qubits (0,
    // 1)
    for (idx z = 0; z <= 1; ++z) {
        for (idx x = 0; x <= 1; ++x) {
            // Locate the 2-element subspace belonging to Bob (qubit 2)
            // given Alice's outcome state |z x>
            idx base_idx = (z * 4) + (x * 2);

            ket psi_B = ket::Zero(2);
            psi_B(0) = input_aAB(base_idx + 0); // |z x 0>
            psi_B(1) = input_aAB(base_idx + 1); // |z x 1>

            // Re-normalize Bob's conditional state (probability of each branch
            // is exactly 1/4)
            psi_B *= 2.0;

            // Perform the correction on Bob's qubit: psi_B = Z^z * X^x * psi_B
            psi_B = powm(gt.Z, z) * powm(gt.X, x) * psi_B;

            // Verification: Norm difference must be practically zero
            EXPECT_NEAR(norm(psi_B - psi_a), 0.0, 1e-5);
        }
    }
}
