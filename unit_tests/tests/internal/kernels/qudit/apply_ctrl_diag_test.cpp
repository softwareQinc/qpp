#include <cmath>
#include <vector>

#include "gtest/gtest.h"

#include "qpp/qpp.hpp"

#include "qpp/internal/kernels/qudit/apply_ctrl_diag.hpp"

using namespace qpp;

TEST(QuditCtrlDiagTest, InPlaceApplication) {
    using namespace qpp::internal::kernels::qudit;

    // 1. Setup system parameters
    const idx d = 3; // Qudit dimension (qutrits)
    const idx n = 3; // Total number of qudits

    // Total state vector dimension: 3^3 = 27
    idx D = 1;
    for (idx i = 0; i < n; ++i) {
        D *= d;
    }

    // Initialize state vector to a uniform superposition
    // |psi> = 1/sqrt(27) * sum_{i=0}^{26} |i>
    ket state = ket::Constant(D, 1.0 / std::sqrt(static_cast<double>(D)));
    ket original_state = state; // Keep a backup to verify changes

    // 2. Define the gate topology
    // Let qudit 0 and qudit 1 be controls, qudit 2 be the target
    std::vector<idx> ctrl = {0, 1};
    std::vector<idx> target = {2};

    // Shifts: shift[0] = 1, shift[1] = 0
    // Effective control value equation: c = ((q0 + 1) % 3) * 3 + ((q1 + 0) % 3)
    std::vector<idx> shift = {1, 0};

    // 3. Define the diagonal gate A for the target (dimension d^k = 3^1 = 3)
    // Let's use phase shifts: A = diag(1, e^{i*pi/4}, e^{i*pi/2})
    using namespace std::complex_literals;
    Eigen::VectorXcd A(d);
    A(0) = 1.0;
    A(1) = std::exp(1i * pi / 4.0);
    A(2) = std::exp(1i * pi / 2.0);

    // 4. Call your function
    apply_ctrl_psi_kq_diag_inplace(state, A, ctrl, target, shift, d, n);

    // 5. Verify the results analytically
    // We iterate over every basis state |q0 q1 q2>
    for (idx q0 = 0; q0 < d; ++q0) {
        for (idx q1 = 0; q1 < d; ++q1) {
            for (idx q2 = 0; q2 < d; ++q2) {

                // Map the multi-qudit coordinates to the global vector index
                // Index structure: q0 * 3^2 + q1 * 3^1 + q2 * 3^0
                idx global_idx = (q0 * 9) + (q1 * 3) + q2;

                // Calculate the modular shift impact on the control power 'c'
                idx eff_q0 = (q0 + shift[0]) % d;
                idx eff_q1 = (q1 + shift[1]) % d;
                idx joint_ctrl_power = (eff_q0 * d) + eff_q1;

                // Target state element phase modification: A(q2)^c
                std::complex<double> expected_phase =
                    std::pow(A(q2), joint_ctrl_power);
                std::complex<double> expected_value =
                    original_state(global_idx) * expected_phase;

                // Test with a floating-point tolerance threshold
                EXPECT_NEAR(state(global_idx).real(), expected_value.real(),
                            1e-12)
                    << "Mismatch at basis state |" << q0 << q1 << q2
                    << "> (Real part)";
                EXPECT_NEAR(state(global_idx).imag(), expected_value.imag(),
                            1e-12)
                    << "Mismatch at basis state |" << q0 << q1 << q2
                    << "> (Imag part)";
            }
        }
    }
}

// Concrete sanity check on a specific shifted configuration
TEST(QuditCtrlDiagTest, SpecificStateVerification) {
    using namespace qpp::internal::kernels::qudit;

    const idx d = 3;
    const idx n = 3;
    ket state = ket::Zero(27);

    // Set a pure computational basis state |0 2 1>
    // index = 0 * 9 + 2 * 3 + 1 * 1 = 7
    state(7) = 1.0;

    std::vector<idx> ctrl = {0, 1};
    std::vector<idx> target = {2};
    std::vector<idx> shift = {1, 0}; // Shifts: q0 shifts by +1, q1 shifts by +0

    using namespace std::complex_literals;
    Eigen::VectorXcd A(d);
    A(0) = 1.0;
    A(1) = 2.5; // Using real scale factors for easier debugging trace visual
                // inspection
    A(2) = 4.0;

    apply_ctrl_psi_kq_diag_inplace(state, A, ctrl, target, shift, d, n);

    // Analytical expectation:
    // q0 = 0 -> effective q0 = (0 + 1) % 3 = 1
    // q1 = 2 -> effective q1 = (2 + 0) % 3 = 2
    // Joint control power 'c' = 1 * 3 + 2 = 5
    // Target state value: q2 = 1 -> A(1) = 2.5
    // Expected scaling factor = (2.5)^5 = 97.65625

    double expected_amplitude = std::pow(2.5, 5); // 97.65625
    EXPECT_NEAR(state(7).real(), expected_amplitude, 1e-9);
    EXPECT_NEAR(state(7).imag(), 0.0, 1e-9);

    // Verify all other elements remain untouched (zero)
    for (idx i = 0; i < 27; ++i) {
        if (i != 7) {
            EXPECT_DOUBLE_EQ(state(i).real(), 0.0);
            EXPECT_DOUBLE_EQ(state(i).imag(), 0.0);
        }
    }
}

TEST(QuditCtrlDensityMatrixTest, SuperpositionDensityMatrix) {
    using namespace qpp::internal::kernels::qudit;

    // 1. Setup system parameters
    const idx d = 3; // Qudit dimension (qutrits)
    const idx n = 3; // Total number of qudits

    // Total state vector dimension: 3^3 = 27 -> Density matrix: 27 x 27
    idx D = 1;
    for (idx i = 0; i < n; ++i) {
        D *= d;
    }

    // Initialize state vector to a uniform superposition: |psi> = 1/sqrt(27) *
    // sum |i>
    ket psi = ket::Constant(D, 1.0 / std::sqrt(static_cast<double>(D)));

    // Construct the initial density matrix: rho = |psi><psi|
    cmat rho = psi * psi.adjoint();
    cmat original_rho = rho; // Backup for verification

    // 2. Define the gate topology
    std::vector<idx> ctrl = {0, 1}; // Qudit 0 is MSD, Qudit 1 is LSD
    std::vector<idx> target = {2};
    std::vector<idx> shift = {1, 0}; // Shift Q0 by +1, Q1 by +0

    // 3. Define the diagonal gate A for the target
    using namespace std::complex_literals;
    Eigen::VectorXcd A(d);
    A(0) = 1.0;
    A(1) = std::exp(1i * pi / 4.0);
    A(2) = std::exp(1i * pi / 2.0);

    // 4. Call your density matrix function
    apply_ctrl_rho_kq_diag_inplace(rho, A, ctrl, target, shift, d, n);

    // 5. Verify row-by-row and column-by-column
    for (idx r0 = 0; r0 < d; ++r0) {
        for (idx r1 = 0; r1 < d; ++r1) {
            for (idx r2 = 0; r2 < d; ++r2) {
                idx row_idx = (r0 * 9) + (r1 * 3) + r2;

                // Calculate control power for the row
                idx eff_r0 = (r0 + shift[0]) % d;
                idx eff_r1 = (r1 + shift[1]) % d;
                idx row_ctrl_power = (eff_r0 * d) + eff_r1;
                std::complex<double> row_phase =
                    std::pow(A(r2), row_ctrl_power);

                for (idx c0 = 0; c0 < d; ++c0) {
                    for (idx c1 = 0; c1 < d; ++c1) {
                        for (idx c2 = 0; c2 < d; ++c2) {
                            idx col_idx = (c0 * 9) + (c1 * 3) + c2;

                            // Calculate control power for the column
                            idx eff_c0 = (c0 + shift[0]) % d;
                            idx eff_c1 = (c1 + shift[1]) % d;
                            idx col_ctrl_power = (eff_c0 * d) + eff_c1;
                            // Column gets complex conjugate (U^\dagger)
                            std::complex<double> col_phase_conj =
                                std::conj(std::pow(A(c2), col_ctrl_power));

                            // Analytical expectation: rho(i,j) * U_row *
                            // U_col^*
                            std::complex<double> expected_value =
                                original_rho(row_idx, col_idx) * row_phase *
                                col_phase_conj;

                            EXPECT_NEAR(rho(row_idx, col_idx).real(),
                                        expected_value.real(), 1e-12)
                                << "Mismatch at rho(|" << r0 << r1 << r2 << "><"
                                << c0 << c1 << c2 << "|) [Real]";
                            EXPECT_NEAR(rho(row_idx, col_idx).imag(),
                                        expected_value.imag(), 1e-12)
                                << "Mismatch at rho(|" << r0 << r1 << r2 << "><"
                                << c0 << c1 << c2 << "|) [Imag]";
                        }
                    }
                }
            }
        }
    }
}

TEST(QuditCtrlDensityMatrixTest, MixedStateSpecificElementVerification) {
    using namespace qpp::internal::kernels::qudit;

    const idx d = 3;
    const idx n = 3;
    cmat rho = cmat::Zero(27, 27);

    // Construct a specific off-diagonal mixed element coherence block
    // Row state: |0 2 1> -> index = 0*9 + 2*3 + 1 = 7
    // Col state: |1 0 2> -> index = 1*9 + 0*3 + 2 = 11
    idx r_idx = 7;
    idx c_idx = 11;
    rho(r_idx, c_idx) = std::complex<double>(0.5, 0.5);

    std::vector<idx> ctrl = {0, 1};
    std::vector<idx> target = {2};
    std::vector<idx> shift = {1, 0};

    using namespace std::complex_literals;
    Eigen::VectorXcd A(d);
    A(0) = 1.0;
    A(1) = 2.0; // Scaled elements for transparent power validation
    A(2) = 3.0;

    apply_ctrl_rho_kq_diag_inplace(rho, A, ctrl, target, shift, d, n);

    // Analytical expectation trace math:
    // Row element |0 2 1>:
    //   r0 = 0 -> eff_r0 = (0 + 1) % 3 = 1
    //   r1 = 2 -> eff_r1 = (2 + 0) % 3 = 2
    //   row_power = 1 * 3 + 2 = 5
    //   r2 = 1 -> A(1) = 2.0 -> U_row = 2^5 = 32.0
    //
    // Col element |1 0 2>:
    //   c0 = 1 -> eff_c0 = (1 + 1) % 3 = 2
    //   c1 = 0 -> eff_c1 = (0 + 0) % 3 = 0
    //   col_power = 2 * 3 + 0 = 6
    //   c2 = 2 -> A(2) = 3.0 -> U_col = 3^6 = 729.0
    //
    // Total Scaling Factor = U_row * conj(U_col) = 32.0 * 729.0 = 23328.0

    std::complex<double> expected_val =
        std::complex<double>(0.5, 0.5) * 23328.0;

    EXPECT_NEAR(rho(r_idx, c_idx).real(), expected_val.real(), 1e-9);
    EXPECT_NEAR(rho(r_idx, c_idx).imag(), expected_val.imag(), 1e-9);

    // Ensure unaffected components remain completely zero
    rho(r_idx, c_idx) = 0.0;
    EXPECT_TRUE(rho.isZero(1e-15));
}
