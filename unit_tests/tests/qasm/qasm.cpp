#include <iostream>

#include "gtest/gtest.h"

#include "qpp/qpp.h"

using namespace qpp;

// Unit testing "qasm.hpp"


/// BEGIN std::unique_ptr<QCircuit> qasm::read_from_file(
///       const std::string& fname)

TEST(qpp_qasm_read_from_file, StdCompliance) {
    // generic circuits
    EXPECT_NO_THROW(qasm::read_from_file(PROJECT_ROOT_DIR
                                         "/qasmtools/qasm/generic/adder.qasm"));
    EXPECT_NO_THROW(qasm::read_from_file(
        PROJECT_ROOT_DIR "/qasmtools/qasm/generic/bigadder.qasm"));
    EXPECT_NO_THROW(qasm::read_from_file(
        PROJECT_ROOT_DIR "/qasmtools/qasm/generic/inverseqft1.qasm"));
    EXPECT_NO_THROW(qasm::read_from_file(
        PROJECT_ROOT_DIR "/qasmtools/qasm/generic/inverseqft2.qasm"));
    EXPECT_NO_THROW(qasm::read_from_file(
        PROJECT_ROOT_DIR "/qasmtools/qasm/generic/ipea_3_pi_8.qasm"));
    EXPECT_NO_THROW(qasm::read_from_file(
        PROJECT_ROOT_DIR "/qasmtools/qasm/generic/pea_3_pi_8.qasm"));
    EXPECT_NO_THROW(qasm::read_from_file(PROJECT_ROOT_DIR
                                         "/qasmtools/qasm/generic/qec.qasm"));
    EXPECT_NO_THROW(qasm::read_from_file(PROJECT_ROOT_DIR
                                         "/qasmtools/qasm/generic/qft.qasm"));
    EXPECT_NO_THROW(qasm::read_from_file(PROJECT_ROOT_DIR
                                         "/qasmtools/qasm/generic/qpt.qasm"));
    EXPECT_NO_THROW(qasm::read_from_file(PROJECT_ROOT_DIR
                                         "/qasmtools/qasm/generic/rb.qasm"));
    EXPECT_NO_THROW(qasm::read_from_file(
        PROJECT_ROOT_DIR "/qasmtools/qasm/generic/teleport.qasm"));
    EXPECT_NO_THROW(qasm::read_from_file(
        PROJECT_ROOT_DIR "/qasmtools/qasm/generic/teleportv2.qasm"));
    EXPECT_NO_THROW(qasm::read_from_file(
        PROJECT_ROOT_DIR "/qasmtools/qasm/generic/W-state.qasm"));

    // ibmqx2 circuits
    EXPECT_NO_THROW(qasm::read_from_file(
        PROJECT_ROOT_DIR "/qasmtools/qasm/ibmqx2/011_3_qubit_grover_50_.qasm"));
    EXPECT_NO_THROW(qasm::read_from_file(
        PROJECT_ROOT_DIR "/qasmtools/qasm/ibmqx2/Deutsch_Algorithm.qasm"));
    EXPECT_NO_THROW(qasm::read_from_file(PROJECT_ROOT_DIR
                                         "/qasmtools/qasm/ibmqx2/iswap.qasm"));
    EXPECT_NO_THROW(qasm::read_from_file(
        PROJECT_ROOT_DIR "/qasmtools/qasm/ibmqx2/qe_qft_3.qasm"));
    EXPECT_NO_THROW(qasm::read_from_file(
        PROJECT_ROOT_DIR "/qasmtools/qasm/ibmqx2/qe_qft_4.qasm"));
    EXPECT_NO_THROW(qasm::read_from_file(
        PROJECT_ROOT_DIR "/qasmtools/qasm/ibmqx2/qe_qft_5.qasm"));
    EXPECT_NO_THROW(qasm::read_from_file(PROJECT_ROOT_DIR
                                         "/qasmtools/qasm/ibmqx2/W3test.qasm"));
}

TEST(qpp_qasm_read_from_file, BuiltinGates) {
    QCircuit q_circuit = qasm::read_from_file(
        PROJECT_ROOT_DIR
        "/unit_tests/tests/qasm/circuits/units/builtingates.qasm");
    QEngine engine{q_circuit};
    engine.execute();

    // Final state
    ket psi1 = engine.get_psi();
    idx c0 = engine.get_dit(0);
    idx c1 = engine.get_dit(1);
    ket m = kron(mket({c0}), mket({c1}));

    // Check U gate against the equational definitions
    ket psi2 = gt.RZ(0.2) * gt.RY(0.1) * gt.RZ(0.3) * 0_ket;
    ket psi3 = gt.RZ(0.2 + 3 * pi) * gt.RX(pi / 2) * gt.RZ(0.1 + pi) *
               gt.RX(pi / 2) * gt.RZ(0.3) * 0_ket;
    EXPECT_NEAR(1, norm(adjoint(psi1) * kron(psi2, m)), 1e-5);
    EXPECT_NEAR(1, norm(adjoint(psi1) * kron(psi3, m)), 1e-5);

    // Check CX gate
    EXPECT_EQ(1, c0);
    EXPECT_EQ(1, c1);
}

TEST(qpp_qasm_read_from_file, Teleportation) {
    QCircuit q_circuit = qasm::read_from_file(
        PROJECT_ROOT_DIR
        "/unit_tests/tests/qasm/circuits/units/teleportation.qasm");
    QEngine engine{q_circuit};
    engine.execute();

    // Final state, should be |0> according to the QASM circuit
    auto rho = ptrace(engine.get_psi(), {0, 1}, {2, 2, 2});

    // Check norm
    EXPECT_NEAR(0, norm(rho - prj(0_ket)), 1e-5);
}

TEST(qpp_qasm_read_from_file, MappedGates) {
    QCircuit q_circuit = qasm::read_from_file(
        PROJECT_ROOT_DIR
        "/unit_tests/tests/qasm/circuits/units/mappedgates.qasm");
    QEngine engine{q_circuit};
    engine.execute();

    // Final state
    ket psi1 = engine.get_psi();

    // Reference state
    ket psi2 = QASMTOOLS_QASM2_SPECS
                   ? kron(gt.H * 0_ket * (-1_i), gt.H * 0_ket * (-1_i))
                   : kron(gt.H * 0_ket, gt.H * 0_ket);

    // Check norm
    EXPECT_NEAR(0, norm(psi1 - psi2), 1e-5);
}

TEST(qpp_qasm_read_from_file, NonDestrMeas) {
    QCircuit q_circuit = qasm::read_from_file(
        PROJECT_ROOT_DIR
        "/unit_tests/tests/qasm/circuits/units/nondestrmeas.qasm");
    QEngine engine{q_circuit};
    engine.execute();

    // Measurement result
    idx res = engine.get_dit(0);

    // Final state
    ket psi = engine.get_psi();

    // Reference state
    ket mres = QASMTOOLS_QASM2_SPECS ? mket({res}) * (-1_i) : mket({res});

    // Check norm
    EXPECT_NEAR(0, norm(psi - mres), 1e-5);
}

TEST(qpp_qasm_read_from_file, Reset) {
    QCircuit q_circuit = qasm::read_from_file(
        PROJECT_ROOT_DIR "/unit_tests/tests/qasm/circuits/units/reset.qasm");
    QEngine engine{q_circuit};
    engine.execute();

    // Final state
    ket psi = engine.get_psi();

    // Reference state
    ket psi2 = QASMTOOLS_QASM2_SPECS ? 0_ket * (-1_i) : 0_ket;

    // Check norm
    EXPECT_NEAR(0, norm(psi - psi2), 1e-5);
}

TEST(qpp_qasm_read_from_file, SciNot) {
    QCircuit q_circuit = qasm::read_from_file(
        PROJECT_ROOT_DIR "/unit_tests/tests/qasm/circuits/units/scinot.qasm");
    QEngine engine{q_circuit};
    engine.execute();

    // Final state
    ket psi = engine.get_psi();

    // Check against C++ scientific notation
    ket phi = gt.H * gt.RZ(1.0E-3) * gt.H * 0_ket;

    // Check norm
    EXPECT_NEAR(1, norm(adjoint(psi) * phi), 1e-5);
}

