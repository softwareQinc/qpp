#ifndef PATH
#define PATH ""
#endif

#include "gtest/gtest.h"
#include "qpp.h"

using namespace qpp;

// Unit testing "qasm.hpp"

/******************************************************************************/
/// BEGIN std::unique_ptr<qpp::QCircuit>
///       qpp::qasm::read_from_file(const std::string& fname)
/******************************************************************************/
TEST(qpp_qasm_read_from_file, StdCompliance) {
    // generic circuits
    EXPECT_NO_THROW(
        qasm::read_from_file(PATH "/tests/qasm/circuits/generic/adder.qasm"));
    EXPECT_NO_THROW(qasm::read_from_file(
        PATH "/tests/qasm/circuits/generic/bigadder.qasm"));
    EXPECT_NO_THROW(qasm::read_from_file(
        PATH "/tests/qasm/circuits/generic/inverseqft1.qasm"));
    EXPECT_NO_THROW(qasm::read_from_file(
        PATH "/tests/qasm/circuits/generic/inverseqft2.qasm"));
    EXPECT_NO_THROW(qasm::read_from_file(
        PATH "/tests/qasm/circuits/generic/ipea_3_pi_8.qasm"));
    EXPECT_NO_THROW(qasm::read_from_file(
        PATH "/tests/qasm/circuits/generic/pea_3_pi_8.qasm"));
    EXPECT_NO_THROW(
        qasm::read_from_file(PATH "/tests/qasm/circuits/generic/qec.qasm"));
    EXPECT_NO_THROW(
        qasm::read_from_file(PATH "/tests/qasm/circuits/generic/qft.qasm"));
    EXPECT_NO_THROW(
        qasm::read_from_file(PATH "/tests/qasm/circuits/generic/qpt.qasm"));
    EXPECT_NO_THROW(
        qasm::read_from_file(PATH "/tests/qasm/circuits/generic/rb.qasm"));
    EXPECT_NO_THROW(qasm::read_from_file(
        PATH "/tests/qasm/circuits/generic/teleport.qasm"));
    EXPECT_NO_THROW(qasm::read_from_file(
        PATH "/tests/qasm/circuits/generic/teleportv2.qasm"));
    EXPECT_NO_THROW(
        qasm::read_from_file(PATH "/tests/qasm/circuits/generic/W-state.qasm"));

    // ibmqx2 circuits
    EXPECT_NO_THROW(qasm::read_from_file(
        PATH "/tests/qasm/circuits/ibmqx2/011_3_qubit_grover_50_.qasm"));
    EXPECT_NO_THROW(qasm::read_from_file(
        PATH "/tests/qasm/circuits/ibmqx2/Deutsch_Algorithm.qasm"));
    EXPECT_NO_THROW(
        qasm::read_from_file(PATH "/tests/qasm/circuits/ibmqx2/iswap.qasm"));
    EXPECT_NO_THROW(
        qasm::read_from_file(PATH "/tests/qasm/circuits/ibmqx2/qe_qft_3.qasm"));
    EXPECT_NO_THROW(
        qasm::read_from_file(PATH "/tests/qasm/circuits/ibmqx2/qe_qft_4.qasm"));
    EXPECT_NO_THROW(
        qasm::read_from_file(PATH "/tests/qasm/circuits/ibmqx2/qe_qft_5.qasm"));
    EXPECT_NO_THROW(
        qasm::read_from_file(PATH "/tests/qasm/circuits/ibmqx2/W3test.qasm"));

    // invalid circuits
    EXPECT_THROW(qasm::read_from_file(
                     PATH "/tests/qasm/circuits/invalid/gate_no_found.qasm"),
                 exception::Undeclared);
    EXPECT_THROW(qasm::read_from_file(
                     PATH
                     "/tests/qasm/circuits/invalid/missing_semicolon.qasm"),
                 exception::ParseError);
}
/******************************************************************************/
TEST(qpp_qasm_read_from_file, BuiltinGates) {
    qpp::QCircuit q_circuit = qasm::read_from_file(
        PATH "/tests/qasm/circuits/units/builtingates.qasm");
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
    EXPECT_NEAR(1, norm(adjoint(psi1) * kron(psi2, m)), 1e-7);
    EXPECT_NEAR(1, norm(adjoint(psi1) * kron(psi3, m)), 1e-7);

    // Check CX gate
    EXPECT_EQ(1, c0);
    EXPECT_EQ(1, c1);
}
/******************************************************************************/
TEST(qpp_qasm_read_from_file, Teleportation) {
    qpp::QCircuit q_circuit = qasm::read_from_file(
        PATH "/tests/qasm/circuits/units/teleportation.qasm");
    QEngine engine{q_circuit};
    engine.execute();

    // Final state, should be |0> according to the QASM circuit
    auto rho = ptrace(engine.get_psi(), {0, 1}, {2, 2, 2});

    // Check norm
    EXPECT_NEAR(0, norm(rho - qpp::prj(0_ket)), 1e-7);
}
/******************************************************************************/
TEST(qpp_qasm_read_from_file, MappedGates) {
    qpp::QCircuit q_circuit = qasm::read_from_file(
        PATH "/tests/qasm/circuits/units/mappedgates.qasm");
    QEngine engine{q_circuit};
    engine.execute();

    // Final state
    ket psi1 = engine.get_psi();

    // Reference state
    ket psi2 = kron(gt.H * 0_ket, gt.H * 0_ket);

    // Check norm
    EXPECT_NEAR(0, norm(psi1 - psi2), 1e-7);
}
/******************************************************************************/
TEST(qpp_qasm_read_from_file, NonDestrMeas) {
    qpp::QCircuit q_circuit = qasm::read_from_file(
        PATH "/tests/qasm/circuits/units/nondestrmeas.qasm");
    QEngine engine{q_circuit};
    engine.execute();

    // Measurement result
    idx res = engine.get_dit(0);

    // Final state
    ket psi = engine.get_psi();

    // Check norm
    EXPECT_NEAR(0, norm(psi - mket({res})), 1e-7);
}
/******************************************************************************/
TEST(qpp_qasm_read_from_file, Reset) {
    qpp::QCircuit q_circuit =
        qasm::read_from_file(PATH "/tests/qasm/circuits/units/reset.qasm");
    QEngine engine{q_circuit};
    engine.execute();

    // Final state
    ket psi = engine.get_psi();

    // Check norm
    EXPECT_NEAR(0, norm(psi - 0_ket), 1e-7);
}
/******************************************************************************/
