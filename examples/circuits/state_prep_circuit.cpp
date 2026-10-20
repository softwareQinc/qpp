// Source: ./examples/circuits/state_prep_circuit.cpp
//
// Arbitrary state preparation quantum circuit simulator, see
// https://github.com/softwareQinc/qpp/issues/76

#include <cmath>
#include <iostream>

#include <qpp/qpp.hpp>

int main() {
    using namespace qpp;

    std::cout << ">> Arbitrary state preparation quantum circuit simulation\n";

    // 4-qubit state from https://github.com/softwareQinc/qpp/issues/76
    ket psi(16);
    psi << 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1;
    psi /= std::sqrt(static_cast<realT>(6));
    std::cout << ">> Target state:\n";
    std::cout << disp(dirac(psi)) << "\n\n";

    // quantum circuit that maps |0000> to psi
    QCircuit qc = state_prep_circuit(psi);

    // display the quantum circuit and its corresponding resources
    std::cout << qc << "\n\n" << qc.get_resources() << "\n\n";
    std::cout << ">> CNOT count: " << qc.get_gate_count(gt.CNOT) << "\n\n";

    std::cout << ">> Running...\n";
    QEngine engine{qc};
    engine.execute();
    ket result = engine.get_state();
    std::cout << ">> Prepared state:\n";
    std::cout << disp(dirac(result)) << '\n';
    std::cout << ">> Norm difference: " << norm(result - psi) << "\n\n";

    // a random 3-qutrit state, prepared with joint (block-diagonal) gates
    idx d = 3;
    ket phi = randket(27);
    QCircuit qc_qutrits = state_prep_circuit(phi, d);
    std::cout << qc_qutrits << "\n\n";
    ket result_qutrits = QEngine{qc_qutrits}.execute().get_state();
    std::cout << ">> Random 3-qutrit state, norm difference: "
              << norm(result_qutrits - phi) << '\n';
}
