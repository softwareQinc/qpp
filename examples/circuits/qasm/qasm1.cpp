// Minimal QASM example
// Source: ./examples/circuits/qasm/qasm1.cpp
#include <iostream>

#include "qpp.h"

int main() {
    using namespace qpp;

    // create a qpp::QCircuit from a QASM file
    QCircuit qc = qasm::read_from_file(
        PATH "/examples/circuits/qasm/teleport_minimal.qasm");
    // note that QASM measurements are non-destructive, so the final state after
    // this step when executed on an engine will be a state of 3 qubits; that is
    // why we discard the first two qubits in the next line
    qc.discard({0, 1});

    std::cout << ">> BEGIN CIRCUIT\n";
    std::cout << qc << '\n'; // displays the circuit
    std::cout << ">> END CIRCUIT\n\n";

    QEngine q_engine{qc}; // create an engine out of a quantum circuit
    q_engine.execute();   // execute the circuit
    std::cout << ">> BEGIN ENGINE STATISTICS\n";
    std::cout << q_engine << '\n'; // display the measurement statistics
    std::cout << ">> END ENGINE STATISTICS\n\n";

    // displays the final output state
    std::cout << ">> Final state:\n";
    std::cout << disp(q_engine.get_psi()) << '\n';
}
