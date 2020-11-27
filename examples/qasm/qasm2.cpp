// Another QASM example
// Source: ./examples/qasm/qasm2.cpp

#include <iostream>

#include "qpp.h"

int main() {
    using namespace qpp;

    // read the circuit from the input stream
    auto qc = qasm::read(std::cin);

    // initialize the quantum engine with a circuit
    QEngine q_engine{qc};

    // display the quantum circuit
    std::cout << ">> BEGIN CIRCUIT\n";
    std::cout << q_engine.get_circuit() << '\n';
    std::cout << ">> END CIRCUIT\n\n";

    // execute the quantum circuit
    q_engine.execute();

    // display the measurement statistics
    std::cout << ">> BEGIN ENGINE STATISTICS\n";
    std::cout << q_engine << '\n';
    std::cout << ">> END ENGINE STATISTICS\n\n";

    // display the final state (its transpose to save console space)
    ket psi_final = q_engine.get_psi();
    std::cout << ">> Final state:\n";
    std::cout << disp(transpose(psi_final)) << '\n';
}
