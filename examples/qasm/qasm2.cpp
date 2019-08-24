// Another QASM example
// Source: ./examples/circuits/qasm/qasm2.cpp

#include <iostream>

#include "qpp.h"

int main(int argc, char** argv) {
    using namespace qpp;
    if (argc < 2) {
        std::cerr << "Input file not specified!\n";
        return EXIT_FAILURE;
    }

    // read the circuit
    auto qc = qasm::read_from_file(argv[1]);

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
