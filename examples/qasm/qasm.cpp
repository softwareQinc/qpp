// OpenQASM example, executes an OpenQASM circuit read from the input stream
// (repeatedly if the number of repetitions is passed as an argument)
// Source: ./examples/qasm/qasm.cpp
#include <iostream>

#include "qpp.h"

int main(int argc, char** argv) {
    using namespace qpp;

    QCircuit qc;
    // read the circuit from the input stream
    qc = qasm::read(std::cin);

    idx reps = 1;
    if (argc > 1)
        reps = std::stoi(argv[1]);

    // initialize the quantum engine with a circuit
    QEngine q_engine{qc};

    // display the quantum circuit and its corresponding resources
    std::cout << qc << "\n\n" << qc.get_resources() << "\n\n";

    // execute the quantum circuit
    q_engine.execute(reps);

    // display the measurement statistics
    std::cout << q_engine << '\n';

    // uncomment the following lines if you want to display the final state
    /*
    ket psi_final = q_engine.get_psi();
    std::cout << ">> Final state (transpose):\n"
              << disp(transpose(psi_final)) << '\n';
    */
}
