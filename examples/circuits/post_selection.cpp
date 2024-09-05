// Qubit post-selection circuit simulator
// Source: ./examples/circuits/post_selection.cpp

#include <iostream>

#include "qpp/qpp.hpp"

int main() {
    using namespace qpp;

    std::cout << ">> Post-selection quantum circuit simulator\n\n";

    // quantum circuit with 2 qubits and 2 classical bits
    QCircuit qc{2, 2};
    // prepare equal superposition
    qc.gate_fan(gt.H);
    // post-select non-destructively all qubits on {1, 1}, and write the results
    // starting from the classical bit 0
    qc.post_select({0, 1}, {1, 1}, 0, false);
    // display the quantum circuit and its corresponding resources
    std::cout << qc << "\n\n" << qc.get_resources() << "\n\n";

    // initialize the quantum engine with the circuit
    QEngine engine{qc};

    // un-comment the line below to enforce post-selection, i.e., repeat
    // post-selection steps until success
    // engine.set_ensure_post_selection(true);

    // execute the entire circuit a few times
    engine.execute(1000); // we expect ~250 successful runs

    // display the measurement statistics
    std::cout << engine << "\n\n";
    // display the output quantum state
    std::cout << ">> Final state:\n";
    std::cout << disp(dirac(engine.get_state())) << '\n';
}
