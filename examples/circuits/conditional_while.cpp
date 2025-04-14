// Conditional WHILE quantum circuit simulator
// Source: ./examples/circuits/conditional_while.cpp

#include <iostream>

#include "qpp/qpp.hpp"

int main() {
    using namespace qpp;

    std::cout << ">> Conditional WHILE quantum circuit simulator\n\n";

    // quantum circuit with 2 qubits and 2 classical bits
    QCircuit qc{3, 3};

    // define a boolean predicate of the required form std::vector<idx> -> bool
    auto pred = [](std::vector<idx> dits) {
        // returns true as long as the first two classical dits ARE NOT 1, 1
        return !(dits[0] == 1 && dits[1] == 1);
    };

    // conditional WHILE statement
    // keep preparing the first measuring the first 2 qubits non-destructively
    // until we obtain the result 1, 1
    qc.cond_while(pred);
    // curly braces are optional, used to force code indenting
    {
        qc.reset({0, 1}); // resets the first two qubits to the |00> state
        qc.gate_fan(gt.H, {0, 1});    // next, prepare them in the |++> state
        qc.measure({0, 1}, 0, false); // finally, measure them non-destructively
    } // keep repeating until both measurement results are 1, 1
    qc.cond_end();

    // the WHILE statement finished, flip the state of the third qubit
    qc.gate(gt.X, 2);        // the final state will be |111>
    qc.measure(2, 2, false); // measure the third qubit non-destructively

    // display the quantum circuit and its corresponding resources
    std::cout << qc << "\n\n" << qc.get_resources() << "\n\n";

    // initialize the quantum engine with the circuit
    QEngine engine{qc};

    // execute the entire circuit once
    engine.execute();

    // display the measurement statistics
    std::cout << engine << "\n\n";
    // display the output quantum state
    std::cout << ">> Final state:\n";
    std::cout << disp(dirac(engine.get_state())) << '\n';
}
