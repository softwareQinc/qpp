// Source: ./examples/circuits/conditional_if.cpp
//
// Conditional IF quantum circuit simulator

#include <iostream>

#include "qpp/qpp.hpp"

int main() {
    using namespace qpp;

    std::cout << ">> Conditional IF quantum circuit simulator\n\n";

    // quantum circuit with 2 qubits and 2 classical bits
    QCircuit qc{3, 3};
    // prepare the first qubit in the |+> state
    qc.gate(gt.H, 0);
    // measure the first qubit non-destructively
    qc.measure(0, 0, false);

    // define a boolean predicate of the required form
    // qpp::internal::const_proxy_to_engine_dits_t -> bool
    auto pred = [](auto dits) {
        // returns true when the first dit is 1 at runtime (when run by a
        // quantum engine); in our case, this corresponds to the result of the
        // measurement result of the first qubit
        return dits[0] == 1;
    };
    // conditional IF statement
    // flips the second qubit when the predicate above was true, otherwise
    // flips the third qubit
    qc.cond_if(pred);
    // curly braces are optional, used to force code indenting
    {
        qc.gate(gt.X, 1); // the final state will be |110>
    }
    qc.cond_else();
    {
        qc.gate(gt.X, 2); // the final state will be |001>
    }
    qc.cond_end();

    // measure the second and third qubits non-destructively
    qc.measure(1, 1, false);
    qc.measure(2, 2, false);

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
