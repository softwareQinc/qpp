// Source: ./examples/circuits/runtime_set_dits.cpp
//
// Overwriting quantum engine dits at runtime quantum circuit simulator

#include <iostream>

#include "qpp/qpp.hpp"

int main() {
    using namespace qpp;

    std::cout << ">> Overwriting quantum engine dits at runtime quantum "
                 "circuit simulator\n\n";

    // quantum circuit with 2 qubits and 2 classical bits
    QCircuit qc{1, 1};
    // prepare the first qubit in the |+> state
    qc.gate(gt.H, 0);
    // measure the first qubit non-destructively
    qc.measure(0, 0, false);

    // if the measurement result is 1
    qc.cond_if([](auto dits) { return dits[0] == 1; });
    // curly braces are optional, used to force code indenting
    {
        // define a functor of the required form proxy_to_engine_dits_t -> void
        auto functor = [](auto dits) { dits[0] = 100; };

        // runtime set dit statement
        // set the corresponding dit to 100, why not?
        qc.set_dits_runtime(functor);
    }
    qc.cond_end();

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
