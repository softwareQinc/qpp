// Samples from the output of an OpenQASM program read from the input stream,
// repeatedly if the number of repetitions is passed as the first argument. If
// there is more than one argument (i.e., argc > 2), then the corresponding
// arguments are interpreted as a list of qubits that will be sampled from.

// Note: make sure you remove (comment) the measurements that correspond to the
// final result from the OpenQASM program, as otherwise the sampling will
// produce trivial results (the result of the measurement repeated multiple
// times)

// Source: ./examples/qasm/qpp_qasm_sample.cpp
#include <iostream>

#include "qpp.h"

int main(int argc, char** argv) {
    using namespace qpp;

    // read the circuit from the input stream
    QCircuit qc = qasm::read(std::cin);

    // initialize the quantum engine with a circuit
    QEngine q_engine{qc};

    // display the quantum circuit and its corresponding resources
    std::cout << qc << "\n\n" << qc.get_resources() << "\n\n";

    // execute the quantum circuit
    idx reps = argc > 1 ? std::stoi(argv[1]) : 1; // repetitions

    std::vector<idx> target;
    if (argc > 2) {
        for (idx i = 2; i < static_cast<idx>(argc); ++i) {
            target.emplace_back(std::stoi(argv[i]));
        }
    }

    // execute the quantum circuit and sample from its output state repeatedly
    QEngine::Statistics stats = target.empty()
                                    ? q_engine.execute_sample(reps)
                                    : q_engine.execute_sample(target, reps);

    // display the statistics
    std::cout << stats << '\n';
}
