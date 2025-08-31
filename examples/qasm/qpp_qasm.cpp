// Executes an OpenQASM program read from the input stream, repeatedly if the
// number of repetitions is passed as the first argument. If there is a second
// argument (i.e., argc == 3), then the final quantum state is displayed. If
// there are three or more arguments (i.e., argc  > 3), then the projector onto
// the final state is displayed.
// Source: ./examples/qasm/qpp_qasm.cpp

#include <cstdlib>
#include <iostream>

#include "qpp/qpp.hpp"

int main(int argc, char** argv) {
    using namespace qpp;

    // read the circuit from the input stream
    QCircuit qc = qasm::read(std::cin);

    // initialize the quantum engine with the circuit
    QEngine q_engine{qc};

    // display the quantum circuit and its corresponding resources
    std::cout << qc << "\n\n" << qc.get_resources() << "\n";

    // execute the quantum circuit
    idx reps = argc > 1 ? std::stoi(argv[1]) : 1; // repetitions
    // if no reps, exit
    if (reps == 0) {
        std::exit(EXIT_SUCCESS);
    }

    q_engine.execute(reps);

    // display the measurement statistics
    std::cout << '\n' << q_engine << '\n';

    // display the final state on demand
    if (argc == 3) {
        std::cout << "\n>> Final state:\n";
        std::cout << disp(dirac(q_engine.get_state())) << '\n';
    } else if (argc > 3) {
        std::cout << "\n>> Final density operator:\n";
        std::cout << disp(dirac(prj(q_engine.get_state()))) << '\n';
    }
}
