// Used for testing, do not use it as an example
#include <iostream>

#include "qpp.h"

int main() {
    /////////// testing ///////////
    using namespace qpp;

    QCircuit q_circuit{2, 2};
    q_circuit.nop();
    q_circuit.measureZ(0, 0);
    q_circuit.measureZ(1, 1);
    std::cout << q_circuit << "\n\n";

    QNoisyEngine<QubitBitFlipNoise> q_noisy_engine{q_circuit,
                                                   QubitBitFlipNoise{0.1}};
    q_noisy_engine.execute(1000);
    std::cout << q_noisy_engine;
}
