// Used for testing, do not use it as an example
#include <iostream>

#include "qpp.h"

int main() {
    using namespace qpp;

    /////////// testing ///////////

    std::cout << "testing...\n";

    const Noise<StateDependentNoise>& noise = QubitPhaseDampingNoise(0.4);

    ket psi = 0.8 * 0_ket + 0.6 * 1_ket;

    cmat psi_out = noise(psi, 0);
    std::cout << disp(psi_out) << "\n\n";

    std::cout << disp(noise.get_probs(), " ") << "\n";
    std::cout << noise.get_last_idx() << "\n";
    std::cout << noise.get_last_p() << "\n";
    std::cout << disp(noise.get_last_K()) << "\n";
    std::cout << norm(psi_out) << "\n";
}
