// Used for testing, do not use it as an example
#include <iostream>

#include "qpp.h"

int main() {
    using namespace qpp;

    /////////// testing ///////////

    std::cout << "testing...\n";
    auto qubit_dep_noise = QubitDephasingNoise(0.3);
    std::cout << disp<cmat>(qubit_dep_noise) << '\n';

    ket psi = st.mes(2);

    cmat x = QubitAmplitudeDampingNoise(0.3, psi, 1);
    std::cout << disp(x) << "\n\n";

    cmat y = QubitPhaseDampingNoise(0.4, psi, 1);
    std::cout << disp(y) << "\n\n";

    //std::cout << disp(psi) << "\n\n" << norm(psi) << "\n\n";
}
