// Used for testing, do not use it as an example
#include <iostream>

#include "qpp.h"

int main() {
    using namespace qpp;

    /////////// testing ///////////
    // we simulate the action fully depolarizing noise applied multiple times
    ket psi = 0.8 * 0_ket + 0.6 * 1_ket; // 0.8|0> + 0.6|1>
    auto noise = QubitDepolarizingNoise(3. / 4);
    idx N = 10000;
    cmat result = cmat::Zero(2, 2);
    for (idx i = 0; i < N; ++i) {
        ket out = noise(psi, 0);
        result = result + prj(out);
    }
    result = 1. / N * result;
    std::cout << disp(result) << '\n'; // we expect to be close to identity
}