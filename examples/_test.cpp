// Used for testing, do not use it as an example
#include <iostream>

#include "qpp.h"

int main() {
    using namespace qpp;

    /////////// testing ///////////
    // we simulate the action fully depolarizing noise applied multiple times
    ket psi = randket(2);
    auto noise = QubitDepolarizingNoise(3. / 4);
    idx N = 10000;
    cmat result = cmat::Zero(2, 2);
    for (idx i = 0; i < N; ++i) {
        ket out = noise(psi, 0);
        result = result + prj(out);
    }
    result /= N;
    std::cout << disp(result) << '\n'; // we expect to be close to identity
}