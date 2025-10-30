// Source: ./examples/randomness.cpp
//
// Randomness

#include <iostream>
#include <vector>

#include "qpp/qpp.hpp"

int main() {
    using namespace qpp;

    std::cout << ">> Generating a random ket on D = 5\n";
    ket rket = randket(5);
    std::cout << disp(rket) << '\n';

    std::vector<realT> probs = abssq(rket);
    std::cout << ">> Probabilities: "
              << disp(probs, IOManipContainerOpts{}.set_sep(", ")) << '\n';

    std::cout << ">> Sum of the probabilities: ";
    std::cout << sum(probs.begin(), probs.end()) << '\n';
}
