// Randomness
// Source: ./examples/randomness.cpp
#include <iostream>
#include <vector>

#include "qpp.h"

int main() {
    using namespace qpp;
    std::cout << ">> Generating a random ket on D = 5\n";
    ket rket = randket(5);
    std::cout << disp(rket) << '\n';

    std::vector<double> probs = abssq(rket);
    std::cout << ">> Probabilities: " << disp(probs, ", ") << '\n';

    std::cout << ">> Sum of the probabilities: ";
    std::cout << sum(probs.begin(), probs.end()) << '\n';
}
