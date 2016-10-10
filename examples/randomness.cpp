// Randomness
// Source: ./examples/randomness.cpp
#include <iostream>
#include <vector>
#include "qpp.h"

using namespace qpp;

int main()
{
    std::cout << ">> Generating a random ket on D = 5" << std::endl;
    ket rket = randket(5);
    std::cout << disp(rket) << std::endl;

    std::vector<double> probs = abssq(rket);
    std::cout << ">> Probabilities: " << disp(probs, ", ") << std::endl;

    std::cout << ">> Sum of the probabilities: ";
    std::cout << sum(probs.begin(), probs.end()) << std::endl;
}
