// Exceptions
// Source: ./examples/exceptions.cpp
#include <exception>
#include <iostream>
#include "qpp.h"

using namespace qpp;

int main()
{
    cmat rho = randrho(16); // 4 qubits (subsystems)
    try
    {
        double mInfo = qmutualinfo(rho, {0}, {4}); // throws qpp::Exception
        std::cout << ">> Mutual information between first and last subsystem: ";
        std::cout << mInfo << '\n';
    }
    catch (const std::exception& e)
    {
        std::cout << ">> Exception caught: " << e.what() << '\n';
    }
}
