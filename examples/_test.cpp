// Used for testing, do not use it as an example
#include <iostream>
#include <fstream>
#include "qpp.h"
#include "experimental/experimental.h"

using namespace qpp;

int main()
{
    std::cout << "Testing...\n";

    experimental::Dynamic_bitset b(9);
    b.rand();
    std::cout << b << '\n';
}
