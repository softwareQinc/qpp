// Used for testing, do not use it as an example
#include <iostream>
#include "qpp.h"
#include "experimental/experimental.h"

using namespace qpp;

int main()
{
    std::cout << "Testing...\n";
    PRINTLN("Testing debug messages");
    ERRORLN("Testing debug error messages");

    std::cout << negativity(prj(st.mes(7)), {7, 7}); // (d - 1) / 2
}
