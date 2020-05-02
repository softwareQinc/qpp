// Used for testing, do not use it as an example
#include <iostream>

#include "qpp.h"
#include "experimental/experimental.h"

int main() {
    /////////// testing ///////////
    using namespace qpp;
    try {
    } catch (std::exception& e) {
        std::cout << e.what() << std::endl;
    }
}
