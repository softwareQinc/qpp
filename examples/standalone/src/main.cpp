#include <iostream>

#include <qpp/qpp.h>

int main() {
    using namespace qpp;
    std::cout << "The |0> state is:\n" << disp(0_ket) << '\n';
}
