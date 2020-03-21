// Used for testing, do not use it as an example
#include <iostream>

#include "qpp.h"
#include "experimental/experimental.h"

int main() {
    /////////// testing ///////////
    using namespace qpp;

    ket in = st.plus();
    qram data{0, 1};
    ket out = experimental::qRAM(in, data);
    std::cout << disp(out) << '\n';
}
