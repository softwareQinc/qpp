// Used for testing, do not use it as an example
#include <iostream>

#include "qpp.h"
#include "experimental/experimental.h"

int main() {
    /////////// testing ///////////
    using namespace qpp;

    ket in = st.plus();
    qram data{1, 2, 3};
    ket out = qRAM(in, data);
    std::cout << disp(out) << '\n';
}
