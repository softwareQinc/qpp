// Used for testing, do not use it as an example
#include <iostream>

#include "qpp.h"
#include "experimental/experimental.h"

int main() {
    using namespace qpp;
    idx d = 3;
    idx n = 4;
    idx D = std::pow(d, n);

    ket input = randket(D);
    ket result = QFT(input, d);
    std::cout << norm(gt.Fd(D) * input - result) << '\n';

    ket a = randket(D);
    ket b = randket(D);
    std::cout << norm(kron(a, b) - gt.SWAPd(D) * kron(b, a)) << "\n\n";

    ket psi = randket(D);
    idx k = 3;
    auto x = applyQFT(prj(psi), {3, 1, 2}, d);
    auto y = apply(prj(psi), gt.Fd(std::pow(d, k)), {3, 1, 2}, d);
//    std::cout << disp(y) << "\n\n";
    std::cout << norm(x - y) << "\n";
}
