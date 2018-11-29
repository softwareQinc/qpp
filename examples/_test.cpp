// Used for testing, do not use it as an example
#include <iostream>

#include "qpp.h"
#include "experimental/experimental.h"

int main() {
    using namespace qpp;
//    idx d = 3;
//    idx n = 4;
//    idx D = std::pow(d, n);
//
//    ket input = randket(D);
//    ket result = QFT(input, d);
//    std::cout << norm(gt.Fd(D) * input - result) << '\n';
//
//    ket a = randket(D);
//    ket b = randket(D);
//    std::cout << norm(kron(a, b) - gt.SWAPd(D) * kron(b, a)) << "\n\n";
//
//    ket psi = randket(D);
//    idx k = 3;
//    auto x = applyQFT(prj(psi), {3, 1, 2}, d);
//    auto y = apply(prj(psi), gt.Fd(std::pow(d, k)), {3, 1, 2}, d);
////    std::cout << disp(y) << "\n\n";
//    std::cout << norm(x - y) << "\n";
//
//
//    ket xx = randket(D);
//    ket yy = applyQFT(xx,{2,0}, d);
//    ket zz = applyINVQFT(yy, {2,0}, d);
//
//    std::cout << norm(xx - zz) << "\n";

//    cmat test = gt.MODMUL(4, 7);
//    std::cout << disp(test) << "\n\n";
//    std::cout << disp(test*adjoint(test)) << "\n\n";

    ket psi = 00_ket;
    ket phi = st.one(2);
    std::cout << disp(phi) << std::endl << std::endl;
    std::cout << disp(01_ket) << std::endl << std::endl;
    //std::cout << norm(phi - 01_ket) << "\n";

    auto res = x2contfrac(3.1415927, 4);
    std::cout << disp(res, " ") << "\n";
    auto conv = convergents(res);
    for(auto&& elem: conv)
        std::cout << "(" << elem.first << ", " << elem.second << ")" << "\n";

    auto conv1 = convergents(x2contfrac(6.32, 10));
    for(auto&& elem: conv1)
        std::cout << "(" << elem.first << ", " << elem.second << ")" << "\n";
}
