// Used for testing, do not use it as an example
#include <iostream>

#include "qpp.h"

int main() {
    using namespace qpp;
    std::cout << "Testing\n";

    const idx bits = 70; // number of bits
    Bit_circuit b{bits};

    const idx trials = 20; // number of trials

    b.rand(); // randomize the vector
    // b.set();
    auto c = b;

    std::random_device rd;
    std::mt19937 gen{rd()};
    std::vector<std::vector<idx>> indices(trials);
    // generate the indices
    for (idx i = 0; i < trials; ++i) {
        std::vector<idx> v(bits);
        std::iota(v.begin(), v.end(), 0);
        std::shuffle(v.begin(), v.end(), gen);
        std::vector<idx> tof(v.data(), v.data() + 3);
        indices[i] = tof;
    }

    // apply the Toffoli gate to random places
    for (idx i = 0; i < trials; ++i) {
        std::cout << "First: ";
        for (auto&& elem : indices[i])
            std::cout << elem << " ";
        std::cout << '\n';
        b.TOF(indices[i]);
    }

    // do it again in reverse order
    for (idx i = trials; i-- > 0;) {
        std::cout << "Second: ";
        for (auto&& elem : indices[i])
            std::cout << elem << " ";
        std::cout << '\n';
        b.TOF(indices[i]);
    }

    std::cout << "Initial: " << b << '\n';
    std::cout << "Final:   " << c << '\n';
    std::cout << "Hamming weight: " << b.count() << '\n';

    std::cout << b.gate_count.NOT << " " << b.gate_count.X << " "
              << b.gate_count.TOF << '\n';

    std::cout << (b == c) << '\n';
    std::cout << (b != c) << '\n';

    Dynamic_bitset bb(9);
    bb.set(1).set(3).set(8);
    std::cout << bb << '\n';
    // std::cout << (32 & ~(!false << 5));

    std::cout << "here\n";
    std::cout << bb.to_string('o', 'X') << '\n';

    Dynamic_bitset vlad(20);
    std::cout << vlad << '\n';

    std::vector<unsigned int> vv(20);
    for (auto& elem : vv) {
        std::cout << elem;
    }
    std::cout << '\n';

    ket x = (10_ket + 01_ket) / std::sqrt(2);
    std::cout << disp(x) << '\n';

    bra y = (10_bra + 01_bra) / std::sqrt(2);
    std::cout << disp(x) << '\n';

    cmat z = 110_prj;
    std::cout << disp(z) << '\n';
}
