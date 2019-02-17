// Used for testing, do not use it as an example
#include <iostream>

#include "qpp.h"

int main() {
    /////////// testing ///////////
    using namespace qpp;

    const idx n = 1;
    const idx d = 2;
    QCircuit qc{n, 1, d, "testing"};
    //    for (idx i = 0; i < n; ++i)
    //        qc.gate(randU(d), i, "randU_" + std::to_string(i));

    qc.gate(gt.Z, 0);
    qc.nop();
    qc.gate(gt.Z, 0);
    qc.nop();
    qc.gate(gt.Z, 0);
    qc.measureZ(0, 0);

    QNoisyEngine<QubitDepolarizingNoise> engine{qc,
                                                QubitDepolarizingNoise{0.75}};
    idx i = 0;
    ket psi_initial;
    engine.set_psi(1_ket);

    for (auto&& step : qc) {
        if (i++ == 0) {
            psi_initial = engine.get_psi();
        }
        engine.execute(step);
    }

    std::cout << "---- CIRCUIT ----\n"
              << qc << std::endl
              << "---- END CIRCUIT ----\n";
    std::cout << "---- ENGINE ----\n"
              << engine << std::endl
              << "---- END ENGINE ----\n";

    std::cout << engine.to_JSON() << '\n';
    std::cout << qc.to_JSON() << '\n';

    std::cout << qc.get_nop_count() << '\n';
}
