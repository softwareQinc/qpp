// Used for testing, do not use it as an example
#include <iostream>

#include "qpp.h"

int main() {
    /////////// testing ///////////
    using namespace qpp;

    const idx n = 3;
    const idx d = 2;
    QCircuit qc{n, 0, d, "testing QFT/TFQ"};
    for (idx i = 0; i < n; ++i)
        qc.gate(randU(d), i, "randU_" + std::to_string(i));
    qc.QFT();
    qc.TFQ();

    QNoisyEngine<QubitAmplitudeDampingNoise> engine{
        qc, QubitAmplitudeDampingNoise{0.99}};
    idx i = 0;
    ket psi_initial;
    for (auto&& step : qc) {
        engine.execute(step);
        if (++i == n) {
            psi_initial = engine.get_psi();
        }
    }

    std::cout << "---- CIRCUIT ----\n"
              << qc << std::endl
              << "---- END CIRCUIT ----\n";
    std::cout << "---- ENGINE ----\n"
              << engine << std::endl
              << "---- END ENGINE ----\n";
    std::cout << "Norm difference: " << norm(psi_initial - engine.get_psi());
}
