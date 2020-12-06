// Used for testing, do not use it as an example
#include <iostream>

#include "qpp.h"
#include "experimental/experimental.hpp"

int main(int argc, char** argv) {
    /////////// testing ///////////
    using namespace qpp;

    assert(argc > 2);
    std::size_t reps = std::stoi(argv[1]);
    bool noisy = std::stoi(argv[2]);

    std::size_t nq = 5, nc = nq;
    QCircuit qCircuit{nq, nc};

    qCircuit.gate(gt.H, 0);
    qCircuit.CTRL(gt.X, 0, 1);
    for (std::size_t i = 0; i < nc; ++i)
        qCircuit.measureZ(i, i);
    std::cout << qCircuit << "\n\n";

    std::unique_ptr<QEngine> qEngine;
    if (noisy)
        qEngine.reset(new QNoisyEngine<QubitBitFlipNoise>{
            qCircuit, QubitBitFlipNoise{0.01}});
    else
        qEngine.reset(new QEngine{qCircuit});

    Timer<> t;
    qEngine->execute(reps);
    std::cout << *qEngine << '\n';
    std::cout << "Took:" << t.toc() << " seconds\n";
}
