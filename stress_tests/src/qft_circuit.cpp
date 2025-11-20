#include <iostream>

#include <qpp/qpp.hpp>

int main(int argc, char** argv) {
    using namespace qpp;
    idx n = std::stoi(argv[1]);
    QCircuit qc(n, n);
    qc.QFT();

    QEngine qe{qc};
    Timer t;
    qe.execute();
    std::cout << "Took " << t.toc() << " seconds\n";
}
