// Standalone example, assumes quantum++ is installed in a system-wide visible
// directory

#include <iostream>

#include <qpp/qpp.hpp>

int main() {
    using namespace qpp;

    QCircuit qc{1, 1, 2, "coin flip"};
    qc.gate(gt.H, 0);
    qc.measure_all();
    std::cout << qc << "\n\n" << qc.get_resources() << "\n\n";

    QEngine qe{qc};
    std::cout << qe.execute(100) << '\n';
}
