#include <iostream>

#include <qpp/qpp.h>

int main() {
    using namespace qpp;
    QCircuit qc{1, 1, 2, "coin flip"};
    qc.gate(gt.H, 0);
    qc.measureZ(0, 0);

    std::cout << qc << "\n\n" << qc.get_resources() << "\n\n";
    std::cout << QEngine{qc}.execute(100) << "\n";
}
