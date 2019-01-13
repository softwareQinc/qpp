// Used for testing, do not use it as an example
#include <iostream>

#include "qpp.h"
#include "experimental/experimental.h"

int main() {
    using namespace qpp;
    using namespace qpp::experimental;

    /*
        QCircuitDescription circ{10, 10};
        // circ.gate(gt.X, 7, "X gate on qubit 7");
        circ.measureV(gt.X, 4, 0, "Initial measurement on qubit 4");
        circ.gate(gt.X, 2, "X gate on qubit 2");
        circ.gate(gt.CNOT, 1, 3, "CNOT on 1 3");
        circ.measureZ(2, 0, "Z meas on qubit 2");
        circ.gate(gt.TOF, 5, 3, 6, "TOF on 5 3 6");
        circ.measureZ(6, 0);
        circ.measureZ(9, 0);
        circ.CTRL(gt.X, 6, 7, "Control X 6->7");
        circ.measureV(gt.H, 0, 0, "Measure in X basis");
        circ.gate(gt.Y, 8);
        circ.measureKs({std::sqrt(2) * gt.X, std::sqrt(2) * gt.Z}, 7, 0,
       "Kraus"); circ.measureZ(8, 0); circ.gate(gt.Z, {1, 2, 3}, "Z on 1 2 3");
        circ.measureZ(5, 9);
        circ.QFT({1, 2, 3});
        circ.TFQ({1, 2, 3});
        // circ.gate(gt.TOF, 4, 6, 7, "TOF 4 6 7");

        std::cout << circ << std::endl;

        QCircuitDescription test{10, 1, 2, "some circuit"};

        test.measureZ(1, 0);
        test.measureZ(3, 0);

        std::cout << test << '\n';
    */
    // TODO consider putting d after name for default fn arguments

    QCircuitDescription teleport_qubit{3, 2, 2, "qubit teleportation"};
    // set the qubit 0 to a random state
    cmat U = randU(2);
    teleport_qubit.gate(U, 0);

    // set the MES between qudits 1 and 2
    teleport_qubit.gate(gt.H, 1);
    teleport_qubit.gate(gt.CNOT, 1, 2);

    // perform the Bell measurement between qudits 0 and 1
    teleport_qubit.gate(gt.CNOT, 0, 1);
    teleport_qubit.gate(gt.H, 0);

    // perform the measurements
    teleport_qubit.measureZ(0, 0);
    teleport_qubit.measureZ(1, 1);

    // apply the classical controls
    teleport_qubit.cCTRL(gt.X, 1, 2);
    teleport_qubit.cCTRL(gt.Z, 0, 2);

    QCircuit qc_teleport_qubit{teleport_qubit};
    qc_teleport_qubit.run();
    std::cout << qc_teleport_qubit << '\n';

    ket psi_initial = U * 0_ket;
    ket psi_final = qc_teleport_qubit.get_psi();
    std::cout << disp(psi_final) << '\n';
    std::cout << "norm difference: " << norm(psi_final - psi_initial) << '\n';
}
