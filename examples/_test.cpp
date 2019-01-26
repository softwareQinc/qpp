// Used for testing, do not use it as an example
#include <iostream>

#include "qpp.h"

int main() {
    using namespace qpp;

    /////////// testing ///////////

    std::cout << "testing...\n";
    QCircuitDescription qcd{3, 3};
    qcd.gate_fan(gt.H).measureV(gt.Id(2), 0, 0);
    qcd.gate(gt.X, 1);
    qcd.gate(gt.X, 1);
    qcd.gate(gt.X, 1);
    qcd.gate(gt.X, 1);
    qcd.gate(gt.X, 1);
    qcd.gate(gt.X, 1);
    qcd.gate(gt.X, 1);
    qcd.measureZ(1, 1);
    qcd.measureV(gt.H, 2, 2);

    QCircuit qc{qcd};
    std::cout << ">> BEGIN CIRCUIT DESCRIPTION\n";
    std::cout << qc.get_circuit_description() << '\n';
    std::cout << ">> END CIRCUIT DESCRIPTION\n\n";

    std::cout << ">> BEGIN RUN\n";
    qc.reset();
    qc.run(true);
    std::cout << qc << '\n';
    std::cout << "psi:\n";
    std::cout << disp(qc.get_psi()) << '\n';
    std::cout << "m_ip_: " << qc.get_m_ip() << ", ";
    std::cout << "q_ip_: " << qc.get_q_ip() << "\n";
    std::cout << ">> END RUN\n";

    std::cout << ">> RESET\n";
    qc.reset();
    std::cout << "m_ip_: " << qc.get_m_ip() << ", ";
    std::cout << "q_ip_: " << qc.get_q_ip() << "\n";

    std::cout << ">> RUN 2\n";
    qc.run(true, 2);

    std::cout << ">> RESET\n";
    qc.reset();

    std::cout << ">> RUN 0\n";
    qc.run(true, 1);
    std::cout << ">> RUN 1\n";
    qc.run(true, 1);
    std::cout << ">> RUN 2\n";
    qc.run(true, 2);
    std::cout << ">> RUN 2, non-verbose\n";
    qc.run(false, 2);
    std::cout << ">> RUN UNTIL END, verbose\n";
    qc.run(true);
    std::cout << ">> END RUN\n";

    std::cout << qc << '\n';
    std::cout << "psi:\n";
    std::cout << disp(qc.get_psi()) << '\n';
    std::cout << "m_ip_: " << qc.get_m_ip() << ", ";
    std::cout << "q_ip_: " << qc.get_q_ip() << ", ";
    std::cout << "ip_: " << qc.get_ip() << "\n\n";

    // std::cout << *qc.get_iter(); // end, throws exception::InvalidIterator()

    std::cout << "measurement count: "
              << qc.get_circuit_description().get_measurement_count() << ", ";
    std::cout << "gate_count: " << qc.get_circuit_description().get_gate_count()
              << ", ";
    std::cout << "total count: "
              << qc.get_circuit_description().get_steps_count() << '\n';

    for (auto&& elem : qcd)
        std::cout << elem << '\n';

    QCircuitDescription a(2, 2);
    a.gate(gt.H, 0);
    std::cout << a << '\n';

    QCircuit q{a};
    q.run(true);
    std::cout << q << '\n';

    std::cout << std::boolalpha << is_iterable<QCircuitDescription>::value;
    std::cout << std::noboolalpha << '\n';
    std::cout << std::boolalpha << is_iterable<IQCircuit>::value;
    std::cout << std::noboolalpha << '\n';
    std::cout << std::boolalpha << is_iterable<QCircuit>::value;
    std::cout << std::noboolalpha << '\n';
}
