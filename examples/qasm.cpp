#include <iostream>

#include "qpp.h"

int main(int argc, char** argv) {
  using namespace qpp;
  if (argc < 2) {
    std::cerr << "Input file not specified.\n";
  }

  // Read the circuit
  auto qc = qasm::read_from_file(argv[1]);

  // initialize the quantum engine with a circuit
  QEngine engine{*qc};

  // display the quantum circuit
  std::cout << ">> BEGIN CIRCUIT\n";
  std::cout << engine.get_circuit() << '\n';
  std::cout << ">> END CIRCUIT\n\n";

  // execute the quantum circuit
  engine.execute();

  // display the measurement statistics
  std::cout << ">> BEGIN ENGINE STATISTICS\n";
  std::cout << engine << '\n';
  std::cout << ">> END ENGINE STATISTICS\n\n";

  // display the final state
  ket psi_final = engine.get_psi();
  std::cout << ">> Final state:\n";
  std::cout << disp(psi_final) << '\n';

  return 1;
}
