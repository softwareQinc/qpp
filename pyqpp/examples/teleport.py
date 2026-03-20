import pyqpp.gates as gt
import pyqpp.states as st
from pyqpp import QCircuit, QEngine, dirac, norm, randU


def main() -> None:
    print("Qubit teleportation quantum circuit simulation\n")
    # Create a circuit with 3 qubits and 2 classical bits
    qc = QCircuit(3, 2)

    # Prepare qubit 0 in a random state
    random_unitary = randU(2)
    qc.gate(random_unitary, 0, "randU")

    # Create a Bell pair between qubits 1 and 2
    qc.gate(gt.H, 1)
    qc.CTRL(gt.X, 1, 2)

    # Perform a Bell-basis measurement on qubits 0 and 1
    qc.CTRL(gt.X, 0, 1)
    qc.gate(gt.H, 0)
    qc.measure([0, 1])

    # Apply classically controlled correction gates to qubit 2
    qc.cCTRL(gt.X, 1, 2)
    qc.cCTRL(gt.Z, 0, 2)

    # Create the execution engine
    qe = QEngine(qc)

    # Show the circuit and resource counts
    print(qc)
    print()
    print(qc.get_resources())
    print()

    # Run the circuit
    qe.execute()

    # Show execution / measurement results
    print(qe)
    print()

    # Check that teleportation succeeded
    input_state = random_unitary @ st.z0
    output_state = qe.get_state()

    print("Teleported state:")
    print(dirac(output_state))
    print("Norm difference:")
    print(norm(output_state - input_state))


if __name__ == "__main__":
    main()
