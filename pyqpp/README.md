pyqpp is a Python 3 wrapper for Quantum++.

# Installation

pyqpp can be installed using `pip`:

```
pip install git+https://github.com/softwareQinc/qpp
```

**Important**: If the installation fails due to your system being unable to
detect the location of the Eigen3 matrix library, set the environment variable
`EIGEN3_INSTALL_DIR` to point to the location of the Eigen3 library
(include the `include/eigen3` part of the path).

# Overview

pyqpp includes `Bit_circuit`, `Dynamic_bitset`, `QCircuit`, `QEngine`, and
several `QNoisyEngine` classes from Quantum++. Additionally, pyqpp provides
commonly used quantum `gates` and `states`, and some basic Eigen operations.
***
Example:

```python3
import numpy as np
from pyqpp import *

print("Qubit teleportation quantum circuit simulation\n")

# quantum circuit with 3 qubits and 2 classical bits
qc = QCircuit(3, 2)
# set the qubit 0 to a random state
U = randU(2)
# apply the gate U with name randU to qubit 0
qc.gate(U, 0, "randU")

# set the MES between qubits 1 and 2
qc.gate(gates.H, 1)
qc.CTRL(gates.X, 1, 2)

# perform the Bell measurement between qubits 0 and 1
qc.CTRL(gates.X, 0, 1)
qc.gate(gates.H, 0)
qc.measure([0, 1])

# apply the classical controls
qc.cCTRL(gates.X, 1, 2)
qc.cCTRL(gates.Z, 0, 2)

# initialize the quantum engine with a circuit
engine = QEngine(qc)

# display the quantum circuit and its corresponding resources
print(qc)
print()
print(qc.get_resources())
print()

# execute the entire circuit
engine.execute()

# display the measurement statistics
print(engine)
print()

# verify that the teleportation was successful
psi_in = np.matmul(U, states.z0)
psi_out = engine.get_psi()
print("Teleported state:\n", psi_out)
print("Norm difference:\n", norm(psi_out - psi_in))
```

## OpenQASM circuits

Use `pyqpp.qasm.read_from_file` to obtain the `QCircuit` representation of an
OpenQASM 2.0 file.

# Custom Bindings

pyqpp was created using pybind11 (
see [`pyqpp/qpp_wrapper.cpp`](https://github.com/softwareQinc/qpp/blob/main/pyqpp/qpp_wrapper.cpp))
. To wrap a custom function, use `pybind11::module::def`.

```C++
template<typename Func, typename ...Extra>
module &def(const char *name_, Func &&f, const Extra&... extra)
```

`Func` can be a plain C++ function, a function pointer, or a lambda function.
***
For example, consider the `qpp::randU` method

```C++
cmat randU(idx D = 2);
```

which is wrapped as

```C++
PYBIND11_MODULE(pyqpp, m) {
    ...

    m.def("randU", &qpp::randU, "Generates a random unitary matrix",
          py::arg("D") = 2);

    ...
}
```

## Template methods

We cannot wrap templated functions; instead, we must explicitly instantiate
them. For example, consider the `qpp::norm` method

```C++
template <typename Derived>
double norm(const Eigen::MatrixBase<Derived>& A);
```

One way to wrap this is

```C++
PYBIND11_MODULE(pyqpp, m) {
    ...

    m.def("norm", [](const cmat& A) { return qpp::norm(A); }, "Frobenius norm");
    m.def("norm", [](const ket& A) { return qpp::norm(A); }, "Frobenius norm");

    ...
}
```

This creates the overloaded `pyqpp.norm` function, which can accept `cmat`
or `ket` types. To avoid repetition of boilerplate code, we can templatize the
binding:

```C++
template<typename T>
void def_norm(pybind11::module &m) {
    m.def("norm", [](const T& A) { return qpp::norm(A); }, "Frobenius norm");
}

PYBIND11_MODULE(pyqpp, m) {
    ...

    def_norm<cmat>(m);
    def_norm<ket>(m);

    ...
}
```
