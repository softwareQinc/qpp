# Python3 bindings for Quantum++

[**pyqpp**](https://github.com/softwareQinc/qpp/blob/main/pyqpp) provides a
modern Python 3 interface to the Quantum++ C++ library. It exposes the core
quantum computing primitives, high-performance circuit simulation engines, and
utility functions of Quantum++ to Python users through a clean and intuitive
API built with [pybind11](https://github.com/pybind/pybind11).

---

## Installation instructions

**pyqpp** requires the same dependencies as Quantum++,
and can be installed using `pip`

```shell
pip install git+https://github.com/softwareQinc/qpp
```

---

## Overview

**pyqpp** provides access to core Quantum++ components such as `Bit_circuit`,
`Dynamic_bitset`, `QCircuitT`, `QEngineT`, `QNoisyEngineT`, and several
additional engine sub-classes. It also provides common quantum gates,
predefined states, and basic Eigen-based operations.

Example:

```python3
import numpy as np
from pyqpp import *

print("Qubit teleportation quantum circuit simulation\n")

# Quantum circuit with 3 qubits and 2 classical bits
qc = QCircuit(3, 2)
# Set the qubit 0 to a random state
U = randU(2)
# Apply the gate U named randU to qubit 0
qc.gate(U, 0, "randU")

# Establish a maximally entangled state between qubits 1 and 2
qc.gate(gates.H, 1)
qc.CTRL(gates.X, 1, 2)

# Perform a Bell measurement between qubits 0 and 1
qc.CTRL(gates.X, 0, 1)
qc.gate(gates.H, 0)
qc.measure([0, 1])

# Apply the classical controls
qc.cCTRL(gates.X, 1, 2)
qc.cCTRL(gates.Z, 0, 2)

# Initialize the quantum engine with a circuit
engine = QEngine(qc)

# Display the quantum circuit and its corresponding resources
print(qc)
print()
print(qc.get_resources())
print()

# Execute the entire circuit
engine.execute()

# Display the measurement statistics
print(engine)
print()

# Verify that the teleportation was successful
psi_in = np.matmul(U, states.z0)
psi_out = engine.get_state()
print("Teleported state:")
print(dirac(psi_out))
print("Norm difference:\n", norm(psi_out - psi_in))
```

---

## OpenQASM circuits

Use `pyqpp.qasm.read_from_file` to obtain the `QCircuit` representation of an
OpenQASM 2.0 file.

---

## Custom Bindings

**pyqpp** was created using [pybind11](https://github.com/pybind/pybind11), see
[pyqpp/qpp_wrapper.cpp](https://github.com/softwareQinc/qpp/blob/main/pyqpp/qpp_wrapper.cpp).
To wrap a custom function, use `pybind11::module::def`.

```cpp
template<typename Func, typename ...Extra>
module &def(const char *name_, Func &&f, const Extra&... extra)
```

`Func` can be a plain C++ function, a function pointer, or a lambda function.
For example, consider the `qpp::randU` method

```cpp
cmat randU(idx D = 2);
```

which is wrapped as

```cpp
PYBIND11_MODULE(pyqpp, m) {
    ...

    m.def("randU", &qpp::randU, "Generates a random unitary matrix",
          py::arg("D") = 2);

    ...
}
```

---

## Template methods

We cannot wrap templated functions; instead, we must explicitly instantiate
them. For example, consider the `qpp::norm` method

```cpp
template <typename Derived>
double norm(const Eigen::MatrixBase<Derived>& A);
```

One way to wrap this is

```cpp
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

```cpp
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

---

## Creating Python stubs for IDE auto-completion and static type checking

In case autocompletion (or static type checking via
[mypy](https://www.mypy-lang.org/)) does not work properly in your editor/IDE,
you may need to create Python stubs for the package. To do this, execute

```shell
mkdir ~/python_stubs
export MYPATH=$MYPATH:~/python_subs # put this in your .profile or .bashrc
. ~/venv/bin/activate
stubgen -p pyqpp -o ~/python_stubs
ln -s ~/python_stubs/pyqpp ~/venv/lib/python3.11/site-packages
```

In the above, we assumed that your platform is UNIX/UNIX-like, and that you
have **pyqpp** installed in a virtual environment under `~/venv`. Please modify
accordingly for your system.
