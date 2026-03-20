# pyqpp - Python 3 bindings for Quantum++

---

## Installation instructions

**pyqpp** shares the same dependencies as **Quantum++** and can be installed
via `pip`

```shell
pip install git+https://github.com/softwareQinc/qpp
```

---

## Overview

**pyqpp** is a modern Python 3 interface to the **Quantum++** C++ library. It
provides Python users with complete access to the full **Quantum++** API,
including quantum computing primitives, circuit construction and manipulation
tools, high-performance circuit simulation engines, and the full suite of
utility functions. The interface is implemented using
[pybind11](https://github.com/pybind/pybind11), offering a clean, intuitive,
and high-performance Pythonic experience while preserving the power and
efficiency of the underlying C++ implementation.

The bindings are easy to extend; see the
[Custom Bindings](#custom-bindings) section for guidance.

Illustrative examples are provided in [`examples`](examples).

---

## OpenQASM circuits

Use `pyqpp.qasm.read_from_file()` to obtain the `qpp::QCircuit` representation
of an OpenQASM 2.0 file.

---

## Custom bindings

To wrap a custom function, use `pybind11::module::def`, such as

```cpp
template<typename Func, typename ...Extra>
module &def(const char *name_, Func &&f, const Extra&... extra)
```

Here `Func` can be a plain C++ function, a function pointer, or a lambda
function. For example, consider the `qpp::randU()` function

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

## Template functions

We cannot wrap templated functions; instead, we must explicitly instantiate
them. For example, consider the `qpp::norm()` function

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

This creates the overloaded `pyqpp.norm()` function, which can accept
`qpp::ket` or `qpp::cmat` types. To avoid repetition of boilerplate code, we
can templatize the binding

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
