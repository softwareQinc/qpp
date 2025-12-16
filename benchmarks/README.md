# Benchmarking suite using [Catch2](https://github.com/catchorg/Catch2)

**Quantum++** includes a benchmarking suite built using
[Catch2](https://github.com/catchorg/Catch2).

---

## Configuring the system

**Quantum++ is required** to build and run the benchmarks. If it is not
installed, the benchmarks can be built by first configuring the main project
from this directory using

```shell
cmake -S .. -B build/qpp
cmake -B build -DCMAKE_PREFIX_PATH=build/qpp
```

If **Quantum++** is installed on your system, replace the commands above with

```shell
cmake -B build
```

---

## Building and running the benchmarks

Once the system is configured, build all benchmarks by executing from this
directory

```shell
cmake --build build --parallel 8
```

This will compile all benchmark executables into `./build/benchmarks`. To run
a specific benchmark, e.g., `qft_bench`, execute

```shell
./build/benchmarks/qft_bench
```

For help and additional options, run any benchmark with the `--help` flag.

> **Note:**
> All **Quantum++**
> [CMake flags and arguments](https://github.com/softwareQinc/qpp/blob/main/INSTALL.md#cmake-optional-flags-and-arguments)
> **do not propagate** to the benchmarks. They must be defined independently if
> they are needed.
