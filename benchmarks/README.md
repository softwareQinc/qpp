# Benchmarking suite using [Catch2](https://github.com/catchorg/Catch2)

Quantum++ includes a benchmarking suite built with
[Catch2](https://github.com/catchorg/Catch2).

To build all benchmarks, execute from this directory

```shell
cmake -B build && cmake --build build --parallel 8
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
