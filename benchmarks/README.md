# Benchmarking suite using [Catch2](https://github.com/catchorg/Catch2)

## Pre-requisites

- [Quantum++](https://github.com/softwareQinc/qpp), follow the [installation
  instructions](https://github.com/softwareQinc/qpp/blob/main/INSTALL.md)
- [CMake](https://cmake.org/)
- C++17 compliant compiler

## Set up and running

Execute

```shell
cmake -B build
cmake --build build --parallel 4
./build/bin/benchmark_executable_name # e.g., ./build/bin/qft
```

For help, execute

```shell
./build/bin/benchmark_executable_name --help # e.g., ./build/bin/qft --help
```
