# Benchmarking suite using [Catch2](https://github.com/catchorg/Catch2)

## Pre-requisites

- [Quantum++](https://github.com/softwareQinc/qpp), follow the [installation
  instructions](https://github.com/softwareQinc/qpp/blob/main/INSTALL.md)
- [CMake](https://cmake.org/)
- C++17 compliant compiler

## Set up and running

Execute from the root of this directory

```bash
cmake -B build
cmake --build build --parallel 4
./build/bench
```
