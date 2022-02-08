# Stress test suite (for POSIX compliant systems)

Tests some quantum circuit/operation as a function of the number of cores of the
machine and number of qubits (hard-coded in the script `run.sh`). Records the
time for each such run in `results_$DATE.csv`, where `$DATE` represents the
current date and time of the run. Requires [OpenMP](https://www.openmp.org/).

To compile, I recommend using [CMake](https://cmake.org), then perform an
out-of-source build. Assuming you are at the root of this project, i.e., inside
`qpp/stress_test`, type

```bash
mkdir build
cmake ..
make -j4
```

To run, change the directory back to `qpp/stress_tests`

```bash
cd ..
```

then type

```bash
bash run.sh <path_to_stress_test_executable> <results_output_directory>
```

e.g.,

```bash
bash run.sh build/qft results
```

The results will be written in `results_output_directory/results_$DATE.csv`. If
the directory `results_output_directory` does not exist, it will be created.

## Python stress tests

We wrote some [Qiskit](https://qiskit.org/) and [QuTiP](https://qutip.org/)
stress tests scripts, all located in the `python` directory. The output format
they produce is exactly the same as the one produced by the C++ version, and
the `run.sh` script works with them with no modification, e.g.,

```bash
bash run.sh python/qft_qutip.py results
```

The command above assumes that [Qiskit](https://qiskit.org/) and
[QuTiP](https://qutip.org/) were successfully installed on your system. For
installation details, please refer to their corresponding documentation.

## Windows (MSVC) support

To build the C++ stress tests with MSVC on Windows, follow the same instructions
as above but replace `make -j4` with `msbuild stress_tests.sln` and adapt the
`run.sh` script accordingly.
