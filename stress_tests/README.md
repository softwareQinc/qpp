# Stress test suite (for POSIX-compliant systems) 

Tests some quantum circuit/operation as a function of the number of cores of
the machine and number of qubits. Records the time for each such run in
`results_$DATE.csv`, where `$DATE` represents the current date of the run.
Requires [OpenMP](http://openmp.org/).

To compile, I recommend using [CMake](http://www.cmake.org), then perform an
out-of-source build. Assuming you are
at the root of the project, i.e. inside `qpp`, type

```bash
cd ./stress_tests
mkdir ./build
cmake ..
make
```

To run, change the directory back to `stress_tests` 

```bash
cd ..
```

then type

```bash
bash run.sh
```

To modify the stress test configuration, edit the `run.sh` file, it should be
self-explanatory. In case you want to switch or to add a stress test, do not forget
to edit the `CMakeLists.txt` and modify the line 

```bash
SET(SOURCE_FILES
        src/your_stress_test.cpp)
```        
accordingly.