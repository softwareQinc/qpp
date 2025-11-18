# Version 6.1.0 - xx November 2025

- Massive performance improvements (10x to 1000x) for qubit state manipulation
  and qubit circuit execution through extensive optimizations of `qpp::apply()`
  and `qpp::applyCTRL()`
- Switched back to MAJOR.MINOR.PATCH release versioning:
  - MAJOR -- introduces significant changes that break backward compatibility
  - MINOR -- may include limited or minor compatibility breaks
  - PATCH -- fully backward compatible fixes and improvements
- New feature: support for setting dits at runtime in
  ["qpp/classes/qcircuit.hpp"], implemented
  `qpp::QCircuit::set_dits_runtime(mutable_dits_functor_t dits)`. Here
  `mutable_dits_functor_t` is an alias to an internal mutable functor defined
  in ["qpp/internal/classes/qcircuit_runtime_step.hpp"]
- Changed `cond_func_t` type alias to `cond_pred_t` for
  Boolean predicates in conditional statements and moved the type alias from
  ["qpp/types.hpp"] to
  ["qpp/internal/classes/qcircuit_conditional_step.hpp"]; alias to an internal
  type.
- Renamed the class `internal::QCircuitConditionalStep` to
  `internal::QCircuitRuntimeStep`, and the corresponding
  ["qpp/internal/classes/qcircuit_conditional_step.hpp"] to
  ["qpp/internal/classes/qcircuit_runtime_step.hpp"]
- Renamed `qpp::QCircuit::has_conditionals()` to
  `qpp::QCircuit::has_runtime_steps()`
- Added
  [["examples/circuits/runtime_set_dits.cpp"](https://github.com/softwareQinc/qpp/blob/main/examples/circuits/runtime_set_dits.cpp)]
- Conditional statements are now indented with tabs when displaying
  qpp::QCircuit instances
- Added benchmarking suite using [Catch2](https://github.com/catchorg/Catch2)
  in ["benchmarks"](https://github.com/softwareQinc/qpp/blob/main/benchmarks)
- Bumped Eigen3 minimum required version to 5.0.0
- Bumped CMake minimum required version to 3.20
- Modernized CMake configuration files
- Removed `EIGEN3_INSTALL_DIR` CMake flag

# Version 6.0 - 14 April 2025

- Breaking change: renamed ["qpp/qpp.h"] to ["qpp/qpp.hpp"]
- New feature: support for post-selection in
  ["qpp/classes/qengine.hpp"] and ["qpp/classes/qcircuit.hpp"]
- Implemented `qpp::QCircuit::`
  - `QCircuit& post_select()` - destructive/non-destructive post-selection of
    single/multiple qudits in the Z basis
  - `QCircuit& post_selectV()` - destructive/non-destructive post-selection of
    single/multiple qudits in an arbitrary orthonormal basis
- Implemented `qpp::QEngineT<>::`
  - `bool get_ensure_post_selection() const` - True if post-selection is
    enforced, false otherwise
  - `QEngineT<>& set_ensure_post_selection()` - Enforces post-selection (i.e.,
    post-selection steps are repeated until success)
  - `idx get_max_post_selection_reps() const` - Maximum number of executions of
    a circuit post-selection step until success
  - `QEngineT<>& set_max_post_selection_reps()` - Sets the maximum number of
    executions of a circuit post-selection step until success
  - `bool post_select_ok() const` - True if post-selection was successful (or
    absent), false otherwise
- Added
  [["examples/circuits/post_selection.cpp"](https://github.com/softwareQinc/qpp/blob/main/examples/circuits/post_selection.cpp)]
  example
- New feature: implemented support for conditional statements in
  ["qpp/classes/qcircuit.hpp"]
- Implemented `qpp::QCircuit::`
  - `QCircuit& cond_if()` - conditional IF statement
  - `QCircuit& cond_else()` - conditional ELSE statement
  - `QCircuit& cond_while()` - conditional WHILE statement
  - `QCircuit& cond_end()` - conditional END block delimiter
  - `bool has_conditionals() const noexcept` - true if and only if the circuit contains
    conditional statements
  - `bool validate_conditionals() const` - true if an only if the conditional
    statements are valid (e.g., matching `cond_end()` to `cond_if()` etc.
- Added
  [["examples/circuits/conditional_if.cpp"](https://github.com/softwareQinc/qpp/blob/main/examples/circuits/conditional_if.cpp)]
  and
  [["examples/circuits/conditional_while.cpp"](htps://github.com/softwareQinc/qpp/blob/main/examples/circuits/conditional_while.cpp)]
  examples
- Refactored `qpp::QCircuit::GateStep/MeasurementStep/NOPStep` into separate
  files ["qpp/internal/classes/qcircuit_gate_step.hpp"],
  ["qpp/internal/classes/qcircuit_measurement_step.hpp"], and
  ["qpp/internal/classes/qcircuit_nop_step.hpp"], respectively
- Refactored `qpp::QCircuit::Resources` into an independent class in a separate
  file ["qpp/internal/classes/qcircuit_resources.hpp"]
- Refactored qpp::QCircuit::iterator class into an independent class, defined
  outside qpp::QCircuit in ["qpp/classes/qcircuit.hpp"]
- Refactored `qpp::internal::QEngineState` and
  `qpp::internal::QEngineStatistics`
  in separate files, ["qpp/internal/classes/qengine_state.hpp"] and
  ["qpp/internal/classes/qengine_statistics.hpp"], respectively
- API changes in ["qpp/classes/qcircuit.hpp"]
  - `qpp::QCircuit::get_measured()` -> `qpp::QCircuit::get_measured_d()`
  - `qpp::QCircuit::get_non_measured()` -> `qpp::QCircuit::get_non_measured_d()`
  - `qpp::QCircuit::was_measured()` -> `qpp::QCircuit::was_measured_d()`
- API changes in ["qpp/classes/qengine.hpp"]
  - `qpp::QEngineT<>::get_measured()` -> `qpp::QEngineT<>::get_measured_d()`
  - `qpp::QEngineT<>::get_non_measured()` ->
    `qpp::QEngineT<>::get_non_measured_d()`
  - `qpp::QEngineT<>::was_measured()` -> `qpp::QEngineT<>::was_measured_d()`
- Bugfix in `qpp::internal::canonical_form()`, the re-ordering is now stable,
  so `qpp::QCircuit` measurement probabilities are not displayed in reversed
  order w.r.t. target
- Simplified MATLAB detection via CMake `find_package()` function. Users should
  only use `-DQPP_MATLAB=ON` when building with MATLAB support, all other
  MATLAB-related CMake flags have been removed.
- Bugfix in `qpp::adjoint(QCircuit)`
- Added `cond_func_t` type alias in ["qpp/types.hpp"] for Boolean predicates of
  the form `std::vector<idx> -> bool`
- Added `qpp::read_from_string()` to ["qpp/qasm/qasm.hpp"] and an associated
  pyqpp wrapper function

# Version 5.1 - 1 March 2024

- Replaced ["CHANGES"] by ["CHANGES.md"],
  as we now use Markdown format to keep track of changes in new releases
- Removed Eigen3, pybind11, and GoogleTest dependencies; if not detected,
  they are installed automatically as build dependencies by CMake
- Bumped GoogleTest version to HEAD latest, as
  [recommended by Google](https://github.com/google/googletest?tab=readme-ov-file#live-at-head)
- Removed `-DWITH_EXAMPLES` and `-DWITH_UNIT_TESTS` CMake flags. Now both
  `examples` and `unit_tests` CMake targets are enabled.
- Renamed ["qpp/classes/circuits/circuits.hpp"] to
  ["qpp/classes/qcircuit.hpp"]
- Introduced ["qpp/classes/qcircuit_traits.hpp"] that implement
  circuit traits at compile-time
- Renamed ["qpp/classes/circuits/engines.hpp"] to
  ["qpp/classes/qengine.hpp"], and refactored the latter into
  - ["qpp/classes/qbase_engine.hpp"] - base class for all engines
  - ["qpp/classes/qengine.hpp"] - ideal quantum engines
  - ["qpp/classes/qnoisy_engine.hpp"] - noisy quantum engines
- Introduced ["qpp/classes/qengine_traits.hpp"] that implement
  engine traits at run-time. All engines are now deriving from it.
  The traits implements `qpp::IQEngineTraits::`
  - `std::string traits_get_name() const` - Engine's name
  - `bool traits_is_noisy() const` - Simulates noisy execution
  - `bool traits_is_pure() const` - Operates on pure states
- API changes in ["qpp/classes/qengine.hpp"] and
  ["qpp/classes/qnoisy_engine.hpp"]
  - Enabled mixed-state engines by refactoring
    - `qpp::QEngine` -> `template<typename T> qpp::QEngineT<T>`
    - `qpp::QNoisyEngine` -> `template<typename T> qpp::QNoisyEngineT<T>`
      The template argument T is restricted to `qpp::ket` (pure state
      engines) and `qpp::cmat` (mixed states engines)
  - The following additional engines are now available
    - `qpp::QEngine` - pure state ideal engine, backwards
      compatibility
    - `qpp::QKetEngine` - same as `qpp::QEngine`
    - `qpp::QDensityEngine` - mixed state ideal engine
    - `qpp::QNoisyEngine` - pure state noisy engine, backwards
      compatibility
    - `qpp::QKetNoisyEngine` - same as `qpp::QNoisyEngine`
    - `qpp::QDensityNoisyEngine` - mixed state noisy engine
  - Renamed `qpp::QEngineT<>::get_psi()` -> `qpp::QEngineT<>::get_state()`
  - Removed `qpp::QEngineT<>::is_noisy()`
  - Added the new engines to **pyqpp**, which now defines the following
    factory functions for instantiating engines
    - `pyqpp.QEngine()` - pure state ideal engine, backwards
      compatibility
    - `pyqpp.QKetEngine()` - same as `pyqpp.QEngine()`
    - `pyqpp.QDensityEngine()` - mixed state ideal engine
    - `pyqpp.QNoisyEngine()` - pure state noisy engine, backwards
      compatibility
    - `pyqpp.QKetNoisyEngine()` - same as `pyqpp.QNoisyEngine()`
    - `pyqpp.QDensityNoisyEngine()` - mixed state noisy engine
  - Removed the default argument `bool try_sampling = true` in
    `qpp::QEngineT<>::execute(idx reps = 1, bool try_sampling = true)` ->
    `qpp::QEngineT<>::execute(idx reps = 1)`
    so now `qpp::QEngineT<>` will always try to sample from the output when
    the circuit is executed multiple times (i.e., when `reps > 1`)
- Introduced no-op (dummy) quantum engines in ["qpp/classes/qdummy_engine.hpp"]
  that provides
  - `template<typename T, typename QCT> qpp::QDummyEngine<T, QCT>`
  - `qpp::QKetDummyEngine` - specialization of `qpp::QDummyEngine<>` with
    `T = qpp::ket`, `QCT = qpp::QCircuit`
  - `qpp::QDensityDummyEngine` - specialization of `qpp::QDummyEngine<>` with
    `T=qpp::cmat`, `QCT=qpp::QCircuit`
- Introduced corresponding `pyqpp.QKetDummyEngine` and
  `pyqpp.QDensityDummyEngine` in **pyqpp**

# Version 5.0 - 10 January 2024

- All header files are moved into ["include/qpp"], so to include the header
  one now must `#include "qpp/qpp.h"` (and similarly for other header files
  in the internal implementation). This change was made for the sake of
  making the include statements look uniform in both non-installed
  (compiling without having **Quantum++** installed) and installed (compiling
  with **Quantum++** installed headers) modes.
- Refactored option-related types into a new header file ["qpp/options.hpp"]
- `qpp::disp()` refactoring, formatting options are now passed via
  formatting option manipulator structures, all defined in
  ["qpp/options.hpp"]. See the ["examples"] directory for usage examples.
- Implemented `qpp::dirac()` ["qpp/functions.hpp"] for converting states and
  matrices to Dirac notation, and the corresponding `std::ostream` manipulator
  `qpp::disp()` ["qpp/input_output.hpp"]. The type returned by `qpp::dirac()`
  is type-defed as `qpp::dirac_t` in ["qpp/types.hpp"]. Options are passed
  via the formatting option manipulator structure `qpp::IOManipDiracOpts`.
- Implemented `pyqpp.dirac()` in **pyqpp**, see above.
- Significant speedup in **pyqpp**'s compilation time
- Added support for two qubit rotations by pi/2 (Molmer-Sorensen two qubit
  gates).
- API changes in ["qpp/classes/circuits/circuits.hpp"]
  - [qpp::QCircuit::]
    - `add_circuit()` -> `compose_circuit()`
    - `match_circuit_left()` -> `couple_circuit_left()`
    - `match_circuit_right()` -> `couple_circuit_right()`
- Implemented ["qpp/classes/circuits/circuits.hpp"] control on another
  quantum circuit,
  - `qpp::QCircuit::compose_CTRL_circuit()`
    see [GitHub Issue \#165](https://github.com/softwareQinc/qpp/issues/165),
    and the corresponding standalone versions
  - `qpp::compose_CTRL_circuit()`
- `qpp::measure_seq()` now returns results in the same order as the order in
  which the qubits are being measured
- API changes in ["qpp/functions.hpp"], `qpp::zket2dits()` now returns
  `std::optional<std::vector<idx>>`
- Bugfix in `qpp::QCircuit::compose_circuit()`

# Version 4.3.4 - 14 August 2023

- Docker update, see the ["docker"] directory
- Fix in ["types.hpp"] for defaulting types when not using CMake

# Version 4.3.3 - 8 August 2023

- Minor bugfix in **pyqpp** ["setup.py"] that prevented pip install from remote
- Bugfix in `qpp::QEngine::execute()` that prevented setting the initial
  state of the engine to a custom state

# Version 4.3.2 - 12 June 2023

- This is a maintenance release
- Compiling errors fixed on SunOS/OpenIndiana

# Version 4.3.1 - 5 June 2023

- CMake dependent flag name changes. These flags can be used in standalone
  projects when configuring with CMake.
  - `USE_OPENQASM2_SPECS` -> `QASMTOOLS_QASM2_SPECS`
  - `WITH_MATLAB` -> `QPP_MATLAB`
  - `WITH_OPENMP` -> `QPP_OPENMP`
  - `TYPE_BIGINT` -> `QPP_BIGINT`
  - `TYPE_FP` -> `QPP_FP`
  - `TYPE_IDX` -> `QPP_IDX`
- Compile-time error fix when using `FP_TYPE="long double"` with
  `QASMTOOLS_QASM2_SPECS=ON`

# Version 4.3 - 26 May 2023

- This is a maintenance release
- Fixed **pyqpp** installation on Windows under MSVC
- Migrated all continuous integration to GitHub actions
- Minor updates in qasmtools

# Version 4.2 - 13 May 2023

- **Quantum++** is now available on Homebrew (macOS/Linux), and can be
  installed with `brew install quantum++`
- New types in ["types.hpp"]
  - `realT` - floating point type, replaces `double` in the source code
  - `ubigint` - unsigned big integer, defined in terms of bigint as
    `using ubigint = std::make_unsigned<bigint>::type;`
- The fundamental types in ["types.hpp"] are now `qpp::idx`, `qpp::bigint`,
  and `qpp::realT`. They can be changed at compile time (see below). All
  other types are dependent on those, please do not change.
- Added support for changing the underlying index type via the CMake
  property `TYPE_IDX`. Choices are `default`, `short`, `int`, `long`,
  `long long`, `unsigned short`, `unsigned int`, `unsigned long`, and
  `unsigned long long`. The `default` type is set in ["types.hpp"]. To
  change the underlying type with CMake, pass the variable `-DTYPE_IDX` to
  CMake, e.g., `cmake -B build -DTYPE_IDX="unsigned long long"`.
- Added support for changing the underlying big integer type via the CMake
  property `TYPE_IDX`. Choices are `default`, `short`, `int`, `long`,
  `long long`. The `default` type is set in ["types.hpp"]. To change the
  underlying type with CMake, pass the variable `-DTYPE_BIGINT` to CMake,
  e.g., `cmake -B build -DTYPE_BIGINT="long long"`.
- Added support for changing the underlying floating point type via the
  CMake property `TYPE_FP`. Choices are `default`, `float`, `double`, and
  `long double`. The `default` type is set in ["types.hpp"]. To change the
  underlying type with CMake, pass the variable `-DTYPE_FP` to CMake, e.g.,
  `cmake -B build -DTYPE_FP="long double"`.
- Refactored the two instances of `qpp::rand(a, b)` ["random.hpp"] to a single
  template function
  ```cpp
  template<typename T, std::enable_if...> qpp::rand(T a, T b)
  ```
  enabled only on arithmetic types by
  `typename std::enable_if_t<std::is_arithmetic_v<T>>* = nullptr`
- Refactored qpp::randn(mean, sigma) ["random.hpp"] to a template function
  ```cpp
   template<typename T, std::enable_if...> qpp::random(T a, T b)
  ```
  enabled only on arithmetic types by
  `typename std::enable_if_t<std::is_arithmetic_v<T>>* = nullptr`
- CMake minimum required version bumped to 3.15

# Version 4.1 - 2 May 2023

- Implemented ["classes/circuits/circuits.hpp"]
  - `qpp::qpe_circuit()` - Quantum phase estimation circuit
  - `qpp::QCircuit::has_measurements()` - True if the quantum circuit
    description contains measurement steps, false otherwise
  - `qpp::QCircuit::removes_qudits()` - True if the quantum circuit
    description contains any operations/measurements that remove qudits,
    false otherwise
- Added new exception class ["classes/exception.hpp"]:
  - `qpp::exception::NotFound` - Element not found
- Implemented ["pyqpp/qpp_wrapper.cpp"] (**pyqpp** Python wrapper)
  - `set_prng_seed()` - sets the seed of the PRNG to a specific value, or,
    if invoked with no argument, to a random value
- Migrated **pyqpp** installation method to ["pyproject.toml"]
- Extensive code refactoring, transitioned `qpp::QCircuit` logic to
  `std::visit()` over `std::variant`
- API changes in ["classes/circuits/circuits.hpp"], all functions with
  `name` and `shift` defaulted empty/zeroed arguments are now wrapped in
  `std::optional<>`, as well as some `qpp::QCircuit` getters, i.e., have the
  signature `qpp::QCircuit::`:
  - `get_gate_count(std::optional<cmat> U = std::nullopt)`
  - `get_gate_depth(std::optional<cmat> U = std::nullopt)`
  - `get_measurement_count(std::optional<cmat> V = std::nullopt)`
  - `get_measurement_depth(std::optional<cmat> V = std::nullopt)`
    In addition, `std::optional<>` is also added as arguments to the free
    functions:
  - `qpp::random_circuit_count()`
  - `qpp::random_circuit_depth()`
- API changes in ["classes/circuits/circuits.hpp"], more `std::optional<>`,
  updated the following `qpp::QCircuit` member functions signatures to
  `qpp::QCircuit::`:
  - `add_circuit(QCircuit other, bigint pos_qudit, std::optional<idx> pos_dit = std::nullopt)`
  - `match_circuit_left(QCircuit other, const std::vector<idx>& target, std::optional<idx> pos_dit = std::nullopt)`
  - `match_circuit_right(QCircuit other, const std::vector<idx>& target, std::optional<idx> pos_dit = std::nullopt)`
  - `add_circuit(QCircuit qc1, const QCircuit& qc2, bigint pos_qudit, std::optional<idx> pos_dit = std::nullopt)`
  - `match_circuit_left(QCircuit qc1, const QCircuit& qc2, const std::vector<idx>& target, std::optional<idx> pos_dit = std::nullopt)`
  - `match_circuit_right(QCircuit qc1, const QCircuit& qc2, const std::vector<idx>& target, std::optional<idx> pos_dit = std::nullopt)`
- API changes in ["classes/circuits/circuits.hpp"], renamed `qpp::QCircuit::`:
  - `was_measured(idx i) -> was_measured(idx i)`
  - `get_measured_nd(idx i) -> was_measured_nd(idx i)`
- API changes in ["classes/circuits/engines.hpp"], renamed `qpp::QEngine::`:
  - `was_measured(idx i) -> was_measured(idx i)`
  - `get_measured_nd(idx i) -> was_measured_nd(idx i)`
- API changes in ["classes/reversible.hpp"], `qpp::BitCircuit::count/depth`
  getters now wrap their gate name argument in `std::optional<>`, i.e.,
  have the signature `qpp::BitCircuit::`:
  - `get_gate_count(std::optional<std::string> name = std::nullopt)`
  - `get_gate_depth(std::optional<std::string> name = std::nullopt)`
- API changes in ["functions.hpp"], `qpp::n2multiidx()` and `qpp::multiidx2n()`
  now take their argument underlying type (number, array) as template
  (must be an integral type)
- Removed `qpp::QCircuit::is_non_CTRL()`
- Refactored **pyqpp**, organized the source code in separate
  ["src"]/["include"] directory structure that mirrors the qpp source code
  structure

# Version 4.0.1 - 2 April 2023

- This is a maintenance release with minor changes
- Updated pybind11 to version 2.10.4
- Bugfix in `qpp::QFT()`/`TFQ()` when applied to density matrices

# Version 4.0 - 20 March 2023

- Major performance improvement in `qpp::QEngine::execute()`, which now takes
  a second default argument `bool try_sampling = true`; if possible, the
  engine will sample from the output of the circuit instead of executing
  the measurement steps.
- API changes in ["classes/circuits/circuits.hpp"], renamed:
  - `qpp::QCircuit::gate_joint()` -> `qpp::QCircuit::gate()`
  - `qpp::QCircuit::CTRL_joint()` -> `qpp::QCircuit::CTRL()`
  - `qpp::QCircuit::cCTRL_joint()` -> `qpp::QCircuit::cCTRL()`
  - `qpp::QCircuit::CTRL()` -> `qpp::QCircuit::CTRL_fan()`
  - `qpp::QCircuit::cCTRL()` -> `qpp::QCircuit::cCTRL_fan()`
  - `qpp::QCircuit::measureZ()` -> `qpp::QCircuit::measure()`
- API changes in `qpp::QCircuit::measure()` for multiple qudits, now the
  results are stored sequentially (instead as a decimal representation
  of the multi-index)
- API changes in `qpp::measure_seq()`, now, as part of its return tuple,
  returns the full list of outcome probabilities
- API changes for
  - `qpp::QCircuit::get_gate_count()`
  - `qpp::QCircuit::get_gate_depth()`
  - `qpp::QCircuit::get_measurement_count()`
  - `qpp::QCircuit::get_measurement_depth()`
    so now the old `std::string` overloads take instead a `qpp::cmat` as an
    argument
- API changes for
  - `qpp::Q[Noisy]Circuit::execute()`
    removed `bool reset_stats = true` parameter as it was confusing, suffices
    to invoke `Q[Noisy]Engine::reset(bool reset_stats = true)` for resetting
    the measurement statistics; added `bool try_sampling = true` argument, see
    above
- Noise is now added before measurements as well when running with a
  `qpp::QNoisyEngine<>`
- Implemented ["classes/circuits/circuits.hpp"]:
  - `qpp::QCircuit::measure_all()` - measures all remaining available qudits
- Implemented the random `qpp::QCircuit` free function generators in
  ["classes/circuits/circuits.hpp"]:
  - `qpp::random_circuit_count()` - random `QCircuit` w.r.t. gate count
  - `qpp::random_circuit_depth()` - random `QCircuit` w.r.t. gate depth
- Implemented ["classes/gates.hpp"]:
  - `qpp::Gates::GATE()` - constructs the matrix representation of a
    multi-partite gate that acts on a subsystem
- Implemented ["operations.hpp"]:
  - `qpp::applyCTRL_fan()` - applies CTRL-CTRL-...-CTRL-U-U-...-U
- Implemented ["instruments.hpp"]:
  - `qpp::sample()` - samples from a quantum state. Use this function
    whenever you are not interested in the output quantum
    state, but interested only in statistics; it is
    significantly faster than its `qpp::measure()`-type
    functions.
- Added copy/deepcopy support for all **pyqpp** classes:
  - `QCircuit`, `QEngine`, `QNoisyEngine<>`, `Dynamic_Bitset`, `Bit_circuit`
- Added `context` to exception throwing, so now it is obvious which argument
  triggers an exception
- `EIGEN3_INSTALL_DIR` can now be set as an OS environment variable (in
  addition to a CMake variable). If both are set, then the CMake variable
  trumps. If none are set, the system tries to detect the location of the
  Eigen3 library automatically.
- Marked high-performance functions with the custom attribute
  `[[qpp::critical]]` and functions that use parallelization with the custom
  attribute `[[qpp::parallel]]`. Since C++17, unknown attributes are supposed
  to be ignored by the compiler (and one can use this technique to define
  custom attributes). Warnings are explicitly disabled for
  Clang/GCC/MSVC/Intel. If your compiler emits a warning, please disable it
  with the corresponding `#pragma` directive at the beginning of ["qpp.h"]
- Added from-string constructor for `qpp::Bit_circuit` and
  `qpp::Dynamic_bitset`
- When building with CMake, the new `QPP_VERSION_NUM` (numeric) and
  `QPP_VERSION_STR` (string) preprocessor definitions are automatically
  injected in the code
- Added initializer lists overloads for `qpp::prod()` and `qpp::sum()`
  in ["functions.h"]:
  - `qpp::prod(const std::initializer_list<T>&)`
  - `qpp::sum(const std::initializer_list<T>&)`
- Enhanced `qpp::prod()` and `qpp::sum()` to allow list of matrices
- Bugfix in `qpp::measure()` [https://github.com/softwareQinc/qpp/issues/132]
- Bugfix in `QCircuit::CTRL` [https://github.com/softwareQinc/qpp/issues/130]
- Bumped
  [GoogleTest version to 1.12.1](https://github.com/google/googletest/commit/58d77fa8070e8cec2dc1ed015d66b454c8d78850)

# Version 3.1 - 11 January 2022

- Added **pyqpp**, a Python wrapper around **Quantum++**. See
  [**pyqpp** documentation](https://github.com/softwareQinc/qpp/wiki/8.-pyqpp)
  for more details.
- Minor update of Shor's algorithm example ["examples/shor.cpp"]
- Renamed ["examples/qasm/qasm.cpp"] to ["examples/qasm/qpp_qasm.cpp"]
- Due to phase discrepancies, by default the parser now uses Qiskit
  definitions (which are also the usual ones used in QC textbooks). To
  switch to standard OpenQASM 2.0 gate definitions, configure the project
  with `cmake -DUSE_OPENQASM2_SPECS=ON`.
- Implemented `qpp::QCircuit` methods:
  - `get_dirty_dits()` - vector of used classical dits in the circuit
    (by either cCTRL or measurements)
  - `get_dirty_qudits()` - vector of used qudits in the circuit
  - `was_measured_nd()` - vector of qudits that were measured
    non-destructively
  - `get_measurement_dits()` - vector of classical dits that participated in
    measurements
  - `is_measurement_dit()` - whether a classical dit participated in a
    measurement
- Changed the ordering of bits in `qpp::Bit_circuit::display()`, so that the
  zero-th bit (top bit) is now the leftmost bit in the textual
  representation. That is, `qpp::Bit_circuit` assumes big-endian order.
- Bugfixes

# Version 3.0 - 5 October 2021

- Major release, bumped up C++ standard to C++17
- Decoupled the OpenQASM parser from the main codebase. A hard copy of
  qasmtools
  ([https://github.com/softwareQinc/qasmtools](https://github.com/softwareQinc/qasmtools))
  is now common to both **staq** and **Quantum++**, and by default uses
  standard OpenQASM 2.0 gate definitions. To switch to Qiskit definitions
  (which are also the usual ones used in QC textbooks), configure the project
  with
  `cmake -DUSE_QISKIT_SPECS=ON`.
- Refactored `qpp::internal::Singleton` ["internal/classes/singleton.hpp"]
  so that `qpp::internal::Singleton::get_instance()` returns a thread local
  instance whenever the compiler supports thread_local. For more
  fine-grained control, one can use
  - `qpp::internal::Singleton::get_no_thread_local_instance()`
  - `qpp::internal::Singleton::get_thread_local_instance()`
- Added `qpp::schmidt()` ["entanglement.hpp"] and improved the performance
  of singular value decomposition, thanks @antoine-bussy
  ([https://github.com/antoine-bussy](https://github.com/antoine-bussy))

# Version 2.7 - 1 August 2021

- Added installation support and auto package detection via CMake
  `find_package(qpp REQUIRED)`
- Simplified installation instructions ["INSTALL.md"]
- Added an almost complete BB84 example ["examples/bb84.cpp"], does not
  include error correction and privacy amplification
- Added CircleCI continuous integration
- CMake minimum required version bumped to 3.12 for automatic unit tests
  detection by CMake and other CMake-related "good practices"
- Unit testing and examples are now separate CMake targets, one needs to
  explicitly type `make unit_testing` and/or `make examples` to build the
  unit testing suite and/or the examples, respectively
- Simplified unit testing, now one can run tests with `ctest` or
  `make test` (after explicitly built with `make unit_testing`). Use
  `GTEST_COLOR=1 ARGS="-V" make test` or `GTEST_COLOR=1 ctest -V` for
  coloured verbose testing output.
- Fixed minor inadvertence in OpenQASM `crz` implementation, see
  [GitHub Issue \#99](https://github.com/softwareQinc/qpp/issues/99)
- Changed the signature of `qpp::cwise()` ["functions.hpp"] so it takes the
  scalar argument by value (and not constant reference as before). This
  avoids an internal compiler error (ICE) triggered when compiling
  ["examples/functor.cpp"] in MinGW.
- Removed `<iostream>` implicit dependencies ["qpp.h"]
- Added optional context to `qpp::Exception` ["classes/exception.hpp"]
- Made `qpp::QCircuit` default-constructible ["classes/circuits/circuits.hpp"]
  so it constructs by default a 1-qubit circuit with no classical bits
- `qpp::QEngine::reset()` now accepts a default argument
  `bool reset_stats = true` and cleans the engine's statistics by default
  ["classes/circuits/engines.hpp"]
- Renamed ["examples/qasm/qasm1.cpp"] and ["examples/qasm/qasm2.cpp"]
  to ["examples/qasm/qasm_teleport_minimal.cpp"] and
  ["examples/qasm/qasm.cpp"]; the latter executes an arbitrary OpenQASM
  program read from the standard input or from a file (if specified)
- Added the aggregate `qpp::QCircuit::Resources`
  ["classes/circuits/circuits.hpp"] class for quantum circuit resource
  estimations. It derives from `qpp::IDisplay` and `qpp::IJSON`, so one can
  display instances of it as text or in JSON format.
- Implemented `qpp::QCircuit::get_resources()` which returns an
  instance of `qpp::QCircuit::Resources`
- Renamed ["examples/circuits/quantum_phase_estimation.cpp"] to
  ["examples/circuits/qpe_circuit.cpp"]
- Added ["examples/qasm/coin_flip.qasm"] coin flipping OpenQASM example
  that can be run, e.g., with `qasm < coin_flip.qasm [number of reps]`
- `qpp::measure()` deducts the input state type (`qpp::ket` or `qpp::cmat`)
  and returns accordingly a set of kets or density matrices in its return
  tuple, thanks [@anoine-bussy](https://github.com/antoine-bussy)
- Implemented `qpp::QCircuit::match_circuit_left()` and
  `qpp::QCircuit::match_circuit_right()`, that match circuit qudits in
  arbitrary order
- Implemented `qpp::zket2dits` ["functions.h"] which extracts the dits from a
  normalized multi-partite pure state in the computational basis.
  Behaves like the inverse of `qpp::mket()`.
- Removed `cmath_cygwin.patch` as the problem was fixed since 2016.
- Minor refactorings, code simplifications and bugfixes
- Removed Travis CI

# Version 2.6 - 9 January 2021

- Added Quantum Phase Estimation low-level API example in
  ["examples/qpe.cpp"], courtesy of [@ryanhill1](https://github.com/ryanhill1),
  see
  [GitHub Pull Request \#91](https://github.com/softwareQinc/qpp/pull/91)
- Added Quantum Phase Estimation high-level API example in
  ["examples/circuits/quantum_phase_estimation.cpp"], thanks @DevelopDaily for
  the suggestion,
  [GitHub Issue \#96](https://github.com/softwareQinc/qpp/issues/96)
- `qpp::load()` and `qpp::save()` ["input_output.hpp"] now use C++ I/O streams
  instead of file names, and load/save in text format (instead of
  binary format), while preserving the required precision. The old load/save
  functions that use binary files for I/O are now moved to the `qpp::obsolete`
  namespace and are deprecated.
- Minor API changes:
  ["MATLAB/matlab.hpp"]
  `qpp::saveMATLAB()` -> `qpp::save_MATLAB()`
  `qpp::loadMATLAB()` -> `qpp::load_MATLAB()`
  ["classes/circuits/circuits.hpp"]
  `qpp::cCTRL_custom()` -> `qpp::cCTRL_joint()`
  `qpp::CTRL_custom()` -> `qpp::CTRL_joint()`
  `qpp::gate_custom()` -> `qpp::gate_joint()`
- Documentation improvements in ["classes/circuits/circuits.hpp"]
- Updated Travis CI configuration

# Version 2.5 - 28 November 2020

- Significant improvements in `qpp::QEngine::execute()` for multiple runs
- Changed all header extensions from `.h` to `.hpp` (except for `qpp.h`)
- Automatic OpenMP detection by CMake
  - CMake minimum required version bumped to 3.9 (3.12 for macOS) for
    automatic OpenMP detection
  - OpenQASM bugfixes
  - Renamed `master` branch to `main`

# Version 2.4 - 13 May 2020

- Enabled documentation searching by setting `SEARCHENGINE = ON` in
  `Doxyfile`
- Implemented the following static member functions `qpp::QCircuit`::
  - `is_CTRL()` - gate step is a controlled gate or not
  - `is_cCTRL()` - gate step is a classically-controlled gate or not
  - `is_non_CTRL()` - gate step is a non-controlled gate or not
- Changed `qpp::QCircuit`::
  - `is_clean()` -> `is_clean_qudit()`
  - `get_clean()` -> `get_clean_qudits()`
- Split `qpp::QCircuit::remove_clean()` into `qpp::QCircuit`::
  - `remove_clean_qudit()` - removes single clean qudit from quantum
    circuit description
  - `remove_clean_qudits()` - removes list of clean qudits from quantum
    circuit description
- Implemented `qpp::QCircuit`::
  - `is_clean_dit()` - true if the classical dit is clean (not used),
    false otherwise
  - `get_clean_dits()` - list of clean classical dits
  - `remove_clean_dit()` - removes single clean classical dit from
    quantum circuit description
  - `remove_clean_dits()` - removes list of clean classical dits from
    quantum circuit description
- Modified the signature of `qpp::QCircuit::compress()` to
  `qpp::QCircuit::compress(bool compress_dits = false)` so one can compress
  the clean classical dits as well
- Removed the measured/non-measured qudits from `qpp::QCircuit::display()` and
  `qpp::QEngine::display()` since they make the output very large for quantum
  circuit descriptions with many qudits; kept them in the corresponding
  `qpp::QCircuit::to_JSON()` and `qpp::QEngine::to_JSON()` member functions
- Implemented `qpp::QEngine::set_dits()` which sets all engine's classical
  dits at once

# Version 2.3 - 2 May 2020

- Added the following examples in the ["examples"] directory
  - `qram.cpp` - quantumly-accessible Random Access Memory over
    classical data
  - `layouts.cpp` - various physical qubit layouts
- Implemented `qpp::QCircuit`::
  - `bool operator==() const noexcept`
  - `bool operator!=() const noexcept`
    for quantum circuits simple (in)equality testing
- Implemented
  - `qpp::super2kraus()` - extracts a set of orthogonal (in the
    Hilbert-Schmidt norm) Kraus operators from a
    superoperator matrix ["operations.h"]
  - `qpp::bernoulli()` - generates Bernoulli-p random Booleans
    ["random.h"]
- Kraus-related functions on full systems (not subsystems) now support
  input and output spaces of different dimensions, i.e., can now use
  non-square Kraus operators. All other Kraus-related functions for which
  the Kraus operators act on subsystems continue to require square Kraus
  operators.
- Bugfix in `qpp::QCircuit::add_circuit()`
- Fixed all Doxygen warnings (except spurious ones)

# Version 2.2 - 14 April 2020

- Updated Google Test to version 1.10.0
- Removed the Visual Studio solution (since it can be automatically
  generated by CMake), as we prefer to use CMake uniformly across all
  platforms
- Implemented
  - `qpp::qRAM()` - quantumly-accessible Random Access Memory over
    classical data (qRAM) ["operations.h"]
- Added the new type `qpp::qram` for representing qRAM data ["types.h"]
- Thanks to @DevelopDaily (https://github.com/softwareQinc/qpp/issues/72),
  implemented `qpp::QCircuit`::
  - `is_clean()` - true if the qudit is clean (not used), false
    otherwise
  - `get_clean()` - list of clean qudits
  - `remove_clean()` - removes single/list of clean qudits from quantum
    circuit description
  - `compress()` - removes unused qudits from quantum circuit
    description
- More robust MATLAB integration, detects MATLAB automatically when the
  CMake flag `-DWITH_MATLAB=ON` is specified; if CMake cannot detect it
  automatically, pass the additional CMake flag
  `-DMATLAB_INSTALL_DIR=/path/to/MATLAB` to manually specify MATLAB's
  location
- Updated MATLAB mx-API, detects `MX_HAS_INTERLEAVED_COMPLEX` for MATLAB
  versions >= R2018a, see
  [MATLAB Support for Interleaved Complex](https://www.mathworks.com/help/matlab/matlab_external/matlab-support-for-interleaved-complex.html),
  but also continues to work with older MATLAB versions < R2018a
- Changed the CMake flag `EIGEN3_INCLUDE_DIR` to `EIGEN3_INSTALL_DIR` so it is
  consistent with the `MATLAB_INSTALL_DIR` CMake flag

# Version 2.1 - 14 March 2020 (3.14.2020 Pi day release)

- Migrated the repository to
  [https://github.com/softwareQinc/qpp](https://github.com/softwareQinc/qpp)
- Changed gate definitions in ["qasm/preprocessor.h"] so they agree with
  Qiskit, see
  [https://github.com/softwareQinc/qpp/issues/65](https://github.com/softwareQinc/qpp/issues/65)
  and
  [https://github.com/softwareQinc/qpp/issues/70](https://github.com/softwareQinc/qpp/issues/70)
- Minor API change in the `QEngine::execute()`, `reset()`, `reset_stats()`, now
  they all return a reference to `*this`
- Added `Bit_circuit::display()` and `Bit_circuit::to_JSON()` in
  ["reversible.h"]
- Added the option for shifted control gates in:
  - `qpp::applyCTRL()`
  - `qpp::Gates::CTRL()`
  - `qpp::QCircuit::CTRL()`
- Implemented
  - `qpp::States::j() const` - |j⟩ computational basis state of a single
    qudit

# Version 2.0 - 24 August 2019

- Added support for OpenQASM via the interface ["qasm/qasm.h"] containing:
  - `qpp::QCircuit qpp::qasm::read_from_file(const std::string& fname)` - reads
    an OpenQASM file and returns the qpp::QCircuit representation
  - `qpp::QCircuit qpp::qasm::read(std::ifstream& stream)` - reads an input
    stream and returns the qpp::QCircuit representation
- Added corresponding OpenQASM examples in the directory
  ["examples/circuits/qasm"]
- Implemented qpp::QCircuit::
  - `kron()` - Kronecker product with another quantum circuit description
  - `measureZ()` - new overload that allows multiple targets at once
  - `adjoint()` - in place adjoint of quantum circuit description
  - `kron()` - in place Kronecker product with another quantum circuit
    description
- Implemented the corresponding free functions as friends of qpp::QCircuit:
  - `qpp::adjoint()` - adjoint of quantum circuit description
  - `qpp::kron()` - Kronecker product between two quantum circuit descriptions
  - `qpp::replicate()` - replicates (repeats) a quantum circuit description
  - `qpp::add_circuit()` - adds two quantum circuit descriptions
- Implemented qudit resetting and discarding functions qpp::QCircuit::
  - `reset()` - resets qudits in a quantum circuit description by measuring
    them non-destructively in the computational basis followed by shifting them
    back to the |0> state
  - `discard()` - discards qudits in a quantum circuit description by measuring
    them destructively and ignoring the result of the measurement
- Implemented the corresponding free functions in ["instruments.h"]:
  - `qpp::reset()` - resets qudits in a multipartite state by measuring them
    non-destructively in the computational basis followed by shifting them back
    to the |0> state
  - `qpp::discard()` - discards qudits by measuring them destructively and
    ignoring the result of the measurement
- Implemented `qpp::QCircuit::set_name()` - set the QCircuit name after
  construction
- API change in `qpp::QCircuit::cCTRL*()`, now accepts an additional shift on
  the control dits
- Split `qpp::QCircuit::get_depth()` into `qpp::QCircuit::get_gate_depth()` -
  gate depth and `qpp::QCircuit::get_measurement_depth()` - measurement depth
  for finer control over depth computations
- Added support for non-destructive measurements in ["instruments.h"] by
  modifying the API for:
  - `qpp::measure()`
  - `qpp::measure_seq()`
- Added ["classes/layouts.h"] header file for representing different physical
  qudit layouts. Implemented the class(es):
  - `qpp::Lattice` - qudits on a multi-dimensional orthogonal grid
  - `qpp::PeriodicBoundaryLattice` - qudits on a multi-dimensional orthogonal
    grid with periodic boundary conditions
- Removed the documentation directory (["doc"]). To generate the documentation
  in HTML/LaTeX formats, run `doxygen` from the root of the project.

# Version 1.3 - 25 July 2019

- Minimum required CMake version is now 3.2
- Added MSVC support for CMake
- Updated Google Test to version 1.8.1
- Added `enum {RES, PROB, ST}` in ["constants.h"] for using with `std::get<>`
  on the result of `qpp::measure()` etc.
- Split the ["classes/circuits.h"] into ["classes/circuits/circuits.h"] and
  ["classes/circuits/engines.h"]
- Added noisy quantum engine for non-correlated noise
- Added noisy teleportation example (that uses a noisy quantum engine) in
  ["examples/circuits/noisy_teleport_qubit_circuit.cpp"]
- Implemented `qpp::QCircuit::`:
  - `QFT()` - Quantum Fourier Transform
  - `TFQ()` - inverse Quantum Fourier Transform
  - `get_gate_depth()` - circuit gate depth
  - `add_qudit()` - adds qudit(s) on the fly
  - `add_dit()` - adds classical dit(s) on the fly
  - `add_circuit()` - adds a quantum circuit description to the current one
  - `replicate()` - replicates (repeats) the current quantum circuit
    description
- Implemented `qpp::Bit_circuit::`:
  - `get_gate_count()` - classical reversible circuit gate count
  - `get_depth()` - classical reversible circuit gate depth

# Version 1.2 - 10 February 2019

- Added new function in ["functions.h"]:
  - `qpp::normalize()` - normalizes a state vector (ket or bra) or density
    matrix
- Added new functions in ["number_theory.h"]:
  - `qpp::convergents()` - overloads that compute the convergents from a
    continued fraction or a real number, see
    https://mathworld.wolfram.com/Convergent.html
- Added new functions in ["operations.h"]:
  - `qpp::TFQ()` - inverse qudit quantum Fourier transform (includes qubit QFT
    as a special case) on an entire quantum state (state vector or density
    matrix)
  - `qpp::applyTFQ()` - apply the inverse qudit quantum Fourier transform on a
    subset of qudits of a state vector or density matrix
- Added new member functions for `qpp::Gates`:
  - `qpp::Gates::MODMUL()` - quantum modular multiplication
- Added Shor's algorithm example in ["examples/shor.cpp"]
- Introduced the largest unsigned index `qpp::idx_inf` in ["constants.h"]
- Added new exception classes in ["classes/exception.h"]:
  - `qpp::exception::QuditAlreadyMeasured` - self-explanatory
  - `qpp::exception::Duplicates` - thrown when a system such as a std::vector
    contains duplicates
  - `qpp::exception::NotImplemented` - code not yet implemented
  - `qpp::exception::InvalidIterator` - invalid iterator (not initialized,
    past end, etc.)
- Added `std::string qpp::Gates::get_name(const cmat&)` function which returns
  the names of the most common (one, two and three) qubit gates
- Added new header file ["classes/noise.h"] and corresponding classes for
  simulating quantum noise (e.g. depolarizing, dephasing etc.) and a
  corresponding simple example in ["examples/noise.cpp"]
- Added new abstract class `qpp::IJSON` in ["classes/idisplay.h"] for very
  basic JSON serialization support, and corresponding
  `qpp::QCircuit::to_JSON()` and `qpp::QEngine::to_JSON()`
- Added new classes in the new header file ["classes/circuits.h"] for
  simulating qudit quantum circuits:
  - `qpp::QCircuit` - Quantum circuit description, does not effectively
    perform any operations
  - `qpp::QEngine` - Quantum engine for simulation of quantum circuits
- Added new circuit simulator examples in ["examples/circuits"]:
  - `teleport_qubit_circuit.cpp` - qubit teleportation circuit example
  - `teleport_qudit_circuit.cpp` - qudit teleportation circuit example
- Removed `qpp::eps` constant from ["constants.h"], as one should simply use
  direct comparison with 0 when dealing with small doubles
- Added support for hashing Eigen matrices/vectors/expressions via
  `qpp::hash_eigen()` in ["functions.h"]; code based on
  `boost::hash_combine()`, see
  https://www.boost.org/doc/libs/1_69_0/doc/html/hash/reference.html#boost.hash_combine

# Version 1.1 - 26 November 2018

- Added new functions in ["operations.h"]:
  - `qpp::QFT()` - qudit quantum Fourier transform (includes qubit QFT as a
    special case) on an entire quantum state (state vector or density matrix)
  - `qpp::applyQFT()` - apply the qudit quantum Fourier transform on a subset
    of qudits of a state vector or density matrix
- Added Quantum Fourier transform example in ["examples/qft.cpp"]
- Added new member functions for `qpp::Gates`:
  - `qpp::Gates::RX()` - rotation around X
  - `qpp::Gates::RY()` - rotation around Y
  - `qpp::Gates::RZ()` - rotation around Z
  - `qpp::Gates.SWAPd()` - qudit SWAP gate
- Added a suite of stress tests in ["stress_tests"] (only for POSIX systems)
- Added Qiskit and QuTiP stress tests in ["stress_tests/python"]

# Version 1.0 - 3 July 2018

- Added full support for reversible classical circuits, added
  ["classes/reversible.h"] header file that contains the following classes:
  - `qpp::Dynamic_bitset` - bitsets of variable length
  - `qpp::Bit_circuit` - classical reversible circuits
- Added ["examples/reversible.cpp"] example
- The path to MATLAB can now be specified in the CMake command line (not hard
  coded anymore), such as
  `cmake .. -DWITH_MATLAB="/Applications/MATLAB_R2017b.app"`
- `qpp::ket operator "" \_ket()`, `qpp::ket operator "" _bra()` and
  `qpp::ket operator "" _prj()` are now inside the newly introduced inline
  namespace `qpp::literals`
- Minor bugfix in `qpp::Exception::what()`, thanks Nick Lewycky for spotting it
- Updated the Wiki

# Version 1.0-rc4 - Release Candidate 4, 24 January 2018

- Changed `qpp::ket operator "" \_q()` to `qpp::ket operator "" \_ket()`
- Added `qpp::bra operator "" \_bra()` and `qpp::cmat operator "" \_prj()`
- Updated documentation

# Version 1.0-rc3 - Release Candidate 3, 21 January 2018

- License change from GPL to MIT
- Added Visual Studio 2017 solution
- Added support for AppVeyor continuous integration
- Added support for reversible classical circuits, via the following classes in
  ["experimental.h"]:
  - `qpp::experimental::Dynamic_bitset` - bitsets of variable length
  - `qpp::experimental::Bit_circuit` - classical reversible circuits
- Added `qpp::ket operator "" \_q()` helper in ["functions.h"] for constructing
  multi-partite qubit kets
- Added comprehensive Wiki instead of the quick starting guide

# Version 1.0-rc2 - Release Candidate 2, 6 September 2017

- Added serialization capabilities for the PRNG in `qpp::RandomDevices` from
  ["classes/random_devices.h"]:
  - `qpp::RandomDevices::load()` - loads the state of PRNG from a stream
  - `qpp::RandomDevices::save()` - saves the state of PRNG to a stream
- Added a getter for `qpp::RandomDevices::prng_` (now with private access):
  - `qpp::RandomDevices::get_prng()`
- Eigen3 is detected automatically by CMake (if Eigen3 was installed via a
  package manager). If not, it can be specified manually by passing the
  `-DEIGEN3_INCLUDE_DIR=/path/to/eigen3` flag to the CMake command line.

# Version 1.0-rc1 - Release Candidate 1, 11 November 2016

- Slight performance improvements (lambda workers do not capture by value
  anymore)
- API change in ["classes/exception.h"]: now all **Quantum++** exceptions have
  their own type, being derived from the base exception class
  `qpp::exception::Exception`. All exception classes are now in their separate
  namespace `qpp::exception`.

# Version 1.0.0-beta4 - 2 November 2016

- Added new functions in ["classes/states.h"]:
  - `qpp::States::zero()` - zero state of n qudits
  - `qpp::States::one()` - one state of n qudits
  - `qpp::States::jn()` - |jj...j> state of n qudits
  - `qpp::States::plus()` - plus state of n qubits
  - `qpp::States::minus()` - minus state of n qubits
- Corrected `qpp::is_matrix_expression<>` type trait so now it works fine
  Thanks to @davidhigh https://stackoverflow.com/a/40293333/3093378
- Bugfix in `qpp::cor()`
- Renamed `qpp::loadMATLABmatrix()` and `saveMATLABmatrix()` to
  `qpp::loadMATLAB()` and `qpp::saveMATLAB()`, respectively. Slightly changed
  their implementation so now one can load/save any Eigen type, similarly
  to `qpp::load()` and `qpp::save()`.
- Almost full code coverage via unit testing

# Version 1.0.0-beta3 - 22 October 2016

- Added support for Travis CI continuous integration
- Added new functions in ["classes/states.h"]:
  - `qpp::States::mes()` - maximally entangled state of 2 qudits
- Added new functions in ["random.h"]:
  - `qpp::randprob()` - random probability vector uniformly distributed
    over the probability simplex

# Version 1.0.0-beta2 - 13 October 2016

- Updated to GoogleTest 1.8.0 for unit testing
- Added d = 2 overloads for all functions that take a list of dimensions
  as one of their arguments
- Additional code sanitizing and unit testing

# Version 1.0.0-beta1 - 10 October 2016

- GoogleMock 1.7.0 is now included in the project for portable unit testing
  (see the README.md file for more details)
- Removed the `qpp::ubigint` type so there is no more danger of mixing signed
  and unsigned types. All related functions now use `qpp::bigint`.
- Fixed overflow bug in `qpp::powm()` when numbers were very large
- `qpp::isprime()` is now based on the Miller-Rabin primality test
- `qpp::rand(Type a, Type b)` family of functions now throw if b > a for
  integer types or when a == b for floating-point types
- Added new functions in ["number_theory.h"]:
  - `qpp::randprime()` - random prime number generator
  - `qpp::modmul()` - modular multiplication that avoids overflows
- Additional unit testing, added unit testing skeleton for all codebase

# Version 0.8.8.2 - 3 October 2016

- This is a bugfix release
- Compile-time bugfix in `qpp::measure_seq()` due to trying to modify a const
  copy, which should have been non-const

# Version 0.8.8 - 30 September 2016

- Added OpenMP support for clang >= 3.7 in `CmakeLists.txt`
- Required minimum CMake version is now 3.0
- Added more unit testing, ["number_theory.h"] is now fully tested
- Split unit testing into separate files according to the tested header file
- Fixed typo in `CMakeLists.txt` which in effect disabled OpenMP
  for version 0.8.6. PLEASE use the `CMakeLists.txt` from the current version!
- SFINAE for all functions that have STL-containers as parameters
- Improved performance via removing unnecessary copy, see Issue \#31, thanks
  @vasekp
- Corrected code typo in `qpp::internal::Singleton<T>`, see Issue \#35, thanks
  @titaschanda
- Critical bugfix in `qpp::lcm(const std::vector<bigint>&)`
- Critical bugfix in `qpp::x2contfrac()`
- Added new functions in ["number_theory.h"]:
  - `qpp::egcd()` - extended greatest common divisor algorithm
  - `qpp::modinv()` - modular inverse of a mod p

# Version 0.8.6 - 1 November 2015

- Critical bugfix in `qpp::applyCTRL()`, which was failing before for
  control gates with 2 or more controls. Added unit test in
  ["unit_tests/testing.cpp"].
- Removed ["macros.h"] preprocessor debug macros header added in v0.8.4

# Version 0.8.4 - 14 October 2015

- Added ["macros.h"] preprocessor debug macros header, containing:
  - `PRINT(x)` - `std::cout << (x)`
  - `PRINTLN(x)` - `std::cout << (x) << std::endl`
  - `ERROR(x)` - `std::cerr << (x)`
  - `ERRORLN(x)` - `std::cerr << (x) << std::endl`
- Added new functions in ["operations.h"]:
  - `qpp::ip()` - Generalized inner product
- Modified `qpp::measure()` so that the output consists of pure states
  whenever the input is a pure state and the measurement operators are
  rank-one POVMs (technically rank-one Kraus operators). Includes basis
  measurements as a particular case.
- Improved `qpp::Timer<>`, now both the duration and the clock's type can be
  specified as template parameters, defaulted to
  `std::chrono::duration<double>` and `std::chrono::steady_clock`, respectively

# Version 0.8.2 - 20 May 2015

- Added ["traits.h"] type traits header file with the following type traits:
  - `qpp::is_iterable<>` - checks for iterable STL-like containers
  - `qpp::is_matrix_expression<>` - checks for `Eigen::MatrixBase<>`
    expressions
  - `qpp::is_complex<>` - checks whether the type is `std::complex<>`
- Added ["statistics.h"] header file containing the following new functions:
  - `qpp::uniform()` - uniform probability distribution vector
  - `qpp::marginalX()` - marginal distribution
  - `qpp::marginalY()` - marginal distribution
  - `qpp::avg()` - average
  - `qpp::sigma()` - standard deviation
  - `qpp::cov()` - covariance
  - `qpp::cor()` - correlation
  - `qpp::var()` - variance
- Added new examples in the directory ["examples"]

# Version 0.8 - 8 May 2015

- Added new data type aliases in ["types.h"]:
  - `bigint` - Big integer, alias for `long long int`
  - `ubigint` - Non-negative big integer, alias for `unsigned long long int`
- Added new overloads in ["random.h"]:
  - `qpp::rand()` - overloads for `bigint` and `ubigint`
- Added new functions in ["number_theory.h"]:
  - `qpp::modpow()` - computes a^n mod p for non-negative integers
- Added `IDisplay` as an abstract base class for enforcing the override of
  its member function:
  - `std::ostream& IDisplay::display(std::ostream&) const`
    in all derived classes. This latter function performs all of the stream
    extraction operator work, and it is invoked automatically by the
    `friend std::ostream& operator<<(std::ostream&, const qpp::IDisplay&)`
    defined inline in `qpp::IDisplay`. In other words, the extraction operator
    delegates the work to the virtual `display()` function. Finally, all
    classes that support stream extraction (`qpp::IOManipEigen`,
    `qpp::IOManipPointer`, `qpp::IOManipRange`, and `qpp::Timer`) now inherit
    from `qpp::IDisplay`.
- Replaced `std::basic_ostream<charT, traits>` by `std::ostream` in all code,
  as using traits over-complicated the design and it was not worth it.
- Marked final all classes that are not intended to be used as base classes.
- Modified the `CMakeLists.txt`, which now:
  - supports only GCC and LLVM/Apple clang; otherwise, emits an error during
    the CMake running phase
  - adds `-D_NO_THREAD_LOCAL` definition for conditional compiling the
    source code without support for `std::thread_local` (if using e.g.,
    LLVM/Apple clang with libc++, as libc++ does not yet support
    `std::thread_local`)
  - uses OpenMP only if the compiler is detected as GCC

# Version 0.7 - 22 April 2015

- Marked all non-template functions defined in headers as `inline` and all
  singletons as `static` since otherwise including ["qpp.h"] in separate
  compilation units led to linker errors due to duplicated symbols.

# Version 0.6 - 16 April 2015

- Reverted changes to `qpp::internal::Singleton<T>` (no need anymore for
  `static void operator delete(void*)`), but clients have to explicitly mark
  both constructor and destructor private.
- Minor cosmetic changes.

# Version 0.5 - 19 March 2015

- Bugfix in class `qpp::internal::Singleton<T>` which allowed deleting the
  instance via a pointer. Now `static void operator delete(void*)` is marked
  private.

# Version 0.4 - 8 March 2015

- Added `noexcept` specification where suitable
- Added the following exception types to the `qpp::Exception::Type` strong
  enumeration: `NOT_QUBIT_MATRIX` (replaces previous `NOT_QUBIT_GATE`),
  `NOT_QUBIT_CVECTOR`, `NOT_QUBIT_RVECTOR` and `NOT_QUBIT_VECTOR`
- Added overloads of `qpp::ptrace1()`, `qpp::ptrace2()`, `qpp::ptrace()` and
  `qpp::ptranspose()` that allow a pure state input, as it was the case
  before with `qpp::syspermute()`
- Added overloads of `qpp::prod()` and `qpp::sum()` for STL-like containers
- Added new functions in ["functions.h"]:
  - `qpp::complement()` - computes the complement of a subsystem vector
  - `qpp::rho2bloch()` - qubit density matrix to Bloch vector
  - `qpp::bloch2rho()` - Bloch vector to qubit density matrix
- Added new functions in ["instruments.h"]:
- `qpp::measure_seq()` - measures subsystems sequentially
- Removed support for the C-style random number engine, as it is not thread
  safe. Hence, `Eigen::Matrix::Random()` should not be used, use instead
  `qpp::rand()` functions.
- Added support for unit testing via Google Mock/Google Test, see
  ["unit_tests/testing.cpp"]

# Version 0.3 - 15 February 2015

- Added `qpp::internal::Singleton<T>::get_thread_local_instance()`

# Version 0.2 - 26 January 2015

- Minor fixes
- Added prime-number related functions in ["number_theory.h"]:
  - `qpp::factors()` - prime factors of positive integer
  - `qpp::isprime()` - primality test

# Version 0.1 - 21 December 2014

- First stable version, initial public release
