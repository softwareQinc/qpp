"""
Python 3 bindings for Quantum++ (https://github.com/softwareQinc/qpp)
"""
from __future__ import annotations
import collections.abc
import numpy
import numpy.typing
import typing
from . import codes
from . import gates
from . import qasm
from . import random_devices
from . import states
__all__: list[str] = ['BitCircuit', 'CondWhile', 'DynamicBitset', 'QCircuit', 'QDensityDummyEngine', 'QDensityEngine', 'QDensityNoisyEngine', 'QEngine', 'QFT', 'QKetDummyEngine', 'QKetEngine', 'QKetNoisyEngine', 'QNoisyEngine', 'QubitAmplitudeDampingNoise', 'QubitBitFlipNoise', 'QubitBitPhaseFlipNoise', 'QubitDepolarizingNoise', 'QubitPhaseDampingNoise', 'QubitPhaseFlipNoise', 'QuditDepolarizingNoise', 'TFQ', 'absm', 'abssq', 'adjoint', 'anticomm', 'apply', 'applyCTRL', 'applyCTRL_fan', 'applyQFT', 'applyTFQ', 'avg', 'bernoulli', 'bloch2rho', 'choi2kraus', 'choi2super', 'codes', 'comm', 'complement', 'compose_CTRL_circuit', 'compose_circuit', 'compperm', 'concurrence', 'conjugate', 'const_proxy_to_engine_dits', 'contfrac2x', 'convergents', 'cor', 'cosm', 'couple_circuit_left', 'couple_circuit_right', 'cov', 'det', 'dirac', 'dirac_t', 'dirsum', 'dirsumpow', 'discard', 'ee', 'egcd', 'eig', 'entanglement', 'entropy', 'evals', 'evects', 'expm', 'factors', 'funm', 'gates', 'gcd', 'gconcurrence', 'grams', 'hash_eigen', 'heig', 'hevals', 'hevects', 'infty', 'inverse', 'invperm', 'ip', 'isprime', 'kraus2choi', 'kraus2super', 'kron', 'kronpow', 'lcm', 'load_cmat', 'load_rmat', 'logdet', 'logm', 'lognegativity', 'marginalX', 'marginalY', 'measure', 'measure_seq', 'mket', 'modinv', 'modmul', 'modpow', 'mprj', 'multiidx2n', 'n2multiidx', 'negativity', 'norm', 'normalize', 'omega', 'pi', 'powm', 'prj', 'prod', 'proxy_to_engine_dits', 'ptrace', 'ptrace1', 'ptrace2', 'ptranspose', 'qRAM', 'qasm', 'qmutualinfo', 'qpe_circuit', 'rand', 'randH', 'randU', 'randV', 'randidx', 'randket', 'randkraus', 'randn', 'random_circuit_count', 'random_circuit_depth', 'random_devices', 'randperm', 'randprime', 'randprob', 'randrho', 'renyi', 'replicate', 'reset', 'reshape', 'rho2bloch', 'rho2pure', 'sample', 'save', 'schatten', 'schmidt', 'schmidtA', 'schmidtB', 'schmidtcoeffs', 'schmidtprobs', 'set_prng_seed', 'sigma', 'sinm', 'spectralpowm', 'sqrtm', 'states', 'sum', 'super2choi', 'super2kraus', 'svals', 'svd', 'svdU', 'svdV', 'syspermute', 'trace', 'transpose', 'tsallis', 'uniform', 'var', 'x2contfrac', 'zket2dits']
class BitCircuit(DynamicBitset):
    __hash__: typing.ClassVar[None] = None
    def CNOT(self, ctrl: typing.SupportsInt | typing.SupportsIndex, target: typing.SupportsInt | typing.SupportsIndex) -> BitCircuit:
        """
        Controlled-NOT gate
        """
    def FRED(self, i: typing.SupportsInt | typing.SupportsIndex, j: typing.SupportsInt | typing.SupportsIndex, k: typing.SupportsInt | typing.SupportsIndex) -> BitCircuit:
        """
         Fredkin gate (Controlled-SWAP)
        """
    def NOT(self, i: typing.SupportsInt | typing.SupportsIndex) -> BitCircuit:
        """
        NOT gate (bit flip)
        """
    def SWAP(self, i: typing.SupportsInt | typing.SupportsIndex, j: typing.SupportsInt | typing.SupportsIndex) -> BitCircuit:
        """
        Swap gate
        """
    def TOF(self, i: typing.SupportsInt | typing.SupportsIndex, j: typing.SupportsInt | typing.SupportsIndex, k: typing.SupportsInt | typing.SupportsIndex) -> BitCircuit:
        """
        Toffoli gate
        """
    def X(self, i: typing.SupportsInt | typing.SupportsIndex) -> BitCircuit:
        """
        NOT gate (bit flip)
        """
    def __copy__(self) -> BitCircuit:
        ...
    def __deepcopy__(self, arg0: dict) -> BitCircuit:
        ...
    def __eq__(self, arg0: BitCircuit) -> bool:
        ...
    @typing.overload
    def __init__(self, n: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @typing.overload
    def __init__(self, str: str, zero: str = '0', one: str = '1') -> None:
        ...
    @typing.overload
    def __init__(self, arg0: DynamicBitset) -> None:
        ...
    def __ne__(self, arg0: BitCircuit) -> bool:
        ...
    def __repr__(self) -> str:
        ...
    def get_gate_count(self, name: str | None = None) -> int:
        """
        (Total) Bit circuit gate count. Possible names are NOT (X), CNOT, SWAP, TOF, FRED
        """
    def get_gate_depth(self, name: str | None = None) -> int:
        """
        (Total) Bit circuit gate depth. Possible names are NOT (X), CNOT, SWAP, TOF, FRED
        """
    def reset(self) -> BitCircuit:
        """
        Resets the circuit to all-zero, clears all gates
        """
    def to_JSON(self, enclosed_in_curly_brackets: bool = True) -> str:
        """
        Displays the bit circuit in JSON format
        """
    def to_string(self, zero: str = '0', one: str = '1') -> str:
        """
        String representation
        """
class CondWhile:
    def __enter__(self) -> typing.Any:
        ...
    def __exit__(self, arg0: type | None, arg1: typing.Any | None, arg2: typing.Any | None) -> bool:
        ...
class DynamicBitset:
    __hash__: typing.ClassVar[None] = None
    def __copy__(self) -> DynamicBitset:
        ...
    def __deepcopy__(self, arg0: dict) -> DynamicBitset:
        ...
    def __eq__(self, arg0: DynamicBitset) -> bool:
        ...
    @typing.overload
    def __init__(self, n: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @typing.overload
    def __init__(self, str: str, zero: str = '0', one: str = '1') -> None:
        ...
    def __ne__(self, arg0: DynamicBitset) -> bool:
        ...
    def __repr__(self) -> str:
        ...
    def __sub__(self, arg0: DynamicBitset) -> int:
        """
        Number of places the two bitsets differ (Hamming distance)
        """
    def all(self) -> bool:
        """
        True if all of the bits are set
        """
    def any(self) -> bool:
        """
        True if any of the bits is set
        """
    def count(self) -> int:
        """
        Number of bits set to one in the bitset (Hamming weight)
        """
    @typing.overload
    def flip(self) -> DynamicBitset:
        """
        Flips all bits
        """
    @typing.overload
    def flip(self, pos: typing.SupportsInt | typing.SupportsIndex) -> DynamicBitset:
        """
        Flips the bit at position pos
        """
    def get(self, pos: typing.SupportsInt | typing.SupportsIndex) -> bool:
        """
        The value of the bit at position pos
        """
    def none(self) -> bool:
        """
        True if none of the bits are set
        """
    @typing.overload
    def rand(self, p: typing.SupportsFloat | typing.SupportsIndex = 0.5) -> DynamicBitset:
        """
        Sets all bits according to a Bernoulli(p) distribution
        """
    @typing.overload
    def rand(self, pos: typing.SupportsInt | typing.SupportsIndex, p: typing.SupportsFloat | typing.SupportsIndex = 0.5) -> DynamicBitset:
        """
        Sets the bit at position pos according to a Bernoulli(p) distribution
        """
    @typing.overload
    def reset(self) -> DynamicBitset:
        """
        Sets all bits to false
        """
    @typing.overload
    def reset(self, pos: typing.SupportsInt | typing.SupportsIndex) -> DynamicBitset:
        """
        Sets the bit at position pos to false
        """
    @typing.overload
    def set(self) -> DynamicBitset:
        """
        Sets all bits to true
        """
    @typing.overload
    def set(self, pos: typing.SupportsInt | typing.SupportsIndex, value: bool = True) -> DynamicBitset:
        """
        Sets the bit at position pos
        """
    def size(self) -> int:
        """
        Number of bits stored in the bitset
        """
    def storage_size(self) -> int:
        """
        Size of the underlying storage space (in units of qpp::DynamicBitset::value_type, unsigned int by default)
        """
    def to_string(self, zero: str = '0', one: str = '1') -> str:
        """
        String representation
        """
class QCircuit:
    class Resources:
        def __repr__(self) -> str:
            ...
        @property
        def d(self) -> int:
            ...
        @property
        def gate_count(self) -> int:
            ...
        @property
        def gate_depth(self) -> int:
            ...
        @property
        def measurement_count(self) -> int:
            ...
        @property
        def measurement_depth(self) -> int:
            ...
        @property
        def name(self) -> str | None:
            ...
        @property
        def nc(self) -> int:
            ...
        @property
        def nq(self) -> int:
            ...
        @property
        def step_count(self) -> int:
            ...
        @property
        def total_depth(self) -> int:
            ...
    __hash__: typing.ClassVar[None] = None
    @staticmethod
    def is_CTRL(gate_step: qpp::internal::QCircuitGateStep) -> bool:
        """
        True if the gate step is a controlled gate, false otherwise
        """
    @staticmethod
    def is_cCTRL(gate_step: qpp::internal::QCircuitGateStep) -> bool:
        """
        True if the gate step is a classically-controlled gate, false otherwise
        """
    @typing.overload
    def CTRL(self, U: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], ctrl: typing.SupportsInt | typing.SupportsIndex, target: typing.SupportsInt | typing.SupportsIndex, shift: typing.SupportsInt | typing.SupportsIndex | None = None, name: str | None = None) -> QCircuit:
        """
        Applies the single qudit controlled gate U with control qudit ctrl and target qudit target, i.e., CTRL-U
        """
    @typing.overload
    def CTRL(self, U: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], ctrl: typing.SupportsInt | typing.SupportsIndex, target: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], shift: typing.SupportsInt | typing.SupportsIndex | None = None, name: str | None = None) -> QCircuit:
        """
        Jointly applies the single qudit controlled gate U with control qudit ctrl on the qudit indexes specified by target, i.e., CTRL-U_{joint}
        """
    @typing.overload
    def CTRL(self, U: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], ctrl: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], target: typing.SupportsInt | typing.SupportsIndex, shift: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex] | None = None, name: str | None = None) -> QCircuit:
        """
        Applies the multiple-qudit controlled gate U with multiple control qudits listed in ctrl on the target qudit specified by target, i.e., CTRL-CTRL-...-CTRL-U
        """
    @typing.overload
    def CTRL(self, U: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], ctrl: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], target: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], shift: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex] | None = None, name: str | None = None) -> QCircuit:
        """
        Jointly applies the multiple-qudit controlled gate U with multiple control qudits listed in ctrl on the qudit indexes specified by target, i.e., CTRL-CTRL-...-CTRL-U_{joint}
        """
    @typing.overload
    def CTRL_fan(self, U: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], ctrl: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], target: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], shift: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex] | None = None, name: str | None = None) -> QCircuit:
        """
        Applies the single qudit controlled gate U with multiple control qudits listed in ctrl on every qudit listed in target, i.e., CTRL-CTRL-...-CTRL-U-U-...-U
        """
    @typing.overload
    def CTRL_fan(self, U: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], ctrl: typing.SupportsInt | typing.SupportsIndex, target: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], shift: typing.SupportsInt | typing.SupportsIndex | None = None, name: str | None = None) -> QCircuit:
        """
        Applies the single qudit controlled gate U with control qudit ctrl on every qudit listed in target, i.e., CTRL-U-U-...-U
        """
    @typing.overload
    def QFT(self, target: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], swap: bool = True) -> QCircuit:
        """
        Applies the quantum Fourier transform on the qudit indexes specified by target
        """
    @typing.overload
    def QFT(self, swap: bool = True) -> QCircuit:
        """
        Applies the quantum Fourier transform on all of remaining non-measured qudits
        """
    @typing.overload
    def TFQ(self, swap: bool = True) -> QCircuit:
        """
        Applies the inverse quantum Fourier transform on all of remaining non-measured qudits
        """
    @typing.overload
    def TFQ(self, target: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], swap: bool = True) -> QCircuit:
        """
        Applies the inverse quantum Fourier transform on the qudit indexes specified by target
        """
    def __copy__(self) -> QCircuit:
        ...
    def __deepcopy__(self, arg0: dict) -> QCircuit:
        ...
    def __eq__(self, arg0: QCircuit) -> bool:
        ...
    @typing.overload
    def __init__(self, nq: typing.SupportsInt | typing.SupportsIndex = 1, nc: typing.SupportsInt | typing.SupportsIndex = 0, d: typing.SupportsInt | typing.SupportsIndex = 2, name: str | None = None) -> None:
        ...
    @typing.overload
    def __init__(self, arg0: QCircuit) -> None:
        ...
    def __ne__(self, arg0: QCircuit) -> bool:
        ...
    def __repr__(self) -> str:
        ...
    @typing.overload
    def add_dit(self, n: typing.SupportsInt | typing.SupportsIndex = 1) -> QCircuit:
        """
        Adds n additional classical dits after the last qudit
        """
    @typing.overload
    def add_dit(self, n: typing.SupportsInt | typing.SupportsIndex, pos: typing.SupportsInt | typing.SupportsIndex) -> QCircuit:
        """
        Adds n additional classical dits before qudit pos
        """
    @typing.overload
    def add_qudit(self, n: typing.SupportsInt | typing.SupportsIndex = 1) -> QCircuit:
        """
        Adds n additional qudits after the last qudit
        """
    @typing.overload
    def add_qudit(self, n: typing.SupportsInt | typing.SupportsIndex, pos: typing.SupportsInt | typing.SupportsIndex) -> QCircuit:
        """
        Adds n additional qudits before qudit pos
        """
    def adjoint(self) -> QCircuit:
        """
        Adjoint quantum circuit description, in place
        """
    @typing.overload
    def cCTRL(self, U: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], ctrl_dit: typing.SupportsInt | typing.SupportsIndex, target: typing.SupportsInt | typing.SupportsIndex, shift: typing.SupportsInt | typing.SupportsIndex | None = None, name: str | None = None) -> QCircuit:
        """
        Applies the single qudit controlled gate U with classical control dit ctrl and target qudit target, i.e., cCTRL-U
        """
    @typing.overload
    def cCTRL(self, U: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], ctrl_dit: typing.SupportsInt | typing.SupportsIndex, target: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], shift: typing.SupportsInt | typing.SupportsIndex | None = None, name: str | None = None) -> QCircuit:
        """
        Jointly applies the single qudit controlled gate U with classical control dit ctrl on the qudit indexes specified by target, i.e., cCTRL-U_{joint}
        """
    @typing.overload
    def cCTRL(self, U: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], ctrl_dits: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], target: typing.SupportsInt | typing.SupportsIndex, shift: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex] | None = None, name: str | None = None) -> QCircuit:
        """
        Applies the single qudit controlled gate U with multiple classical control dits listed in ctrl on the target qudit target, i.e., cCTRL-cCTRL-...-cCTRL-U
        """
    @typing.overload
    def cCTRL(self, U: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], ctrl_dits: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], target: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], shift: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex] | None = None, name: str | None = None) -> QCircuit:
        """
        Jointly applies the multiple-qudit controlled gate U with multiple classical control dits listed in ctrl on the qudit indexes specified by target, i.e., cCTRL-cCTRL-...-cCTRL-U_{joint}
        """
    @typing.overload
    def cCTRL_fan(self, U: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], ctrl_dits: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], target: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], shift: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex] | None = None, name: str | None = None) -> QCircuit:
        """
        Applies the single qudit controlled gate U with multiple classical control dits listed in ctrl on every qudit listed in target, i.e., cCTRL-cCTRL-...-cCTRL-U-U-...-U
        """
    @typing.overload
    def cCTRL_fan(self, U: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], ctrl_dit: typing.SupportsInt | typing.SupportsIndex, target: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], shift: typing.SupportsInt | typing.SupportsIndex | None = None, name: str | None = None) -> QCircuit:
        """
        Applies the single qudit controlled gate U with classical control dit ctrl on every qudit listed in target, i.e.,cCTRL-U-U-...-U
        """
    def compose_CTRL_circuit(self, ctrl: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], qc_target: QCircuit, pos_qudit: typing.SupportsInt | typing.SupportsIndex, shift: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex] | None = None, pos_dit: typing.SupportsInt | typing.SupportsIndex | None = None) -> QCircuit:
        """
        Composes (appends) a controlled quantum circuit description to the end of the current one, with the current instance acting as the control
        """
    def compose_circuit(self, other: QCircuit, pos_qudit: typing.SupportsInt | typing.SupportsIndex, pos_dit: typing.SupportsInt | typing.SupportsIndex | None = None) -> QCircuit:
        """
        Composes (appends) a quantum circuit description to the end of the current one
        """
    def compress(self, compress_dits: bool = False) -> QCircuit:
        """
        Removes all clean qudits and relabels the rest of the qudits accordingly
        """
    def cond_else(self) -> QCircuit:
        """
        Adds conditional else
        """
    def cond_end(self) -> QCircuit:
        """
        Adds conditional end
        """
    def cond_if(self, pred: collections.abc.Callable[[const_proxy_to_engine_dits], bool]) -> QCircuit:
        """
        Adds conditional if
        """
    def cond_while(self, pred: collections.abc.Callable[[const_proxy_to_engine_dits], bool]) -> QCircuit:
        """
        Adds conditional while
        """
    def cond_while_ctx(self, pred: collections.abc.Callable[[const_proxy_to_engine_dits], bool]) -> CondWhile:
        """
        Adds conditional while context
        """
    def couple_circuit_left(self, other: QCircuit, target: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], pos_dit: typing.SupportsInt | typing.SupportsIndex | None = None) -> QCircuit:
        """
        Couples (in place) a quantum circuit description to the current one, placed at the left (beginning) of the current one
        """
    def couple_circuit_right(self, other: QCircuit, target: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], pos_dit: typing.SupportsInt | typing.SupportsIndex | None = None) -> QCircuit:
        """
        Couples (in place) a quantum circuit description to the current one, placed at the right (end) of the current one
        """
    @typing.overload
    def discard(self, target: typing.SupportsInt | typing.SupportsIndex, name: str | None = 'discard') -> QCircuit:
        """
        Discards single qudit by measuring it destructively in the Z-basis
        """
    @typing.overload
    def discard(self, target: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], name: str | None = 'discard') -> QCircuit:
        """
        Discards multiple qudits by measuring them destructively in the Z-basis
        """
    @typing.overload
    def gate(self, U: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], i: typing.SupportsInt | typing.SupportsIndex, name: str | None = None) -> QCircuit:
        """
        Applies the single qudit gate U on single qudit i
        """
    @typing.overload
    def gate(self, U: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], i: typing.SupportsInt | typing.SupportsIndex, j: typing.SupportsInt | typing.SupportsIndex, name: str | None = None) -> QCircuit:
        """
        Applies the two qudit gate U on qudits i and j
        """
    @typing.overload
    def gate(self, U: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], i: typing.SupportsInt | typing.SupportsIndex, j: typing.SupportsInt | typing.SupportsIndex, k: typing.SupportsInt | typing.SupportsIndex, name: str | None = None) -> QCircuit:
        """
        Applies the three qudit gate U on qudits i, j and k
        """
    @typing.overload
    def gate(self, U: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], target: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], name: str | None = None) -> QCircuit:
        """
        Jointly applies the multiple-qudit gate U on the qudit indexes specified by target
        """
    @typing.overload
    def gate_fan(self, U: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], target: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], name: str | None = None) -> QCircuit:
        """
        Applies the single qudit gate U on every qudit listed in target
        """
    @typing.overload
    def gate_fan(self, U: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], name: str | None = None) -> QCircuit:
        """
        Applies the single qudit gate U on all of the remaining non-measured qudits
        """
    def get_clean_dits(self) -> list[int]:
        """
        Vector of clean classical dits
        """
    def get_clean_qudits(self) -> list[int]:
        """
        Vector of clean qudits
        """
    def get_d(self) -> int:
        """
        Qudit dimension
        """
    def get_depth(self) -> int:
        """
        Quantum circuit description total depth
        """
    def get_dirty_dits(self) -> list[int]:
        """
        Vector of dirty classical dits
        """
    def get_dirty_qudits(self) -> list[int]:
        """
        Vector of dirty qudits
        """
    def get_gate_count(self, U: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"] | None = None) -> int:
        """
        (Total) Gate count
        """
    def get_gate_depth(self, U: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"] | None = None) -> int:
        """
        (Total) Gate depth
        """
    def get_measured_d(self) -> list[int]:
        """
        Vector of already measured qudit indexes
        """
    def get_measured_nd(self) -> list[int]:
        """
        Vector of already measured non-destructively qudit indexes
        """
    def get_measurement_count(self, V: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"] | None = None) -> int:
        """
        (Total) Measurement count
        """
    def get_measurement_depth(self, V: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"] | None = None) -> int:
        """
        (Total) Measurement depth
        """
    def get_measurement_dits(self) -> list[int]:
        """
        Vector of classical dits that were used to store results of measurements (either destructive or non-destructive)
        """
    def get_name(self) -> str | None:
        """
        Description name
        """
    def get_nc(self) -> int:
        """
        Number of classical dits
        """
    def get_non_measured_d(self) -> list[int]:
        """
        Non-measured qudit indexes
        """
    def get_nop_count(self) -> int:
        """
        No-op count
        """
    def get_nq(self) -> int:
        """
        Number of qudits
        """
    def get_resources(self) -> qpp::internal::QCircuitResources:
        """
        Quantum circuit resources
        """
    def get_step_count(self) -> int:
        """
        Total (gates + measurements) count
        """
    def has_measurements(self) -> bool:
        """
        True if the quantum circuit description contains measurements, false otherwise
        """
    def has_runtime_steps(self) -> bool:
        """
        True if the quantum circuit description contains runtime steps, false otherwise
        """
    def is_clean_dit(self, i: typing.SupportsInt | typing.SupportsIndex) -> bool:
        """
        Whether classical dit i in the quantum circuit description was used before or not
        """
    def is_clean_qudit(self, i: typing.SupportsInt | typing.SupportsIndex) -> bool:
        """
        Whether qudit i in the quantum circuit description was used before or not
        """
    def is_measurement_dit(self, i: typing.SupportsInt | typing.SupportsIndex) -> bool:
        """
        Whether classical dit i in the quantum circuit description was used to store the result of a measurement (either destructive or non-destructive)
        """
    def kron(self, qc: QCircuit) -> QCircuit:
        """
        Kronecker product with another quantum circuit description, in place
        """
    @typing.overload
    def measure(self, target: typing.SupportsInt | typing.SupportsIndex, c_reg: typing.SupportsInt | typing.SupportsIndex, destructive: bool = True, name: str | None = 'mZ') -> QCircuit:
        """
        Z measurement of single qudit
        """
    @typing.overload
    def measure(self, target: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], c_reg: typing.SupportsInt | typing.SupportsIndex = 0, destructive: bool = True, name: str | None = 'mZ') -> QCircuit:
        """
        Z measurement of multiple qudits
        """
    @typing.overload
    def measureV(self, V: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], target: typing.SupportsInt | typing.SupportsIndex, c_reg: typing.SupportsInt | typing.SupportsIndex, destructive: bool = True, name: str | None = None) -> QCircuit:
        """
        Measurement of single qudit in the orthonormal basis specified by the columns of matrix V
        """
    @typing.overload
    def measureV(self, V: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], target: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], c_reg: typing.SupportsInt | typing.SupportsIndex, destructive: bool = True, name: str | None = None) -> QCircuit:
        """
        Measurement of multiple qudits in the orthonormal basis specified by the columns of matrix V
        """
    def measure_all(self, c_reg: typing.SupportsInt | typing.SupportsIndex = 0, destructive: bool = True, name: str | None = 'mZ') -> QCircuit:
        """
        Z measurement of all qudits
        """
    def nop(self) -> QCircuit:
        """
        No operation (no-op)
        """
    @typing.overload
    def post_select(self, target: typing.SupportsInt | typing.SupportsIndex, ps_val: typing.SupportsInt | typing.SupportsIndex, c_reg: typing.SupportsInt | typing.SupportsIndex, destructive: bool = True, name: str | None = None) -> QCircuit:
        """
        Z post-selection of single qudit
        """
    @typing.overload
    def post_select(self, target: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], ps_vals: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], c_reg: typing.SupportsInt | typing.SupportsIndex, destructive: bool = True, name: str | None = None) -> QCircuit:
        """
        Z post-selection of multiple qudits
        """
    @typing.overload
    def post_selectV(self, V: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], target: typing.SupportsInt | typing.SupportsIndex, ps_val: typing.SupportsInt | typing.SupportsIndex, c_reg: typing.SupportsInt | typing.SupportsIndex, destructive: bool = True, name: str | None = None) -> QCircuit:
        """
        Post-selection of single qudit in the orthonormal basis specified by the columns of matrix V
        """
    @typing.overload
    def post_selectV(self, V: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], target: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], ps_val: typing.SupportsInt | typing.SupportsIndex, c_reg: typing.SupportsInt | typing.SupportsIndex, destructive: bool = True, name: str | None = None) -> QCircuit:
        """
        Post-selection of multiple qudits in the orthonormal basis specified by the columns of matrix V
        """
    def remove_clean_dit(self, target: typing.SupportsInt | typing.SupportsIndex) -> QCircuit:
        """
        Removes clean classical dit and relabels the rest of the classical dits accordingly
        """
    def remove_clean_dits(self, target: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> QCircuit:
        """
        Removes clean classical dits and relabels the rest of the classical dits accordingly
        """
    def remove_clean_qudit(self, target: typing.SupportsInt | typing.SupportsIndex) -> QCircuit:
        """
        Removes clean qudit and relabels the rest of the qudits accordingly
        """
    def remove_clean_qudits(self, target: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> QCircuit:
        """
        Removes clean qudits and relabels the rest of the qudits accordingly
        """
    def removes_qudits(self) -> bool:
        """
        True if the quantum circuit description contains measurements that remove qudits, false otherwise
        """
    def replicate(self, n: typing.SupportsInt | typing.SupportsIndex) -> QCircuit:
        """
        Replicates the quantum circuit description, in place
        """
    @typing.overload
    def reset(self, target: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], name: str | None = 'reset') -> QCircuit:
        """
        Reset multiple qudits
        """
    @typing.overload
    def reset(self, target: typing.SupportsInt | typing.SupportsIndex, name: str | None = 'reset') -> QCircuit:
        """
        Reset single qudit
        """
    def set_dits_runtime(self, functor: collections.abc.Callable[[proxy_to_engine_dits], None]) -> QCircuit:
        """
        Set dits at runtime
        """
    def set_name(self, name: str) -> QCircuit:
        """
        Sets name
        """
    def to_JSON(self, enclosed_in_curly_brackets: bool = True) -> str:
        """
        Displays the quantum circuit description in JSON format
        """
    def validate_conditionals(self) -> bool:
        """
        True if valid conditionals, false otherwise
        """
    def was_measured_d(self, i: typing.SupportsInt | typing.SupportsIndex) -> bool:
        """
        Whether qudit i was already measured
        """
    def was_measured_nd(self, i: typing.SupportsInt | typing.SupportsIndex) -> bool:
        """
        Whether qudit i was already measured non-destructively
        """
class QubitAmplitudeDampingNoise:
    def __init__(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    def get_Ks(self) -> list[typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]]:
        """
        Vector of noise operators
        """
    def get_d(self: qpp::NoiseBase<qpp::NoiseType::StateDependent>) -> int:
        """
        Qudit dimension
        """
    def get_last_K(self) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]:
        """
        Last occurring noise element
        """
    def get_last_idx(self) -> int:
        """
        Index of the last occurring noise element
        """
    def get_last_p(self) -> float:
        """
        Probability of the last occurring noise element
        """
    def get_probs(self) -> list[float]:
        """
        Vector of probabilities corresponding to each noise operator
        """
class QubitBitFlipNoise:
    def __init__(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    def get_Ks(self) -> list[typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]]:
        """
        Vector of noise operators
        """
    def get_d(self: qpp::NoiseBase<qpp::NoiseType::StateIndependent>) -> int:
        """
        Qudit dimension
        """
    def get_last_K(self) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]:
        """
        Last occurring noise element
        """
    def get_last_idx(self) -> int:
        """
        Index of the last occurring noise element
        """
    def get_last_p(self) -> float:
        """
        Probability of the last occurring noise element
        """
    def get_probs(self) -> list[float]:
        """
        Vector of probabilities corresponding to each noise operator
        """
class QubitBitPhaseFlipNoise:
    def __init__(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    def get_Ks(self) -> list[typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]]:
        """
        Vector of noise operators
        """
    def get_d(self: qpp::NoiseBase<qpp::NoiseType::StateIndependent>) -> int:
        """
        Qudit dimension
        """
    def get_last_K(self) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]:
        """
        Last occurring noise element
        """
    def get_last_idx(self) -> int:
        """
        Index of the last occurring noise element
        """
    def get_last_p(self) -> float:
        """
        Probability of the last occurring noise element
        """
    def get_probs(self) -> list[float]:
        """
        Vector of probabilities corresponding to each noise operator
        """
class QubitDepolarizingNoise:
    def __init__(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    def get_Ks(self) -> list[typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]]:
        """
        Vector of noise operators
        """
    def get_d(self: qpp::NoiseBase<qpp::NoiseType::StateIndependent>) -> int:
        """
        Qudit dimension
        """
    def get_last_K(self) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]:
        """
        Last occurring noise element
        """
    def get_last_idx(self) -> int:
        """
        Index of the last occurring noise element
        """
    def get_last_p(self) -> float:
        """
        Probability of the last occurring noise element
        """
    def get_probs(self) -> list[float]:
        """
        Vector of probabilities corresponding to each noise operator
        """
class QubitPhaseDampingNoise:
    def __init__(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    def get_Ks(self) -> list[typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]]:
        """
        Vector of noise operators
        """
    def get_d(self: qpp::NoiseBase<qpp::NoiseType::StateDependent>) -> int:
        """
        Qudit dimension
        """
    def get_last_K(self) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]:
        """
        Last occurring noise element
        """
    def get_last_idx(self) -> int:
        """
        Index of the last occurring noise element
        """
    def get_last_p(self) -> float:
        """
        Probability of the last occurring noise element
        """
    def get_probs(self) -> list[float]:
        """
        Vector of probabilities corresponding to each noise operator
        """
class QubitPhaseFlipNoise:
    def __init__(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    def get_Ks(self) -> list[typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]]:
        """
        Vector of noise operators
        """
    def get_d(self: qpp::NoiseBase<qpp::NoiseType::StateIndependent>) -> int:
        """
        Qudit dimension
        """
    def get_last_K(self) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]:
        """
        Last occurring noise element
        """
    def get_last_idx(self) -> int:
        """
        Index of the last occurring noise element
        """
    def get_last_p(self) -> float:
        """
        Probability of the last occurring noise element
        """
    def get_probs(self) -> list[float]:
        """
        Vector of probabilities corresponding to each noise operator
        """
class QuditDepolarizingNoise:
    def __init__(self, arg0: typing.SupportsFloat | typing.SupportsIndex, arg1: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    def get_Ks(self) -> list[typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]]:
        """
        Vector of noise operators
        """
    def get_d(self: qpp::NoiseBase<qpp::NoiseType::StateIndependent>) -> int:
        """
        Qudit dimension
        """
    def get_last_K(self) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]:
        """
        Last occurring noise element
        """
    def get_last_idx(self) -> int:
        """
        Index of the last occurring noise element
        """
    def get_last_p(self) -> float:
        """
        Probability of the last occurring noise element
        """
    def get_probs(self) -> list[float]:
        """
        Vector of probabilities corresponding to each noise operator
        """
class _QDensityDummyEngine:
    def __copy__(self) -> _QDensityDummyEngine:
        ...
    def __deepcopy__(self, arg0: dict) -> _QDensityDummyEngine:
        ...
    def __init__(self, arg0: QCircuit) -> None:
        ...
    def __repr__(self) -> str:
        ...
    def execute(self, reps: typing.SupportsInt | typing.SupportsIndex = 1) -> qpp::QBaseEngine<Eigen::Matrix<std::__1::complex<double>, -1, -1, 0, -1, -1>, qpp::QCircuit>:
        """
        Executes the entire quantum circuit description
        """
    def get_circuit(self) -> QCircuit:
        """
        Underlying quantum circuit description
        """
    def get_state(self) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]:
        """
        Underlying quantum state
        """
    def set_state(self, state: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"]) -> qpp::QBaseEngine<Eigen::Matrix<std::__1::complex<double>, -1, -1, 0, -1, -1>, qpp::QCircuit>:
        """
        Sets the underlying quantum state
        """
    def to_JSON(self, enclosed_in_curly_brackets: bool = True) -> str:
        """
        State of the engine in JSON format
        """
    def traits_get_name(self) -> str:
        """
        Engine name
        """
    def traits_is_noisy(self) -> bool:
        """
        Noisy engine?
        """
    def traits_is_pure(self) -> bool:
        """
        Pure state engine?
        """
class _QDensityEngine:
    def __copy__(self) -> _QDensityEngine:
        ...
    def __deepcopy__(self, arg0: dict) -> _QDensityEngine:
        ...
    def __init__(self, arg0: QCircuit) -> None:
        ...
    def __repr__(self) -> str:
        ...
    def execute(self, reps: typing.SupportsInt | typing.SupportsIndex = 1) -> _QDensityEngine:
        """
        Executes the entire quantum circuit description
        """
    def get_circuit(self) -> QCircuit:
        """
        Underlying quantum circuit description
        """
    def get_dit(self, i: typing.SupportsInt | typing.SupportsIndex) -> int:
        """
        Underlying classical dit at position i
        """
    def get_dits(self) -> list[int]:
        """
        Underlying classical dits
        """
    def get_ensure_post_selection(self) -> bool:
        """
        True if post-selection is enforced (must succeed), false otherwise
        """
    def get_max_post_selection_reps(self) -> int:
        """
        Maximum number of repetitions of a cirucit post-selection step until success
        """
    def get_measured_d(self) -> list[int]:
        """
        Vector of already destructively measured qudit indexes
        """
    def get_non_measured_d(self) -> list[int]:
        """
        Vector of qudit indexes that were not measured destructively
        """
    def get_probs(self) -> list[float]:
        """
        Underlying measurement outcome probabilities
        """
    def get_state(self) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]:
        """
        Underlying quantum state
        """
    def get_stats(self) -> dict[str, int]:
        """
        Measurement statistics for multiple runs
        """
    def get_stats_to_JSON(self) -> str:
        """
        Measurement statistics for multiple runs, JSON format
        """
    def post_select_ok(self) -> bool:
        """
        True if post-selection was successful (or absent), false otherwise
        """
    def reset(self, qstate: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"] | None = None, reset_stats: bool = True) -> _QDensityEngine:
        """
        Resets the engine
        """
    def reset_stats(self) -> _QDensityEngine:
        """
        Resets the collected measurement statistics hash table
        """
    def set_dit(self, i: typing.SupportsInt | typing.SupportsIndex, value: typing.SupportsInt | typing.SupportsIndex) -> _QDensityEngine:
        """
        Sets the classical dit at position i
        """
    def set_dits(self, dits: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> _QDensityEngine:
        """
        Set the classical dits
        """
    def set_ensure_post_selection(self, val: bool) -> _QDensityEngine:
        """
        Enforces post-selection (must succeed)
        """
    def set_max_post_selection_reps(self, max_post_selection_reps: typing.SupportsInt | typing.SupportsIndex) -> _QDensityEngine:
        """
        Sets the maximum number of repetitions of a circuit post-selection step until success
        """
    def set_state(self, state: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"]) -> _QDensityEngine:
        """
        Sets the underlying quantum state
        """
    def to_JSON(self, enclosed_in_curly_brackets: bool = True) -> str:
        """
        State of the engine in JSON format
        """
    def traits_get_name(self) -> str:
        """
        Engine name
        """
    def traits_is_noisy(self) -> bool:
        """
        Noisy engine?
        """
    def traits_is_pure(self) -> bool:
        """
        Pure state engine?
        """
    def was_measured_d(self, i: typing.SupportsInt | typing.SupportsIndex) -> bool:
        """
        Whether qudit i was already measured destructively
        """
class _QDensityNoisyEngine_QubitAmplitudeDampingNoise(_QDensityEngine):
    def __copy__(self) -> _QDensityNoisyEngine_QubitAmplitudeDampingNoise:
        ...
    def __deepcopy__(self, arg0: dict) -> _QDensityNoisyEngine_QubitAmplitudeDampingNoise:
        ...
    def __init__(self, arg0: QCircuit, arg1: QubitAmplitudeDampingNoise) -> None:
        ...
    def get_noise_results(self) -> list[list[int]]:
        """
        Vector of noise results obtained before every step in the circuit
        """
class _QDensityNoisyEngine_QubitBitFlipNoise(_QDensityEngine):
    def __copy__(self) -> _QDensityNoisyEngine_QubitBitFlipNoise:
        ...
    def __deepcopy__(self, arg0: dict) -> _QDensityNoisyEngine_QubitBitFlipNoise:
        ...
    def __init__(self, arg0: QCircuit, arg1: QubitBitFlipNoise) -> None:
        ...
    def get_noise_results(self) -> list[list[int]]:
        """
        Vector of noise results obtained before every step in the circuit
        """
class _QDensityNoisyEngine_QubitBitPhaseFlipNoise(_QDensityEngine):
    def __copy__(self) -> _QDensityNoisyEngine_QubitBitPhaseFlipNoise:
        ...
    def __deepcopy__(self, arg0: dict) -> _QDensityNoisyEngine_QubitBitPhaseFlipNoise:
        ...
    def __init__(self, arg0: QCircuit, arg1: QubitBitPhaseFlipNoise) -> None:
        ...
    def get_noise_results(self) -> list[list[int]]:
        """
        Vector of noise results obtained before every step in the circuit
        """
class _QDensityNoisyEngine_QubitDepolarizingNoise(_QDensityEngine):
    def __copy__(self) -> _QDensityNoisyEngine_QubitDepolarizingNoise:
        ...
    def __deepcopy__(self, arg0: dict) -> _QDensityNoisyEngine_QubitDepolarizingNoise:
        ...
    def __init__(self, arg0: QCircuit, arg1: QubitDepolarizingNoise) -> None:
        ...
    def get_noise_results(self) -> list[list[int]]:
        """
        Vector of noise results obtained before every step in the circuit
        """
class _QDensityNoisyEngine_QubitPhaseDampingNoise(_QDensityEngine):
    def __copy__(self) -> _QDensityNoisyEngine_QubitPhaseDampingNoise:
        ...
    def __deepcopy__(self, arg0: dict) -> _QDensityNoisyEngine_QubitPhaseDampingNoise:
        ...
    def __init__(self, arg0: QCircuit, arg1: QubitPhaseDampingNoise) -> None:
        ...
    def get_noise_results(self) -> list[list[int]]:
        """
        Vector of noise results obtained before every step in the circuit
        """
class _QDensityNoisyEngine_QubitPhaseFlipNoise(_QDensityEngine):
    def __copy__(self) -> _QDensityNoisyEngine_QubitPhaseFlipNoise:
        ...
    def __deepcopy__(self, arg0: dict) -> _QDensityNoisyEngine_QubitPhaseFlipNoise:
        ...
    def __init__(self, arg0: QCircuit, arg1: QubitPhaseFlipNoise) -> None:
        ...
    def get_noise_results(self) -> list[list[int]]:
        """
        Vector of noise results obtained before every step in the circuit
        """
class _QDensityNoisyEngine_QuditDepolarizingNoise(_QDensityEngine):
    def __copy__(self) -> _QDensityNoisyEngine_QuditDepolarizingNoise:
        ...
    def __deepcopy__(self, arg0: dict) -> _QDensityNoisyEngine_QuditDepolarizingNoise:
        ...
    def __init__(self, arg0: QCircuit, arg1: QuditDepolarizingNoise) -> None:
        ...
    def get_noise_results(self) -> list[list[int]]:
        """
        Vector of noise results obtained before every step in the circuit
        """
class _QKetDummyEngine:
    def __copy__(self) -> _QKetDummyEngine:
        ...
    def __deepcopy__(self, arg0: dict) -> _QKetDummyEngine:
        ...
    def __init__(self, arg0: QCircuit) -> None:
        ...
    def __repr__(self) -> str:
        ...
    def execute(self, reps: typing.SupportsInt | typing.SupportsIndex = 1) -> qpp::QBaseEngine<Eigen::Matrix<std::__1::complex<double>, -1, 1, 0, -1, 1>, qpp::QCircuit>:
        """
        Executes the entire quantum circuit description
        """
    def get_circuit(self) -> QCircuit:
        """
        Underlying quantum circuit description
        """
    def get_state(self) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, 1]"]:
        """
        Underlying quantum state
        """
    def set_state(self, state: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, 1]"]) -> qpp::QBaseEngine<Eigen::Matrix<std::__1::complex<double>, -1, 1, 0, -1, 1>, qpp::QCircuit>:
        """
        Sets the underlying quantum state
        """
    def to_JSON(self, enclosed_in_curly_brackets: bool = True) -> str:
        """
        State of the engine in JSON format
        """
    def traits_get_name(self) -> str:
        """
        Engine name
        """
    def traits_is_noisy(self) -> bool:
        """
        Noisy engine?
        """
    def traits_is_pure(self) -> bool:
        """
        Pure state engine?
        """
class _QKetEngine:
    def __copy__(self) -> _QKetEngine:
        ...
    def __deepcopy__(self, arg0: dict) -> _QKetEngine:
        ...
    def __init__(self, arg0: QCircuit) -> None:
        ...
    def __repr__(self) -> str:
        ...
    def execute(self, reps: typing.SupportsInt | typing.SupportsIndex = 1) -> _QKetEngine:
        """
        Executes the entire quantum circuit description
        """
    def get_circuit(self) -> QCircuit:
        """
        Underlying quantum circuit description
        """
    def get_dit(self, i: typing.SupportsInt | typing.SupportsIndex) -> int:
        """
        Underlying classical dit at position i
        """
    def get_dits(self) -> list[int]:
        """
        Underlying classical dits
        """
    def get_ensure_post_selection(self) -> bool:
        """
        True if post-selection is enforced (must succeed), false otherwise
        """
    def get_max_post_selection_reps(self) -> int:
        """
        Maximum number of repetitions of a cirucit post-selection step until success
        """
    def get_measured_d(self) -> list[int]:
        """
        Vector of already destructively measured qudit indexes
        """
    def get_non_measured_d(self) -> list[int]:
        """
        Vector of qudit indexes that were not measured destructively
        """
    def get_probs(self) -> list[float]:
        """
        Underlying measurement outcome probabilities
        """
    def get_state(self) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]:
        """
        Underlying quantum state
        """
    def get_stats(self) -> dict[str, int]:
        """
        Measurement statistics for multiple runs
        """
    def get_stats_to_JSON(self) -> str:
        """
        Measurement statistics for multiple runs, JSON format
        """
    def post_select_ok(self) -> bool:
        """
        True if post-selection was successful (or absent), false otherwise
        """
    def reset(self, qstate: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, 1]"] | None = None, reset_stats: bool = True) -> _QKetEngine:
        """
        Resets the engine
        """
    def reset_stats(self) -> _QKetEngine:
        """
        Resets the collected measurement statistics hash table
        """
    def set_dit(self, i: typing.SupportsInt | typing.SupportsIndex, value: typing.SupportsInt | typing.SupportsIndex) -> _QKetEngine:
        """
        Sets the classical dit at position i
        """
    def set_dits(self, dits: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> _QKetEngine:
        """
        Set the classical dits
        """
    def set_ensure_post_selection(self, val: bool) -> _QKetEngine:
        """
        Enforces post-selection (must succeed)
        """
    def set_max_post_selection_reps(self, max_post_selection_reps: typing.SupportsInt | typing.SupportsIndex) -> _QKetEngine:
        """
        Sets the maximum number of repetitions of a circuit post-selection step until success
        """
    def set_state(self, state: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, 1]"]) -> _QKetEngine:
        """
        Sets the underlying quantum state
        """
    def to_JSON(self, enclosed_in_curly_brackets: bool = True) -> str:
        """
        State of the engine in JSON format
        """
    def traits_get_name(self) -> str:
        """
        Engine name
        """
    def traits_is_noisy(self) -> bool:
        """
        Noisy engine?
        """
    def traits_is_pure(self) -> bool:
        """
        Pure state engine?
        """
    def was_measured_d(self, i: typing.SupportsInt | typing.SupportsIndex) -> bool:
        """
        Whether qudit i was already measured destructively
        """
class _QKetNoisyEngine_QubitAmplitudeDampingNoise(_QKetEngine):
    def __copy__(self) -> _QKetNoisyEngine_QubitAmplitudeDampingNoise:
        ...
    def __deepcopy__(self, arg0: dict) -> _QKetNoisyEngine_QubitAmplitudeDampingNoise:
        ...
    def __init__(self, arg0: QCircuit, arg1: QubitAmplitudeDampingNoise) -> None:
        ...
    def get_noise_results(self) -> list[list[int]]:
        """
        Vector of noise results obtained before every step in the circuit
        """
class _QKetNoisyEngine_QubitBitFlipNoise(_QKetEngine):
    def __copy__(self) -> _QKetNoisyEngine_QubitBitFlipNoise:
        ...
    def __deepcopy__(self, arg0: dict) -> _QKetNoisyEngine_QubitBitFlipNoise:
        ...
    def __init__(self, arg0: QCircuit, arg1: QubitBitFlipNoise) -> None:
        ...
    def get_noise_results(self) -> list[list[int]]:
        """
        Vector of noise results obtained before every step in the circuit
        """
class _QKetNoisyEngine_QubitBitPhaseFlipNoise(_QKetEngine):
    def __copy__(self) -> _QKetNoisyEngine_QubitBitPhaseFlipNoise:
        ...
    def __deepcopy__(self, arg0: dict) -> _QKetNoisyEngine_QubitBitPhaseFlipNoise:
        ...
    def __init__(self, arg0: QCircuit, arg1: QubitBitPhaseFlipNoise) -> None:
        ...
    def get_noise_results(self) -> list[list[int]]:
        """
        Vector of noise results obtained before every step in the circuit
        """
class _QKetNoisyEngine_QubitDepolarizingNoise(_QKetEngine):
    def __copy__(self) -> _QKetNoisyEngine_QubitDepolarizingNoise:
        ...
    def __deepcopy__(self, arg0: dict) -> _QKetNoisyEngine_QubitDepolarizingNoise:
        ...
    def __init__(self, arg0: QCircuit, arg1: QubitDepolarizingNoise) -> None:
        ...
    def get_noise_results(self) -> list[list[int]]:
        """
        Vector of noise results obtained before every step in the circuit
        """
class _QKetNoisyEngine_QubitPhaseDampingNoise(_QKetEngine):
    def __copy__(self) -> _QKetNoisyEngine_QubitPhaseDampingNoise:
        ...
    def __deepcopy__(self, arg0: dict) -> _QKetNoisyEngine_QubitPhaseDampingNoise:
        ...
    def __init__(self, arg0: QCircuit, arg1: QubitPhaseDampingNoise) -> None:
        ...
    def get_noise_results(self) -> list[list[int]]:
        """
        Vector of noise results obtained before every step in the circuit
        """
class _QKetNoisyEngine_QubitPhaseFlipNoise(_QKetEngine):
    def __copy__(self) -> _QKetNoisyEngine_QubitPhaseFlipNoise:
        ...
    def __deepcopy__(self, arg0: dict) -> _QKetNoisyEngine_QubitPhaseFlipNoise:
        ...
    def __init__(self, arg0: QCircuit, arg1: QubitPhaseFlipNoise) -> None:
        ...
    def get_noise_results(self) -> list[list[int]]:
        """
        Vector of noise results obtained before every step in the circuit
        """
class _QKetNoisyEngine_QuditDepolarizingNoise(_QKetEngine):
    def __copy__(self) -> _QKetNoisyEngine_QuditDepolarizingNoise:
        ...
    def __deepcopy__(self, arg0: dict) -> _QKetNoisyEngine_QuditDepolarizingNoise:
        ...
    def __init__(self, arg0: QCircuit, arg1: QuditDepolarizingNoise) -> None:
        ...
    def get_noise_results(self) -> list[list[int]]:
        """
        Vector of noise results obtained before every step in the circuit
        """
class const_proxy_to_engine_dits:
    def __getitem__(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> int:
        ...
class dirac_t:
    __hash__: typing.ClassVar[None] = None
    def __copy__(self) -> dirac_t:
        ...
    def __deepcopy__(self, arg0: dict) -> dirac_t:
        ...
    def __eq__(self, arg0: dirac_t) -> bool:
        ...
    def __ne__(self, arg0: dirac_t) -> bool:
        ...
    def __repr__(self) -> str:
        ...
class proxy_to_engine_dits:
    def __getitem__(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> int:
        ...
    def __setitem__(self, arg0: typing.SupportsInt | typing.SupportsIndex, arg1: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
def QDensityDummyEngine(arg0: QCircuit) -> _QDensityDummyEngine:
    ...
def QDensityEngine(arg0: QCircuit) -> _QDensityEngine:
    ...
@typing.overload
def QDensityNoisyEngine(arg0: QCircuit, arg1: QubitBitFlipNoise) -> _QDensityNoisyEngine_QubitBitFlipNoise:
    ...
@typing.overload
def QDensityNoisyEngine(arg0: QCircuit, arg1: QubitBitPhaseFlipNoise) -> _QDensityNoisyEngine_QubitBitPhaseFlipNoise:
    ...
@typing.overload
def QDensityNoisyEngine(arg0: QCircuit, arg1: QubitDepolarizingNoise) -> _QDensityNoisyEngine_QubitDepolarizingNoise:
    ...
@typing.overload
def QDensityNoisyEngine(arg0: QCircuit, arg1: QubitPhaseFlipNoise) -> _QDensityNoisyEngine_QubitPhaseFlipNoise:
    ...
@typing.overload
def QDensityNoisyEngine(arg0: QCircuit, arg1: QubitAmplitudeDampingNoise) -> _QDensityNoisyEngine_QubitAmplitudeDampingNoise:
    ...
@typing.overload
def QDensityNoisyEngine(arg0: QCircuit, arg1: QubitPhaseDampingNoise) -> _QDensityNoisyEngine_QubitPhaseDampingNoise:
    ...
@typing.overload
def QDensityNoisyEngine(arg0: QCircuit, arg1: QuditDepolarizingNoise) -> _QDensityNoisyEngine_QuditDepolarizingNoise:
    ...
def QFT(A: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], d: typing.SupportsInt | typing.SupportsIndex = 2, swap: bool = True) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]:
    """
    Qudit quantum Fourier transform
    """
def QKetDummyEngine(arg0: QCircuit) -> _QKetDummyEngine:
    ...
def QKetEngine(arg0: QCircuit) -> _QKetEngine:
    ...
@typing.overload
def QKetNoisyEngine(arg0: QCircuit, arg1: QubitBitFlipNoise) -> _QKetNoisyEngine_QubitBitFlipNoise:
    ...
@typing.overload
def QKetNoisyEngine(arg0: QCircuit, arg1: QubitBitPhaseFlipNoise) -> _QKetNoisyEngine_QubitBitPhaseFlipNoise:
    ...
@typing.overload
def QKetNoisyEngine(arg0: QCircuit, arg1: QubitDepolarizingNoise) -> _QKetNoisyEngine_QubitDepolarizingNoise:
    ...
@typing.overload
def QKetNoisyEngine(arg0: QCircuit, arg1: QubitPhaseFlipNoise) -> _QKetNoisyEngine_QubitPhaseFlipNoise:
    ...
@typing.overload
def QKetNoisyEngine(arg0: QCircuit, arg1: QubitAmplitudeDampingNoise) -> _QKetNoisyEngine_QubitAmplitudeDampingNoise:
    ...
@typing.overload
def QKetNoisyEngine(arg0: QCircuit, arg1: QubitPhaseDampingNoise) -> _QKetNoisyEngine_QubitPhaseDampingNoise:
    ...
@typing.overload
def QKetNoisyEngine(arg0: QCircuit, arg1: QuditDepolarizingNoise) -> _QKetNoisyEngine_QuditDepolarizingNoise:
    ...
def TFQ(A: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], d: typing.SupportsInt | typing.SupportsIndex = 2, swap: bool = True) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]:
    """
    Inverse (adjoint) qudit quantum Fourier transform
    """
def absm(arg0: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"]) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]:
    """
    Matrix absolute value
    """
@typing.overload
def abssq(arg0: collections.abc.Sequence[typing.SupportsComplex | typing.SupportsFloat | typing.SupportsIndex]) -> list[float]:
    """
    Absolute values squared of vector
    """
@typing.overload
def abssq(arg0: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"]) -> list[float]:
    """
    Absolute values squared of matrix elements
    """
@typing.overload
def adjoint(arg0: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"]) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]:
    """
    Adjoint (Hermitian conjugate)
    """
@typing.overload
def adjoint(qc: QCircuit, name: str | None = None) -> QCircuit:
    """
    Adjoint quantum circuit description
    """
def anticomm(arg0: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], arg1: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"]) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]:
    """
    Anticommutator
    """
@typing.overload
def apply(state: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], A: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], target: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], dims: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]:
    """
    Applies the gate A to the part target of the multi-partite state vector or density matrix state
    """
@typing.overload
def apply(state: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], A: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], target: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], d: typing.SupportsInt | typing.SupportsIndex = 2) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]:
    """
    Applies the gate A to the part target of the multi-partite state vector or density matrix state
    """
@typing.overload
def apply(A: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], Ks: collections.abc.Sequence[typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"]]) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]:
    """
    Applies the channel specified by the set of Kraus operators Ks to the density matrix A
    """
@typing.overload
def apply(A: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], Ks: collections.abc.Sequence[typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"]], target: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], dims: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]:
    """
    Applies the channel specified by the set of Kraus operators Ks to the part target of the multi-partite density matrix A
    """
@typing.overload
def apply(A: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], Ks: collections.abc.Sequence[typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"]], target: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], d: typing.SupportsInt | typing.SupportsIndex = 2) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]:
    """
    Applies the channel specified by the set of Kraus operators Ks to the part target of the multi-partite density matrix A
    """
@typing.overload
def applyCTRL(state: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], A: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], ctrl: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], target: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], dims: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], shift: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex] | None = None) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]:
    """
    Applies the controlled-gate A to the part target of the multi-partite state vector or density matrix state
    """
@typing.overload
def applyCTRL(state: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], A: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], ctrl: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], target: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], d: typing.SupportsInt | typing.SupportsIndex = 2, shift: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex] | None = None) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]:
    """
    Applies the controlled-gate A to the part target of the multi-partite state vector or density matrix state
    """
@typing.overload
def applyCTRL_fan(state: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], A: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], ctrl: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], target: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], dims: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], shift: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex] | None = None) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]:
    """
    Applies the single qudit controlled-gate A with multiple control qudits listed in ctrl on every qudit listed in target
    """
@typing.overload
def applyCTRL_fan(state: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], A: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], ctrl: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], target: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], d: typing.SupportsInt | typing.SupportsIndex = 2, shift: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex] | None = None) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]:
    """
    Applies the single qudit controlled-gate A with multiple control qudits listed in ctrl on every qudit listed in target
    """
def applyQFT(A: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], target: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], d: typing.SupportsInt | typing.SupportsIndex = 2, swap: bool = True) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]:
    """
    Applies the qudit quantum Fourier transform to the part target of the multi-partite state vector or density matrix A
    """
def applyTFQ(A: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], target: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], d: typing.SupportsInt | typing.SupportsIndex = 2, swap: bool = True) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]:
    """
    Applies the inverse (adjoint) qudit quantum Fourier transform to the part target of the multi-partite state vector or density matrix A
    """
def avg(prob: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex], X: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]) -> float:
    """
    Computes the average E[X] = sum(p_i * x_i).
    """
def bernoulli(p: typing.SupportsFloat | typing.SupportsIndex = 0.5) -> bool:
    """
    Generates a random boolean from a Bernoulli-p distribution.
    """
def bloch2rho(r: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]:
    """
    Computes the density matrix corresponding to the 3-dimensional real Bloch vector r
    """
@typing.overload
def choi2kraus(A: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], Din: typing.SupportsInt | typing.SupportsIndex, Dout: typing.SupportsInt | typing.SupportsIndex) -> list[typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]]:
    """
    Orthogonal Kraus operators from Choi matrix. Extracts a set of orthogonal Kraus operators from the Choi matrix A
    """
@typing.overload
def choi2kraus(A: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"]) -> list[typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]]:
    """
    Orthogonal Kraus operators from Choi matrix. Extracts a set of orthogonal Kraus operators from the Choi matrix A, assuming square Kraus operators
    """
@typing.overload
def choi2super(A: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], Din: typing.SupportsInt | typing.SupportsIndex, Dout: typing.SupportsInt | typing.SupportsIndex) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]:
    """
    Converts Choi matrix to superoperator matrix
    """
@typing.overload
def choi2super(A: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"]) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]:
    """
    Converts Choi matrix to superoperator matrix, assuming square Kraus operators
    """
def comm(arg0: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], arg1: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"]) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]:
    """
    Commutator
    """
def complement(arg0: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], arg1: typing.SupportsInt | typing.SupportsIndex) -> list[int]:
    """
    Complement of a subsystem
    """
def compose_CTRL_circuit(qc_ctrl: QCircuit, ctrl: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], qc_target: QCircuit, pos_qudit: typing.SupportsInt | typing.SupportsIndex, shift: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex] | None = None, pos_dit: typing.SupportsInt | typing.SupportsIndex | None = None, name: str | None = None) -> QCircuit:
    """
    Composes (appends) the qc_target controlled quantum circuit description to the end of the qc_ctrl quantum circuit description
    """
def compose_circuit(qc1: QCircuit, qc2: QCircuit, pos_qudit: typing.SupportsInt | typing.SupportsIndex, name: str | None = None, pos_dit: typing.SupportsInt | typing.SupportsIndex | None = None) -> QCircuit:
    """
    Composes (appends) the second quantum circuit description to the end of the first one; qc_ctrl controls the qc_target.
    """
def compperm(arg0: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], arg1: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> list[int]:
    """
    Composition of two permutations
    """
def concurrence(A: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"]) -> float:
    """
    Wootters concurrence of the bi-partite qubit mixed state A.
    """
def conjugate(arg0: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"]) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]:
    """
    Complex conjugate
    """
def contfrac2x(cf: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], N: typing.SupportsInt | typing.SupportsIndex = 18446744073709551615) -> float:
    """
    Real representation of a simple continued fraction
    """
@typing.overload
def convergents(arg0: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> list[tuple[int, int]]:
    """
    Vector of convergent pairs (a_k, b_k) from a continued fraction
    """
@typing.overload
def convergents(arg0: typing.SupportsFloat | typing.SupportsIndex, arg1: typing.SupportsInt | typing.SupportsIndex) -> list[tuple[int, int]]:
    """
    Vector of convergent pairs (a_k, b_k) approximating real x
    """
def cor(probXY: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, n]"], X: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex], Y: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]) -> float:
    """
    Computes the correlation coefficient of random variables X and Y.
    """
def cosm(arg0: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"]) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]:
    """
    Matrix cosine
    """
def couple_circuit_left(qc1: QCircuit, qc2: QCircuit, target: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], name: str | None = None, pos_dit: typing.SupportsInt | typing.SupportsIndex | None = None) -> QCircuit:
    """
    Couples (in place) the second quantum circuit description to the left (beginning) of the first one
    """
def couple_circuit_right(qc1: QCircuit, qc2: QCircuit, target: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], name: str | None = None, pos_dit: typing.SupportsInt | typing.SupportsIndex | None = None) -> QCircuit:
    """
    Couples (in place) the second quantum circuit description to the right (end) of the first one
    """
def cov(probXY: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, n]"], X: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex], Y: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]) -> float:
    """
    Computes the covariance of random variables X and Y.
    """
def det(arg0: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"]) -> complex:
    """
    Determinant
    """
@typing.overload
def dirac(A: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], dims_rows: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], dims_cols: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> qpp::dirac_t<std::__1::complex<double>>:
    """
    Dirac notation
    """
@typing.overload
def dirac(A: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], d: typing.SupportsInt | typing.SupportsIndex = 2) -> qpp::dirac_t<std::__1::complex<double>>:
    """
    Dirac notation
    """
def dirsum(arg0: collections.abc.Sequence[typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"]]) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]:
    """
    Direct sum of multiple matrices
    """
def dirsumpow(arg0: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], arg1: typing.SupportsInt | typing.SupportsIndex) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]:
    """
    Direct sum power
    """
@typing.overload
def discard(A: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], target: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], dims: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]:
    """
    Discards qudits from the multi-partite state A by performing a destructive measurement and discarding the results. Returns the reduced density matrix.
    """
@typing.overload
def discard(A: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], target: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], dims: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]:
    """
    Discards qudits from the multi-partite state A by performing a destructive measurement and discarding the results. Returns the reduced density matrix.
    """
def egcd(arg0: typing.SupportsInt | typing.SupportsIndex, arg1: typing.SupportsInt | typing.SupportsIndex) -> tuple[int, int, int]:
    """
    Extended GCD. Returns (m, n, gcd) such that ma + nb = gcd(a, b).
    """
def eig(arg0: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"]) -> tuple[typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, 1]"], typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]]:
    """
    Full eigen decomposition
    """
def entanglement(A: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], dims: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> float:
    """
    von-Neumann entanglement entropy of the bi-partite pure state A.
    """
@typing.overload
def entropy(A: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"]) -> float:
    """
    Computes the von-Neumann entropy of a density matrix (base 2).
    """
@typing.overload
def entropy(prob: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]) -> float:
    """
    Computes the Shannon entropy of a probability distribution (base 2).
    """
def evals(arg0: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"]) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, 1]"]:
    """
    Eigenvalues
    """
def evects(arg0: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"]) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]:
    """
    Eigenvectors
    """
def expm(arg0: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"]) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]:
    """
    Matrix exponential
    """
def factors(arg0: typing.SupportsInt | typing.SupportsIndex) -> list[int]:
    """
    Prime factor decomposition of an integer
    """
def funm(A: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], f: std::__1::complex<double> (std::__1::complex<double> const&)) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]:
    """
    Applies the scalar function f to the eigenvalues of matrix A
    """
@typing.overload
def gcd(arg0: typing.SupportsInt | typing.SupportsIndex, arg1: typing.SupportsInt | typing.SupportsIndex) -> int:
    """
    Greatest common divisor of two integers
    """
@typing.overload
def gcd(arg0: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> int:
    """
    GCD of a list of integers
    """
def gconcurrence(A: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"]) -> float:
    """
    G-concurrence of the bi-partite pure state A.
    """
@typing.overload
def grams(arg0: collections.abc.Sequence[typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"]]) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]:
    """
    Gram-Schmidt orthogonalization (from list)
    """
@typing.overload
def grams(arg0: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"]) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]:
    """
    Gram-Schmidt orthogonalization (from columns)
    """
def hash_eigen(A: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], seed: typing.SupportsInt | typing.SupportsIndex = 0) -> int:
    """
    Hash an Eigen expression
    """
def heig(arg0: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"]) -> tuple[typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"], typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]]:
    """
    Full eigen decomposition of Hermitian matrix
    """
def hevals(arg0: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"]) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
    """
    Hermitian eigenvalues
    """
def hevects(arg0: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"]) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]:
    """
    Hermitian eigenvectors
    """
def inverse(arg0: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"]) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]:
    """
    Inverse
    """
def invperm(arg0: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> list[int]:
    """
    Inverse of a permutation
    """
@typing.overload
def ip(phi: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, 1]"], psi: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, 1]"], subsys: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], dims: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, 1]"]:
    """
    Generalized inner product. Computes the projection of the multi-partite state psi onto the subsystem phi.
    """
@typing.overload
def ip(phi: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, 1]"], psi: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, 1]"], subsys: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], d: typing.SupportsInt | typing.SupportsIndex = 2) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, 1]"]:
    """
    Generalized inner product. Computes the projection of the multi-partite state psi onto the subsystem phi, assuming uniform subsystem dimensions.
    """
def isprime(p: typing.SupportsInt | typing.SupportsIndex, k: typing.SupportsInt | typing.SupportsIndex = 80) -> bool:
    """
    Miller-Rabin primality test
    """
def kraus2choi(Ks: collections.abc.Sequence[typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"]]) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]:
    """
    Choi matrix. Constructs the Choi matrix of the channel specified by the set of Kraus operators Ks
    """
def kraus2super(Ks: collections.abc.Sequence[typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"]]) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]:
    """
    Superoperator matrix. Constructs the superoperator matrix of the channel specified by the set of Kraus operators Ks
    """
@typing.overload
def kron(arg0: collections.abc.Sequence[typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"]]) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]:
    """
    Kronecker product of multiple matrices
    """
@typing.overload
def kron(qc1: QCircuit, qc2: QCircuit, name: str | None = None) -> QCircuit:
    """
    Kronecker product between two quantum circuit descriptions
    """
def kronpow(arg0: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], arg1: typing.SupportsInt | typing.SupportsIndex) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]:
    """
    Kronecker power
    """
@typing.overload
def lcm(arg0: typing.SupportsInt | typing.SupportsIndex, arg1: typing.SupportsInt | typing.SupportsIndex) -> int:
    """
    Least common multiple of two integers
    """
@typing.overload
def lcm(arg0: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> int:
    """
    LCM of a list of integers
    """
def load_cmat(filename: str) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]:
    """
    Loads a complex matrix from a text file
    """
def load_rmat(filename: str) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, n]"]:
    """
    Loads a real matrix from a text file
    """
def logdet(arg0: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"]) -> complex:
    """
    Logarithm of the determinant
    """
def logm(arg0: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"]) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]:
    """
    Matrix logarithm
    """
def lognegativity(A: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], dims: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> float:
    """
    Logarithmic negativity of the bi-partite mixed state A.
    """
def marginalX(probXY: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, n]"]) -> list[float]:
    """
    Computes the marginal distribution of X from a joint probability matrix.
    """
def marginalY(probXY: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, n]"]) -> list[float]:
    """
    Computes the marginal distribution of Y from a joint probability matrix.
    """
@typing.overload
def measure(A: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], Ks: collections.abc.Sequence[typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"]]) -> tuple[int, list[float], list[typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]]]:
    """
    Measures the state vector or density operator A using the set of Kraus operators Ks. Returns a tuple of: 1. Result of the measurement (index), 2. Vector of outcome probabilities, and 3. Vector of post-measurement normalized states.
    """
@typing.overload
def measure(A: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], U: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"]) -> tuple[int, list[float], list[typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]]]:
    """
    Measures the state vector or density matrix A in the orthonormal basis specified by the unitary matrix U. Returns a tuple of: 1. Result of the measurement, 2. Vector of outcome probabilities, and 3. Vector of post-measurement normalized states.
    """
@typing.overload
def measure(A: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], Ks: collections.abc.Sequence[typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"]], target: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], dims: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], destructive: bool = True) -> tuple[int, list[float], list[typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]]]:
    """
    Measures the part 'target' of the multi-partite state vector or density matrix A using the set of Kraus operators Ks. If 'destructive' is True, the measured subsystems are traced away.
    """
@typing.overload
def measure(A: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], Ks: collections.abc.Sequence[typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"]], target: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], d: typing.SupportsInt | typing.SupportsIndex = 2, destructive: bool = True) -> tuple[int, list[float], list[typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]]]:
    """
    Measures the part 'target' of the multi-partite state vector or density matrix A using the set of Kraus operators Ks, assuming uniform subsystem dimensions. If 'destructive' is True, the measured subsystems are traced away.
    """
@typing.overload
def measure(A: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], V: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], target: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], dims: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], destructive: bool = True) -> tuple[int, list[float], list[typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]]]:
    """
    Measures the part 'target' of the multi-partite state vector or density matrix A in the orthonormal basis specified by the columns of matrix V. If 'destructive' is True, the measured subsystems are traced away.
    """
@typing.overload
def measure(A: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], V: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], target: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], d: typing.SupportsInt | typing.SupportsIndex = 2, destructive: bool = True) -> tuple[int, list[float], list[typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]]]:
    """
    Measures the part 'target' of the multi-partite state vector or density matrix A in the orthonormal basis specified by the columns of matrix V, assuming uniform subsystem dimensions. If 'destructive' is True, the measured subsystems are traced away.
    """
@typing.overload
def measure_seq(A: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], target: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], dims: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], destructive: bool = True) -> tuple[list[int], list[float], typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]]:
    """
    Sequentially measures the subsystems specified by 'target' in the computational basis. Returns a tuple of: 1. Vector of measurement outcomes, 2. Vector of probabilities for each outcome, and 3. The resulting post-measurement normalized state.
    """
@typing.overload
def measure_seq(A: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], target: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], d: typing.SupportsInt | typing.SupportsIndex = 2, destructive: bool = True) -> tuple[list[int], list[float], typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]]:
    """
    Sequentially measures the subsystems specified by 'target' in the computational basis, assuming uniform subsystem dimensions. Returns a tuple of: 1. Vector of measurement outcomes, 2. Outcome probabilities, and 3. The final post-measurement normalized state.
    """
@typing.overload
def mket(states: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], dims: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, 1]"]:
    """
    Creates a multi-partite ket from indices and dimensions
    """
@typing.overload
def mket(states: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], d: typing.SupportsInt | typing.SupportsIndex = 2) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, 1]"]:
    """
    Creates a multi-partite ket from indices with uniform dimension d
    """
def modinv(arg0: typing.SupportsInt | typing.SupportsIndex, arg1: typing.SupportsInt | typing.SupportsIndex) -> int:
    """
    Modular inverse of a mod p (a and p must be co-prime)
    """
def modmul(arg0: typing.SupportsInt | typing.SupportsIndex, arg1: typing.SupportsInt | typing.SupportsIndex, arg2: typing.SupportsInt | typing.SupportsIndex) -> int:
    """
    Modular multiplication (a*b % p) without overflow
    """
def modpow(arg0: typing.SupportsInt | typing.SupportsIndex, arg1: typing.SupportsInt | typing.SupportsIndex, arg2: typing.SupportsInt | typing.SupportsIndex) -> int:
    """
    Fast modular exponentiation (a^n % p)
    """
@typing.overload
def mprj(arg0: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], arg1: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]:
    """
    Projector onto multi-partite qudit state
    """
@typing.overload
def mprj(mask: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], d: typing.SupportsInt | typing.SupportsIndex = 2) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]:
    """
    Projector onto multi-partite qudit state (uniform d)
    """
def multiidx2n(arg0: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], arg1: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> int:
    """
    Multi-index to integer conversion
    """
def n2multiidx(arg0: typing.SupportsInt | typing.SupportsIndex, arg1: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> list[int]:
    """
    Integer to multi-index conversion
    """
def negativity(A: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], dims: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> float:
    """
    Negativity of the bi-partite mixed state A.
    """
def norm(arg0: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"]) -> float:
    """
    Frobenius norm
    """
def normalize(arg0: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"]) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]:
    """
    Normalize state vector or density matrix
    """
def omega(D: typing.SupportsInt | typing.SupportsIndex) -> complex:
    """
    D-th root of unity
    """
def powm(arg0: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], arg1: typing.SupportsInt | typing.SupportsIndex) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]:
    """
    Fast matrix power (square-and-multiply)
    """
def prj(arg0: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"]) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]:
    """
    Projector onto state vector
    """
@typing.overload
def prod(arg0: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"]) -> complex:
    """
    Element-wise product of matrix
    """
@typing.overload
def prod(arg0: collections.abc.Sequence[typing.SupportsComplex | typing.SupportsFloat | typing.SupportsIndex]) -> complex:
    """
    Product of vector elements
    """
@typing.overload
def ptrace(A: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], target: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], dims: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]:
    """
    Partial trace of the multi-partite state vector or density matrix over the list target of subsystems
    """
@typing.overload
def ptrace(A: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], target: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], d: typing.SupportsInt | typing.SupportsIndex = 2) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]:
    """
    Partial trace of the multi-partite state vector or density matrix over the list target of subsystems, assuming uniform subsystem dimensions
    """
@typing.overload
def ptrace1(A: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], dims: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]:
    """
    Partial trace over the first subsystem of a bi-partite state vector or density matrix
    """
@typing.overload
def ptrace1(A: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], d: typing.SupportsInt | typing.SupportsIndex = 2) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]:
    """
    Partial trace over the first subsystem of a bi-partite state vector or density matrix, assuming uniform subsystem dimensions
    """
@typing.overload
def ptrace2(A: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], dims: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]:
    """
    Partial trace over the second subsystem of a bi-partite state vector or density matrix
    """
@typing.overload
def ptrace2(A: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], d: typing.SupportsInt | typing.SupportsIndex = 2) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]:
    """
    Partial trace over the second subsystem of a bi-partite state vector or density matrix, assuming uniform subsystem dimensions
    """
@typing.overload
def ptranspose(A: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], target: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], dims: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]:
    """
    Partial transpose of the multi-partite state vector or density matrix over the list target of subsystems
    """
@typing.overload
def ptranspose(A: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], target: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], d: typing.SupportsInt | typing.SupportsIndex = 2) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]:
    """
    Partial transpose of the multi-partite state vector or density matrix over the list target of subsystems, assuming uniform subsystem dimensions
    """
@typing.overload
def qRAM(psi: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, 1]"], data: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], DqRAM: typing.SupportsInt | typing.SupportsIndex) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, 1]"]:
    """
    Quantumly-accessible Random Access Memory (qRAM) over classical data. Implements the mapping sum_j alpha_j |j> -> sum_j alpha_j |j>|m_j>.
    """
@typing.overload
def qRAM(psi: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, 1]"], data: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, 1]"]:
    """
    Quantumly-accessible Random Access Memory (qRAM) over classical data. The qRAM subsystem dimension is automatically set to 1 + maximum value stored in the data.
    """
@typing.overload
def qmutualinfo(A: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], subsysA: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], subsysB: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], dims: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> float:
    """
    Computes the quantum mutual information between two subsystems A and B with given dimensions.
    """
@typing.overload
def qmutualinfo(A: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], subsysA: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], subsysB: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], d: typing.SupportsInt | typing.SupportsIndex = 2) -> float:
    """
    Computes the quantum mutual information between two subsystems A and B assuming uniform dimensions d.
    """
def qpe_circuit(U: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], n: typing.SupportsInt | typing.SupportsIndex, omit_measurements: bool = True, d: typing.SupportsInt | typing.SupportsIndex = 2, name: str | None = 'qpe') -> QCircuit:
    """
    Quantum phase estimation circuit with n bits of precision
    """
@typing.overload
def rand(a: typing.SupportsFloat | typing.SupportsIndex, b: typing.SupportsFloat | typing.SupportsIndex) -> float:
    """
    Generates a random real number in [a, b).
    """
@typing.overload
def rand(rows: typing.SupportsInt | typing.SupportsIndex, cols: typing.SupportsInt | typing.SupportsIndex, a: typing.SupportsFloat | typing.SupportsIndex = 0, b: typing.SupportsFloat | typing.SupportsIndex = 1) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]:
    """
    Generates a random complex matrix with entries uniformly distributed in [a, b).
    """
def randH(D: typing.SupportsInt | typing.SupportsIndex = 2) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]:
    """
    Generates a random Hermitian matrix.
    """
def randU(D: typing.SupportsInt | typing.SupportsIndex = 2) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]:
    """
    Generates a random unitary matrix (Haar measure).
    """
def randV(Din: typing.SupportsInt | typing.SupportsIndex, Dout: typing.SupportsInt | typing.SupportsIndex) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]:
    """
    Generates a random isometry matrix.
    """
def randidx(a: typing.SupportsInt | typing.SupportsIndex = 0, b: typing.SupportsInt | typing.SupportsIndex = 18446744073709551615) -> int:
    """
    Generates a random index uniformly distributed in [a, b].
    """
def randket(D: typing.SupportsInt | typing.SupportsIndex = 2) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, 1]"]:
    """
    Generates a random normalized pure state (ket).
    """
@typing.overload
def randkraus(N: typing.SupportsInt | typing.SupportsIndex, Din: typing.SupportsInt | typing.SupportsIndex, Dout: typing.SupportsInt | typing.SupportsIndex) -> list[typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]]:
    """
    Generates a set of random Kraus operators (rectangular).
    """
@typing.overload
def randkraus(N: typing.SupportsInt | typing.SupportsIndex, D: typing.SupportsInt | typing.SupportsIndex = 2) -> list[typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]]:
    """
    Generates a set of random Kraus operators (square).
    """
@typing.overload
def randn(mean: typing.SupportsFloat | typing.SupportsIndex = 0, sigma: typing.SupportsFloat | typing.SupportsIndex = 1) -> float:
    """
    Generates a random real number normally distributed in N(mean, sigma).
    """
@typing.overload
def randn(rows: typing.SupportsInt | typing.SupportsIndex, cols: typing.SupportsInt | typing.SupportsIndex, mean: typing.SupportsFloat | typing.SupportsIndex = 0, sigma: typing.SupportsFloat | typing.SupportsIndex = 1) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]:
    """
    Generates a random complex matrix with entries normally distributed in N(mean, sigma).
    """
def random_circuit_count(nq: typing.SupportsInt | typing.SupportsIndex, d: typing.SupportsInt | typing.SupportsIndex, gate_count: typing.SupportsInt | typing.SupportsIndex, p_two: typing.SupportsFloat | typing.SupportsIndex | None = None, with_respect_to_gate: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"] | None = None, one_qudit_gate_set: collections.abc.Sequence[typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"]] | None = None, two_qudit_gate_set: collections.abc.Sequence[typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"]] | None = None, one_qudit_gate_names: collections.abc.Sequence[str] | None = None, two_qudit_gate_names: collections.abc.Sequence[str] | None = None) -> QCircuit:
    """
    Random quantum circuit description generator for fixed gate count
    """
def random_circuit_depth(nq: typing.SupportsInt | typing.SupportsIndex, d: typing.SupportsInt | typing.SupportsIndex, gate_depth: typing.SupportsInt | typing.SupportsIndex, p_two: typing.SupportsFloat | typing.SupportsIndex | None = None, with_respect_to_gate: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"] | None = None, one_qudit_gate_set: collections.abc.Sequence[typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"]] | None = None, two_qudit_gate_set: collections.abc.Sequence[typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"]] | None = None, one_qudit_gate_names: collections.abc.Sequence[str] | None = None, two_qudit_gate_names: collections.abc.Sequence[str] | None = None) -> QCircuit:
    """
    Random quantum circuit description generator for fixed gate depth
    """
def randperm(N: typing.SupportsInt | typing.SupportsIndex) -> list[int]:
    """
    Generates a random permutation of [0, ..., N-1].
    """
def randprime(a: typing.SupportsInt | typing.SupportsIndex, b: typing.SupportsInt | typing.SupportsIndex, N: typing.SupportsInt | typing.SupportsIndex = 1000) -> int:
    """
    Generates a random prime in the interval [a, b]
    """
def randprob(N: typing.SupportsInt | typing.SupportsIndex) -> list[float]:
    """
    Generates a random probability vector uniformly over the simplex.
    """
def randrho(D: typing.SupportsInt | typing.SupportsIndex = 2) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]:
    """
    Generates a random density matrix.
    """
@typing.overload
def renyi(A: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], alpha: typing.SupportsFloat | typing.SupportsIndex) -> float:
    """
    Computes the Renyi-alpha entropy of a density matrix.
    """
@typing.overload
def renyi(prob: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex], alpha: typing.SupportsFloat | typing.SupportsIndex) -> float:
    """
    Computes the Renyi-alpha entropy of a probability distribution.
    """
def replicate(qc: QCircuit, n: typing.SupportsInt | typing.SupportsIndex, name: str | None = None) -> QCircuit:
    """
    Replicates a quantum circuit description
    """
@typing.overload
def reset(A: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], target: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], dims: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]:
    """
    Resets qudits in the multi-partite state A by performing a non-destructive computational basis measurement on the 'target' qudits and shifting them back to the |0...0> state.
    """
@typing.overload
def reset(A: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], target: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], d: typing.SupportsInt | typing.SupportsIndex = 2) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]:
    """
    Resets qudits in the multi-partite state A, assuming uniform subsystem dimensions. Performs a non-destructive computational basis measurement and shifts target qudits back to the |0...0> state.
    """
def reshape(arg0: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], arg1: typing.SupportsInt | typing.SupportsIndex, arg2: typing.SupportsInt | typing.SupportsIndex) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]:
    """
    Reshape matrix
    """
def rho2bloch(A: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"]) -> list[float]:
    """
    Computes the 3-dimensional real Bloch vector corresponding to the qubit density matrix A
    """
def rho2pure(arg0: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"]) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, 1]"]:
    """
    Finds pure state representation of rank-1 matrix
    """
@typing.overload
def sample(A: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], target: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], dims: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> tuple[list[int], float]:
    """
    Samples from a quantum state in the computational basis (Z-basis)
    """
@typing.overload
def sample(A: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], target: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], d: typing.SupportsInt | typing.SupportsIndex = 2) -> tuple[list[int], float]:
    """
    Samples from a quantum state in the computational basis (Z-basis)
    """
@typing.overload
def sample(num_samples: typing.SupportsInt | typing.SupportsIndex, A: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], target: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], dims: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> dict[str, int]:
    """
    Samples repeatedly from a quantum state in the computational basis (Z-basis)
    """
@typing.overload
def sample(num_samples: typing.SupportsInt | typing.SupportsIndex, A: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], target: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], d: typing.SupportsInt | typing.SupportsIndex = 2) -> dict[str, int]:
    """
    Samples repeatedly from a quantum state in the computational basis (Z-basis)
    """
@typing.overload
def save(A: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], filename: str) -> None:
    """
    Saves a complex matrix to a text file
    """
@typing.overload
def save(A: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, n]"], filename: str) -> None:
    """
    Saves a real matrix to a text file
    """
def schatten(arg0: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], arg1: typing.SupportsFloat | typing.SupportsIndex) -> float:
    """
    Schatten matrix norm
    """
def schmidt(A: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], dims: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> tuple[typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"], typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"], typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"], typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]]:
    """
    Full Schmidt decomposition: returns tuple(U, V, coeffs, probs).
    """
def schmidtA(A: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], dims: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]:
    """
    Schmidt basis on Alice side.
    """
def schmidtB(A: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], dims: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]:
    """
    Schmidt basis on Bob side.
    """
@typing.overload
def schmidtcoeffs(A: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], dims: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
    """
    Schmidt coefficients of the bi-partite pure state A.
    """
@typing.overload
def schmidtcoeffs(A: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], d: typing.SupportsInt | typing.SupportsIndex = 2) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
    """
    Schmidt coefficients of the bi-partite pure state A (uniform dimensions).
    """
def schmidtprobs(A: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], dims: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> list[float]:
    """
    Schmidt probabilities (squares of coefficients) of the bi-partite pure state A.
    """
@typing.overload
def set_prng_seed(seed: typing.SupportsInt | typing.SupportsIndex) -> None:
    """
    Sets the prng seed to a specific value
    """
@typing.overload
def set_prng_seed() -> None:
    """
    Sets the prng seed to a random value
    """
def sigma(prob: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex], X: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]) -> float:
    """
    Computes the standard deviation of random variable X.
    """
def sinm(arg0: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"]) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]:
    """
    Matrix sine
    """
def spectralpowm(arg0: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], arg1: typing.SupportsComplex | typing.SupportsFloat | typing.SupportsIndex) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]:
    """
    Matrix power via spectral decomposition
    """
def sqrtm(arg0: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"]) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]:
    """
    Matrix square root
    """
@typing.overload
def sum(arg0: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"]) -> complex:
    """
    Element-wise sum of matrix
    """
@typing.overload
def sum(arg0: collections.abc.Sequence[typing.SupportsComplex | typing.SupportsFloat | typing.SupportsIndex]) -> complex:
    """
    Sum of vector elements
    """
def super2choi(A: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"]) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]:
    """
    Converts superoperator matrix to Choi matrix
    """
def super2kraus(A: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"]) -> list[typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]]:
    """
    Orthogonal Kraus operators from superoperator matrix. Extracts a set of orthogonal Kraus operators from the superoperator matrix A
    """
def svals(arg0: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"]) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
    """
    Singular values
    """
def svd(arg0: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"]) -> tuple[typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"], typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"], typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]]:
    """
    Full singular value decomposition
    """
def svdU(arg0: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"]) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]:
    """
    Left singular vectors
    """
def svdV(arg0: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"]) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]:
    """
    Right singular vectors
    """
@typing.overload
def syspermute(A: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], perm: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], dims: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]:
    """
    Subsystem permutation. Permutes the subsystems of a state vector or density matrix. The subsystem perm[i] is moved to the location i.
    """
@typing.overload
def syspermute(A: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], perm: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], d: typing.SupportsInt | typing.SupportsIndex = 2) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]:
    """
    Subsystem permutation. Permutes the subsystems of a state vector or density matrix, assuming uniform subsystem dimensions.
    """
def trace(arg0: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"]) -> complex:
    """
    Trace
    """
def transpose(arg0: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"]) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]:
    """
    Transpose
    """
@typing.overload
def tsallis(A: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], q: typing.SupportsFloat | typing.SupportsIndex) -> float:
    """
    Computes the Tsallis-q entropy of a density matrix.
    """
@typing.overload
def tsallis(prob: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex], q: typing.SupportsFloat | typing.SupportsIndex) -> float:
    """
    Computes the Tsallis-q entropy of a probability distribution.
    """
def uniform(N: typing.SupportsInt | typing.SupportsIndex) -> list[float]:
    """
    Generates a uniform probability distribution vector of size N.
    """
def var(prob: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex], X: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]) -> float:
    """
    Computes the variance of random variable X.
    """
def x2contfrac(x: typing.SupportsFloat | typing.SupportsIndex, N: typing.SupportsInt | typing.SupportsIndex, cut: typing.SupportsInt | typing.SupportsIndex = 10000) -> list[int]:
    """
    Simple continued fraction expansion of x
    """
@typing.overload
def zket2dits(psi: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], dims: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], precision: typing.SupportsFloat | typing.SupportsIndex = 1e-12) -> list[int] | None:
    """
    Extracts dits from a computational basis state given specific subsystem dimensions.
    Returns None if psi is not a computational basis state.
    """
@typing.overload
def zket2dits(psi: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"], d: typing.SupportsInt | typing.SupportsIndex = 2, precision: typing.SupportsFloat | typing.SupportsIndex = 1e-12) -> list[int] | None:
    """
    Extracts dits from a computational basis state assuming uniform subsystem dimension d.
    Returns None if psi is not a computational basis state.
    """
ee: float = 2.718281828459045
infty: float  # value = inf
pi: float = 3.141592653589793
QEngine = QKetEngine
QNoisyEngine = QKetNoisyEngine
