from __future__ import annotations
import collections.abc
import numpy
import numpy.typing
import typing
__all__: list[str] = ['CNOT', 'CNOTba', 'CZ', 'FRED', 'Fd', 'H', 'Id', 'Id2', 'MODMUL', 'RX', 'RXX', 'RY', 'RYY', 'RZ', 'Rn', 'S', 'SWAP', 'SWAPd', 'T', 'TOF', 'X', 'Xd', 'Y', 'Z', 'Zd', 'get_name']
def Fd(D: typing.SupportsInt | typing.SupportsIndex = 2) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]:
    """
    Quantum Fourier transform gate for qudits
    """
def Id(D: typing.SupportsInt | typing.SupportsIndex = 2) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]:
    """
    Identity gate
    """
def MODMUL(a: typing.SupportsInt | typing.SupportsIndex, N: typing.SupportsInt | typing.SupportsIndex, n: typing.SupportsInt | typing.SupportsIndex) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]:
    """
    Modular multiplication gate for qubits
    """
def RX(theta: typing.SupportsFloat | typing.SupportsIndex) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]:
    """
    Qubit rotation of theta about the X axis
    """
def RY(theta: typing.SupportsFloat | typing.SupportsIndex) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]:
    """
    Qubit rotation of theta about the Y axis
    """
def RZ(theta: typing.SupportsFloat | typing.SupportsIndex) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]:
    """
    Qubit rotation of theta about the Z axis
    """
def Rn(theta: typing.SupportsFloat | typing.SupportsIndex, n: typing.Annotated[collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex], "FixedSize(3)"]) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]:
    """
    Qubit rotation of theta about the 3-dimensional real (unit) vector n
    """
def SWAPd(D: typing.SupportsInt | typing.SupportsIndex = 2) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]:
    """
    SWAP gate for qudits
    """
def Xd(D: typing.SupportsInt | typing.SupportsIndex = 2) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]:
    """
    Generalized X gate for qudits
    """
def Zd(D: typing.SupportsInt | typing.SupportsIndex = 2) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]:
    """
    Generalized Z gate for qudits
    """
def get_name(U: typing.Annotated[numpy.typing.ArrayLike, numpy.complex128, "[m, n]"]) -> str | None:
    """
    Get the name of the most common qubit gates
    """
CNOT: numpy.ndarray  # value = array([[1.+0.j, 0.+0.j, 0.+0.j, 0.+0.j],...
CNOTba: numpy.ndarray  # value = array([[1.+0.j, 0.+0.j, 0.+0.j, 0.+0.j],...
CZ: numpy.ndarray  # value = array([[ 1.+0.j,  0.+0.j,  0.+0.j,  0.+0.j],...
FRED: numpy.ndarray  # value = array([[1.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j],...
H: numpy.ndarray  # value = array([[ 0.70710678+0.j,  0.70710678+0.j],...
Id2: numpy.ndarray  # value = array([[1.+0.j, 0.+0.j],...
RXX: numpy.ndarray  # value = array([[0.70710678+0.j        , 0.        +0.j        ,...
RYY: numpy.ndarray  # value = array([[ 0.70710678+0.j        ,  0.        +0.j        ,...
S: numpy.ndarray  # value = array([[1.+0.j, 0.+0.j],...
SWAP: numpy.ndarray  # value = array([[1.+0.j, 0.+0.j, 0.+0.j, 0.+0.j],...
T: numpy.ndarray  # value = array([[1.        +0.j        , 0.        +0.j        ],...
TOF: numpy.ndarray  # value = array([[1.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j],...
X: numpy.ndarray  # value = array([[0.+0.j, 1.+0.j],...
Y: numpy.ndarray  # value = array([[ 0.+0.j, -0.-1.j],...
Z: numpy.ndarray  # value = array([[ 1.+0.j,  0.+0.j],...
