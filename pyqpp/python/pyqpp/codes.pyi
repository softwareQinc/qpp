from __future__ import annotations
import numpy
import numpy.typing
import typing
__all__: list[str] = ['FIVE_QUBIT', 'SHOR_NINE_QUBIT', 'STEANE_SEVEN_QUBIT', 'Type', 'codeword']
class Type:
    """
    Members:
    
      FIVE_QUBIT
    
      STEANE_SEVEN_QUBIT
    
      SHOR_NINE_QUBIT
    """
    FIVE_QUBIT: typing.ClassVar[Type]  # value = <Type.FIVE_QUBIT: 0>
    SHOR_NINE_QUBIT: typing.ClassVar[Type]  # value = <Type.SHOR_NINE_QUBIT: 2>
    STEANE_SEVEN_QUBIT: typing.ClassVar[Type]  # value = <Type.STEANE_SEVEN_QUBIT: 1>
    __members__: typing.ClassVar[dict[str, Type]]  # value = {'FIVE_QUBIT': <Type.FIVE_QUBIT: 0>, 'STEANE_SEVEN_QUBIT': <Type.STEANE_SEVEN_QUBIT: 1>, 'SHOR_NINE_QUBIT': <Type.SHOR_NINE_QUBIT: 2>}
    def __eq__(self, other: typing.Any) -> bool:
        ...
    def __getstate__(self) -> int:
        ...
    def __hash__(self) -> int:
        ...
    def __index__(self) -> int:
        ...
    def __init__(self, value: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    def __int__(self) -> int:
        ...
    def __ne__(self, other: typing.Any) -> bool:
        ...
    def __repr__(self) -> str:
        ...
    def __setstate__(self, state: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    def __str__(self) -> str:
        ...
    @property
    def name(self) -> str:
        ...
    @property
    def value(self) -> int:
        ...
def codeword(type: Type, i: typing.SupportsInt | typing.SupportsIndex) -> typing.Annotated[numpy.typing.NDArray[numpy.complex128], "[m, n]"]:
    """
    Returns the i-th codeword of the specified quantum code type
    """
FIVE_QUBIT: Type  # value = <Type.FIVE_QUBIT: 0>
SHOR_NINE_QUBIT: Type  # value = <Type.SHOR_NINE_QUBIT: 2>
STEANE_SEVEN_QUBIT: Type  # value = <Type.STEANE_SEVEN_QUBIT: 1>
