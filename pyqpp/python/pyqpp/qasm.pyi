from __future__ import annotations
import pyqpp._pyqpp
__all__: list[str] = ['read_from_file', 'read_from_string']
def read_from_file(arg0: str) -> pyqpp._pyqpp.QCircuit:
    """
    Reads an OpenQASM circuit from a file and returns a qpp::QCircuit object
    """
def read_from_string(arg0: str) -> pyqpp._pyqpp.QCircuit:
    """
    Reads an OpenQASM circuit from a string and returns a qpp::QCircuit object
    """
