from . import _pyqpp as _core


def __getattr__(name: str):
    return getattr(_core.qasm, name)


def __dir__():
    return dir(_core.qasm)
