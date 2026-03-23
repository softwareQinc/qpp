from ._pyqpp import qasm
import sys

# Explicitly re-export everything from the C++ submodule
# This allows IDEs to follow the path: pyqpp -> qasm -> _pyqpp
__all__ = [name for name in dir(qasm) if not name.startswith("_")]

# Map the internal names so they are accessible directly in this module
globals().update({name: getattr(qasm, name) for name in __all__})

# This "becomes" the C++ submodule in the eyes of the Python interpreter
sys.modules[__name__] = qasm
