from ._pyqpp import random_devices
import sys

# Explicitly re-export everything from the C++ submodule
# This allows IDEs to follow the path: pyqpp -> random_devices -> _pyqpp
__all__ = [name for name in dir(random_devices) if not name.startswith("_")]

# Map the internal names so they are accessible directly in this module
globals().update({name: getattr(random_devices, name) for name in __all__})

# This "becomes" the C++ submodule in the eyes of the Python interpreter
sys.modules[__name__] = random_devices
