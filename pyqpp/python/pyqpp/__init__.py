from ._pyqpp import *

from . import codes as codes
from . import gates as gates
from . import qasm as qasm
from . import random_devices as random_devices
from . import states as states

__version__ = "7.0.1"


def _version_to_number(version: str) -> int:
    """
    Convert a semantic version string to a numeric value.

    Example:
        "7.0.0" -> 70000
        "7.1.3" -> 70103
    """
    parts = version.split(".")
    if len(parts) != 3:
        raise ValueError(f"Expected version 'MAJOR.MINOR.PATCH', got: {version}")

    major, minor, patch = map(int, parts)
    return major * 10000 + minor * 100 + patch


QPP_VERSION_STR = __version__
QPP_VERSION_NUM = _version_to_number(__version__)

# -----------------------------------------------------------------------------
# Public API
# -----------------------------------------------------------------------------
try:
    from ._pyqpp import __all__ as _extension_exports
except ImportError:
    _extension_exports = []

__all__ = list(_extension_exports) + [
    "codes",
    "gates",
    "qasm",
    "random_devices",
    "states",
    "__version__",
    "QPP_VERSION_STR",
    "QPP_VERSION_NUM",
]
