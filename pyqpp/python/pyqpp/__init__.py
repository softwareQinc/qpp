from ._pyqpp import *

from . import gates as gates
from . import states as states
from . import qasm as qasm

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
    # Prefer explicit export list from the extension if available
    from ._pyqpp import __all__ as _core_all
except ImportError:
    # Fallback: export everything that does not start with "_"
    _core_all = [name for name in dir() if not name.startswith("_")]

__all__ = list(_core_all) + [
    "gates",
    "states",
    "qasm",
    "__version__",
    "QPP_VERSION_STR",
    "QPP_VERSION_NUM",
]
