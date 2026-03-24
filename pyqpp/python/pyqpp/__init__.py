"""
pyqpp: A Python wrapper for the Quantum++ library.
"""

# Metadata and versioning logic
__version__ = "7.0.3"


def _version_to_number(version: str) -> int:
    try:
        # Split by '-' to handle pre-releases like '7.0.2-beta'
        parts = version.split("-")[0].split(".")
        major, minor, patch = map(int, parts)
        return major * 10000 + minor * 100 + patch
    except (ValueError, IndexError):
        return 0


QPP_VERSION_STR = __version__
QPP_VERSION_NUM = _version_to_number(__version__)

# Binary extension loading
try:
    from . import _pyqpp
except ImportError as e:
    raise ImportError(f"Failed to load pyqpp binary: {e}") from e

# Handle exports from the binary
_ext_exports = getattr(
    _pyqpp, "__all__", [k for k in dir(_pyqpp) if not k.startswith("_")]
)
globals().update({k: getattr(_pyqpp, k) for k in _ext_exports})

# Submodule imports
from . import codes, gates, qasm, random_devices, states

# Define __all__
_submodules = ["codes", "gates", "qasm", "random_devices", "states"]
_metadata = ["__version__", "QPP_VERSION_STR", "QPP_VERSION_NUM"]
__all__ = list(dict.fromkeys([*_submodules, *_ext_exports, *_metadata]))
