from .pyqpp import *

__version__ = "7.0.0"


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
