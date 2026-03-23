from __future__ import annotations
import typing
__all__: list[str] = ['get_prng', 'load', 'mt19937', 'save', 'seed']
class mt19937:
    def __call__(self) -> int:
        """
        Generate a random number
        """
    def seed(self, s: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
def get_prng() -> mt19937:
    """
    Returns a reference to the internal Mersenne Twister (mt19937) object
    """
def load(state: str) -> None:
    """
    Loads the PRNG state from a string
    """
def save() -> str:
    """
    Saves the current state of the PRNG to a string
    """
def seed(seed_val: typing.SupportsInt | typing.SupportsIndex) -> None:
    """
    Seeds the internal Mersenne Twister PRNG
    """
