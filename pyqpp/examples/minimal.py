import pyqpp.states as st
from pyqpp import dirac


def main() -> None:
    print("Hello Quantum++!\nThis is the |0> state:\n")
    print(dirac(st.z0))


if __name__ == "__main__":
    main()
