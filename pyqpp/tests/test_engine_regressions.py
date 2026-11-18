"""Regression tests for pyqpp engine lifetime and execution."""

import copy
import gc
import subprocess
import sys
from pathlib import Path

import pytest

from pyqpp import (
    QCircuit,
    QEngine,
    QDensityEngine,
    QNoisyEngine,
    QDensityNoisyEngine,
    QubitDepolarizingNoise,
    gates,
)

ENGINE_TYPES = (
    "ket",
    "density",
    "noisy_ket",
    "noisy_density",
)


def make_circuit():
    qc = QCircuit(2, 2)
    qc.gate(gates.H, 0)
    qc.CTRL(gates.X, 0, 1)
    qc.measure_all()
    return qc


def make_engine(kind, qc):
    if kind == "ket":
        return QEngine(qc)
    if kind == "density":
        return QDensityEngine(qc)
    if kind == "noisy_ket":
        return QNoisyEngine(qc, QubitDepolarizingNoise(0.01))
    if kind == "noisy_density":
        return QDensityNoisyEngine(qc, QubitDepolarizingNoise(0.01))
    raise ValueError(kind)


def total_stats(engine):
    return sum(engine.get_stats().values())


def case_constructor_lifetime(kind):
    qc = make_circuit()
    engine = make_engine(kind, qc)

    del qc
    gc.collect()

    engine.execute(100)
    assert total_stats(engine) == 100


def case_copy_lifetime(kind):
    qc = make_circuit()
    engine = make_engine(kind, qc)

    shallow = copy.copy(engine)
    deep = copy.deepcopy(engine)

    del qc, engine
    gc.collect()

    shallow.execute(100)
    deep.execute(100)

    assert total_stats(shallow) == 100
    assert total_stats(deep) == 100


def case_chained_copy(kind):
    qc = make_circuit()
    engine = make_engine(kind, qc)

    first = copy.copy(engine)
    second = copy.copy(first)

    del qc, engine, first
    gc.collect()

    second.execute(100)
    assert total_stats(second) == 100


def case_copy_after_execution(kind):
    qc = make_circuit()
    engine = make_engine(kind, qc)

    engine.execute(100)

    shallow = copy.copy(engine)
    deep = copy.deepcopy(engine)

    assert shallow.get_stats() == engine.get_stats()
    assert deep.get_stats() == engine.get_stats()

    del qc, engine
    gc.collect()

    for copied in (shallow, deep):
        copied.reset()
        copied.execute(100)
        assert total_stats(copied) == 100


def case_reset_and_reexecute(kind):
    qc = make_circuit()
    engine = make_engine(kind, qc)

    engine.execute(100)
    assert total_stats(engine) == 100

    engine.reset(None, False)
    engine.execute(100)
    assert total_stats(engine) == 200

    engine.reset()
    assert total_stats(engine) == 0

    engine.execute(100)
    assert total_stats(engine) == 100


def case_empty_circuit(kind):
    qc = QCircuit(2, 2)
    engine = make_engine(kind, qc)

    engine.execute(1)
    engine.execute(100)

    assert total_stats(engine) == 0


def case_repetitions(kind):
    qc = make_circuit()
    engine = make_engine(kind, qc)

    engine.execute(1000)

    assert total_stats(engine) == 1000

    if kind in ("ket", "density"):
        stats = engine.get_stats()

        assert set(stats).issubset({"00", "11"})
        assert 400 < stats.get("00", 0) < 600
        assert 400 < stats.get("11", 0) < 600


CASES = {
    "constructor_lifetime": case_constructor_lifetime,
    "copy_lifetime": case_copy_lifetime,
    "chained_copy": case_chained_copy,
    "copy_after_execution": case_copy_after_execution,
    "reset_and_reexecute": case_reset_and_reexecute,
    "empty_circuit": case_empty_circuit,
    "repetitions": case_repetitions,
}


@pytest.mark.parametrize("kind", ENGINE_TYPES)
@pytest.mark.parametrize("case", CASES)
def test_engine_regression(case, kind):
    # Run in a subprocess so native crashes are isolated.
    result = subprocess.run(
        [
            sys.executable,
            "-X",
            "faulthandler",
            str(Path(__file__).resolve()),
            case,
            kind,
        ],
        capture_output=True,
        text=True,
        timeout=60,
    )

    assert result.returncode == 0, (
        f"{case}/{kind} failed "
        f"(exit code {result.returncode})\n"
        f"{result.stdout}\n{result.stderr}"
    )


if __name__ == "__main__":
    case, kind = sys.argv[1:]
    CASES[case](kind)
