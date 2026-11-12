from pathlib import Path

import pytest
from ase import Atoms
from ase.io import read

STRUCTURES = Path(__file__).parent.parent / "structures"


@pytest.fixture(scope="session")
def structure():
    def load(name: str) -> Atoms:
        return read(STRUCTURES / name)

    return load
