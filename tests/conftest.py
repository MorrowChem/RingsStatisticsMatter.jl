from importlib.util import find_spec
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


# Importing juliacall here would start Julia before the backend sets its signal handling.
needs_julia = pytest.mark.skipif(find_spec("juliacall") is None, reason="juliacall not installed")


@pytest.fixture(params=["rust", pytest.param("julia", marks=needs_julia)])
def backend(request):
    return request.param
