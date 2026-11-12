from itertools import pairwise

import numpy as np
import pytest
from ase import Atoms
from ase.neighborlist import neighbor_list

from rings_statistics_matter.neighbours import max_pair_cutoff, neighbour_csr
from rings_statistics_matter.rings import _supercell_repeats


def _sorted_rows(offsets, indices):
    return [sorted(indices[a:b]) for a, b in pairwise(offsets)]


def _ase_rows(atoms, cutoff):
    i, j = neighbor_list("ij", atoms, cutoff)
    return [sorted(j[i == k]) for k in range(len(atoms))]


def test_float_cutoff_matches_ase(structure):
    atoms = structure("aSi_500atom_test.xyz")
    assert _sorted_rows(*neighbour_csr(atoms, 2.85)) == _ase_rows(atoms, 2.85)


def test_per_atom_cutoff_matches_ase(structure):
    atoms = structure("GST_5k_amorphous.xyz")
    radii = np.where(atoms.numbers == 32, 1.4, 1.6)
    assert _sorted_rows(*neighbour_csr(atoms, radii)) == _ase_rows(atoms, list(radii))


def test_pair_cutoff_matches_ase(structure):
    atoms = structure("GST_5k_amorphous.xyz")
    cutoff = {("Ge", "Ge"): 0.0, (32, 51): 3.4, ("Ge", "Te"): 3.4, (51, 52): 3.0}
    assert _sorted_rows(*neighbour_csr(atoms, cutoff)) == _ase_rows(atoms, cutoff)


@pytest.mark.parametrize(
    ("cutoff", "expected"), [(2, 2.0), (np.float64(2.5), 2.5), ([1.0, 1.5], 3.0), ({(1, 1): 1.2}, 1.2)]
)
def test_max_pair_cutoff(cutoff, expected):
    assert max_pair_cutoff(cutoff) == expected


def test_supercell_repeats_use_cell_widths():
    sheared = Atoms("Si", cell=[[10, 0, 0], [9, 3, 0], [0, 0, 10]], pbc=[True, True, False])
    # widths are 3, 3 and 10 although the vector lengths are 10, 9.5 and 10
    np.testing.assert_array_equal(_supercell_repeats(sheared, span=12.0), [4, 4, 1])
