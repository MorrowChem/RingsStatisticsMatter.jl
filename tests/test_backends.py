import numpy as np
import pytest

from rings_statistics_matter import ring_statistics
from tests.conftest import needs_julia


def _canonical(ring):
    ring = list(ring)
    return min(tuple(d[s:] + d[:s]) for d in (ring, ring[::-1]) for s in range(len(ring)))


def _ring_sets(rings):
    return [sorted(map(_canonical, by_size)) for by_size in rings]


@needs_julia
@pytest.mark.parametrize("kwargs", [{"maxlvl": 8}, {"maxlvl": 9}, {"supercell": False}])
def test_backends_agree(structure, kwargs):
    atoms = structure("aSi_500atom_test.xyz")
    rust_counts, rust_rings = ring_statistics(atoms, cutoff=2.85, backend="rust", **kwargs)
    julia_counts, julia_rings = ring_statistics(atoms, cutoff=2.85, backend="julia", **kwargs)
    np.testing.assert_array_equal(rust_counts, julia_counts)
    assert _ring_sets(rust_rings) == _ring_sets(julia_rings)


@pytest.mark.parametrize("refnodes", [[], [0], list(range(0, 4096, 7))])
def test_refnodes_do_not_change_result(structure, refnodes):
    atoms = structure("aSi_500atom_test.xyz")
    auto_counts, auto_rings = ring_statistics(atoms, cutoff=2.85)
    counts, rings = ring_statistics(atoms, cutoff=2.85, refnodes=refnodes)
    np.testing.assert_array_equal(counts, auto_counts)
    assert all(np.array_equal(a, b) for a, b in zip(rings, auto_rings, strict=True))


def test_rings_are_unique_and_bonded(structure):
    atoms = structure("aSi_500atom_test.xyz")
    counts, rings = ring_statistics(atoms, cutoff=2.85)
    for size, by_size in enumerate(rings, start=1):
        assert by_size.shape == (counts[size - 1], size)
        assert len(set(map(_canonical, by_size))) == len(by_size)
    distances = atoms.get_all_distances(mic=True)
    six = rings[5]
    assert (distances[six, np.roll(six, 1, axis=1)] < 2.85).all()


def test_without_deduplication_each_ring_appears_once_per_node(structure):
    atoms = structure("GST_5k_amorphous.xyz")
    cutoff = {(32, 32): 3.0, (32, 51): 3.4, (32, 52): 3.4, (51, 52): 3.0}
    counts, rings = ring_statistics(atoms, cutoff=cutoff, deduplicate=False, supercell=False)
    for size, by_size in enumerate(rings, start=1):
        assert len(by_size) == round(size * counts[size - 1])


def test_thread_count_does_not_change_result(structure):
    atoms = structure("aSi_500atom_test.xyz")
    one, rings_one = ring_statistics(atoms, cutoff=2.85, num_threads=1)
    many, rings_many = ring_statistics(atoms, cutoff=2.85, num_threads=4)
    np.testing.assert_array_equal(one, many)
    assert all(np.array_equal(a, b) for a, b in zip(rings_one, rings_many, strict=True))


def test_unknown_backend(structure):
    with pytest.raises(ValueError, match="Unknown backend"):
        ring_statistics(structure("POSCAR.mp-149_Si"), cutoff=2.85, backend="fortran")
