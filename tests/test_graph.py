import numpy as np
import pytest

from rings_statistics_matter import ring_statistics, ring_statistics_graph
from rings_statistics_matter.neighbours import neighbour_csr

# One node per cell
PCU = ([[0, 0]] * 3, [[1, 0, 0], [0, 1, 0], [0, 0, 1]])
# Two nodes per primitive fcc cell
DIA = ([[0, 1]] * 4, [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])


def per_cell(counts):
    return {size + 1: count for size, count in enumerate(counts) if count}


def test_finite_graph():
    hexagon_with_tail = [[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 0], [0, 6]]
    counts, rings = ring_statistics_graph(hexagon_with_tail)
    assert per_cell(counts) == {6: 1}
    assert sorted(rings[5][0]) == [0, 1, 2, 3, 4, 5]


def test_fused_squares_have_no_six_ring():
    counts, _ = ring_statistics_graph([[0, 1], [1, 2], [2, 3], [3, 0], [1, 4], [4, 5], [5, 2]])
    assert per_cell(counts) == {4: 2}


@pytest.mark.parametrize(("net", "expected"), [(PCU, {4: 3, 6: 4}), (DIA, {6: 4})])
def test_periodic_nets(net, expected):
    edges, images = net
    counts, rings = ring_statistics_graph(edges, images=images, maxlvl=8)
    assert per_cell(counts) == expected
    assert all(r.max() <= np.max(edges) for r in rings if r.size)


def test_longer_images_give_the_same_net():
    edges, images = PCU
    counts, _ = ring_statistics_graph(edges, images=2 * np.array(images), maxlvl=8)
    assert per_cell(counts) == {4: 3, 6: 4}


def test_matches_atoms_route(structure):
    atoms = structure("POSCAR-SiO2-a-quartz")
    reference, _ = ring_statistics(atoms, cutoff=1.8, maxlvl=14)

    repeats = (4, 4, 4)
    big = atoms.repeat(repeats)
    offsets, indices = neighbour_csr(big, 1.8)
    source = np.repeat(np.arange(len(big)), np.diff(offsets))
    home = source < len(atoms)
    cell = indices[home] // len(atoms)
    images = np.stack(np.unravel_index(cell, repeats), axis=1)
    images = np.where(images > 2, images - 4, images)
    edges = np.stack([source[home], indices[home] % len(atoms)], axis=1)

    counts, _ = ring_statistics_graph(edges, n_nodes=len(atoms), images=images, maxlvl=14)
    assert per_cell(counts) == per_cell(reference)


def test_rejects_bad_input():
    with pytest.raises(ValueError, match="node indices"):
        ring_statistics_graph([[0, 3]], n_nodes=2)
    with pytest.raises(ValueError, match="one row per edge"):
        ring_statistics_graph([[0, 1]], images=[[0, 0, 0], [1, 0, 0]])
