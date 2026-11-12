import numpy as np

from rings_statistics_matter import ring_statistics


def test_diamond_si(structure):
    atoms = structure("POSCAR.mp-149_Si")
    counts, rings = ring_statistics(atoms, cutoff=2.85)
    assert counts.sum() == 16 and counts[5] == 16
    assert rings[5].shape == (16, 6)
    assert rings[5].max() < len(atoms)


def test_alpha_quartz(structure):
    counts, _ = ring_statistics(structure("POSCAR-SiO2-a-quartz"), cutoff=1.8, maxlvl=20)
    assert counts.sum() == 18 and counts[11] == 3 and counts[15] == 15


def test_per_atom_cutoff_with_supercell(structure):
    atoms = structure("POSCAR.mp-149_Si")
    counts, _ = ring_statistics(atoms, cutoff=np.full(len(atoms), 1.4))
    assert counts.sum() == 16 and counts[5] == 16


def test_default_cutoff(structure):
    counts, _ = ring_statistics(structure("POSCAR.mp-149_Si"))
    assert counts[5] == 16


def test_refnodes_do_not_change_result(structure):
    atoms = structure("aSi_500atom_test.xyz")
    auto, _ = ring_statistics(atoms, cutoff=2.85, maxlvl=8)
    single, _ = ring_statistics(atoms, cutoff=2.85, maxlvl=8, refnodes=[0])
    np.testing.assert_array_equal(auto, single)
