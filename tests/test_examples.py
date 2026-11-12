from rings_statistics_matter import ring_statistics


def test_gst_homopolar(structure):
    cutoff = {(32, 32): 3.0, (32, 51): 3.4, (32, 52): 3.4, (51, 52): 3.0}
    counts, _ = ring_statistics(structure("GST_5k_amorphous.xyz"), cutoff=cutoff)
    assert counts.sum() == 987 and counts[6] == 37


def test_gst_no_ge_ge_bonds(structure):
    cutoff = {(32, 32): 0.0, (32, 51): 3.4, (32, 52): 3.4, (51, 52): 3.0}
    counts, _ = ring_statistics(structure("GST_5k_amorphous.xyz"), cutoff=cutoff)
    assert counts.sum() == 944 and counts[6] == 29
