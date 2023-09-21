import unittest
import os
from ase.io import read
from julia_rings.rings import ring_statistics

TEST_DIR = os.path.dirname(os.path.abspath(__file__))

class TestCrystalsRings(unittest.TestCase):
    
    def setUp(self):
        self.output = os.path.join(TEST_DIR, "tmp_rings_stats.npy")
    
    def tearDown(self):
        os.remove(self.output)

    def test_Si(self):
        ats = read(os.path.join(TEST_DIR, "../structures/POSCAR.mp-149_Si"))
        rs, _ = ring_statistics(ats, cutoff=2.85, outfile=self.output)
        assert(sum(rs)==16. and rs[5]==16.)
        
    def test_SiO2(self):
        ats = read(os.path.join(TEST_DIR, "../structures/POSCAR-SiO2-a-quartz"))
        rs, _ = ring_statistics(ats, cutoff=2.2, maxlvl=24, outfile=self.output)
        assert(sum(rs)==18. and rs[11]==3. and rs[15]==15.)
        
if __name__ == '__main__':
    unittest.main()