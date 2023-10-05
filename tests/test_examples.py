import unittest
import os
from ase.io import read
from julia_rings.rings import ring_statistics

TEST_DIR = os.path.dirname(os.path.abspath(__file__))

class TestExamples(unittest.TestCase):
    
    def setUp(self):
        self.ats = read(os.path.join(TEST_DIR, "../structures/GST_5k_amorphous.xyz"))
        self.output = [os.path.join(TEST_DIR, "tmp_rings_stats.npy")]
        
        
    def tearDown(self):
        try:
            for i in self.output:
                os.remove(i)
        except FileNotFoundError:
            pass
    
    def test_hompolar(self):
        cutoff = {(32, 32): 3.0, (32, 51): 3.4,
                  (32, 52): 3.4, (51,52): 3.0}
        
        rs, _ = ring_statistics(self.ats, cutoff=cutoff)
        assert(sum(rs)==987.0 and rs[6]==37.)
        
    
    def test_nonhomopolar(self):
        cutoff = {(32, 32): 0.0, (32, 51): 3.4,
                  (32, 52): 3.4, (51,52): 3.0}
        
        rs, _ = ring_statistics(self.ats, cutoff=cutoff)
        print(sum(rs), rs[6])
        assert(sum(rs)==944. and rs[6]==29.)