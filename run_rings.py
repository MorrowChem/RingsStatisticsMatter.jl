import argparse
from julia_rings.rings import ring_statistics
from ase.io import read
import time
import numpy as np

parser = argparse.ArgumentParser(description="Ring Statistics Script")
parser.add_argument("input_file", help="Input file in ASE-supported format")
args = parser.parse_args()

ats = read(args.input_file, "-1")

st = time.time()
rs, rings = ring_statistics(ats)
np.savetxt('rings_stats.npy', rs)
et = time.time()

elapsed = et - st
print(f"Elapsed time: {elapsed:.2f} seconds")
