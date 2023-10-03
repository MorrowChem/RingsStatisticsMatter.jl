import argparse
from julia_rings.rings import ring_statistics
from ase.io import read
import time
import numpy as np

parser = argparse.ArgumentParser(description="Ring Statistics Script")
parser.add_argument("input_file", help="Input file in ASE-supported format")
parser.add_argument("-hp", "--homopolar", help="Whether to include homopolar bonds",
                    choices=['y', 'n'], default='y')
args = parser.parse_args()

if args.homopolar == 'y':
    homopolar = True
else:
    homopolar = False

ats = read(args.input_file, "-1")

if homopolar:
    cutoff = {(32, 32): 3.0, (32, 51): 3.4, (32, 52): 3.4, (51,52): 3.0}
else:
    cutoff = {(32, 32): 0.0, (32, 51): 3.4, (32, 52): 3.4, (51,52): 3.0}

st = time.time()
rs, rings = ring_statistics(ats, cutoff=cutoff)
np.savetxt('rings_stats.npy', rs)
rings_array = np.array(rings, dtype=object)
np.save('rings.npy', rings_array)
et = time.time()

elapsed = et - st
print(f"Elapsed time: {elapsed:.2f} seconds")
