import argparse
from julia_rings.rings import ring_statistics
from ase.io import read
import time
import numpy as np

parser = argparse.ArgumentParser(description="Ring Statistics Script")
parser.add_argument("input_file", help="Input file in ASE-supported format")
parser.add_argument("-v", "--verbosity", help="Verbosity level 0, 1", default=0)
parser.add_argument("-o", "--output", help="Output file name", default="rings_lists.json")
parser.add_argument("--maxpths", help="Maximum number of paths to consider at each node", default=1000)
parser.add_argument("--maxlvl", help="Rings of size up to 2*maxlvl can be found", default=12)
parser.add_argument("-c", "--cutoff", help="Cutoff for neighbour search. "
                    "This accepts a global float. Adjust run_rings.py for more control",
                    default=2.2)
args = parser.parse_args()

v = int(args.verbosity)

ats = read(args.input_file, "-1")

st = time.time()
rs, rings = ring_statistics(ats, verbosity=v, mxpths=int(args.maxpths),
                            outfile=args.output, maxlvl=int(args.maxlvl),
                            cutoff=float(args.cutoff))
np.savetxt('rings_stats.npy', rs)
et = time.time()

elapsed = et - st
print(f"Elapsed time: {elapsed:.2f} seconds")
