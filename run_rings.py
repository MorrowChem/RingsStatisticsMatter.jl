import argparse
import json
import time

import numpy as np
from ase.io import read

from rings_statistics_matter import ring_statistics

parser = argparse.ArgumentParser(description="Ring statistics for the last frame of a structure file")
parser.add_argument("input_file", help="Input file in ASE-supported format")
parser.add_argument("-o", "--output", help="JSON file for ring node lists", default="rings_lists.json")
parser.add_argument("--mxpths", help="Maximum number of paths to consider at each node", type=int, default=1000)
parser.add_argument("--maxlvl", help="Search depth: rings of up to maxlvl - 1 nodes are found", type=int, default=12)
parser.add_argument("-c", "--cutoff", help="Global bond cutoff (default: 1.3 x natural radii)", type=float)
parser.add_argument("--no-supercell", help="Do not repeat small cells", action="store_true")
args = parser.parse_args()

atoms = read(args.input_file, "-1")

start = time.time()
counts, rings = ring_statistics(
    atoms, cutoff=args.cutoff, maxlvl=args.maxlvl, supercell=not args.no_supercell, mxpths=args.mxpths
)
print(f"Elapsed time: {time.time() - start:.2f} seconds")

for size, count in enumerate(counts, start=1):
    if count:
        print(f"{size:9d} | {count:10.1f}")

np.savetxt("rings_stats.npy", counts)
with open(args.output, "w") as f:
    json.dump([r.tolist() for r in rings], f)
