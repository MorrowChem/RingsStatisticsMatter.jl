import argparse
import json
import time

import numpy as np
from ase.io import read

from rings_statistics_matter import ring_statistics

parser = argparse.ArgumentParser(description="Ring statistics for Ge-Sb-Te with element-pair cutoffs")
parser.add_argument("input_file", help="Input file in ASE-supported format")
parser.add_argument("-hp", "--homopolar", help="Whether to include Ge-Ge bonds", choices=["y", "n"], default="y")
args = parser.parse_args()

atoms = read(args.input_file, "-1")
cutoff = {(32, 32): 3.0 if args.homopolar == "y" else 0.0, (32, 51): 3.4, (32, 52): 3.4, (51, 52): 3.0}

start = time.time()
counts, rings = ring_statistics(atoms, cutoff=cutoff)
print(f"Elapsed time: {time.time() - start:.2f} seconds")

np.savetxt("rings_stats.npy", counts)
with open("rings_lists.json", "w") as f:
    json.dump([r.tolist() for r in rings], f)
