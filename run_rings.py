import argparse
from rings import ring_statistics
from ase.io import read

parser = argparse.ArgumentParser(description="Ring Statistics Script")
parser.add_argument("input_file", help="Input file in ASE-supported format")
args = parser.parse_args()

ats = read(args.input_file, "-1")

rs, rings = ring_statistics(ats)