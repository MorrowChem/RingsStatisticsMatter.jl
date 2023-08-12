import numpy as np
from ase.io import read, write
from ase.neighborlist import NeighborList, NewPrimitiveNeighborList
import json
from os.path import abspath, dirname, join
import julia
from julia_rings.utilities import halton_sequence, place_points_in_cube

# Initialize Julia
# julia.install() # you may have to run this the first time

# Create a Julia instance
jl = julia.Julia(compiled_modules=False)
import julia.Main as Main
Main.include(join(abspath(dirname(__file__)), "rings.jl"))


def ring_statistics(ats, refnodes='auto', index='-1', cutoff=2.85):
    
    if type(ats) == str:
        ats = read(ats, index=index)

    nl = NeighborList(cutoffs=cutoff, self_interaction=False, bothways=True,
                  primitive=NewPrimitiveNeighborList)
    nl.update(ats)
    # + 1 for Julia indexing
    neighs = [nl.get_neighbors(i)[0]+1 for i in range(len(ats))]

    if refnodes == 'auto':
        # Define the number of points you want to place
        vol = ats.get_volume()
        num_points = int(np.ceil((vol/(20**3)))) + 2
        print('Num refnodes requested: ', num_points)

        # Place points maximally separated in the periodic unit cube
        points = place_points_in_cube(num_points, ats)

        spos = ats.get_scaled_positions()
        refnodes = []

        for v in points:
            diff = spos - v
            ref = np.argmin(np.linalg.norm(diff, axis=1))
            refnodes.append(ref+1) # Julia indexing
        
        refnodes.append(1) # just a convenience placeholder, not used
        refnodes = np.array(refnodes)
        print('Refnodes chosen: ', refnodes[:-1])
    
    else:
        refnodes = np.array(refnodes)

    results = Main.rings.ring_statistics(len(ats), neighs, refnodes)
    rs, ngf, rings = results

    return rs, rings