import numpy as np
from ase.io import read, write
from ase.neighborlist import NeighborList, NewPrimitiveNeighborList, neighbor_list, natural_cutoffs
from os.path import abspath, dirname, join
import julia
from julia_rings.utilities import halton_sequence, place_points_in_cube

# Initialize Julia
# julia.install() # you may have to run this the first time

# Create a Julia instance
jl = julia.Julia(compiled_modules=False)
import julia.Main as Main
Main.include(join(abspath(dirname(__file__)), "rings.jl"))


def ring_statistics(ats, refnodes='auto', index='-1', cutoff=None,
                    maxlevel=12, no_supercell=False, verbosity=0, **kwargs):
    """Python wrapper for calling Julia ring statistics function

    Args:
        ats ase.atoms.Atoms: Atoms object or ASE-readable filename representing a single structure
        refnodes (str, optional or list): 'auto' or list of reference nodes for distance map. 
                                          'auto' places a number of points scaling linearly with len(ats)
                                          fairly evenly-spaced on non-special positions
        index (str, optional): index of structure if filename specified for ats. Default to last structure
        cutoff(float or dict, optional):
            Cutoff for neighbour search (ASE NeighborList). It can be (from ASE docs)
            * A single float: This is a global bond cutoff for all elements.
            * A dictionary: This specifies cutoff values for element
              pairs. Specification accepts element numbers or symbols.
              Example: {(1, 6): 1.1, (1, 1): 1.0, ('C', 'C'): 1.85}
            * A list/array with a per atom value: This specifies the radius of
              an atomic sphere for each atoms. If spheres overlap, atoms are
              within each others neighborhood, i.e. this should be about half the
              value of the other cutoff methods above.
        maxlevel (Int, optional): Max search depth when finding rings. Rings of size <= 2*maxlevel will be found.
                                  Exponentially more expensive with increasing maxlevel 
    kwargs: addtional keyword arguments for fine-control of ring statistics. See Julia code for details
        
    Returns:
        rs (np.ndarray): Array of ring sizes from 1-N
        rings (list of NxM np.ndarray): Array of node indices for every ring found, separated into sublists by size
    """

    if cutoff is None:
        cutoff = [i*1.3 for i in natural_cutoffs(ats)]
        print("Cutoffs auto-set to: ", np.unique(cutoff))
        
    if type(cutoff) == dict:
        maxcutoff = max(cutoff.values())
    elif type(cutoff) == float:
        maxcutoff = cutoff
    elif type(cutoff) in (np.ndarray, list):
        maxcutoff = max(cutoff)
    else:
        raise TypeError(f'Cutoff type {type(cutoff)} not supported.')
    
    if type(ats) == str:
        ats = read(ats, index=index)
    
    rsf = None
    if np.min(ats.cell.diagonal()) < maxlevel*maxcutoff and not no_supercell:
        rsf = int(np.floor(maxlevel*maxcutoff/np.min(ats.cell.diagonal())))
        ats = ats * [rsf for i in range(3)]
        print(f"WARNING: Cell too small for maxlevel. Scaling cell by {rsf} x {rsf} x {rsf}")
        print(f"Now calculating for {len(ats)} atoms")
        print("Ring lists will include nodes from supercells.")
        kwargs['rsf'] = rsf**3 # for scaling rs back to original cell size

    if type(cutoff) in (float, np.ndarray, list):       
        nl = NeighborList(cutoffs=cutoff, skin=0.0, self_interaction=False, bothways=True,
                      primitive=NewPrimitiveNeighborList)
        nl.update(ats)
        # + 1 for Julia indexing
        neighs = [nl.get_neighbors(i)[0]+1 for i in range(len(ats))]
        
    elif type(cutoff) == dict:
        nl = neighbor_list('ij', ats, cutoff)
        neighs = [[] for i in range(len(ats))]
        for ct in range(len(nl[0])):
            # + 1 for Julia indexing
            neighs[nl[0][ct]].append(nl[1][ct]+1)
        neighs = [np.array(i, dtype=np.int64) for i in neighs]
        
    else:
        raise NotImplementedError(f'Cutoff type {type(cutoff)} not supported.'
                                  'Try a float, array of floats, or dict')
        

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
            refnode = np.argmin(np.linalg.norm(diff, axis=1))
            refnodes.append(refnode+1) # Julia indexing
        
        refnodes.append(1) # just a convenience padding element, not used
        refnodes = np.array(refnodes)
        if verbosity > 0:
            print('Refnodes chosen: ', refnodes[:-1])
    
    else:
        refnodes = np.array(refnodes)

    results = Main.rings.ring_statistics(len(ats), neighs, refnodes, **kwargs)
    rs, ngf, rings = results
    if not rsf is None:
        rs = rs/rsf**3 # undo supercell scaling

    return rs, rings