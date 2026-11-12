import numpy as np
from ase import Atoms
from scipy.spatial import cKDTree

REFNODE_SPACING = 20.0  # Å, one reference node per cube of this edge


def _halton(num_points: int, offset: float) -> np.ndarray:
    """Quasi-uniform points in the unit cube, shifted off special positions."""
    points = np.full((num_points, 3), offset)
    for dim, base in enumerate((2, 3, 4)):
        for n in range(num_points):
            f, i = 1.0 / base, n
            while i > 0:
                points[n, dim] += f * (i % base)
                i //= base
                f /= base
    return points


def auto_refnodes(atoms: Atoms) -> np.ndarray:
    """Indices of atoms spread evenly through the cell, for the referential distance maps."""
    cubes = atoms.get_volume() / REFNODE_SPACING**3
    points = _halton(int(np.ceil(cubes)) + 2, offset=0.2 * cubes ** (-1 / 3))
    _, nearest = cKDTree(atoms.get_scaled_positions()).query(points)
    return nearest.astype(np.int64)
