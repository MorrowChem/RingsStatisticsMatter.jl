from collections.abc import Mapping, Sequence

import numpy as np
from ase import Atoms
from ase.data import atomic_numbers
from vesin import NeighborList

Cutoff = float | Sequence[float] | np.ndarray | Mapping[tuple[int | str, int | str], float]


def max_pair_cutoff(cutoff: Cutoff) -> float:
    """Largest bond length the cutoff specification can accept."""
    if isinstance(cutoff, Mapping):
        return float(max(cutoff.values()))
    if np.ndim(cutoff) == 0:
        return float(cutoff)
    return 2.0 * float(np.max(cutoff))


def _pair_cutoff_table(cutoff: Mapping[tuple[int | str, int | str], float]) -> np.ndarray:
    numbers = {pair: tuple(atomic_numbers[x] if isinstance(x, str) else int(x) for x in pair) for pair in cutoff}
    table = np.zeros((max(max(z) for z in numbers.values()) + 1,) * 2)
    for pair, (z1, z2) in numbers.items():
        table[z1, z2] = table[z2, z1] = cutoff[pair]
    return table


def neighbour_csr(atoms: Atoms, cutoff: Cutoff) -> tuple[np.ndarray, np.ndarray]:
    """Bonded neighbours of every atom in CSR form.

    Args:
        atoms: Structure.
        cutoff: A float is a global bond cutoff. A per-atom sequence gives sphere radii,
            atoms bond when spheres overlap. A dict maps element pairs (numbers or symbols)
            to bond cutoffs, unlisted pairs do not bond.

    Returns:
        offsets (n_atoms + 1) and indices, both int64: the neighbours of atom i are
        indices[offsets[i]:offsets[i + 1]].
    """
    i, j, d = NeighborList(cutoff=max_pair_cutoff(cutoff), full_list=True).compute(
        points=atoms.positions, box=atoms.cell.array, periodic=atoms.pbc, quantities="ijd"
    )
    if isinstance(cutoff, Mapping):
        table = _pair_cutoff_table(cutoff)
        z = np.minimum(atoms.numbers, len(table) - 1)
        listed = atoms.numbers < len(table)
        bonded = (d < table[z[i], z[j]]) & listed[i] & listed[j]
    elif np.ndim(cutoff) == 0:
        bonded = d < cutoff
    else:
        radii = np.asarray(cutoff, dtype=float)
        if radii.shape != (len(atoms),):
            raise ValueError(f"Per-atom cutoff has shape {radii.shape}, expected ({len(atoms)},)")
        bonded = d < radii[i] + radii[j]

    i, j = i[bonded], j[bonded]
    order = np.argsort(i, kind="stable")
    offsets = np.zeros(len(atoms) + 1, dtype=np.int64)
    np.cumsum(np.bincount(i, minlength=len(atoms)), out=offsets[1:])
    return offsets, j[order].astype(np.int64)
