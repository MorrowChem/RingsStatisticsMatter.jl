import logging
import warnings
from collections.abc import Mapping, Sequence
from typing import Literal

import numpy as np
from ase import Atoms
from ase.neighborlist import natural_cutoffs

from rings_statistics_matter import _kernel, julia_backend
from rings_statistics_matter.neighbours import Cutoff, max_pair_cutoff, neighbour_csr
from rings_statistics_matter.refnodes import auto_refnodes

logger = logging.getLogger(__name__)

NATURAL_CUTOFF_SCALE = 1.3


def _supercell_repeats(atoms: Atoms, span: float) -> np.ndarray:
    """Repeats along each periodic cell vector so that every cell width is at least span."""
    cell = atoms.cell.array
    normals = np.cross(cell[[1, 2, 0]], cell[[2, 0, 1]])
    widths = atoms.get_volume() / np.linalg.norm(normals, axis=1)
    return np.where(atoms.pbc, np.ceil(span / widths), 1).astype(int)


def ring_statistics(
    atoms: Atoms,
    cutoff: Cutoff | None = None,
    maxlvl: int = 12,
    refnodes: Sequence[int] | Literal["auto"] = "auto",
    supercell: bool = True,
    deduplicate: bool = True,
    backend: Literal["rust", "julia"] = "rust",
    num_threads: int | None = None,
    mxpths: int = 1000,
) -> tuple[np.ndarray, list[np.ndarray]]:
    """Primitive ring statistics (Yuan & Cormack 2002).

    Args:
        atoms: Structure.
        cutoff: Bond cutoff, see `neighbour_csr`. Defaults to 1.3 x the ASE natural radii.
        maxlvl: Search depth. Rings of up to maxlvl - 1 nodes are found. Cost grows
            exponentially with it.
        refnodes: Atom indices for the referential distance maps, or "auto".
        supercell: Repeat cells narrower than 2 * maxlvl * cutoff, so a ring cannot meet
            its own image. Counts and node indices are mapped back to the input cell.
        deduplicate: Return each ring once. Otherwise once per node it contains.
        backend: Kernel implementation. "julia" needs the julia extra.
        num_threads: Rust backend only, defaults to all cores. Julia threads are set by
            JULIA_NUM_THREADS before first use.
        mxpths: Julia backend only. Buffer for shortest paths from one source to one node.

    Returns:
        counts: counts[n - 1] is the number of n-membered rings in the input cell.
        rings: rings[n - 1] is an (n_rings, n) array of atom indices.
    """
    if cutoff is None:
        cutoff = NATURAL_CUTOFF_SCALE * np.asarray(natural_cutoffs(atoms))
    per_atom = not isinstance(cutoff, Mapping) and np.ndim(cutoff) == 1

    replicas = 1
    if supercell:
        repeats = _supercell_repeats(atoms, span=2 * maxlvl * max_pair_cutoff(cutoff))
        replicas = int(np.prod(repeats))
        if replicas > 1:
            logger.info("Cell too small for maxlvl=%d, repeating %s: %d atoms", maxlvl, repeats, len(atoms) * replicas)
            atoms = atoms.repeat(repeats)
            if per_atom:
                cutoff = np.tile(cutoff, replicas)

    offsets, indices = neighbour_csr(atoms, cutoff)
    refnodes = auto_refnodes(atoms) if isinstance(refnodes, str) else np.asarray(refnodes, dtype=np.int64)

    if backend == "rust":
        counts, rings = _kernel.ring_statistics_csr(
            offsets,
            indices,
            refnodes,
            replicas=replicas,
            maxlvl=maxlvl,
            deduplicate=deduplicate,
            num_threads=num_threads,
        )
    elif backend == "julia":
        counts, rings = julia_backend.ring_statistics_csr(
            offsets, indices, refnodes, replicas=replicas, maxlvl=maxlvl, deduplicate=deduplicate, mxpths=mxpths
        )
    else:
        raise ValueError(f"Unknown backend {backend!r}")
    if not np.allclose(counts, np.round(counts)):
        warnings.warn(
            "Non-integer ring counts: rings were missed or double counted. Use a larger cell or supercell=True.",
            stacklevel=2,
        )
    return counts, rings
