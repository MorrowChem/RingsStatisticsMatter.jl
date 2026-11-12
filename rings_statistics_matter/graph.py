from collections.abc import Sequence

import numpy as np

from rings_statistics_matter import _kernel

NUM_REFNODES = 8


def _periodic_edges(edges: np.ndarray, images: np.ndarray, n_nodes: int, repeats: np.ndarray) -> np.ndarray:
    """Edges of the graph repeated along each lattice direction, node index = cell * n_nodes + node."""
    cells = np.stack(np.meshgrid(*[np.arange(r) for r in repeats], indexing="ij"), axis=-1).reshape(-1, 3)
    strides = np.array([repeats[1] * repeats[2], repeats[2], 1])
    source = (cells @ strides)[:, None] * n_nodes + edges[:, 0]
    target = (((cells[:, None, :] + images) % repeats) @ strides) * n_nodes + edges[:, 1]
    return np.stack([source.ravel(), target.ravel()], axis=1)


def ring_statistics_graph(
    edges: Sequence[Sequence[int]],
    n_nodes: int | None = None,
    images: Sequence[Sequence[int]] | None = None,
    maxlvl: int = 12,
    deduplicate: bool = True,
    num_threads: int | None = None,
) -> tuple[np.ndarray, list[np.ndarray]]:
    """Primitive ring statistics of a graph given by its edges.

    Args:
        edges: (n_edges, 2) node indices. Direction and duplicates do not matter.
        n_nodes: Number of nodes, defaults to the largest index + 1.
        images: (n_edges, 3) lattice offsets for a periodic graph: edge k joins edges[k, 0] in
            the home cell to edges[k, 1] in cell images[k]. None for a finite graph.
        maxlvl: Search depth. Rings of up to maxlvl - 1 nodes are found.
        deduplicate: Return each ring once. Otherwise once per node it contains.
        num_threads: Defaults to all cores.

    Returns:
        counts: counts[n - 1] is the number of n-membered rings per cell.
        rings: rings[n - 1] is an (n_rings, n) array of node indices.
    """
    edges = np.asarray(edges, dtype=np.int64).reshape(-1, 2)
    if n_nodes is None:
        n_nodes = int(edges.max()) + 1 if len(edges) else 0
    if len(edges) and (edges.min() < 0 or edges.max() >= n_nodes):
        raise ValueError(f"edges must hold node indices below {n_nodes}")

    replicas = 1
    if images is not None:
        images = np.asarray(images, dtype=np.int64).reshape(-1, 3)
        if len(images) != len(edges):
            raise ValueError("images must have one row per edge")
        # A closed walk around the repeated cell then has at least maxlvl steps
        repeats = np.maximum(1, maxlvl * np.abs(images).max(axis=0, initial=0))
        replicas = int(np.prod(repeats))
        edges = _periodic_edges(edges, images, n_nodes, repeats)

    total = n_nodes * replicas
    edges = edges[edges[:, 0] != edges[:, 1]]
    directed = np.unique(np.concatenate([edges, edges[:, ::-1]]), axis=0)
    offsets = np.concatenate([[0], np.cumsum(np.bincount(directed[:, 0], minlength=total))]).astype(np.int64)
    refnodes = np.unique(np.linspace(0, max(total - 1, 0), NUM_REFNODES).astype(np.int64))
    return _kernel.ring_statistics_csr(
        offsets,
        np.ascontiguousarray(directed[:, 1]),
        refnodes,
        replicas=replicas,
        maxlvl=maxlvl,
        deduplicate=deduplicate,
        num_threads=num_threads,
    )
