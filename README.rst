RingsStatisticsMatter
=====================

This package provides a tool for finding primitive rings efficiently for large models of amorphous materials.

Designed for large structural models, this can handle 100,000-atoms in under 10 minutes multithreaded on a 16GB M2 Mac laptop.

This package implements an algorithm based on the work described in:

| Yuan, X. & Cormack, A. N. 
| `"Efficient Algorithm for Primitive Ring Statistics in Topological Networks." <https://doi.org/10.1016/S0927-0256(01)00256-7>`_.
| Computational Materials Science, vol. **24**, pp. 343-360, 2002.  
|  

.. image:: https://raw.githubusercontent.com/MorrowChem/RingsStatisticsMatter/main/images/100k_9rings.jpg
   :width: 200
   :alt: 9-membered primitive rings in a 100k-atom model of a-Si

.. image:: https://raw.githubusercontent.com/MorrowChem/RingsStatisticsMatter/main/images/ring_stats_1M.png
   :width: 200
   :alt: Ring statistics of a 1M-atom model

100k-atom model of a-Si with the 9-membered primitive rings highlighted

.. image:: https://raw.githubusercontent.com/MorrowChem/RingsStatisticsMatter/main/images/8k_567rings.png
   :width: 400
   :alt: 5-, 6- and 7-membered primitive rings in an 8k-atom model of a-C

8k-atom model of a-C with the 5-, 6-, and 7-membered primitive rings highlighted

Features
--------

- Python wrappers for easy plotting and incorporation into more complex analysis frameworks.
- Returns node indices for each ring in your structure
- Flexible multi-element support (ignore or allow homopolar bonds, per-atom radii or element-pair cutoffs). Neighbours from vesin.
- Shared-memory parallelism in both backends (rayon in Rust, native threads in Julia).

Installation
------------

.. code-block:: bash

   pip install rings-statistics-matter

Wheels are built for Linux, macOS and Windows. Installing from source needs a `Rust toolchain <https://rustup.rs>`_.

Backends
--------

The same algorithm is implemented twice. Tests check that both return the same rings.

- ``backend="rust"`` (default): ``src/kernel.rs``. No runtime dependencies. Threads via ``num_threads``.
- ``backend="julia"``: ``rings_statistics_matter/rings.jl``. Install with the ``julia`` extra;
  `juliacall <https://juliapy.github.io/PythonCall.jl/stable/juliacall/>`_ installs Julia on first use.
  Threads via ``JULIA_NUM_THREADS``.

Usage
-----

.. code-block:: python

   from ase.io import read
   from rings_statistics_matter import ring_statistics

   counts, rings = ring_statistics(read("structures/aSi_500atom_test.xyz"), cutoff=2.85)

``counts[n - 1]`` is the number of n-membered rings, ``rings[n - 1]`` an ``(n_rings, n)`` array of atom indices.

``cutoff`` is a global bond length, per-atom radii, or a dict of element-pair bond lengths
(a pair set to 0 never bonds). The default is 1.3 x the ASE natural radii.

Command line:

.. code-block:: bash

   python run_rings.py -c 2.85 structures/aSi_500atom_test.xyz

Multi-component systems, such as phase-change Ge-Sb-Te without Ge-Ge bonds:

.. code-block:: bash

   python examples/run_rings_GeSbTe.py --homopolar n structures/GST_5k_amorphous.xyz

Tests
-----

.. code-block:: bash

   uv sync && uv run pytest

Disclaimer
==============
If you see non-integer entries in the ring statistics, or the nodes don't seem right, this is most likely from using small unit cells.   

We try to automatically create supercells to fix this, but there may be edge cases where this isn't enough. For large enough cells, this shouldn't be an issue.


Extra Arguments
===============

- ``maxlvl``: search depth. Rings of up to ``maxlvl - 1`` nodes are found.

- ``mxpths`` (Julia backend): buffer for shortest paths between a source and one node. Raise it if the kernel reports more paths than ``mxpths``.


Supported Platforms
-------------------

This package has been tested on macOS and Linux.

Releasing
---------

Bump ``version`` in ``pyproject.toml`` and ``Cargo.toml``, then push a ``v*`` tag. The Release workflow builds
the wheels and publishes them to PyPI by trusted publishing.

Contributing
------------

Contributions are welcome! Feel free to open issues or pull requests on the GitHub repository:

https://github.com/MorrowChem/RingsStatisticsMatter

License
-------

This project is licensed under the MIT License - see the LICENSE file for details.
