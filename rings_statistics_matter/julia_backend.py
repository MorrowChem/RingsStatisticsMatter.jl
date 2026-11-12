import os
import sys
import warnings
from functools import cache
from pathlib import Path

import numpy as np

KERNEL = Path(__file__).with_name("rings.jl")


@cache
def _kernel():
    """Start Julia and load the kernel. Thread count is fixed here, on first use."""
    # Julia must handle signals itself when it runs multithreaded under Python
    if "juliacall" in sys.modules:
        if os.environ.get("PYTHON_JULIACALL_HANDLE_SIGNALS") != "yes":
            warnings.warn(
                "juliacall was imported without PYTHON_JULIACALL_HANDLE_SIGNALS=yes: "
                "the Julia backend may crash when it runs on several threads.",
                stacklevel=3,
            )
    else:
        os.environ.setdefault("PYTHON_JULIACALL_HANDLE_SIGNALS", "yes")
        if "JULIA_NUM_THREADS" not in os.environ:
            os.environ.setdefault("PYTHON_JULIACALL_THREADS", "auto")
    from juliacall import Main as jl  # importing juliacall boots Julia

    jl.include(str(KERNEL))
    return jl.rings


def ring_statistics_csr(
    offsets: np.ndarray,
    indices: np.ndarray,
    refnodes: np.ndarray,
    *,
    replicas: int,
    maxlvl: int,
    deduplicate: bool,
    mxpths: int,
) -> tuple[np.ndarray, list[np.ndarray]]:
    """Ring counts per original cell and ring node matrices from the Julia kernel."""
    counts, _, rings = _kernel().ring_statistics_csr(
        offsets,
        indices,
        refnodes,
        maxlvl=maxlvl,
        mxpths=mxpths,
        rsf=replicas,
        deduplicate_rings=deduplicate,
    )
    return np.asarray(counts), [np.asarray(r) for r in rings]
