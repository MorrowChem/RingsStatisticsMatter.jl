use numpy::ndarray::Array2;
use numpy::{IntoPyArray, PyArray1, PyArray2, PyReadonlyArray1};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use rayon::prelude::*;

mod kernel;

fn to_nodes(values: &[i64], nodes: usize, what: &str) -> PyResult<Vec<u32>> {
    values
        .iter()
        .map(|&v| {
            if v < 0 || v as usize >= nodes {
                Err(PyValueError::new_err(format!(
                    "{what} {v} is not a node index below {nodes}"
                )))
            } else {
                Ok(v as u32)
            }
        })
        .collect()
}

/// Ring statistics for a CSR neighbour list.
///
/// Returns counts per original cell (index n - 1 for n-membered rings) and one
/// (n_rings, n) array of node indices per ring size.
#[pyfunction]
#[pyo3(signature = (offsets, indices, refnodes, *, replicas, maxlvl, deduplicate, num_threads=None))]
#[allow(clippy::too_many_arguments, clippy::type_complexity)]
fn ring_statistics_csr<'py>(
    py: Python<'py>,
    offsets: PyReadonlyArray1<'py, i64>,
    indices: PyReadonlyArray1<'py, i64>,
    refnodes: PyReadonlyArray1<'py, i64>,
    replicas: usize,
    maxlvl: usize,
    deduplicate: bool,
    num_threads: Option<usize>,
) -> PyResult<(Bound<'py, PyArray1<f64>>, Vec<Bound<'py, PyArray2<i64>>>)> {
    let offsets = offsets.as_slice()?;
    let indices = indices.as_slice()?;
    if offsets.is_empty() || offsets[0] != 0 || *offsets.last().unwrap() as usize != indices.len() {
        return Err(PyValueError::new_err(
            "offsets must start at 0 and end at len(indices)",
        ));
    }
    if offsets.windows(2).any(|w| w[0] > w[1]) {
        return Err(PyValueError::new_err("offsets must not decrease"));
    }
    let nodes = offsets.len() - 1;
    if nodes >= u32::MAX as usize {
        return Err(PyValueError::new_err("too many nodes"));
    }
    if !(2..=250).contains(&maxlvl) {
        return Err(PyValueError::new_err("maxlvl must be between 2 and 250"));
    }
    if replicas == 0 || nodes % replicas != 0 {
        return Err(PyValueError::new_err(
            "the number of nodes must be a multiple of replicas",
        ));
    }
    let offsets: Vec<usize> = offsets.iter().map(|&o| o as usize).collect();
    let indices = to_nodes(indices, nodes, "neighbour")?;
    let refnodes = to_nodes(refnodes.as_slice()?, nodes, "reference node")?;
    let maxlvl = maxlvl as u8;
    let n_original = (nodes / replicas) as u32;

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads.unwrap_or(0))
        .build()
        .map_err(|e| PyValueError::new_err(e.to_string()))?;

    let (counts, rings) = py.allow_threads(|| {
        pool.install(|| {
            let graph = kernel::Graph {
                offsets: &offsets,
                indices: &indices,
            };
            let found = kernel::find_rings(&graph, &refnodes, maxlvl);
            let counts: Vec<f64> = found
                .counts
                .iter()
                .enumerate()
                .map(|(i, &count)| count as f64 / ((i + 1) * replicas) as f64)
                .collect();
            let rings: Vec<Vec<u32>> = found
                .rings
                .into_par_iter()
                .enumerate()
                .map(|(i, rings)| {
                    if deduplicate {
                        kernel::deduplicate(&rings, i + 1, n_original)
                    } else {
                        rings
                    }
                })
                .collect();
            (counts, rings)
        })
    });

    let rings = rings
        .into_iter()
        .enumerate()
        .map(|(i, flat)| {
            let size = i + 1;
            let flat: Vec<i64> = flat.into_iter().map(i64::from).collect();
            Array2::from_shape_vec((flat.len() / size, size), flat)
                .unwrap()
                .into_pyarray(py)
        })
        .collect();
    Ok((counts.into_pyarray(py), rings))
}

#[pymodule]
fn _kernel(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(ring_statistics_csr, module)?)
}
