//! Primitive ring search after Yuan & Cormack, Comput. Mater. Sci. 24, 343 (2002).
//!
//! A port of `rings_statistics_matter/rings.jl`. For every source node, each node closer
//! than maxlvl / 2 is tried as the midpoint ("prime") of a ring through the source. A
//! candidate ring made of two shortest paths is primitive when no pair of its nodes is
//! joined by a path shorter than the ring provides.

use rayon::prelude::*;

/// Marks the midpoint slot of a path, which holds no path node.
const NONE: u32 = u32::MAX;

pub struct Graph<'a> {
    pub offsets: &'a [usize],
    pub indices: &'a [u32],
}

impl Graph<'_> {
    pub fn len(&self) -> usize {
        self.offsets.len() - 1
    }

    fn neighbours(&self, node: u32) -> &[u32] {
        &self.indices[self.offsets[node as usize]..self.offsets[node as usize + 1]]
    }
}

/// Breadth-first levels from one node, up to `maxlvl`. Nodes further away hold `maxlvl + 2`.
struct Levels {
    maxlvl: u8,
    dist: Vec<u8>,
    touched: Vec<u32>,
    queue: Vec<u32>,
}

impl Levels {
    fn new(nodes: usize, maxlvl: u8) -> Self {
        Self {
            maxlvl,
            dist: vec![maxlvl + 2; nodes],
            touched: Vec::new(),
            queue: Vec::new(),
        }
    }

    fn search(&mut self, graph: &Graph, source: u32) {
        for &node in &self.touched {
            self.dist[node as usize] = self.maxlvl + 2;
        }
        self.touched.clear();
        self.queue.clear();

        self.dist[source as usize] = 0;
        self.touched.push(source);
        self.queue.push(source);
        let mut head = 0;
        while head < self.queue.len() {
            let node = self.queue[head];
            head += 1;
            let level = self.dist[node as usize] + 1;
            for &next in graph.neighbours(node) {
                if self.dist[next as usize] > level {
                    self.dist[next as usize] = level;
                    self.touched.push(next);
                    if level < self.maxlvl {
                        self.queue.push(next);
                    }
                }
            }
        }
    }
}

/// Levels from the reference nodes, stored per node for the references within `maxlvl`.
/// Entries of a node are ordered by reference.
pub struct ReferenceMaps {
    offsets: Vec<usize>,
    refs: Vec<u32>,
    dists: Vec<u8>,
}

impl ReferenceMaps {
    pub fn new(graph: &Graph, refnodes: &[u32], maxlvl: u8) -> Self {
        let reached: Vec<Vec<(u32, u8)>> = refnodes
            .par_iter()
            .map_init(
                || Levels::new(graph.len(), maxlvl),
                |levels, &refnode| {
                    levels.search(graph, refnode);
                    levels
                        .touched
                        .iter()
                        .map(|&node| (node, levels.dist[node as usize]))
                        .collect()
                },
            )
            .collect();

        let mut offsets = vec![0usize; graph.len() + 1];
        for &(node, _) in reached.iter().flatten() {
            offsets[node as usize + 1] += 1;
        }
        for node in 0..graph.len() {
            offsets[node + 1] += offsets[node];
        }
        let mut cursor = offsets.clone();
        let mut refs = vec![0u32; offsets[graph.len()]];
        let mut dists = vec![0u8; offsets[graph.len()]];
        for (reference, nodes) in reached.iter().enumerate() {
            for &(node, dist) in nodes {
                let slot = &mut cursor[node as usize];
                refs[*slot] = reference as u32;
                dists[*slot] = dist;
                *slot += 1;
            }
        }
        Self {
            offsets,
            refs,
            dists,
        }
    }

    fn of(&self, node: u32) -> (&[u32], &[u8]) {
        let range = self.offsets[node as usize]..self.offsets[node as usize + 1];
        (&self.refs[range.clone()], &self.dists[range])
    }
}

/// Shortest paths from one prime back to the source. Each path holds the nodes at
/// levels 1..lvlprim, nearest the source first.
#[derive(Default)]
struct Paths {
    stride: usize,
    nodes: Vec<u32>,
    count: usize,
}

impl Paths {
    fn path(&self, index: usize) -> &[u32] {
        &self.nodes[index * self.stride..(index + 1) * self.stride]
    }

    /// Node of a path at a 1-based level; `NONE` at the level of the prime.
    fn at(&self, index: usize, level: usize) -> u32 {
        if level <= self.stride {
            self.nodes[index * self.stride + level - 1]
        } else {
            NONE
        }
    }

    fn shares_node(&self, index: usize, other: &Paths, other_index: usize) -> bool {
        // nodes at different levels are distinct, so compare level by level
        self.path(index)
            .iter()
            .zip(other.path(other_index))
            .any(|(a, b)| a == b)
    }
}

struct Candidate {
    path1: usize,
    path2: usize,
    alive: bool,
}

#[derive(Clone)]
pub struct Found {
    /// counts[n - 1]: n-membered rings, each counted once per source that found it.
    pub counts: Vec<u64>,
    /// rings[n - 1]: node lists of n-membered rings, flattened.
    pub rings: Vec<Vec<u32>>,
}

impl Found {
    fn new(maxlvl: u8) -> Self {
        let sizes = 2 * maxlvl as usize;
        Self {
            counts: vec![0; sizes],
            rings: vec![Vec::new(); sizes],
        }
    }

    fn merge(mut self, other: Found) -> Found {
        for (a, b) in self.counts.iter_mut().zip(other.counts) {
            *a += b;
        }
        for (a, b) in self.rings.iter_mut().zip(other.rings) {
            a.extend(b);
        }
        self
    }
}

struct Search<'a> {
    graph: &'a Graph<'a>,
    maps: &'a ReferenceMaps,
    maxlvl: u8,
    levels: Levels,
    primes: Vec<u32>,
    scratch: Vec<u32>,
    paths1: Paths,
    paths2: Paths,
    candidates: Vec<Candidate>,
    found: Found,
}

impl<'a> Search<'a> {
    fn new(graph: &'a Graph<'a>, maps: &'a ReferenceMaps, maxlvl: u8) -> Self {
        Self {
            graph,
            maps,
            maxlvl,
            levels: Levels::new(graph.len(), maxlvl),
            primes: Vec::new(),
            scratch: Vec::new(),
            paths1: Paths::default(),
            paths2: Paths::default(),
            candidates: Vec::new(),
            found: Found::new(maxlvl),
        }
    }

    fn level(&self, node: u32) -> u8 {
        self.levels.dist[node as usize]
    }

    fn record_paths(&mut self, prime: u32, second: bool) {
        let stride = self.level(prime) as usize - 1;
        let mut paths = std::mem::take(if second {
            &mut self.paths2
        } else {
            &mut self.paths1
        });
        paths.stride = stride;
        paths.nodes.clear();
        paths.count = 0;
        self.scratch.clear();
        self.scratch.resize(stride, NONE);
        let mut scratch = std::mem::take(&mut self.scratch);
        self.descend(prime, &mut scratch, &mut paths);
        self.scratch = scratch;
        *(if second {
            &mut self.paths2
        } else {
            &mut self.paths1
        }) = paths;
    }

    fn descend(&self, node: u32, current: &mut [u32], paths: &mut Paths) {
        let level = self.level(node);
        for &next in self.graph.neighbours(node) {
            let next_level = self.level(next);
            if next_level == 0 {
                paths.nodes.extend_from_slice(current);
                paths.count += 1;
            } else if next_level + 1 == level {
                current[next_level as usize - 1] = next;
                self.descend(next, current, paths);
            }
        }
    }

    /// Is `goal` reachable from `node` in fewer steps than the candidate ring provides?
    /// `depth` steps have been taken; the level windows around `goal` shrink with each.
    fn shortcut(&self, node: u32, goal: u32, reach: i32, depth: i32) -> bool {
        if node == goal {
            return true;
        }
        let maxlvl = self.maxlvl;
        let (goal_refs, goal_dists) = self.maps.of(goal);
        let goal_level = self.level(goal) as i32;
        let outside = |value: i32, centre: i32| {
            value >= centre + reach - depth || value <= centre - reach + depth
        };

        'neighbours: for &next in self.graph.neighbours(node) {
            let (refs, dists) = self.maps.of(next);
            let (mut i, mut j) = (0, 0);
            while i < refs.len() && j < goal_refs.len() {
                if refs[i] < goal_refs[j] {
                    i += 1;
                } else if refs[i] > goal_refs[j] {
                    j += 1;
                } else {
                    if dists[i] < maxlvl && outside(dists[i] as i32, goal_dists[j] as i32) {
                        continue 'neighbours;
                    }
                    i += 1;
                    j += 1;
                }
            }
            if outside(self.level(next) as i32, goal_level) {
                continue;
            }
            if self.shortcut(next, goal, reach, depth + 1) {
                return true;
            }
        }
        false
    }

    /// Count the candidates that are primitive rings. `odd` rings have two primes.
    fn primitive_rings(&mut self, odd: bool, source: u32, prime1: u32, prime2: u32) {
        let lvlprim = self.level(prime1) as usize;
        let mut candidates = std::mem::take(&mut self.candidates);
        let paths1 = std::mem::take(&mut self.paths1);
        let paths2 = if odd {
            std::mem::take(&mut self.paths2)
        } else {
            Paths::default()
        };
        let second = if odd { &paths2 } else { &paths1 };

        'candidates: for index in 0..candidates.len() {
            if !candidates[index].alive {
                continue;
            }
            let (path1, path2) = (candidates[index].path1, candidates[index].path2);
            for lvlmax in lvlprim..=lvlprim + odd as usize {
                for lvlchk in 1..lvlmax {
                    let lvlmid = lvlmax - lvlchk;
                    let mut check = paths1.at(path1, lvlchk);
                    let mut mid = second.at(path2, lvlmid);
                    if check == NONE {
                        check = prime1;
                    }
                    if mid == NONE {
                        mid = prime2;
                    }
                    if !self.shortcut(check, mid, lvlprim as i32 - 1, 0) {
                        continue;
                    }
                    // every later candidate through the same two nodes fails the same way
                    for later in &mut candidates[index + 1..] {
                        if !later.alive {
                            continue;
                        }
                        let forward = paths1.at(later.path1, lvlchk) == check
                            && second.at(later.path2, lvlmid) == mid;
                        let backward = !odd
                            && paths1.at(later.path2, lvlchk) == check
                            && paths1.at(later.path1, lvlmid) == mid;
                        if forward || backward {
                            later.alive = false;
                        }
                    }
                    continue 'candidates;
                }
            }

            let size = 2 * lvlprim + odd as usize;
            self.found.counts[size - 1] += 1;
            let ring = &mut self.found.rings[size - 1];
            ring.push(source);
            ring.extend_from_slice(paths1.path(path1));
            ring.push(prime1);
            if odd {
                ring.push(prime2);
            }
            ring.extend(second.path(path2).iter().rev());
        }

        candidates.clear();
        self.candidates = candidates;
        self.paths1 = paths1;
        if odd {
            self.paths2 = paths2;
        }
    }

    fn rings_through(&mut self, source: u32) {
        let graph = self.graph;
        self.levels.search(graph, source);
        let half = self.maxlvl / 2;
        let mut primes = std::mem::take(&mut self.primes);
        primes.clear();
        primes.extend(
            self.levels
                .touched
                .iter()
                .filter(|&&node| (1..half).contains(&self.level(node))),
        );
        primes.sort_unstable();

        for &prime in &primes {
            let level = self.level(prime);

            // even ring: two shortest paths meet at the prime
            if graph
                .neighbours(prime)
                .iter()
                .filter(|&&n| self.level(n) + 1 == level)
                .count()
                > 1
            {
                self.record_paths(prime, false);
                for path1 in 0..self.paths1.count {
                    for path2 in path1..self.paths1.count {
                        let disjoint = if path1 == path2 {
                            self.paths1.stride == 0
                        } else {
                            !self.paths1.shares_node(path1, &self.paths1, path2)
                        };
                        if disjoint {
                            self.candidates.push(Candidate {
                                path1,
                                path2,
                                alive: true,
                            });
                        }
                    }
                }
                self.primitive_rings(false, source, prime, prime);
            }

            // odd ring: shortest paths end on two bonded primes of equal level
            let mut recorded = false;
            for &prime2 in graph.neighbours(prime) {
                if self.level(prime2) != level {
                    continue;
                }
                if !recorded {
                    self.record_paths(prime, false);
                    recorded = true;
                }
                if prime2 <= prime {
                    continue;
                }
                self.record_paths(prime2, true);
                for path1 in 0..self.paths1.count {
                    for path2 in 0..self.paths2.count {
                        if !self.paths1.shares_node(path1, &self.paths2, path2) {
                            self.candidates.push(Candidate {
                                path1,
                                path2,
                                alive: true,
                            });
                        }
                    }
                }
                self.primitive_rings(true, source, prime, prime2);
            }
        }
        self.primes = primes;
    }
}

/// Rings through every node. Each n-membered ring is found once from each of its nodes.
pub fn find_rings(graph: &Graph, refnodes: &[u32], maxlvl: u8) -> Found {
    let maps = ReferenceMaps::new(graph, refnodes, maxlvl);
    let sources: Vec<u32> = (0..graph.len() as u32).collect();
    let chunk = (sources.len() / (8 * rayon::current_num_threads())).max(1);
    sources
        .par_chunks(chunk)
        .map(|chunk| {
            let mut search = Search::new(graph, &maps, maxlvl);
            for &source in chunk {
                search.rings_through(source);
            }
            search.found
        })
        .reduce(|| Found::new(maxlvl), Found::merge)
}

/// Smallest rotation of a ring or its reverse, in lexicographic order.
fn canonical(ring: &[u32]) -> Vec<u32> {
    let n = ring.len();
    let reversed: Vec<u32> = ring.iter().rev().copied().collect();
    let mut best = ring.to_vec();
    let mut candidate = Vec::with_capacity(n);
    for direction in [ring, &reversed] {
        for shift in 0..n {
            candidate.clear();
            candidate.extend_from_slice(&direction[shift..]);
            candidate.extend_from_slice(&direction[..shift]);
            if candidate < best {
                best.clone_from(&candidate);
            }
        }
    }
    best
}

/// Map nodes into the first `n_original` and keep one ring of each set equal up to
/// rotation and reversal. `rings` is flattened, `size` nodes per ring.
pub fn deduplicate(rings: &[u32], size: usize, n_original: u32) -> Vec<u32> {
    let mut unique: Vec<Vec<u32>> = rings
        .par_chunks(size)
        .map(|ring| {
            canonical(
                &ring
                    .iter()
                    .map(|node| node % n_original)
                    .collect::<Vec<_>>(),
            )
        })
        .collect();
    unique.par_sort_unstable();
    unique.dedup();
    unique.concat()
}
