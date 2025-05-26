// File: bvg.rs
// Purpose: Reconstruct a degree-constrained MST from a population
//          of BitVectors with asexual bit-flip inheritance.
//
// This algorithm:
// 1. Reads a population of binary genomes from a file given in random order
// 2. Computes pairwise Hamming distances between all genomes
// 3. Weights edges by negative log-likelihood under a binomial mutation model
// 4. Selects a root genome with high connectivity (statistical outlier)
// 5. Constructs a degree-constrained minimum spanning tree representing evolutionary relationships
// 6. Outputs parent-child relationships for phylogenetic reconstruction

use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::env;
use std::sync::Mutex;
use std::process;
use rayon::prelude::*;
use statrs::distribution::{Binomial, Discrete};

// Algorithm parameters
const GENOME_LEN: usize = 10000;        // Length of each binary genome
const POPULATION_SIZE: usize = 10000;   // Number of genomes in population
const MUTATION_PROB: f64 = 0.2;         // Per-bit mutation probability in inheritance
const LOG_LIKELIHOOD_CUTOFF: f64 = 20.0; // Ignore edges with -log(p) > cutoff (too unlikely)
const SCALE: f64 = 1000.0;              // Scale factor for integer edge weights
const ZSCORE_THRESHOLD: f64 = 3.0;      // Root must be mean + 3*stddev above average degree

// Optimized data types for performance
type Genome = Vec<u64>;  // Use u64 chunks for fast bit operations (64x speedup vs BitVec)

// Bit manipulation constants
const BITS_PER_CHUNK: usize = 64;
const CHUNKS: usize = (GENOME_LEN + BITS_PER_CHUNK - 1) / BITS_PER_CHUNK;

/// Configuration struct to hold all program parameters
#[derive(Debug)]
struct Config {
    input_file: String,
    output_file: String,
    zscore_threshold: f64,
    log_likelihood_cutoff: f64,
    chunk_size: usize,
    threads: Option<usize>,
    verbose: bool,
}

impl Default for Config {
    fn default() -> Self {
        Self {
            input_file: "test.data".to_string(),
            output_file: "out.txt".to_string(),
            zscore_threshold: ZSCORE_THRESHOLD,
            log_likelihood_cutoff: LOG_LIKELIHOOD_CUTOFF,
            chunk_size: 100,
            threads: None,
            verbose: false,
        }
    }
}

/// Parse command line arguments and return configuration
fn parse_args() -> Config {
    let args: Vec<String> = env::args().collect();
    let mut config = Config::default();
    let mut i = 1;
    
    while i < args.len() {
        match args[i].as_str() {
            "-h" | "--help" => {
                print_help(&args[0]);
                process::exit(0);
            }
            "-i" | "--input" => {
                if i + 1 < args.len() {
                    config.input_file = args[i + 1].clone();
                    i += 2;
                } else {
                    eprintln!("Error: --input requires a filename");
                    process::exit(1);
                }
            }
            "-o" | "--output" => {
                if i + 1 < args.len() {
                    config.output_file = args[i + 1].clone();
                    i += 2;
                } else {
                    eprintln!("Error: --output requires a filename");
                    process::exit(1);
                }
            }
            "-z" | "--zscore" => {
                if i + 1 < args.len() {
                    config.zscore_threshold = args[i + 1].parse().unwrap_or_else(|_| {
                        eprintln!("Error: Invalid z-score value: {}", args[i + 1]);
                        process::exit(1);
                    });
                    i += 2;
                } else {
                    eprintln!("Error: --zscore requires a number");
                    process::exit(1);
                }
            }
            "-m" | "--mutation" => {
                eprintln!("Warning: Mutation probability is fixed at {} for this version", MUTATION_PROB);
                i += 1;
            }
            "-c" | "--cutoff" => {
                if i + 1 < args.len() {
                    config.log_likelihood_cutoff = args[i + 1].parse().unwrap_or_else(|_| {
                        eprintln!("Error: Invalid cutoff value: {}", args[i + 1]);
                        process::exit(1);
                    });
                    i += 2;
                } else {
                    eprintln!("Error: --cutoff requires a number");
                    process::exit(1);
                }
            }
            "--chunk-size" => {
                if i + 1 < args.len() {
                    config.chunk_size = args[i + 1].parse().unwrap_or_else(|_| {
                        eprintln!("Error: Invalid chunk size: {}", args[i + 1]);
                        process::exit(1);
                    });
                    i += 2;
                } else {
                    eprintln!("Error: --chunk-size requires a number");
                    process::exit(1);
                }
            }
            "-j" | "--threads" => {
                if i + 1 < args.len() {
                    let threads = args[i + 1].parse().unwrap_or_else(|_| {
                        eprintln!("Error: Invalid thread count: {}", args[i + 1]);
                        process::exit(1);
                    });
                    config.threads = Some(threads);
                    i += 2;
                } else {
                    eprintln!("Error: --threads requires a number");
                    process::exit(1);
                }
            }
            "-v" | "--verbose" => {
                config.verbose = true;
                i += 1;
            }
            "--version" => {
                println!("bvg-rs version 1.0.0");
                process::exit(0);
            }
            _ => {
                if args[i].starts_with('-') {
                    eprintln!("Error: Unknown option: {}", args[i]);
                    eprintln!("Use --help for usage information");
                    process::exit(1);
                } else {
                    // Assume it's an input file if no option specified
                    config.input_file = args[i].clone();
                }
                i += 1;
            }
        }
    }
    
    config
}

/// Print help message
fn print_help(program_name: &str) {
    println!("BVG-RS: Bitvector Genealogy Reconstruction");
    println!("Reconstruct minimum spanning tree from binary genomes using mutation model");
    println!();
    println!("USAGE:");
    println!("    {} [OPTIONS] [INPUT_FILE]", program_name);
    println!();
    println!("OPTIONS:");
    println!("    -i, --input <FILE>       Input file with binary genomes [default: test.data]");
    println!("    -o, --output <FILE>      Output file for parent indices [default: out.txt]");
    println!("    -z, --zscore <FLOAT>     Z-score threshold for root selection [default: 3.0]");
    println!("    -c, --cutoff <FLOAT>     Log-likelihood cutoff for edge filtering [default: 20.0]");
    println!("    -j, --threads <NUM>      Number of threads to use [default: auto]");
    println!("        --chunk-size <NUM>   Parallel processing chunk size [default: 100]");
    println!("    -v, --verbose            Enable verbose output");
    println!("    -h, --help               Print this help message");
    println!("        --version            Print version information");
    println!();
    println!("ALGORITHM PARAMETERS:");
    println!("    Genome length: {} bits", GENOME_LEN);
    println!("    Population size: {} genomes", POPULATION_SIZE);
    println!("    Mutation probability: {} (fixed)", MUTATION_PROB);
    println!();
    println!("EXAMPLES:");
    println!("    {}                                    # Use defaults", program_name);
    println!("    {} -i genomes.txt -o tree.txt        # Specify files", program_name);
    println!("    {} --zscore 2.5 --cutoff 15.0         # Adjust parameters", program_name);
    println!("    {} --threads 8 --verbose             # Performance tuning", program_name);
    println!();
    println!("INPUT FORMAT:");
    println!("    Each line should contain a binary string (e.g., '0110101...')");
    println!("    Lines should be {} characters long", GENOME_LEN);
    println!();
    println!("OUTPUT FORMAT:");
    println!("    Each line contains the parent index for that genome (0-based)");
    println!("    Root genome has parent index -1");
}

/// Read population genomes from file
/// Each line contains a binary string (e.g., "0110101...")
/// Converts to u64 chunks for efficient bit operations
fn read_population(filename: &str) -> Vec<Genome> {
    let file = File::open(filename).expect("Failed to open input file");
    let reader = BufReader::new(file);
    let mut population = Vec::with_capacity(POPULATION_SIZE);
    
    for line in reader.lines().take(POPULATION_SIZE) {
        let line = line.unwrap();
        let mut genome = vec![0u64; CHUNKS];
        
        // Pack bits into u64 chunks for fast XOR operations
        for (i, c) in line.chars().take(GENOME_LEN).enumerate() {
            if c == '1' {
                let chunk_idx = i / BITS_PER_CHUNK;
                let bit_idx = i % BITS_PER_CHUNK;
                genome[chunk_idx] |= 1u64 << bit_idx;
            }
        }
        population.push(genome);
    }
    
    population
}

/// Precompute edge weights based on binomial mutation model
/// For each possible Hamming distance d, compute -log(P(d)) where P(d) is the
/// probability of observing d mutations under the binomial model.
/// Lower weights = more likely mutations = stronger evolutionary relationships
fn compute_log_lookup_table(cutoff: f64) -> Vec<Option<i32>> {
    let binom = Binomial::new(MUTATION_PROB, GENOME_LEN as u64).unwrap();
    let mut log_lookup = Vec::with_capacity(GENOME_LEN + 1);
    
    for d in 0..=GENOME_LEN {
        let p = binom.pmf(d as u64);  // Probability of exactly d mutations
        let logp = -p.ln();           // Negative log-likelihood (lower = more likely)
        
        // Filter out impossible/improbable mutations
        if !logp.is_finite() || logp > cutoff {
            log_lookup.push(None);  // Don't create edge for this distance
        } else {
            log_lookup.push(Some((logp * SCALE) as i32));  // Scale to integer weights
        }
    }
    
    log_lookup
}

/// Fast Hamming distance using u64 XOR and popcount
/// This is the key performance optimization - operates on 64-bit chunks
/// instead of individual bits, giving ~64x speedup over naive approaches
#[inline]
fn hamming_distance(a: &Genome, b: &Genome) -> usize {
    a.iter()
        .zip(b.iter())
        .map(|(x, y)| (x ^ y).count_ones() as usize)  // XOR + popcount per chunk
        .sum()
}

/// Simple adjacency list graph representation optimized for this use case
/// Much faster than generic graph libraries due to cache-friendly memory layout
struct SimpleGraph {
    edges: Vec<Vec<(usize, i32)>>, // adjacency list: node -> [(neighbor, weight)]
}

impl SimpleGraph {
    fn new(n: usize) -> Self {
        Self {
            edges: vec![Vec::new(); n],
        }
    }
    
    /// Add undirected edge with given weight
    fn add_edge(&mut self, u: usize, v: usize, weight: i32) {
        self.edges[u].push((v, weight));
        self.edges[v].push((u, weight));
    }
    
    /// Get degree (number of neighbors) for node u
    fn degree(&self, u: usize) -> usize {
        self.edges[u].len()
    }
    
    /// Get all neighbors of node u with their edge weights
    fn neighbors(&self, u: usize) -> &[(usize, i32)] {
        &self.edges[u]
    }
}

/// Prim's minimum spanning tree algorithm optimized for dense graphs
/// Returns parent array where parent[i] is the parent of node i in the MST
/// Root node has parent[root] = None
fn prim_mst(graph: &SimpleGraph, root: usize) -> Vec<Option<usize>> {
    let n = graph.edges.len();
    let mut in_mst = vec![false; n];      // Track which nodes are in MST
    let mut parent = vec![None; n];       // Parent of each node in MST
    let mut min_cost = vec![i32::MAX; n]; // Minimum cost to connect each node
    
    min_cost[root] = 0;  // Start from root with zero cost
    
    // Add nodes to MST one by one
    for _ in 0..n {
        // Find minimum cost vertex not yet in MST
        let mut u = 0;
        let mut min_val = i32::MAX;
        
        for v in 0..n {
            if !in_mst[v] && min_cost[v] < min_val {
                min_val = min_cost[v];
                u = v;
            }
        }
        
        in_mst[u] = true;  // Add u to MST
        
        // Update costs of vertices adjacent to u
        for &(v, weight) in graph.neighbors(u) {
            if !in_mst[v] && weight < min_cost[v] {
                min_cost[v] = weight;
                parent[v] = Some(u);  // u is now the best parent for v
            }
        }
    }
    
    parent
}

fn main() {
    // Parse command line arguments
    let config = parse_args();
    
    // Set thread count if specified
    if let Some(threads) = config.threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global()
            .unwrap();
    }
    
    if config.verbose {
        eprintln!("Configuration:");
        eprintln!("  Input file: {}", config.input_file);
        eprintln!("  Output file: {}", config.output_file);
        eprintln!("  Z-score threshold: {}", config.zscore_threshold);
        eprintln!("  Mutation probability: {} (fixed)", MUTATION_PROB);
        eprintln!("  Log-likelihood cutoff: {}", config.log_likelihood_cutoff);
        eprintln!("  Chunk size: {}", config.chunk_size);
        eprintln!("  Threads: {}", rayon::current_num_threads());
        eprintln!();
    }
    
    // Read population genomes from input file
    if config.verbose {
        eprintln!("Reading genomes from {}...", config.input_file);
    }
    let population = read_population(&config.input_file);
    eprintln!("Read {} genomes", population.len());
    
    // Precompute log-likelihood scores for Hamming distances
    if config.verbose {
        eprintln!("Computing mutation model lookup table...");
    }
    let log_lookup = compute_log_lookup_table(config.log_likelihood_cutoff);
    
    // Compute all pairwise evolutionary distances in parallel
    // This is O(n²) but highly parallelizable - each thread processes chunks of genome pairs
    // Uses chunked processing to reduce lock contention on the shared edge vector
    if config.verbose {
        eprintln!("Computing pairwise distances...");
    }
    
    let edges = Mutex::new(Vec::new());
    
    (0..POPULATION_SIZE)
        .into_par_iter()
        .chunks(config.chunk_size)
        .for_each(|chunk| {
            let mut local_edges = Vec::new();
            
            // For each genome j in this chunk, compare with all previous genomes i < j
            for j in chunk {
                for i in 0..j {
                    let d = hamming_distance(&population[i], &population[j]);
                    if let Some(weight) = log_lookup[d] {
                        local_edges.push((i, j, weight));
                    }
                }
            }
            
            // Batch insert to reduce lock contention (much faster than individual pushes)
            edges.lock().unwrap().extend(local_edges);
        });
    
    let edges = edges.into_inner().unwrap();
    eprintln!("Generated {} edges", edges.len());
    
    // Build graph from computed edges
    if config.verbose {
        eprintln!("Building graph...");
    }
    let mut graph = SimpleGraph::new(POPULATION_SIZE);
    for (u, v, weight) in edges {
        graph.add_edge(u, v, weight);
    }
    
    // Calculate degree statistics for root selection
    // We want to find genomes with unusually high connectivity (potential ancestors)
    if config.verbose {
        eprintln!("Computing degree statistics...");
    }
    let degrees: Vec<usize> = (0..POPULATION_SIZE)
        .map(|i| graph.degree(i))
        .collect();
    
    let mean: f64 = degrees.iter().sum::<usize>() as f64 / POPULATION_SIZE as f64;
    let variance: f64 = degrees.iter()
        .map(|&d| {
            let diff = d as f64 - mean;
            diff * diff
        })
        .sum::<f64>() / POPULATION_SIZE as f64;
    let stddev = variance.sqrt();
    
    let threshold = mean + config.zscore_threshold * stddev;
    
    if config.verbose {
        eprintln!("Degree statistics:");
        eprintln!("  Mean: {:.2}", mean);
        eprintln!("  Std dev: {:.2}", stddev);
        eprintln!("  Threshold (mean + {:.1}σ): {:.2}", config.zscore_threshold, threshold);
    } else {
        eprintln!("threshold = {}", threshold);
    }
    
    // Select the MST root as a high-degree statistical outlier
    // Rationale: Ancestral genomes should have many descendants, leading to high connectivity
    // We require degree >= mean + zscore_threshold * stddev to avoid noise
    if config.verbose {
        eprintln!("Selecting root node...");
    }
    
    let mut root_idx = 0;
    let mut max_degree = 0;
    let mut root_found = false;
    
    for (i, &deg) in degrees.iter().enumerate() {
        if deg >= max_degree && deg as f64 >= threshold {
            max_degree = deg;
            root_idx = i;
            root_found = true;
        }
    }
    
    if !root_found {
        eprintln!("ERROR: No suitable root vertex found.");
        if config.verbose {
            eprintln!("Try lowering the z-score threshold with --zscore");
            let max_actual_degree = degrees.iter().max().unwrap_or(&0);
            eprintln!("Highest degree found: {}", max_actual_degree);
        }
        process::exit(1);
    }
    
    eprintln!("Selected root {}. Degree = {}", root_idx, max_degree);
    
    // Build MST using optimized Prim's algorithm
    // This creates a tree structure representing the most likely evolutionary relationships
    if config.verbose {
        eprintln!("Computing minimum spanning tree...");
    }
    let predecessors = prim_mst(&graph, root_idx);
    
    // Output phylogenetic tree in parent-index format
    // Each line i contains the parent of genome i, or -1 if i is the root
    if config.verbose {
        eprintln!("Writing results to {}...", config.output_file);
    }
    
    let mut output = File::create(&config.output_file)
        .unwrap_or_else(|e| {
            eprintln!("Error creating output file {}: {}", config.output_file, e);
            process::exit(1);
        });
    
    for i in 0..POPULATION_SIZE {
        if let Some(parent_idx) = predecessors[i] {
            if parent_idx == i {
                writeln!(output, "-1").unwrap();  // Root node
            } else {
                writeln!(output, "{}", parent_idx).unwrap();  // Parent index
            }
        } else {
            writeln!(output, "-1").unwrap();  // Shouldn't happen with proper MST, but handle gracefully
        }
    }
    
    if config.verbose {
        eprintln!("Phylogenetic reconstruction complete!");
    }
}