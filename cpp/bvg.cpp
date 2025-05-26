// File: bvg.cpp
// Purpose: Reconstruct a degree-constrained MST from a population
//          of BitVectors with asexual bit-flip inheritance.

#include <algorithm>
#include <fstream>
#include <iostream>
#include <tuple>
#include <vector>
#include <string>
#include <cmath>

// Boost
#include <boost/dynamic_bitset.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/prim_minimum_spanning_tree.hpp>
#include <boost/math/distributions/binomial.hpp>

// TBB
#include <tbb/global_control.h>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/concurrent_queue.h>

using namespace std;
using namespace boost;
using namespace tbb;

// Constants and Definitions
namespace
{
constexpr int GENOME_LEN = 10000;               // Bitvector genome length
constexpr int POPULATION_SIZE = 10000;          // Population size
constexpr double MUTATION_PROB = 0.2;           // Per-bit flip probability
constexpr double LOG_LIKELIHOOD_CUTOFF = 20.0;  // Ignore improbable mutations
constexpr double SCALE = 1000.0;                // Edge weight scaling
constexpr double ZSCORE_THRESHOLD = 3.0;        // root must be at least mean + 3*stddev


using Genome = dynamic_bitset<>;
using Edge = std::tuple<int, int, int>;         // (u, v, weight)
}

int main(int argc, char* argv[])
{
  double zscore_threshold = ZSCORE_THRESHOLD;

  // Parse optional CLI flags:
  // --zscore Z: Use z-score thresholding with Z = ZSCORE_THRESHOLD
  for (int i = 1; i < argc; ++i) {
    string arg = argv[i];
    if (arg == "--zscore" && i + 1 < argc) {
      zscore_threshold = atof(argv[++i]);
    }
  }


  concurrent_queue<Edge> edge_set;

  constexpr char filename[] = "test.data";
  // constexpr size_t EXPECTED_MAX_DEGREE = 10;

  using Graph = adjacency_list<vecS, vecS, undirectedS,
        property<vertex_distance_t, int>,
        property<edge_weight_t, int>>;
  using Vertex = graph_traits<Graph>::vertex_descriptor;

  Graph graph(POPULATION_SIZE);
  auto weight_map = get(edge_weight, graph);
  vector<Vertex> predecessors(POPULATION_SIZE);  // MST output

  // Read population genomes from test.data
  vector<Genome> population(POPULATION_SIZE, Genome(GENOME_LEN));
  ifstream ifs(filename);
  string line;
  size_t index = 0;
  while (getline(ifs, line) && index < POPULATION_SIZE) {
    for (size_t i = 0; i < GENOME_LEN; ++i) {
      if (line[i] == '1') population[index].set(i);
    }
    ++index;
  }

  // Precompute log-likelihood scores for Hamming distances
  //
  // For each possible Hamming distance d between parent and child,
  // compute -log(P(d)) under the binomial model where each bit has
  // an independent flip probability MUTATION_PROB.
  //
  // This yields a similarity score: lower = more likely = better edge weight.
  // Values above a cutoff are treated as noise and excluded.
  vector<int> log_lookup(GENOME_LEN + 1, -1); // logarithmic lookup table
  boost::math::binomial_distribution<> binom(GENOME_LEN, MUTATION_PROB);
  for (int d = 0; d <= GENOME_LEN; ++d) {
    double logp = -log(boost::math::pdf(binom, d));
    log_lookup[d] = (!isfinite(logp) || logp > LOG_LIKELIHOOD_CUTOFF) ? -1
                    : static_cast<int>(logp * SCALE);
  }

  // edge-weight calculator
  auto calculate_edge_weights = [&](const blocked_range<int>& block) {
    for (int j = block.begin(); j < block.end(); ++j) {
      // Full pairwise comparison: all i < j
      for (int i = 0; i < j; ++i) {
        int d = (population[j] ^ population[i]).count();
        int weight = log_lookup[d];
        if (weight != -1) edge_set.push({i, j, weight});
      }
    }
  };

  // Parallel block to compute candidate edges between all pairs (i, j)
  // Adds edge (i, j) if the mutation likelihood (in -log space) is plausible.
  // This is a full pairwise comparison, so O(n^2) complexity, but in practice
  // no slower than KNN or other heuristics.
  parallel_for(blocked_range<int>(0, POPULATION_SIZE),
               calculate_edge_weights,
               auto_partitioner());

  // Construct Boost Graph from edge set
  for (Edge edge; edge_set.try_pop(edge); ) {
    auto [u, v, w] = edge;
    auto [e, inserted] = add_edge(u, v, graph);
    if (inserted) weight_map[e] = w;
  }


  vector<size_t> degrees(POPULATION_SIZE);
  double sum = 0.0;
  for (size_t i = 0; i < POPULATION_SIZE; ++i) {
    degrees[i] = degree(i, graph);
    sum += degrees[i];
  }
  double mean = sum / POPULATION_SIZE;

  double sq_sum = 0.0;
  for (size_t d : degrees) sq_sum += (d - mean) * (d - mean);
  double stddev = sqrt(sq_sum / POPULATION_SIZE);

  double threshold = mean + zscore_threshold * stddev;
  cerr << "threshold = " << threshold << '\n';
  // auto degree_threshold = EXPECTED_MAX_DEGREE + 3;

  // Select the MST root as the vertex with highest degree
  // among those that are at least zscore threshold stddevs above the mean.
  // This favors informative statistical outliers (e.g., the true ancestor).
  // If no such vertex exists, select the one with highest degree.
  auto [vi_start, vi_end] = vertices(graph);
  auto root = vi_start;
  size_t max_degree = 0;
  bool root_found = false;
  for (auto vi = vi_start; vi != vi_end; ++vi) {
    size_t deg = degree(*vi, graph);
    if (deg >= max_degree && deg >= threshold) {
      max_degree = deg;
      root = vi;
      root_found = true;
    }
  }

  if (!root_found) {
    cerr << "ERROR: No suitable root vertex found.\n";
    return 1;
  }

  cerr << "Selected root. Degree = " << max_degree << "\n";

  // Run Prim's MST algorithm
  auto distance_map = get(vertex_distance, graph);
  auto index_map = get(vertex_index, graph);

  prim_minimum_spanning_tree(graph, *root,
                             make_iterator_property_map(predecessors.begin(), index_map),
                             distance_map, weight_map, index_map,
                             default_dijkstra_visitor());


  // Output result: Parent index per genome (one per line)
  // Genesis gets parent = -1
  for (size_t i = 0; i < POPULATION_SIZE; ++i) {
    if (predecessors[i] != i) {
      cout << predecessors[i] << '\n';
    } else {
      cout << "-1\n";
    }
  }
}
