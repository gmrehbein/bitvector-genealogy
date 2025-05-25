// ---------------------------------------------------------
// File: bvg.cc
// Purpose: Reconstruct a degree-constrained MST from a population
//          of BitVectors with asexual bit-flip inheritance.
// ---------------------------------------------------------

#include <fstream>
#include <iostream>
#include <tuple>
#include <vector>
#include <string>
#include <cmath>
#include <cstdlib>
#include <algorithm>

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

// -------------------------------
// Constants and Definitions
// -------------------------------
namespace
{
constexpr int GENOME_LEN = 10000;               // Bitvector genome length
constexpr int POPULATION_SIZE = 10000;          // Population size
constexpr double MUTATION_PROB = 0.2;           // Per-bit flip probability
constexpr double LOG_LIKELIHOOD_CUTOFF = 20.0;  // Ignore improbable mutations
constexpr double SCALE = 1000.0;                // Edge weight scaling

using Genome = dynamic_bitset<>;
using Edge = std::tuple<int, int, int>;              // (u, v, weight)
}

int main(int argc, char* argv[])
{
  bool use_knn = false;
  int knn_k = 10;

  // Parse optional CLI flag: --knn N
  for (int i = 1; i < argc; ++i) {
    string arg = argv[i];
    if (arg == "--knn" && i + 1 < argc) {
      knn_k = atoi(argv[++i]);
      use_knn = true;
    }
  }

  concurrent_queue<Edge> edge_set;
  // Use max available parallelism
  global_control control(global_control::max_allowed_parallelism, std::thread::hardware_concurrency());

  constexpr char filename[] = "test.data";
  constexpr size_t EXPECTED_MAX_DEGREE = 10;

  using Graph = adjacency_list<vecS, vecS, undirectedS,
        property<vertex_distance_t, int>,
        property<edge_weight_t, int>>;
  using Vertex = graph_traits<Graph>::vertex_descriptor;

  Graph graph(POPULATION_SIZE);
  auto weight_map = get(edge_weight, graph);
  vector<Vertex> predecessors(POPULATION_SIZE);  // MST output

  // ------------------------------------------
  // Read population genomes from test.data
  // ------------------------------------------

  vector<Genome> population(POPULATION_SIZE, Genome(GENOME_LEN));
  ifstream ifs(filename);
  string line;
  size_t index = 0;
  while (getline(ifs, line) && index < POPULATION_SIZE) {
    for (size_t i = 0; i < GENOME_LEN; ++i)
      if (line[i] == '1') population[index].set(i);
    ++index;
  }

  cerr << "Read " << index << " genomes.\n";
  cerr << (use_knn ? "Using KNN sparsification (K = " + to_string(knn_k) + ')'
          : "Using full pairwise edge construction") << '\n';

  // ------------------------------------------
  // Precompute log-likelihoods for all d = 0..10000
  // ------------------------------------------

  vector<int> log_lut(GENOME_LEN + 1, -1); // logarithmic lookup table
  boost::math::binomial_distribution<> binom(GENOME_LEN, MUTATION_PROB);
  for (int d = 0; d <= GENOME_LEN; ++d) {
    double logp = -log(boost::math::pdf(binom, d));
    log_lut[d] = (!isfinite(logp) || logp > LOG_LIKELIHOOD_CUTOFF) ? -1
                 : static_cast<int>(logp * SCALE);
  }

  // edge-weight calculator
  auto calculate_edge_weights = [&](const blocked_range<int>& r) {
    for (int j = r.begin(); j < r.end(); ++j) {
      if (use_knn) {
        // Find top-k most similar ancestors for j
        vector<pair<int, int>> neighbors;
        for (int i = 0; i < j; ++i) {
          int d = (population[j] ^ population[i]).count(); // Hamming distance
          neighbors.emplace_back(i, d);
        }

        // Retain K nearest by Hamming distance
        partial_sort(neighbors.begin(), neighbors.begin() + min(knn_k, (int)neighbors.size()),
        neighbors.end(), [](auto& a, auto& b) {
          return a.second < b.second;
        });

        for (int k = 0; k < min(knn_k, (int)neighbors.size()); ++k) {
          int i = neighbors[k].first;
          int d = neighbors[k].second;
          int weight = log_lut[d]; // precomputed -log P(d)
          if (weight != -1)
            edge_set.push({i, j, weight});
        }
      } else {
        // Full pairwise comparison: all i < j
        for (int i = 0; i < j; ++i) {
          int d = (population[j] ^ population[i]).count();
          int weight = log_lut[d];
          if (weight != -1)
            edge_set.push({i, j, weight});
        }
      }
    }
  };

  // ------------------------------------------
  // Generate candidate edges in parallel
  // ------------------------------------------

  parallel_for(blocked_range<int>(0, POPULATION_SIZE),
               calculate_edge_weights,
               auto_partitioner());

  // ------------------------------------------
  // Construct Boost Graph from edge set
  // ------------------------------------------

  for (Edge triple; edge_set.try_pop(triple); ) {
    auto [u, v, w] = triple;
    auto [e, inserted] = add_edge(u, v, graph);
    if (inserted) weight_map[e] = w;
  }

  // ------------------------------------------
  // Select root vertex with highest viable degree
  // ------------------------------------------

  auto [vi_start, vi_end] = vertices(graph);
  auto v_max = vi_start;
  size_t max_degree = 0;
  for (auto vi = vi_start; vi != vi_end; ++vi) {
    size_t deg = degree(*vi, graph);
    if (deg >= max_degree && deg <= EXPECTED_MAX_DEGREE + 3) {
      max_degree = deg;
      v_max = vi;
    }
  }

  if (v_max == vi_end) {
    cerr << "ERROR: No suitable root vertex found.\n";
    return 1;
  }

  cerr << "Selected root. Degree = " << max_degree << "\n";

  // ------------------------------------------
  // Run Prim's MST algorithm
  // ------------------------------------------

  auto distance_map = get(vertex_distance, graph);
  auto index_map = get(vertex_index, graph);

  prim_minimum_spanning_tree(graph, *v_max,
                             make_iterator_property_map(predecessors.begin(), index_map),
                             distance_map, weight_map, index_map,
                             default_dijkstra_visitor());

  // ------------------------------------------
  // Output: Parent index per genome (one per line)
  // Genesis gets parent = -1
  // ------------------------------------------

  for (size_t i = 0; i < POPULATION_SIZE; ++i) {
    if (predecessors[i] != i) {
        cout << predecessors[i] << '\n';
    } else {
        cout << "-1\n";
    }
  }
}
