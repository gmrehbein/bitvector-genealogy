
// File: bvg.cc (Log-Likelihood Edge Weights)
// Author: Gregory Rehbein

#include <cmath>
#include <fstream>
#include <iostream>
#include <ranges>
#include <string>
#include <tuple>
#include <vector>

#include <boost/dynamic_bitset.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/prim_minimum_spanning_tree.hpp>
#include <boost/math/distributions/binomial.hpp>

#include <tbb/global_control.h>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/concurrent_queue.h>

using namespace std;
using namespace boost;
using namespace tbb;

namespace {
using Genome = dynamic_bitset<>;
using Edge = std::tuple<int, int, int>;

constexpr int GENOME_LEN = 10000;

concurrent_queue<Edge> edge_set;

class CalculateEdgeWeights {
  const vector<Genome>& populatin_;

public:
  CalculateEdgeWeights(const vector<Genome>& population) : populatin_(population) {}

  void operator()(const blocked_range<int>& r) const {
    constexpr double MUTATION_PROB = 0.2;
    constexpr double LOG_LIKELIHOOD_CUTOFF = 20.0; // roughly equivalent to 9σ
    constexpr double SCALE = 1000.0;

    static const boost::math::binomial_distribution<> binom(GENOME_LEN, MUTATION_PROB);
    static vector<int> log_lut(GENOME_LEN + 1, -1);

    for (int j = r.end() - 1; j >= r.begin(); --j) {
      for (int i = j - 1; i >= 0; --i) {
        int d = static_cast<int>((populatin_[j] ^ populatin_[i]).count());

        if (log_lut[d] == -1) {
          double logp = -std::log(boost::math::pdf(binom, d));
          if (!std::isfinite(logp) || logp > LOG_LIKELIHOOD_CUTOFF) {
            log_lut[d] = -1;
            continue;
          }
          log_lut[d] = static_cast<int>(logp * SCALE);
        }

        if (log_lut[d] != -1)
          edge_set.push({i, j, log_lut[d]});
      }
    }
  }
};
}

int main() {
  global_control control(global_control::max_allowed_parallelism, std::thread::hardware_concurrency());

  constexpr char filename[] = "test.data";
  constexpr size_t EXPECTED_MAX_DEGREE = 10;
  constexpr size_t VERTEX_COUNT = GENOME_LEN;

  using Graph = adjacency_list<vecS, vecS, undirectedS,
                               property<vertex_distance_t, int>,
                               property<edge_weight_t, int>>;
  using Vertex = graph_traits<Graph>::vertex_descriptor;

  Graph graph(VERTEX_COUNT);
  auto weight_map = get(edge_weight, graph);
  vector<Vertex> predecessors(VERTEX_COUNT);  // Predecessor map

  // Read population from file
  vector<Genome> population(VERTEX_COUNT, Genome(GENOME_LEN));
  ifstream ifs(filename);
  string line;
  size_t index = 0;

  while (getline(ifs, line) && index < VERTEX_COUNT) {
    for (size_t i = 0; i < GENOME_LEN; ++i)
      if (line[i] == '1') population[index].set(i);
    ++index;
  }

  // Compute edges
  parallel_for(blocked_range<int>(0, static_cast<int>(VERTEX_COUNT)),
               CalculateEdgeWeights(population),
               auto_partitioner());

  for (Edge triple; edge_set.try_pop(triple); ) {
    auto [u, v, w] = triple;
    auto [e, inserted] = add_edge(u, v, graph);
    if (inserted) weight_map[e] = w;
  }


  // Choose root vertex with max degree ≤ threshold
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
    cerr << "Error: No suitable root vertex found.\n";
    return 1;
  }

  cerr << "Selected root with degree " << max_degree << ". Running Prim's...\n";

  auto distance_map = get(vertex_distance, graph);
  auto index_map = get(vertex_index, graph);

  prim_minimum_spanning_tree(graph, *v_max,
                             make_iterator_property_map(predecessors.begin(), index_map),
                             distance_map, weight_map, index_map,
                             default_dijkstra_visitor());

  cerr << "MST complete. Outputting...\n";

  for (size_t i = 0; i < predecessors.size(); ++i) {
    if (predecessors[i] != i) {
        cout << predecessors[i] << '\n';
    } else {
        cout << "-1\n";
    }
  }
}
