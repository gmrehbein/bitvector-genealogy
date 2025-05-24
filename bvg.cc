// ---------------------------------------------------------
// File: bvg.cc (Modernized to C++23)
// Author: Gregory Rehbein
// ---------------------------------------------------------

#include <fstream>
#include <tuple>
#include <vector>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <ranges>
#include <algorithm>
#include <iostream>

// Boost
#include <boost/dynamic_bitset.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/prim_minimum_spanning_tree.hpp>

// TBB
#include <tbb/global_control.h>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/concurrent_queue.h>

using namespace std;
using namespace boost;
using namespace tbb;

namespace
{
using Genome = dynamic_bitset<>;
using Edge = std::tuple<int, int, int>;

concurrent_queue<Edge> edge_set;

class CalculateEdgeWeights
{
  const vector<Genome>& population_;

public:
  CalculateEdgeWeights(const vector<Genome>& population) : population_(population) {}

  void operator()(const blocked_range<int>& r) const
  {
    constexpr int expected_value = 2000;
    constexpr int std_deviation = 40;

    for (int j = r.end() - 1; j >= r.begin(); --j) {
      for (int i = j - 1; i >= 0; --i) {
        Genome hamming = population_[j] ^ population_[i];
        int distance = abs(static_cast<int>(hamming.count()) - expected_value);
        if (distance <= 5 * std_deviation)
          edge_set.push({i, j, distance});
      }
    }
  }
};
}

int main()
{
  global_control control(global_control::max_allowed_parallelism, std::thread::hardware_concurrency());

  constexpr size_t GENOME_SIZE = 10000;
  constexpr size_t EXPECTED_MAX_DEGREE = 10;
  constexpr size_t VERTEX_COUNT = GENOME_SIZE;
  constexpr char input_filename[] = "test.data";

  using Graph = adjacency_list<vecS, vecS, undirectedS,
        property<vertex_distance_t, int>,
        property<edge_weight_t, int>>;

  using Vertex = graph_traits<Graph>::vertex_descriptor;

  Graph graph(VERTEX_COUNT);
  auto weight_map = get(edge_weight, graph);
  vector<Vertex> predecessors(VERTEX_COUNT);  // Predecessor map

  // Read genome data
  vector<Genome> population(VERTEX_COUNT, Genome(GENOME_SIZE));
  ifstream ifs(input_filename);
  string line;
  size_t index = 0;

  while (getline(ifs, line) && index < VERTEX_COUNT) {
    for (size_t i = 0; i < GENOME_SIZE; ++i) {
      if (line[i] == '1') population[index].set(i);
    }
    ++index;
  }

  parallel_for(blocked_range<int>(0, static_cast<int>(VERTEX_COUNT)),
               CalculateEdgeWeights(population),
               auto_partitioner());

  for (Edge triple; edge_set.try_pop(triple); ) {
    auto [u, v, w] = triple;
    auto [e, inserted] = add_edge(u, v, graph);
    if (inserted) weight_map[e] = w;
  }

  // Find root vertex with highest degree (within constraint)
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

  auto distance_map = get(vertex_distance, graph);
  auto index_map = get(vertex_index, graph);

  prim_minimum_spanning_tree(graph, *v_max,
                             make_iterator_property_map(predecessors.begin(), index_map),
                             distance_map, weight_map, index_map,
                             default_dijkstra_visitor());

  for (size_t i = 0; i < predecessors.size(); ++i) {
    printf("%zu\n", (predecessors[i] != i) ? predecessors[i] : static_cast<size_t>(-1));
  }
}
