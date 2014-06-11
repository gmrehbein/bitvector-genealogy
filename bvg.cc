// ---------------------------------------------------------
// File: bvg.cc
// Author: Gregory Rehbein
// Copyright (C) 2014 Gregory Rehbein <gmrehbein@gmail.com>
//-----------------------------------------------------------

// C++
#include <fstream>
#include <tuple>

// C
#include <cstdio>

// Boost
#include <boost/dynamic_bitset.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/prim_minimum_spanning_tree.hpp>

// TBB
#include <tbb/task_scheduler_init.h>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/concurrent_queue.h>

using std::vector;
using std::ifstream;

using namespace boost;
using namespace tbb;

namespace {
concurrent_queue<std::tuple<int,int,int> > edge_set;

class CalculateEdgeWeights
{
  const dynamic_bitset<>* m_population;

public:
  CalculateEdgeWeights(const dynamic_bitset<> population[]): m_population(population) {}
  void operator() (const blocked_range<int>& r) const
  {
    const int expected_value = 2000;
    const int std_deviation = 40;
    for (int j = r.end() - 1; j >= r.begin(); --j) {
      for (int i = j - 1; i > -1; --i) {
        dynamic_bitset<> hamming(m_population[j] ^ m_population[i]);
		int distance = abs(hamming.count() - expected_value);
        if (distance <= 5*std_deviation)
			edge_set.push({i,j,distance});
      }
    }
  }
};
}

int main(int argc, char* argv[])
{
  task_scheduler_init init;

  const size_t kGenomeLen = 10000;
  const char* filename = "test.data";
  const size_t expected_max_degree = 10;

  typedef adjacency_list<vecS, vecS, undirectedS,
          property<vertex_distance_t, int>,
          property<edge_weight_t, int> > Graph;

  typedef graph_traits<Graph>::vertex_descriptor Vertex;

  const size_t VERTEX_COUNT = kGenomeLen;
  Graph graph(VERTEX_COUNT);

  property_map<Graph, edge_weight_t>::type weight_map = get(edge_weight, graph);
  vector<Vertex> p(num_vertices(graph)); // predecessor map

  char genome[kGenomeLen + 1];
  dynamic_bitset<> population[VERTEX_COUNT];

  ifstream ifs(filename, ifstream::in);
  int index = 0;
  while (ifs.getline(genome, kGenomeLen + 1).good()) {
    population[index].resize(kGenomeLen);
    for (size_t i = 0; i < kGenomeLen; ++i) {
      if ('1' == genome[i]) {
        population[index].set(i);
      }
    }
    ++index;
  }
  ifs.close();

  parallel_for(blocked_range<int> (0, VERTEX_COUNT),
               CalculateEdgeWeights(population),
               auto_partitioner());
  
  std::tuple<int,int,int> triple;
  while (edge_set.try_pop(triple)) {
	  graph_traits<Graph>::edge_descriptor e;
	  bool inserted;
      tie(e,inserted) = add_edge(std::get<0>(triple), std::get<1>(triple), graph);
      weight_map[e] = std::get<2>(triple); 	
  }

  // Search for the vertex *v_max with the highest degree.
  // Nominate this as the root vertex in Prim's MST algorithm
  graph_traits<Graph>::vertex_iterator vi_curr, vi_end, v_max;
  graph_traits<Graph>::degree_size_type max_degree = 0;
  graph_traits<Graph>::degree_size_type curr_degree = 0;
  for (tie(vi_curr, vi_end) = vertices(graph); vi_curr != vi_end; ++vi_curr) {
    curr_degree = degree(*vi_curr, graph);
    if (curr_degree >= max_degree && curr_degree <= expected_max_degree + 3) {
      max_degree = curr_degree;
      v_max = vi_curr;
    }
  }

  // Compute the minimum spanning tree with root vertex *v_max
  property_map<Graph, vertex_distance_t>::type distance = get(vertex_distance, graph);
  property_map<Graph, vertex_index_t>::type index_map = get(vertex_index, graph);
  prim_minimum_spanning_tree(graph, *v_max, &p[0], distance,
                             weight_map, index_map, default_dijkstra_visitor());

  // Output result
  for (size_t i = 0; i != p.size(); ++i) {
    if (p[i] != i)
      printf("%zu\n", p[i]);
    else
      printf("-1\n");
  }
}
