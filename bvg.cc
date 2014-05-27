// ---------------------------------------------------------
// File: bvg.cc
// Author: Gregory Rehbein
// Copyright (C) 2014 Gregory Rehbein <gmrehbein@gmail.com>
//-----------------------------------------------------------

// C++
#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>

// Boost
#include <boost/dynamic_bitset.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/prim_minimum_spanning_tree.hpp>

using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::ios_base;
using std::vector;

using namespace boost;

int main(int argc, char* argv[])
{
  const size_t kGenomeLen = 10000;
  const int expected_value = 2000;
  const uint16_t std_deviation = 40;
  const char* filename = "test.data";
  const uint8_t expected_max_degree = 10;

  typedef adjacency_list<vecS, vecS, undirectedS,
          property<vertex_distance_t, int>,
          property<edge_weight_t, int> > Graph;

  typedef graph_traits<Graph>::vertex_descriptor Vertex;

  const size_t VERTEX_COUNT = kGenomeLen;
  Graph graph(VERTEX_COUNT);

  property_map<Graph, edge_weight_t>::type weight_map = get(edge_weight, graph);
  vector<Vertex> p(num_vertices(graph));

  char genome[kGenomeLen + 1];
  dynamic_bitset<> population[kGenomeLen];

  ifstream ifs(filename, ifstream::in);
  int index = 0;
  while (ifs.getline(genome, kGenomeLen + 1).good()) {
    population[index].resize(kGenomeLen);
    for (size_t i = 0; i < kGenomeLen; ++i) {
      if ('1' == genome[i]) {
        population[index].set(i);
      }
    }
    for (int i = index - 1; i > -1; --i) {
      dynamic_bitset<> hamming(population[index] ^ population[i]);
      int distance = abs(static_cast<int>(hamming.count()) - expected_value);
      if (distance <= 5*std_deviation) {
        graph_traits<Graph>::edge_descriptor e;
        bool inserted;
        boost::tie(e, inserted) = add_edge(index, i, graph);
        weight_map[e] = distance;
      }
    }
    ++index;
  }
  ifs.close();

  vector<int> component(num_vertices(graph));
  int num = connected_components(graph, &component[0]);
  cout << "total number of components = " << num << endl;
 
  // Search for the vertex with the highest degree
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
  
 // cout << "max degree = " << max_degree << endl;

  property_map<Graph, vertex_distance_t>::type distance = get(vertex_distance, graph);
  property_map<Graph, vertex_index_t>::type index_map = get(vertex_index, graph);
  prim_minimum_spanning_tree(graph, *v_max, &p[0], distance, 
                             weight_map, index_map, default_dijkstra_visitor());
  
  for (size_t i = 0; i != p.size(); ++i) {
    if (p[i] != i)
      cout << p[i] << endl;
    else
      cout << -1 << endl;
  }
}
