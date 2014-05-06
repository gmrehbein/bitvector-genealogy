// --------------------------------------------------------
// File: generator.cc
// Author: Gregory Rehbein
//
// Copyright (C) 2014 Gregory Rehbein <gmrehbein@gmail.com>
//----------------------------------------------------------

// C
#include <cstdlib>

// C++
#include <vector>
#include <iostream>
#include <fstream>

// Boost
#include <boost/dynamic_bitset.hpp>
#include <boost/random.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/bernoulli_distribution.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/random_device.hpp>
#include <boost/unordered_map.hpp>

using std::vector;
using std::random_shuffle;
using std::find;
using std::ofstream;
using std::ios_base;
using std::cout;
using std::endl;

using namespace boost;

int parent(int index, mt19937& rand_gen)
{
  uniform_int<> distribution(0, index);
  variate_generator<mt19937&, uniform_int<> > parent_generator(rand_gen, distribution);
  return parent_generator();
}

int main()
{
  const int GENOME_LEN = 10000;
  random_device rd;

  mt19937 mutation_gen(rd());
  mt19937 parent_gen(rd());
  srand(rd()); // for random_shuffle

  bernoulli_distribution<> coin_toss(0.2);

  // generate genesis block
  vector<dynamic_bitset<> > population;
  dynamic_bitset<> genesis(GENOME_LEN);
  for (int i = 0; i < GENOME_LEN; ++i) {
    if (1 == coin_toss(mutation_gen)) genesis.set(i);
  }
  population.push_back(genesis);

  unordered_map<int, int> child_parent_map;
  child_parent_map[0] = -1;

  for (int index = 1; index < GENOME_LEN; ++index) {
    int parent_index = parent(index - 1, parent_gen);
    child_parent_map[index] = parent_index;
    dynamic_bitset<> child(population[parent_index]);
    for (int j = 0; j < GENOME_LEN; ++j) {
      if (1 == coin_toss(mutation_gen)) child.flip(j);
    }
    population.push_back(child);
  }

  int population_permutation[GENOME_LEN];
  for (int i = 0; i < GENOME_LEN; ++i) {
    population_permutation[i] =i;
  }

  random_shuffle(population_permutation, population_permutation + GENOME_LEN);

  ofstream output("test.data", ios_base::out | ios_base::trunc);
  for (int i = 0; i < GENOME_LEN; ++i) {
    output << population.at(population_permutation[i]) << endl;
  }
  output.close();

  // locate parent of vertex in permuted order
  ofstream parent_output("parent.txt", ios_base::out | ios_base::trunc);
  for (int i = 0; i < GENOME_LEN; ++i) {
    int* parent_ptr = find(population_permutation, population_permutation + GENOME_LEN,
                           child_parent_map[population_permutation[i]]);
    ptrdiff_t parent_index = parent_ptr - population_permutation;
    parent_index == GENOME_LEN ? parent_output << -1 << endl :
        parent_output << parent_index << endl;
  }

  parent_output.close();
}
