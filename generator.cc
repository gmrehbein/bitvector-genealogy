// --------------------------------------------------------
// File: generator.cc (Modernized to C++23)
// Author: Gregory Rehbein
// --------------------------------------------------------

// C++
#include <fstream>
#include <vector>
#include <numeric>
#include <random>
#include <ranges>

// Boost
#include <boost/dynamic_bitset.hpp>

using namespace std;
using Genome = boost::dynamic_bitset<>;

int main()
{
  constexpr int GENOME_LEN = 10000;
  constexpr int POPULATION_SIZE = 10000;
  std::random_device rd;

  auto mutate = [&](Genome& genome) {
    static mt19937 mutation_gen(rd());
    static bernoulli_distribution coin_toss(0.2);
    for (int i = 0; i < GENOME_LEN; ++i) {
      if (coin_toss(mutation_gen)) genome.flip(i);
    }
  };

  auto create_genesis_block = [&]() -> Genome {
    mt19937 gen(rd());
    bernoulli_distribution coin_toss(0.5);
    Genome genome(GENOME_LEN);
    for (int i = 0; i < GENOME_LEN; ++i) {
      if (coin_toss(gen)) genome.flip(i);
    }

    return genome;
  };

  auto pick_random_parent = [&](int upper_bound_inclusive) -> int {
    static mt19937 parent_gen(rd());
    uniform_int_distribution<> dist(0, upper_bound_inclusive);
    return dist(parent_gen);
  };

  // Generate genesis block
  vector<Genome> population;
  population.reserve(POPULATION_SIZE);
  population.push_back(create_genesis_block());

  // maps children at index i to the index of their parent
  vector<int> child_to_parent(POPULATION_SIZE, -1);

  // generate the population of bitvectors
  for (int index = 1; index < GENOME_LEN; ++index) {
    auto parent_index = pick_random_parent(index - 1);
    child_to_parent[index] = parent_index;
    auto child = population[parent_index];
    mutate(child);
    population.push_back(child);
  }

  // Create permutation of indices
  mt19937 permutation_gen(rd());
  vector<int> population_permutation(GENOME_LEN);
  iota(begin(population_permutation), end(population_permutation), 0);
  ranges::shuffle(population_permutation, permutation_gen);

  {
    // Write population to file
    ofstream output("test.data", ios_base::out | ios_base::trunc);
    for (const auto& index : population_permutation) {
      output << population[index] << '\n';
    }
  }

  {
    // Write parent map in permuted order
    ofstream parent_output("parent.txt", ios_base::out | ios_base::trunc);
    for (const auto& child_index : population_permutation) {
      int parent_id = child_to_parent[child_index];
      auto it = std::ranges::find(population_permutation, parent_id);
      if (it == end(population_permutation)) {
        parent_output << -1 << '\n';
      } else {
        parent_output << std::distance(begin(population_permutation), it) << '\n';
      }
    }
  }
}
