#!/usr/bin/env python3

import numpy as np

GENOME_LEN = 10000
POPULATION_SIZE = 10000
MUTATION_PROB = 0.2
GENESIS_BIT_PROB = 0.5

def create_genesis_block(rng):
    return (rng.random(GENOME_LEN) < GENESIS_BIT_PROB).astype(bool)

def mutate(genome, rng):
    mutation_mask = rng.random(GENOME_LEN) < MUTATION_PROB
    genome ^= mutation_mask.astype(bool)
    return genome


def main():
    rng = np.random.default_rng()

    # Initialize population list and genealogy
    population = []
    population.append(create_genesis_block(rng))
    child_to_parent = [-1] * POPULATION_SIZE

    # Main loop to generate children
    for index in range(1, POPULATION_SIZE):
        parent_index = rng.integers(0, index)
        child_to_parent[index] = parent_index
        parent = population[parent_index];
        child = mutate(parent.copy(), rng)
        population.append(child)

    # Generate a random permutation of population indices
    population_permutation = rng.permutation(POPULATION_SIZE)

    # Write genomes to file
    with open("test.data", "w", newline="\n") as f:
        for i in population_permutation:
            genome = population[i].astype(bool)
            genome_str = ''.join('1' if b else '0' for b in genome)
            assert len(genome_str) == 10000
            f.write(genome_str + "\n")


    # Write permuted parent indices
    with open("parent.txt", "w") as f:
        for i in population_permutation:
            parent_id = child_to_parent[i]
            if parent_id == -1:
                f.write("-1\n")
            else:
                idx = np.where(population_permutation == parent_id)[0]
                f.write(f"{idx[0] if idx.size > 0 else -1}\n")

if __name__ == "__main__":
    main()

