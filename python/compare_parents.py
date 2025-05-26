#!/usr/bin/env python3

import sys

def read_parents(path):
    with open(path) as f:
        return [int(line.strip()) for line in f if line.strip()]

def compare_parents(true_path, reconstructed_path):
    true_parents = read_parents(true_path)
    inferred_parents = read_parents(reconstructed_path)

    if len(true_parents) != len(inferred_parents):
        raise ValueError("Files have different lengths")

    total = len(true_parents)
    mismatches = 0

    for i, (true, inferred) in enumerate(zip(true_parents, inferred_parents)):
        if true == -1:
            continue  # skip genesis
        if true != inferred:
            mismatches += 1

    correct = total - mismatches - 1  # -1 for genesis
    accuracy = correct / (total - 1)
    error_rate = mismatches / (total - 1)

    print(f"Total nodes:         {total}")
    print(f"Non-genesis nodes:   {total - 1}")
    print(f"Correct links:       {correct}")
    print(f"Mismatched links:    {mismatches}")
    print(f"Accuracy:            {accuracy:.5f}")
    print(f"Misclassification:   {error_rate:.5f}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: compare_parents.py parent.txt out.txt")
        sys.exit(1)

    compare_parents(sys.argv[1], sys.argv[2])
